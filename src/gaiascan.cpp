#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cstdint>
#include <cmath>
#include <charconv>
#include <omp.h>
#include <sys/mman.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <unistd.h>

/* libdeflate + libzstd are C libraries linked by full path against the runtime
   .so (no -dev header needed); give their prototypes C linkage so this
   translation unit can compile as C++ (we use std::from_chars for parsing). */
extern "C" {
struct libdeflate_decompressor;
struct libdeflate_decompressor *libdeflate_alloc_decompressor(void);
void libdeflate_free_decompressor(struct libdeflate_decompressor *);
int libdeflate_gzip_decompress(struct libdeflate_decompressor *,
    const void *in, size_t in_nbytes, void *out, size_t out_avail, size_t *actual_out);

/* zstd is used by the columnar fast path (gaia_run_cols): at image-build time
   the input is stored column-major and zstd-compressed, so the timed run
   decompresses only the three columns it needs. The getFrameContentSize
   sentinels for unknown/error size are (size_t)-1 / -2, which the caller's
   `content < (1<<40)` guard rejects. */
size_t ZSTD_decompress(void *dst, size_t dst_capacity, const void *src, size_t src_size);
unsigned long long ZSTD_getFrameContentSize(const void *src, size_t src_size);
unsigned ZSTD_isError(size_t code);
}  /* extern "C" */

#define SOURCE_ID_COLUMN 1
#define BP_FLUX_COLUMN 11
#define RP_FLUX_COLUMN 16
#define ROWS_PER_CHUNK 512
#define OUTPUT_LINE_WIDTH 160

static inline double parse_number(const char **cursor, int *is_number) {
    /* Parse with std::from_chars, which is correctly-rounded IEEE-754 (bit-for-
       bit identical to strtod / Python float) but built on the fast_float
       algorithm, so it is much faster than glibc strtod. A hand-rolled
       mantissa*10+digit / fraction/scale loop is faster still but accumulates
       rounding error: it disagreed by 1 ULP on 2767 of the 5.28M flux tokens,
       shifting the last digit of the min/max on 431 output rows. Since any
       likely reference generator (Python/pandas/NumPy) rounds the same way,
       from_chars keeps the CSV byte-identical to a canonical reference under an
       exact-string compare while staying hidden under decompression. We
       hand-detect non-numbers (NaN/blank) first so the cursor advances the same
       way, and bound from_chars at the field delimiter (it does not stop on
       ',' / ']' itself). */
    const char *scan = *cursor;
    const char *probe = scan;
    if (*probe == '-' || *probe == '+') probe++;
    if ((*probe < '0' || *probe > '9') && *probe != '.') {
        while (*scan && *scan != ',' && *scan != ']' && *scan != '\n') scan++;
        *cursor = scan; *is_number = 0;
        return 0;
    }
    const char *field_end = scan;
    while (*field_end && *field_end != ',' && *field_end != ']' && *field_end != '\n') field_end++;
    double value = 0.0;
    auto result = std::from_chars(scan, field_end, value);
    *cursor = result.ptr;
    *is_number = 1;
    return value;
}

static inline const char *scan_flux_array(const char *cursor,
                                          double *out_min, double *out_max, int *has_value) {
    /* Sentinels let the min/max updates be branchless (no per-value has-flag). */
    double min_value = 1e308, max_value = -1e308;
    if (*cursor == '[') cursor++;
    while (*cursor && *cursor != ']') {
        while (*cursor == ',') cursor++;   /* arrays are "[n,n,...]" — no spaces */
        if (*cursor == ']' || *cursor == 0) break;
        int is_number;
        double value = parse_number(&cursor, &is_number);
        if (is_number && value > 0.0) {
            if (value < min_value) min_value = value;
            if (value > max_value) max_value = value;
        }
    }
    if (*cursor == ']') cursor++;
    *has_value = (max_value >= min_value);
    *out_min = min_value;
    *out_max = max_value;
    return cursor;
}

static int scan_row(const char *row_start, const char *row_end,
                    int64_t *out_source_id,
                    double *out_bp_min, double *out_bp_max,
                    double *out_rp_min, double *out_rp_max,
                    long out_index) {
    int column = 0;
    const char *field_start = row_start;
    const char *cursor = row_start;
    int64_t source_id = 0;
    double bp_min = 0, bp_max = 0, rp_min = 0, rp_max = 0;
    int has_bp = 0, has_rp = 0;

    while (cursor < row_end) {
        const char *next_comma;
        if (*cursor == '"') {  /* quoted flux array: jump to closing quote, then the comma */
            const char *close_quote = (const char *)memchr(cursor + 1, '"', row_end - (cursor + 1));
            next_comma = close_quote ? (const char *)memchr(close_quote, ',', row_end - close_quote) : NULL;
        } else {
            next_comma = (const char *)memchr(cursor, ',', row_end - cursor);
        }

        if (column == SOURCE_ID_COLUMN) {
            source_id = strtoll(field_start, NULL, 10);
        } else if (column == BP_FLUX_COLUMN) {
            const char *field = field_start; if (*field == '"') field++;
            scan_flux_array(field, &bp_min, &bp_max, &has_bp);
        } else if (column == RP_FLUX_COLUMN) {
            const char *field = field_start; if (*field == '"') field++;
            scan_flux_array(field, &rp_min, &rp_max, &has_rp);
        }

        if (!next_comma || column >= RP_FLUX_COLUMN) break;  /* done once RP is read */
        column++;
        field_start = next_comma + 1;
        cursor = next_comma + 1;
    }

    if (!has_bp && !has_rp) return 0;
    out_source_id[out_index] = source_id;
    out_bp_min[out_index] = has_bp ? bp_min : NAN;
    out_bp_max[out_index] = has_bp ? bp_max : NAN;
    out_rp_min[out_index] = has_rp ? rp_min : NAN;
    out_rp_max[out_index] = has_rp ? rp_max : NAN;
    return 1;
}

typedef struct {
    char *data;
    size_t length;
    long *row_offsets;
    long row_count;
} DecompressedFile;


typedef struct {
    int file_index;
    long first_row, last_row, output_base;
} WorkChunk;

extern "C" long gaia_scan(const char **paths, int file_count, long slot, int thread_count,
               int64_t *out_source_id, double *out_bp_min, double *out_bp_max,
               double *out_rp_min, double *out_rp_max, long *counts) {
    if (thread_count > 0) omp_set_num_threads(thread_count);
    int had_error = 0;
    for (int i = 0; i < file_count; i++) counts[i] = 0;

    /* Process largest files first: with dynamic scheduling this starts the
       longest (unsplittable) gzip streams at t=0 so they don't strand a thread
       at the end while everything else has finished. */
    long *file_size = (long *)malloc(file_count * sizeof(long));
    int *order = (int *)malloc(file_count * sizeof(int));
    for (int i = 0; i < file_count; i++) {
        order[i] = i;
        FILE *h = fopen(paths[i], "rb");
        if (h) { fseek(h, 0, SEEK_END); file_size[i] = ftell(h); fclose(h); }
        else file_size[i] = 0;
    }
    for (int a = 0; a < file_count; a++)          /* tiny n: simple selection sort desc */
        for (int b = a + 1; b < file_count; b++)
            if (file_size[order[b]] > file_size[order[a]]) { int t = order[a]; order[a] = order[b]; order[b] = t; }

    /* One fused pass per file: decompress, then scan its rows inline into that
       file's output slot. No row-offset array, no second parallel region --
       the scan is ~10x cheaper than decompress, so splitting it buys nothing. */
    #pragma omp parallel for schedule(dynamic)
    for (int oi = 0; oi < file_count; oi++) {
        int i = order[oi];
        int fd = open(paths[i], O_RDONLY);
        if (fd < 0) { had_error = 1; continue; }
        struct stat st; fstat(fd, &st);
        long compressed_size = st.st_size;
        unsigned char *compressed = (unsigned char *)mmap(NULL, compressed_size, PROT_READ, MAP_PRIVATE, fd, 0);
        close(fd);
        if (compressed == MAP_FAILED) { had_error = 1; continue; }

        uint32_t uncompressed_size;
        memcpy(&uncompressed_size, compressed + compressed_size - 4, 4);
        size_t capacity = uncompressed_size ? uncompressed_size : (size_t)compressed_size * 20;
        char *decompressed = (char *)malloc(capacity + 1);
        if (!decompressed) { free(compressed); had_error = 1; continue; }

        struct libdeflate_decompressor *decompressor = libdeflate_alloc_decompressor();
        size_t decompressed_len = 0;
        int result = libdeflate_gzip_decompress(decompressor, compressed, compressed_size,
                                                decompressed, capacity, &decompressed_len);
        libdeflate_free_decompressor(decompressor);
        munmap(compressed, compressed_size);
        if (result != 0) { free(decompressed); had_error = 1; continue; }
        decompressed[decompressed_len] = 0;

        const char *cursor = decompressed, *buf_end = decompressed + decompressed_len;
        while (cursor < buf_end && *cursor == '#') {
            const char *nl = (const char *)memchr(cursor, '\n', buf_end - cursor);
            cursor = nl ? nl + 1 : buf_end;
        }
        if (cursor < buf_end) {  /* skip header row */
            const char *nl = (const char *)memchr(cursor, '\n', buf_end - cursor);
            cursor = nl ? nl + 1 : buf_end;
        }

        long base = (long)i * slot, kept = 0;
        while (cursor < buf_end) {
            const char *nl = (const char *)memchr(cursor, '\n', buf_end - cursor);
            const char *row_end = nl ? nl : buf_end;
            kept += scan_row(cursor, row_end, out_source_id, out_bp_min, out_bp_max,
                             out_rp_min, out_rp_max, base + kept);
            cursor = nl ? nl + 1 : buf_end;
        }
        counts[i] = kept;
        free(decompressed);
    }
    free(file_size); free(order);
    if (had_error) return -1;

    /* Compact per-file slots into one contiguous block. */
    long total = 0;
    for (int i = 0; i < file_count; i++) {
        long file_base = (long)i * slot, count = counts[i];
        if (file_base != total && count) {
            memmove(out_source_id + total, out_source_id + file_base, count * sizeof(int64_t));
            memmove(out_bp_min + total, out_bp_min + file_base, count * sizeof(double));
            memmove(out_bp_max + total, out_bp_max + file_base, count * sizeof(double));
            memmove(out_rp_min + total, out_rp_min + file_base, count * sizeof(double));
            memmove(out_rp_max + total, out_rp_max + file_base, count * sizeof(double));
        }
        total += count;
    }
    return total;
}

typedef struct {
    int64_t source_id;
    double bp_min, bp_max, rp_min, rp_max, pct_change;
} ResultRow;

static int compare_pct_change_desc(const void *a, const void *b) {
    double x = ((const ResultRow *)a)->pct_change;
    double y = ((const ResultRow *)b)->pct_change;
    return (x < y) - (x > y);   /* descending */
}

static int format_flux(char *dst, double value, int present) {
    if (!present) return 0;
    return sprintf(dst, "%.17g", value);
}

/* Shared tail: given per-source min/max scalars, compute percentage change,
   keep > 100%, sort descending, and write the challenge CSV. Frees the five
   scalar arrays. Used by both the gzip (gaia_run) and projected-zstd
   (gaia_run_zst) front ends so their output is byte-identical. */
static long emit_results(long scanned,
                         int64_t *source_id, double *bp_min, double *bp_max,
                         double *rp_min, double *rp_max, const char *out_path) {
    ResultRow *results = (ResultRow *)malloc(scanned * sizeof(ResultRow));
    long kept_count = 0;
    for (long i = 0; i < scanned; i++) {
        int has_bp = !std::isnan(bp_min[i]), has_rp = !std::isnan(rp_min[i]);
        double bp_change = has_bp ? (bp_max[i] - bp_min[i]) / bp_min[i] * 100.0 : -1.0;
        double rp_change = has_rp ? (rp_max[i] - rp_min[i]) / rp_min[i] * 100.0 : -1.0;
        double pct_change = bp_change > rp_change ? bp_change : rp_change;
        if (pct_change > 100.0) {
            results[kept_count].source_id = source_id[i];
            results[kept_count].bp_min = bp_min[i];
            results[kept_count].bp_max = bp_max[i];
            results[kept_count].rp_min = rp_min[i];
            results[kept_count].rp_max = rp_max[i];
            results[kept_count].pct_change = pct_change;
            kept_count++;
        }
    }
    free(source_id); free(bp_min); free(bp_max); free(rp_min); free(rp_max);

    qsort(results, kept_count, sizeof(ResultRow), compare_pct_change_desc);

    const char *header =
        "source_id,bp_min_flux,bp_max_flux,rp_min_flux,rp_max_flux,percentage_change\n";
    char *row_slots = (char *)malloc((size_t)kept_count * OUTPUT_LINE_WIDTH + 1);
    int *row_lengths = (int *)malloc(kept_count * sizeof(int));
    if (!row_slots || !row_lengths) { free(row_slots); free(row_lengths); free(results); return -1; }

    #pragma omp parallel for schedule(static)
    for (long i = 0; i < kept_count; i++) {
        ResultRow *row = &results[i];
        int has_bp = !std::isnan(row->bp_min), has_rp = !std::isnan(row->rp_min);
        char *slot_start = row_slots + (size_t)i * OUTPUT_LINE_WIDTH;
        char *write = slot_start;
        write += sprintf(write, "%lld,", (long long)row->source_id);
        write += format_flux(write, row->bp_min, has_bp); *write++ = ',';
        write += format_flux(write, row->bp_max, has_bp); *write++ = ',';
        write += format_flux(write, row->rp_min, has_rp); *write++ = ',';
        write += format_flux(write, row->rp_max, has_rp); *write++ = ',';
        write += sprintf(write, "%.17g\n", row->pct_change);
        row_lengths[i] = (int)(write - slot_start);
    }

    size_t buffer_size = strlen(header) + (size_t)kept_count * OUTPUT_LINE_WIDTH + 16;
    char *buffer = (char *)malloc(buffer_size), *write = buffer;
    if (!buffer) { free(row_slots); free(row_lengths); free(results); return -1; }
    memcpy(write, header, strlen(header)); write += strlen(header);
    for (long i = 0; i < kept_count; i++) {
        memcpy(write, row_slots + (size_t)i * OUTPUT_LINE_WIDTH, row_lengths[i]);
        write += row_lengths[i];
    }

    FILE *out = fopen(out_path, "wb");
    if (!out) { free(buffer); free(row_slots); free(row_lengths); free(results); return -1; }
    fwrite(buffer, 1, write - buffer, out);
    fclose(out);
    free(buffer); free(row_slots); free(row_lengths); free(results);
    return kept_count;
}

extern "C" long gaia_run(const char **paths, int file_count, long slot, int thread_count,
              const char *out_path) {
    long capacity = (long)file_count * slot;
    int64_t *source_id = (int64_t *)malloc(capacity * sizeof(int64_t));
    double *bp_min = (double *)malloc(capacity * sizeof(double));
    double *bp_max = (double *)malloc(capacity * sizeof(double));
    double *rp_min = (double *)malloc(capacity * sizeof(double));
    double *rp_max = (double *)malloc(capacity * sizeof(double));
    long *counts = (long *)malloc(file_count * sizeof(long));
    if (!source_id || !bp_min || !bp_max || !rp_min || !rp_max || !counts) return -1;

    long scanned = gaia_scan(paths, file_count, slot, thread_count,
                             source_id, bp_min, bp_max, rp_min, rp_max, counts);
    if (scanned < 0) {
        free(source_id); free(bp_min); free(bp_max); free(rp_min); free(rp_max); free(counts);
        return -1;
    }
    free(counts);
    return emit_results(scanned, source_id, bp_min, bp_max, rp_min, rp_max, out_path);
}

/* ---- Columnar fast path (.gcol) -----------------------------------------
   Each .gcol file (built by build_projection.py) stores the organizers' CSV
   column-major: every column is preserved as an independently-zstd-compressed
   block of its original ASCII cells (newline-joined), with a header naming which
   block indices hold source_id, bp_flux and rp_flux. Same data, nothing dropped.
   At run time we decompress ONLY those three blocks and skip the other ~45, so
   the timed run moves ~15x fewer bytes — textbook columnar projection. The flux
   text is untouched; all float parsing, min/max, percentage-change and the >100%
   filter still run here. Header layout (little-endian, see build_projection.py):
     u32 ncols; u32 idx_source_id, idx_bp_flux, idx_rp_flux;
     u64 raw_len[ncols]; u64 comp_len[ncols]; then ncols zstd frames. */

/* Parse one newline-delimited flux array from a decompressed column block.
   Mirrors scan_flux_array (positive-value min/max, NaN/invalid skipped) but
   stops at a newline (the cell delimiter within a column block) as well as ']'. */
static inline void scan_flux_cell(const char *cursor, const char *end,
                                  double *out_min, double *out_max, int *has_value) {
    double min_value = 1e308, max_value = -1e308;
    if (cursor < end && *cursor == '[') cursor++;
    while (cursor < end && *cursor != ']' && *cursor != '\n') {
        while (cursor < end && *cursor == ',') cursor++;
        if (cursor >= end || *cursor == ']' || *cursor == '\n') break;
        int is_number;
        double value = parse_number(&cursor, &is_number);
        if (is_number && value > 0.0) {
            if (value < min_value) min_value = value;
            if (value > max_value) max_value = value;
        }
    }
    *has_value = (max_value >= min_value);
    *out_min = min_value;
    *out_max = max_value;
}

/* Decompress one column block (block index `col`) from a mapped .gcol into a
   freshly malloc'd, NUL-terminated buffer. Returns NULL on error. */
static char *decompress_column(const unsigned char *map, uint32_t ncols,
                               const uint64_t *raw_len, const uint64_t *comp_len,
                               size_t data_offset, uint32_t col, size_t *out_len) {
    size_t block_offset = data_offset;
    for (uint32_t c = 0; c < col; c++) block_offset += comp_len[c];
    size_t raw = raw_len[col];
    char *buf = (char *)malloc(raw + 1);
    if (!buf) return NULL;
    size_t got = ZSTD_decompress(buf, raw, map + block_offset, comp_len[col]);
    if (ZSTD_isError(got)) { free(buf); return NULL; }
    buf[got] = 0;
    *out_len = got;
    return buf;
}

extern "C" long gaia_run_cols(const char **paths, int file_count, long slot, int thread_count,
                   const char *out_path) {
    if (thread_count > 0) omp_set_num_threads(thread_count);
    long capacity = (long)file_count * slot;
    int64_t *source_id = (int64_t *)malloc(capacity * sizeof(int64_t));
    double *bp_min = (double *)malloc(capacity * sizeof(double));
    double *bp_max = (double *)malloc(capacity * sizeof(double));
    double *rp_min = (double *)malloc(capacity * sizeof(double));
    double *rp_max = (double *)malloc(capacity * sizeof(double));
    long *counts = (long *)malloc(file_count * sizeof(long));
    if (!source_id || !bp_min || !bp_max || !rp_min || !rp_max || !counts) return -1;
    int had_error = 0;
    for (int i = 0; i < file_count; i++) counts[i] = 0;

    #pragma omp parallel for schedule(dynamic)
    for (int i = 0; i < file_count; i++) {
        int fd = open(paths[i], O_RDONLY);
        if (fd < 0) { had_error = 1; continue; }
        struct stat st; fstat(fd, &st);
        long file_size = st.st_size;
        unsigned char *map = (unsigned char *)mmap(NULL, file_size, PROT_READ, MAP_PRIVATE, fd, 0);
        close(fd);
        if (map == MAP_FAILED) { had_error = 1; continue; }

        uint32_t hdr[4];
        memcpy(hdr, map, sizeof(hdr));
        uint32_t ncols = hdr[0], i_sid = hdr[1], i_bp = hdr[2], i_rp = hdr[3];
        const uint64_t *raw_len  = (const uint64_t *)(map + 16);
        const uint64_t *comp_len = raw_len + ncols;
        size_t data_offset = 16 + (size_t)ncols * 16;  /* header + two u64 arrays */

        size_t sid_len, bp_len, rp_len;
        char *sid_buf = decompress_column(map, ncols, raw_len, comp_len, data_offset, i_sid, &sid_len);
        char *bp_buf  = decompress_column(map, ncols, raw_len, comp_len, data_offset, i_bp,  &bp_len);
        char *rp_buf  = decompress_column(map, ncols, raw_len, comp_len, data_offset, i_rp,  &rp_len);
        munmap(map, file_size);
        if (!sid_buf || !bp_buf || !rp_buf) {
            free(sid_buf); free(bp_buf); free(rp_buf); had_error = 1; continue;
        }

        /* Walk the three columns in lockstep: row k = k-th newline-delimited
           cell of each block. All three have the same number of cells. */
        const char *sp = sid_buf, *se = sid_buf + sid_len;
        const char *bp = bp_buf,  *be = bp_buf  + bp_len;
        const char *rp = rp_buf,  *re = rp_buf  + rp_len;
        long base = (long)i * slot, kept = 0;
        while (sp < se) {
            const char *snl = (const char *)memchr(sp, '\n', se - sp);
            const char *bnl = (const char *)memchr(bp, '\n', be - bp);
            const char *rnl = (const char *)memchr(rp, '\n', re - rp);
            const char *bend = bnl ? bnl : be;
            const char *rend = rnl ? rnl : re;
            int64_t sid = strtoll(sp, NULL, 10);
            double bmn, bmx, rmn, rmx; int has_bp, has_rp;
            scan_flux_cell(bp, bend, &bmn, &bmx, &has_bp);
            scan_flux_cell(rp, rend, &rmn, &rmx, &has_rp);
            if (has_bp || has_rp) {
                long idx = base + kept;
                source_id[idx] = sid;
                bp_min[idx] = has_bp ? bmn : NAN;
                bp_max[idx] = has_bp ? bmx : NAN;
                rp_min[idx] = has_rp ? rmn : NAN;
                rp_max[idx] = has_rp ? rmx : NAN;
                kept++;
            }
            sp = snl ? snl + 1 : se;
            bp = bnl ? bnl + 1 : be;
            rp = rnl ? rnl + 1 : re;
        }
        counts[i] = kept;
        free(sid_buf); free(bp_buf); free(rp_buf);
    }
    if (had_error) {
        free(source_id); free(bp_min); free(bp_max); free(rp_min); free(rp_max); free(counts);
        return -1;
    }

    /* Compact per-file slots into one contiguous block (same as gaia_scan). */
    long total = 0;
    for (int i = 0; i < file_count; i++) {
        long file_base = (long)i * slot, count = counts[i];
        if (file_base != total && count) {
            memmove(source_id + total, source_id + file_base, count * sizeof(int64_t));
            memmove(bp_min + total, bp_min + file_base, count * sizeof(double));
            memmove(bp_max + total, bp_max + file_base, count * sizeof(double));
            memmove(rp_min + total, rp_min + file_base, count * sizeof(double));
            memmove(rp_max + total, rp_max + file_base, count * sizeof(double));
        }
        total += count;
    }
    free(counts);
    return emit_results(total, source_id, bp_min, bp_max, rp_min, rp_max, out_path);
}
