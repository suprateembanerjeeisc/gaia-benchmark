#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <math.h>
#include <omp.h>
#include <sys/mman.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <unistd.h>

struct libdeflate_decompressor;
struct libdeflate_decompressor *libdeflate_alloc_decompressor(void);
void libdeflate_free_decompressor(struct libdeflate_decompressor *);
int libdeflate_gzip_decompress(struct libdeflate_decompressor *,
    const void *in, size_t in_nbytes, void *out, size_t out_avail, size_t *actual_out);

#define SOURCE_ID_COLUMN 1
#define BP_FLUX_COLUMN 11
#define RP_FLUX_COLUMN 16
#define ROWS_PER_CHUNK 512
#define OUTPUT_LINE_WIDTH 160

static inline double parse_number(const char **cursor, int *is_number) {
    const char *scan = *cursor;
    int negative = 0;
    if (*scan == '-') { negative = 1; scan++; }
    else if (*scan == '+') scan++;

    if ((*scan < '0' || *scan > '9') && *scan != '.') {
        while (*scan && *scan != ',' && *scan != ']') scan++;
        *cursor = scan; *is_number = 0;
        return 0;
    }

    double mantissa = 0.0;
    while (*scan >= '0' && *scan <= '9') { mantissa = mantissa * 10.0 + (*scan - '0'); scan++; }
    if (*scan == '.') {
        scan++;
        double fraction = 0.0, scale = 1.0;
        while (*scan >= '0' && *scan <= '9') { fraction = fraction * 10.0 + (*scan - '0'); scale *= 10.0; scan++; }
        mantissa += fraction / scale;
    }

    int exponent = 0, exponent_negative = 0;
    if (*scan == 'e' || *scan == 'E') {
        scan++;
        if (*scan == '-') { exponent_negative = 1; scan++; }
        else if (*scan == '+') scan++;
        while (*scan >= '0' && *scan <= '9') { exponent = exponent * 10 + (*scan - '0'); scan++; }
    }

    double value = mantissa;
    if (exponent) {
        double power = 1.0;
        for (int i = 0; i < exponent; i++) power *= 10.0;
        value = exponent_negative ? value / power : value * power;
    }
    *cursor = scan; *is_number = 1;
    return negative ? -value : value;
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

long gaia_scan(const char **paths, int file_count, long slot, int thread_count,
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

long gaia_run(const char **paths, int file_count, long slot, int thread_count,
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

    ResultRow *results = (ResultRow *)malloc(scanned * sizeof(ResultRow));
    long kept_count = 0;
    for (long i = 0; i < scanned; i++) {
        int has_bp = !isnan(bp_min[i]), has_rp = !isnan(rp_min[i]);
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
    free(source_id); free(bp_min); free(bp_max); free(rp_min); free(rp_max); free(counts);

    qsort(results, kept_count, sizeof(ResultRow), compare_pct_change_desc);

    const char *header =
        "source_id,bp_min_flux,bp_max_flux,rp_min_flux,rp_max_flux,percentage_change\n";
    char *row_slots = (char *)malloc((size_t)kept_count * OUTPUT_LINE_WIDTH + 1);
    int *row_lengths = (int *)malloc(kept_count * sizeof(int));
    if (!row_slots || !row_lengths) { free(row_slots); free(row_lengths); free(results); return -1; }

    #pragma omp parallel for schedule(static)
    for (long i = 0; i < kept_count; i++) {
        ResultRow *row = &results[i];
        int has_bp = !isnan(row->bp_min), has_rp = !isnan(row->rp_min);
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
