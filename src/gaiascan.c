#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <math.h>
#include <omp.h>

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
    DecompressedFile *files = (DecompressedFile *)calloc(file_count, sizeof(DecompressedFile));
    int had_error = 0;

    /* ---- Phase 1: decompress every file in parallel, index its rows. ---- */
    #pragma omp parallel for schedule(dynamic)
    for (int i = 0; i < file_count; i++) {
        FILE *handle = fopen(paths[i], "rb");
        if (!handle) { had_error = 1; continue; }
        fseek(handle, 0, SEEK_END);
        long compressed_size = ftell(handle);
        fseek(handle, 0, SEEK_SET);
        unsigned char *compressed = (unsigned char *)malloc(compressed_size);
        if (!compressed || fread(compressed, 1, compressed_size, handle) != (size_t)compressed_size) {
            free(compressed); fclose(handle); had_error = 1; continue;
        }
        fclose(handle);

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
        free(compressed);
        if (result != 0) { free(decompressed); had_error = 1; continue; }
        decompressed[decompressed_len] = 0;

        const char *cursor = decompressed, *buf_end = decompressed + decompressed_len;
        while (cursor < buf_end && *cursor == '#') {
            while (cursor < buf_end && *cursor != '\n') cursor++;
            if (cursor < buf_end) cursor++;
        }
        if (cursor < buf_end) {  /* skip header row */
            while (cursor < buf_end && *cursor != '\n') cursor++;
            if (cursor < buf_end) cursor++;
        }
        long offset_capacity = decompressed_len / 64 + 16;
        long *row_offsets = (long *)malloc(offset_capacity * sizeof(long));
        long row_count = 0;
        while (cursor < buf_end) {
            if (row_count >= offset_capacity) {
                offset_capacity *= 2;
                row_offsets = (long *)realloc(row_offsets, offset_capacity * sizeof(long));
            }
            row_offsets[row_count++] = (long)(cursor - decompressed);
            while (cursor < buf_end && *cursor != '\n') cursor++;
            if (cursor < buf_end) cursor++;
        }
        files[i].data = decompressed;
        files[i].length = decompressed_len;
        files[i].row_offsets = row_offsets;
        files[i].row_count = row_count;
    }
    if (had_error) {
        for (int i = 0; i < file_count; i++) { free(files[i].data); free(files[i].row_offsets); }
        free(files);
        return -1;
    }

    long chunk_total = 0;
    for (int i = 0; i < file_count; i++)
        chunk_total += (files[i].row_count + ROWS_PER_CHUNK - 1) / ROWS_PER_CHUNK;
    WorkChunk *chunks = (WorkChunk *)malloc(chunk_total * sizeof(WorkChunk));
    long chunk_index = 0;
    for (int i = 0; i < file_count; i++) {
        for (long row = 0; row < files[i].row_count; row += ROWS_PER_CHUNK) {
            long last_row = row + ROWS_PER_CHUNK;
            if (last_row > files[i].row_count) last_row = files[i].row_count;
            chunks[chunk_index].file_index = i;
            chunks[chunk_index].first_row = row;
            chunks[chunk_index].last_row = last_row;
            chunks[chunk_index].output_base = (long)i * slot + row;
            chunk_index++;
        }
    }

    long *chunk_kept = (long *)calloc(chunk_total, sizeof(long));

    /* ---- Phase 2: scan all chunks in parallel (the big file is now split up). ---- */
    #pragma omp parallel for schedule(dynamic)
    for (long c = 0; c < chunk_total; c++) {
        DecompressedFile *file = &files[chunks[c].file_index];
        long kept = 0, base = chunks[c].output_base;
        for (long row = chunks[c].first_row; row < chunks[c].last_row; row++) {
            const char *row_start = file->data + file->row_offsets[row];
            const char *row_end = (row + 1 < file->row_count)
                                  ? file->data + file->row_offsets[row + 1]
                                  : file->data + file->length;
            kept += scan_row(row_start, row_end, out_source_id, out_bp_min, out_bp_max,
                             out_rp_min, out_rp_max, base + kept);
        }
        chunk_kept[c] = kept;
    }

    long total = 0;
    for (int i = 0; i < file_count; i++) counts[i] = 0;
    long c = 0;
    for (int i = 0; i < file_count; i++) {
        long file_base = (long)i * slot, write_pos = file_base;
        long chunk_count = (files[i].row_count + ROWS_PER_CHUNK - 1) / ROWS_PER_CHUNK;
        for (long j = 0; j < chunk_count; j++, c++) {
            long read_pos = file_base + chunks[c].first_row, kept = chunk_kept[c];
            if (read_pos != write_pos && kept) {
                memmove(out_source_id + write_pos, out_source_id + read_pos, kept * sizeof(int64_t));
                memmove(out_bp_min + write_pos, out_bp_min + read_pos, kept * sizeof(double));
                memmove(out_bp_max + write_pos, out_bp_max + read_pos, kept * sizeof(double));
                memmove(out_rp_min + write_pos, out_rp_min + read_pos, kept * sizeof(double));
                memmove(out_rp_max + write_pos, out_rp_max + read_pos, kept * sizeof(double));
            }
            write_pos += kept;
        }
        counts[i] = write_pos - file_base;
    }
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

    for (int i = 0; i < file_count; i++) { free(files[i].data); free(files[i].row_offsets); }
    free(files); free(chunks); free(chunk_kept);
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
