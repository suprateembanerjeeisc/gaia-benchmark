/* gaiascan.c — libdeflate gunzip + single-pass BP/RP flux min/max scan.
 *
 * Two-phase, so the largest file never becomes a serial straggler:
 *   Phase 1 (parallel over files): read + libdeflate-decompress each gzip file.
 *   Phase 2 (parallel over CHUNKS): every decompressed buffer is split into row
 *     ranges; all chunks across all files are scanned in one parallel loop, so
 *     the big file's scan spreads across otherwise-idle threads.
 * For each row: extract source_id (col 1) and the min/max of the finite,
 * strictly-positive values in the bp_flux / rp_flux arrays (cols 11 / 16),
 * using a fast custom float parser (strtod is far too slow here).
 *
 * Build (no libdeflate dev header needed; prototypes declared below):
 *   gcc -O3 -fopenmp -fPIC -shared gaiascan.c \
 *       /usr/lib/<triplet>/libdeflate.so.0 -o gaiascan.so
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <math.h>
#include <omp.h>

/* ---- libdeflate (stable ABI; declared here so no -dev package is required) ---- */
struct libdeflate_decompressor;
struct libdeflate_decompressor *libdeflate_alloc_decompressor(void);
void libdeflate_free_decompressor(struct libdeflate_decompressor *);
int libdeflate_gzip_decompress(struct libdeflate_decompressor *,
    const void *in, size_t in_nbytes, void *out, size_t out_avail, size_t *actual_out);

#define SOURCE_COL 1
#define BP_COL 11
#define RP_COL 16
#define CHUNK_ROWS 512   /* rows per phase-2 work chunk */

/* Fast positive-decimal parser for flux tokens ("16144.0028...", "1.22e-14").
 * Advances *pp past the token; *ok=0 for a non-numeric token (null/NaN). */
static inline double parse_num(const char **pp, int *ok) {
    const char *s = *pp;
    int neg = 0;
    if (*s == '-') { neg = 1; s++; } else if (*s == '+') s++;
    if ((*s < '0' || *s > '9') && *s != '.') {
        while (*s && *s != ',' && *s != ']') s++;
        *pp = s; *ok = 0; return 0;
    }
    double mant = 0.0;
    while (*s >= '0' && *s <= '9') { mant = mant * 10.0 + (*s - '0'); s++; }
    if (*s == '.') {
        s++;
        double frac = 0.0, scale = 1.0;
        while (*s >= '0' && *s <= '9') { frac = frac * 10.0 + (*s - '0'); scale *= 10.0; s++; }
        mant += frac / scale;
    }
    int exp = 0, eneg = 0;
    if (*s == 'e' || *s == 'E') {
        s++;
        if (*s == '-') { eneg = 1; s++; } else if (*s == '+') s++;
        while (*s >= '0' && *s <= '9') { exp = exp * 10 + (*s - '0'); s++; }
    }
    double val = mant;
    if (exp) { double p = 1.0; for (int k = 0; k < exp; k++) p *= 10.0; val = eneg ? val / p : val * p; }
    *pp = s; *ok = 1;
    return neg ? -val : val;
}

/* Scan one quoted array "[...]" for min/max of positive values. */
static inline const char *scan_array(const char *s, double *mn, double *mx, int *has) {
    /* Sentinels let the min/max updates be branchless (no per-value has-flag). */
    double lo = 1e308, hi = -1e308;
    if (*s == '[') s++;
    while (*s && *s != ']') {
        while (*s == ',') s++;          /* arrays are "[n,n,...]" — no spaces */
        if (*s == ']' || *s == 0) break;
        int ok;
        double v = parse_num(&s, &ok);
        if (ok && v > 0.0) {
            if (v < lo) lo = v;
            if (v > hi) hi = v;
        }
    }
    if (*s == ']') s++;
    *has = (hi >= lo);
    *mn = lo; *mx = hi;
    return s;
}

/* Scan one CSV row [p, rowend); write result at out index. Returns 1 if kept.
 * Fields are located with memchr (jumping over quoted flux arrays whole) rather
 * than a per-byte scan, and we stop as soon as the last needed column (RP) is
 * read — the columns after it are never touched. */
static int scan_row(const char *p, const char *rowend,
                    int64_t *ids, double *bpmn, double *bpmx,
                    double *rpmn, double *rpmx, long oi) {
    int col = 0;
    const char *field = p, *s = p;
    int64_t sid = 0;
    double bmn = 0, bmx = 0, rmn = 0, rmx = 0;
    int hb = 0, hr = 0;
    while (s < rowend) {
        const char *nc;
        if (*s == '"') {                       /* quoted flux array: jump to closing quote, then comma */
            const char *cq = (const char *)memchr(s + 1, '"', rowend - (s + 1));
            nc = cq ? (const char *)memchr(cq, ',', rowend - cq) : NULL;
        } else {
            nc = (const char *)memchr(s, ',', rowend - s);
        }
        if (col == SOURCE_COL) sid = strtoll(field, NULL, 10);
        else if (col == BP_COL) { const char *q = field; if (*q == '"') q++; scan_array(q, &bmn, &bmx, &hb); }
        else if (col == RP_COL) { const char *q = field; if (*q == '"') q++; scan_array(q, &rmn, &rmx, &hr); }
        if (!nc || col >= RP_COL) break;       /* done once RP is read */
        col++; field = nc + 1; s = nc + 1;
    }
    if (!hb && !hr) return 0;
    ids[oi] = sid;
    bpmn[oi] = hb ? bmn : NAN; bpmx[oi] = hb ? bmx : NAN;
    rpmn[oi] = hr ? rmn : NAN; rpmx[oi] = hr ? rmx : NAN;
    return 1;
}

/* Per-file decompressed buffer + its row (line-start) index. */
typedef struct {
    char *buf;          /* decompressed, NUL-terminated */
    size_t len;
    long *line_off;     /* offsets of each data row start */
    long nrows;         /* number of data rows */
} FileBuf;

/* A phase-2 work chunk: rows [r0, r1) of file f, output starting at out_base. */
typedef struct { int f; long r0, r1, out_base; } Chunk;

long gaia_scan(const char **paths, int nfiles, long slot, int nthreads,
               int64_t *ids, double *bpmn, double *bpmx, double *rpmn, double *rpmx,
               long *counts) {
    if (nthreads > 0) omp_set_num_threads(nthreads);
    FileBuf *fb = (FileBuf *)calloc(nfiles, sizeof(FileBuf));
    int err = 0;

    /* ---- Phase 1: decompress every file in parallel, index its rows. ---- */
    #pragma omp parallel for schedule(dynamic)
    for (int i = 0; i < nfiles; i++) {
        FILE *f = fopen(paths[i], "rb");
        if (!f) { err = 1; continue; }
        fseek(f, 0, SEEK_END); long fsz = ftell(f); fseek(f, 0, SEEK_SET);
        unsigned char *in = (unsigned char *)malloc(fsz);
        if (!in || fread(in, 1, fsz, f) != (size_t)fsz) { free(in); fclose(f); err = 1; continue; }
        fclose(f);
        uint32_t isize; memcpy(&isize, in + fsz - 4, 4);
        size_t cap = isize ? isize : (size_t)fsz * 20;
        char *out = (char *)malloc(cap + 1);
        if (!out) { free(in); err = 1; continue; }
        struct libdeflate_decompressor *d = libdeflate_alloc_decompressor();
        size_t actual = 0;
        int rc = libdeflate_gzip_decompress(d, in, fsz, out, cap, &actual);
        libdeflate_free_decompressor(d); free(in);
        if (rc != 0) { free(out); err = 1; continue; }
        out[actual] = 0;

        /* index data-row starts: skip '#' comment lines + the single header line */
        const char *p = out, *end = out + actual;
        while (p < end && *p == '#') { while (p < end && *p != '\n') p++; if (p < end) p++; }
        if (p < end) { while (p < end && *p != '\n') p++; if (p < end) p++; }  /* header */
        long cap_rows = actual / 64 + 16;
        long *off = (long *)malloc(cap_rows * sizeof(long));
        long nr = 0;
        while (p < end) {
            if (nr >= cap_rows) { cap_rows *= 2; off = (long *)realloc(off, cap_rows * sizeof(long)); }
            off[nr++] = (long)(p - out);
            while (p < end && *p != '\n') p++;
            if (p < end) p++;
        }
        fb[i].buf = out; fb[i].len = actual; fb[i].line_off = off; fb[i].nrows = nr;
    }
    if (err) { for (int i = 0; i < nfiles; i++) { free(fb[i].buf); free(fb[i].line_off); } free(fb); return -1; }

    /* ---- Build the chunk list; each file's rows land in its slot [i*slot, ..). ---- */
    long total_chunks = 0;
    for (int i = 0; i < nfiles; i++) total_chunks += (fb[i].nrows + CHUNK_ROWS - 1) / CHUNK_ROWS;
    Chunk *ch = (Chunk *)malloc(total_chunks * sizeof(Chunk));
    long c = 0;
    for (int i = 0; i < nfiles; i++)
        for (long r = 0; r < fb[i].nrows; r += CHUNK_ROWS) {
            long r1 = r + CHUNK_ROWS; if (r1 > fb[i].nrows) r1 = fb[i].nrows;
            ch[c].f = i; ch[c].r0 = r; ch[c].r1 = r1; ch[c].out_base = (long)i * slot + r;
            c++;
        }

    /* per-chunk kept-count, so we can compact each file's slot afterward */
    long *chunk_kept = (long *)calloc(total_chunks, sizeof(long));

    /* ---- Phase 2: scan all chunks in parallel (big file now split up). ---- */
    #pragma omp parallel for schedule(dynamic)
    for (long k = 0; k < total_chunks; k++) {
        FileBuf *b = &fb[ch[k].f];
        long kept = 0, base = ch[k].out_base;
        for (long r = ch[k].r0; r < ch[k].r1; r++) {
            const char *rp = b->buf + b->line_off[r];
            const char *rend = (r + 1 < b->nrows) ? b->buf + b->line_off[r + 1] : b->buf + b->len;
            kept += scan_row(rp, rend, ids, bpmn, bpmx, rpmn, rpmx, base + kept);
        }
        chunk_kept[k] = kept;
    }

    /* ---- Compact: gather each file's kept rows (they sit at slot start, since
       chunks within a file are contiguous and processed into base+kept, but
       kept<chunk means gaps between chunks) into a packed [0,total). ---- */
    /* First collapse gaps within each file's slot, then across files. */
    long total = 0;
    for (int i = 0; i < nfiles; i++) counts[i] = 0;
    /* recompute per-file kept by summing its chunks, and compact chunk gaps */
    long ci = 0;
    for (int i = 0; i < nfiles; i++) {
        long fbase = (long)i * slot, w = fbase;
        long nchunk = (fb[i].nrows + CHUNK_ROWS - 1) / CHUNK_ROWS;
        for (long j = 0; j < nchunk; j++, ci++) {
            long rb = fbase + ch[ci].r0, kept = chunk_kept[ci];
            if (rb != w && kept) {
                memmove(ids + w, ids + rb, kept * sizeof(int64_t));
                memmove(bpmn + w, bpmn + rb, kept * sizeof(double));
                memmove(bpmx + w, bpmx + rb, kept * sizeof(double));
                memmove(rpmn + w, rpmn + rb, kept * sizeof(double));
                memmove(rpmx + w, rpmx + rb, kept * sizeof(double));
            }
            w += kept;
        }
        counts[i] = w - fbase;
    }
    /* pack files down into one contiguous block */
    for (int i = 0; i < nfiles; i++) {
        long fbase = (long)i * slot, cnt = counts[i];
        if (fbase != total && cnt) {
            memmove(ids + total, ids + fbase, cnt * sizeof(int64_t));
            memmove(bpmn + total, bpmn + fbase, cnt * sizeof(double));
            memmove(bpmx + total, bpmx + fbase, cnt * sizeof(double));
            memmove(rpmn + total, rpmn + fbase, cnt * sizeof(double));
            memmove(rpmx + total, rpmx + fbase, cnt * sizeof(double));
        }
        total += cnt;
    }

    for (int i = 0; i < nfiles; i++) { free(fb[i].buf); free(fb[i].line_off); }
    free(fb); free(ch); free(chunk_kept);
    return total;
}

/* ---- Full pipeline in C: scan + %change + filter>100 + sort + write CSV. ----
 * Avoids the Polars import and any Python-side marshalling. `out_path` receives
 * the challenge CSV. Returns the number of rows written, or -1 on error. */

typedef struct { int64_t sid; double bmn, bmx, rmn, rmx, pc; } Res;

static int cmp_pc_desc(const void *a, const void *b) {
    double x = ((const Res *)a)->pc, y = ((const Res *)b)->pc;
    return (x < y) - (x > y);   /* descending */
}

/* Format a double the way the challenge CSV expects (repr-style, shortest
 * round-trippable). "%.17g" is safe/round-trippable; empty bands print "". */
static int fmt_num(char *dst, double v, int has) {
    if (!has) { return 0; }                 /* empty cell */
    return sprintf(dst, "%.17g", v);
}

long gaia_run(const char **paths, int nfiles, long slot, int nthreads,
              const char *out_path) {
    long cap = (long)nfiles * slot;
    int64_t *ids = (int64_t *)malloc(cap * sizeof(int64_t));
    double *bmn = (double *)malloc(cap * sizeof(double));
    double *bmx = (double *)malloc(cap * sizeof(double));
    double *rmn = (double *)malloc(cap * sizeof(double));
    double *rmx = (double *)malloc(cap * sizeof(double));
    long *counts = (long *)malloc(nfiles * sizeof(long));
    if (!ids || !bmn || !bmx || !rmn || !rmx || !counts) return -1;

    long total = gaia_scan(paths, nfiles, slot, nthreads, ids, bmn, bmx, rmn, rmx, counts);
    if (total < 0) { free(ids); free(bmn); free(bmx); free(rmn); free(rmx); free(counts); return -1; }

    /* compute percentage_change (max of the two bands) and keep > 100 */
    Res *res = (Res *)malloc(total * sizeof(Res));
    long m = 0;
    for (long k = 0; k < total; k++) {
        int hb = !isnan(bmn[k]), hr = !isnan(rmn[k]);
        double bp = hb ? (bmx[k] - bmn[k]) / bmn[k] * 100.0 : -1.0;
        double rp = hr ? (rmx[k] - rmn[k]) / rmn[k] * 100.0 : -1.0;
        double pc = bp > rp ? bp : rp;
        if (pc > 100.0) {
            res[m].sid = ids[k];
            res[m].bmn = bmn[k]; res[m].bmx = bmx[k];
            res[m].rmn = rmn[k]; res[m].rmx = rmx[k];
            res[m].pc = pc; m++;
        }
    }
    free(ids); free(bmn); free(bmx); free(rmn); free(rmx); free(counts);

    qsort(res, m, sizeof(Res), cmp_pc_desc);

    /* Format the CSV in parallel: each row gets a fixed-width slot (SLOT_W bytes),
     * threads format independently, then a serial pass compacts the slots into a
     * contiguous buffer written in one fwrite. Formatting (many sprintf calls)
     * was the dominant tail cost, and it is embarrassingly parallel. */
    const char *hdr = "source_id,bp_min_flux,bp_max_flux,rp_min_flux,rp_max_flux,percentage_change\n";
    const int SLOT_W = 160;
    char *slots = (char *)malloc((size_t)m * SLOT_W + 1);
    int *lens = (int *)malloc(m * sizeof(int));
    if (!slots || !lens) { free(slots); free(lens); free(res); return -1; }
    #pragma omp parallel for schedule(static)
    for (long k = 0; k < m; k++) {
        Res *r = &res[k];
        int hb = !isnan(r->bmn), hr = !isnan(r->rmn);
        char *q = slots + (size_t)k * SLOT_W;
        char *p = q;
        p += sprintf(p, "%lld,", (long long)r->sid);
        p += fmt_num(p, r->bmn, hb); *p++ = ',';
        p += fmt_num(p, r->bmx, hb); *p++ = ',';
        p += fmt_num(p, r->rmn, hr); *p++ = ',';
        p += fmt_num(p, r->rmx, hr); *p++ = ',';
        p += sprintf(p, "%.17g\n", r->pc);
        lens[k] = (int)(p - q);
    }
    size_t bufsz = strlen(hdr) + (size_t)m * SLOT_W + 16;
    char *buf = (char *)malloc(bufsz), *p = buf;
    if (!buf) { free(slots); free(lens); free(res); return -1; }
    memcpy(p, hdr, strlen(hdr)); p += strlen(hdr);
    for (long k = 0; k < m; k++) { memcpy(p, slots + (size_t)k * SLOT_W, lens[k]); p += lens[k]; }
    FILE *f = fopen(out_path, "wb");
    if (!f) { free(buf); free(slots); free(lens); free(res); return -1; }
    fwrite(buf, 1, p - buf, f);
    fclose(f);
    free(buf); free(slots); free(lens); free(res);
    return m;
}
