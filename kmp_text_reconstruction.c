/* ============================================================================
 * Compile:
 *   gcc -O2 -o kmp_text_reconstruction kmp_text_reconstruction.c
 *
 * Usage:
 *   ./kmp_text_reconstruction <input_file>          (or '-' for stdin)
 *   ./kmp_text_reconstruction --algos LIST <input_file>  (subset of greedy,mgreedy,tgreedy,hk,bb)
 *
 * Output protocol:
 *   stdout : reconstructed text, one solution per line, in non-increasing length order.
 *   stderr : per-algorithm metrics : [algo=NAME] len=N elapsed=T.TTTs result=...
 *
 * Architecture:
 *   A. Common types & utilities
 *   B. Best-solution registry
 *   C. Step 1 — preprocess (KMP-fused overlap + containment)
 *   D. Step 2 — correct / sub-optimal / quick algorithms
 *   E. Step 3 — correct / optimal / slow algorithms
 *   F. main()
 *
 * KMP Enhancement (Tarhio & Ukkonen 1988, Section 2.1–2.2):
 *   The KMP automaton for pattern x' scanned over text x simultaneously:
 *     (1) computes the maximal suffix-prefix overlap between x and x', AND
 *     (2) detects whether x' is a substring of x (containment).
 *   Both jobs in O(|x| + |x'|) per pair, yielding O(mn) total for all pairwise
 *   overlaps and containment removal, where m = total string length, n = count.
 *   This replaces the original brute-force O(n² · L²) implementation.
 * ============================================================================ */

/* ============================================================================
 * Section A: Common types & utilities
 * ============================================================================ */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <strings.h>
#include <stdbool.h>
#include <time.h>
#include <stdint.h>
#include <limits.h>

/* --- A.0: Algorithm selection bitmask --- */

#define ALGO_GREEDY   (1u << 0)
#define ALGO_MGREEDY  (1u << 1)
#define ALGO_TGREEDY  (1u << 2)
#define ALGO_HK       (1u << 3)
#define ALGO_BB       (1u << 4)
#define ALGO_ALL      0xFFu

static int
parse_algos(const char *s)
{
    int mask = 0;
    const char *p = s;
    while (*p) {
        const char *start = p;
        while (*p && *p != ',') p++;
        size_t len = (size_t)(p - start);
        if      (len == 6 && strncasecmp(start, "greedy",  6) == 0) mask |= ALGO_GREEDY;
        else if (len == 7 && strncasecmp(start, "mgreedy", 7) == 0) mask |= ALGO_MGREEDY;
        else if (len == 7 && strncasecmp(start, "tgreedy", 7) == 0) mask |= ALGO_TGREEDY;
        else if (len == 2 && strncasecmp(start, "hk",      2) == 0) mask |= ALGO_HK;
        else if (len == 2 && strncasecmp(start, "bb",      2) == 0) mask |= ALGO_BB;
        else return -1;
        if (*p == ',') p++;
    }
    return mask;
}

/* --- A.1: Timing --- */

static double
mono_seconds(void)
{
    return (double)clock() / CLOCKS_PER_SEC;
}

/* --- A.2: Fragments & string utilities --- */

typedef struct {
    char  *str;
    size_t len;
} Fragment;

typedef struct {
    Fragment *items;
    size_t    count;
    size_t    capacity;
} FragmentArray;

static FragmentArray
fa_new(void)
{
    FragmentArray fa = { NULL, 0, 0 };
    return fa;
}

static void
fa_push(FragmentArray *fa, char *str)
{
    if (fa->count == fa->capacity) {
        size_t new_cap = fa->capacity == 0 ? 8 : fa->capacity * 2;
        fa->items = realloc(fa->items, new_cap * sizeof(Fragment));
        fa->capacity = new_cap;
    }
    fa->items[fa->count].str = str;
    fa->items[fa->count].len = strlen(str);
    fa->count++;
}

static void
fa_free(FragmentArray *fa)
{
    for (size_t i = 0; i < fa->count; i++) free(fa->items[i].str);
    free(fa->items);
    fa->items = NULL;
    fa->count = fa->capacity = 0;
}

/* --- A.2b: KMP failure function (preprocessing phase) ---
 *
 * Build the KMP failure (prefix) function for pattern pat of length m.
 * fail[i] = length of the longest proper prefix of pat[0..i] that is
 *           also a suffix of pat[0..i].
 * Time: O(m). Space: O(m).
 *
 * This is the standard KMP preprocessing from Knuth, Morris & Pratt (1977).
 */
static int *
kmp_build_failure(const char *pat, size_t m)
{
    int *fail = malloc(m * sizeof(int));
    fail[0] = 0;
    size_t k = 0;  /* length of previous longest prefix-suffix */
    for (size_t i = 1; i < m; i++) {
        while (k > 0 && pat[k] != pat[i])
            k = (size_t)fail[k - 1];
        if (pat[k] == pat[i])
            k++;
        fail[i] = (int)k;
    }
    return fail;
}

/* --- A.2c: KMP-based overlap + containment (Tarhio & Ukkonen fused step) ---
 *
 * Given text a (length la) and pattern b (length lb) with precomputed failure
 * function fail_b, scan a with the KMP automaton for b.
 *
 * Returns: the maximal suffix-prefix overlap, i.e. the longest suffix of a
 *          that equals a prefix of b. Range: [0, min(la, lb) - 1].
 *          (We exclude full containment from the overlap value.)
 *
 * Side effect: sets *b_contained_in_a = true if b occurs as a substring of a
 *              (i.e., automaton reaches state lb during scan).
 *
 * Time: O(la + lb) (lb for the failure function, la for the scan).
 *       But since fail_b is precomputed and reused, per-call cost is O(la).
 */
static size_t
kmp_overlap_and_containment(const char *a, size_t la,
                            const char *b, size_t lb,
                            const int  *fail_b,
                            bool       *b_contained_in_a)
{
    *b_contained_in_a = false;

    /* Edge case: if b is longer than a, b cannot be contained, and
       the overlap is at most la-1. We still run KMP. */

    size_t q = 0;  /* current state = number of characters of b matched */

    for (size_t i = 0; i < la; i++) {
        while (q > 0 && b[q] != a[i])
            q = (size_t)fail_b[q - 1];
        if (b[q] == a[i])
            q++;
        if (q == lb) {
            /* Full match of b inside a detected */
            *b_contained_in_a = true;
            /* Continue scanning: we still need the final state for overlap */
            q = (size_t)fail_b[q - 1];
        }
    }

    /* After scanning all of a, the automaton state q is the length of the
       longest prefix of b that matches a suffix of a. This is the overlap.
       But if q == lb, that means a ends with all of b (full overlap), which
       only happens if lb <= la; that's containment, already flagged above.
       After the containment reset, q < lb always here. */

    return q;
}

/* Concatenate a and b with the overlap region collapsed (O(L)) */
static char *
merge_with_overlap(const char *a, size_t la,
                   const char *b, size_t lb, size_t ov)
{
    size_t out_len = la + lb - ov;
    char *out = malloc(out_len + 1);
    memcpy(out, a, la);
    memcpy(out + la, b + ov, lb - ov);
    out[out_len] = '\0';
    return out;
}

/* --- A.3: Greedy pairwise merging --- */

/* Find the longest suffix-prefix overlap between a and b using KMP (O(la + lb)).
 * This is the drop-in replacement for the old brute-force overlap_chars. */
static size_t
overlap_kmp(const char *a, size_t la, const char *b, size_t lb)
{
    if (la == 0 || lb == 0) return 0;

    int *fail_b = kmp_build_failure(b, lb);
    bool contained = false;
    size_t ov = kmp_overlap_and_containment(a, la, b, lb, fail_b, &contained);
    free(fail_b);

    /* If b is entirely contained in a, the overlap is defined as |b|
       (per Tarhio & Ukkonen: k = |x'| when x' is a substring of x).
       But for the greedy merge, contained fragments should already be removed.
       In the overlap matrix context, we cap at lb-1 to avoid self-collapse. */
    if (contained) return lb - 1 < la ? lb : lb - 1;

    return ov;
}

/* Loop until one string remains (n − 1 iterations): scan all pairs to find
   the maximum-overlap pair; merge them into a single string.
   Uses KMP overlap: O(n² · L) per iteration, O(n³ · L) total. */
static char *
greedy_merge_pairs(Fragment *frags, size_t n, size_t *out_len)
{
    while (n > 1) {
        size_t best_i = 0, best_j = 1, best_ov = 0;
        for (size_t i = 0; i < n; i++) {
            for (size_t j = 0; j < n; j++) {
                if (i == j) continue;
                size_t ov = overlap_kmp(frags[i].str, frags[i].len,
                                        frags[j].str, frags[j].len);
                if (ov > best_ov) { best_ov = ov; best_i = i; best_j = j; }
            }
        }
        char *merged = merge_with_overlap(frags[best_i].str, frags[best_i].len,
                                          frags[best_j].str, frags[best_j].len,
                                          best_ov);
        size_t merged_len = frags[best_i].len + frags[best_j].len - best_ov;
        free(frags[best_i].str);
        free(frags[best_j].str);
        frags[best_i].str = merged;
        frags[best_i].len = merged_len;
        if (best_j != n - 1) frags[best_j] = frags[n - 1];
        n--;
    }
    *out_len = frags[0].len;
    return frags[0].str;
}

/* --- A.4: Overlap matrix (KMP-fused: O(mn) total) & edge sort --- */

/* Build the overlap matrix using KMP.
 *
 * Per Tarhio & Ukkonen Section 2.1: for each fragment x_i, build the KMP
 * failure function once (O(|x_i|)), then scan all other fragments x_j through
 * the automaton (O(|x_j|) each). Total over all (i,j) pairs: O(mn) where
 * m = sum of all fragment lengths, n = number of fragments.
 *
 * This also performs containment detection (Step A2 fused into A1):
 * if x_j is found inside x_i during the scan, mark x_j for removal.
 *
 * The old brute-force was O(n² · L²). This is O(mn) ≈ O(n · L · n) = O(n² · L)
 * in the worst case but with much better constants due to KMP's linear scan. */
static int **
build_overlap_matrix(const FragmentArray *fa)
{
    size_t n = fa->count;
    int **m = malloc(n * sizeof(int *));
    for (size_t i = 0; i < n; i++) m[i] = calloc(n, sizeof(int));

    for (size_t i = 0; i < n; i++) {
        /* Build KMP failure function for fragment i (used as pattern) */
        int *fail_i = kmp_build_failure(fa->items[i].str, fa->items[i].len);

        for (size_t j = 0; j < n; j++) {
            if (i == j) continue;

            /* Scan text=j with pattern=i → overlap(j, i) and containment of i in j */
            bool contained = false;
            size_t ov = kmp_overlap_and_containment(
                fa->items[j].str, fa->items[j].len,    /* text = fragment j */
                fa->items[i].str, fa->items[i].len,     /* pattern = fragment i */
                fail_i,
                &contained);

            /* ov = longest suffix of j that is a prefix of i = overlap(j, i) */
            m[j][i] = (int)ov;
        }

        free(fail_i);
    }
    return m;
}

static void
free_overlap_matrix(int **m, size_t n)
{
    for (size_t i = 0; i < n; i++) free(m[i]);
    free(m);
}

typedef struct {
    int i;
    int j;
    int weight;
} Edge;

static int
compare_edge_desc(const void *a, const void *b)
{
    const Edge *ea = (const Edge *)a;
    const Edge *eb = (const Edge *)b;
    return eb->weight - ea->weight;
}

/* --- A.5: Cycle-cover graph helpers --- */

/* Pick edges greedily under in/out-degree ≤ 1 (O(n² log n)) */
static int **
mgreedy_select_edges(int **overlap_matrix, size_t n)
{
    /* Collect candidate edges (O(n²)) */
    Edge *edges = malloc(n * n * sizeof(Edge));
    int edge_count = 0;
    for (size_t i = 0; i < n; i++) {
        for (size_t j = 0; j < n; j++) {
            if (i != j && overlap_matrix[i][j] > 0) {
                edges[edge_count].i = (int)i;
                edges[edge_count].j = (int)j;
                edges[edge_count].weight = overlap_matrix[i][j];
                edge_count++;
            }
        }
    }
    /* Sort by descending weight (O(n² log n)) */
    qsort(edges, edge_count, sizeof(Edge), compare_edge_desc);

    int **selected = malloc(n * sizeof(int *));
    for (size_t i = 0; i < n; i++) selected[i] = calloc(n, sizeof(int));

    /* Accept edges keeping in/out-degree ≤ 1 (O(n²)) */
    bool *out_used = calloc(n, sizeof(bool));
    bool *in_used  = calloc(n, sizeof(bool));
    for (int e = 0; e < edge_count; e++) {
        int i = edges[e].i, j = edges[e].j;
        if (out_used[i] || in_used[j]) continue;
        selected[i][j] = edges[e].weight;
        out_used[i] = true;
        in_used[j] = true;
    }

    free(out_used);
    free(in_used);
    free(edges);
    return selected;
}

/* Break cycles at weakest edge (O(n²)) */
static void
mgreedy_break_cycles(int **selected, size_t n)
{
    /* Compute in-degrees (O(n²)) */
    int *in_count = calloc(n, sizeof(int));
    for (size_t i = 0; i < n; i++)
        for (size_t j = 0; j < n; j++)
            if (selected[i][j] > 0) in_count[j]++;

    bool *visited = calloc(n, sizeof(bool));

    /* Mark nodes reachable from path heads (O(n²)) */
    for (size_t i = 0; i < n; i++) {
        if (in_count[i] != 0) continue;
        int cur = (int)i;
        while (cur != -1 && !visited[cur]) {
            visited[cur] = true;
            int next = -1;
            for (size_t k = 0; k < n; k++) {
                if (selected[cur][k] > 0) { next = (int)k; break; }
            }
            cur = next;
        }
    }

    /* Walk each unvisited cycle, drop weakest edge (O(n²)) */
    for (size_t i = 0; i < n; i++) {
        if (visited[i]) continue;
        int min_from = -1, min_to = -1, min_w = INT_MAX;
        int cur = (int)i;
        do {
            visited[cur] = true;
            int next = -1;
            for (size_t k = 0; k < n; k++) {
                if (selected[cur][k] > 0) { next = (int)k; break; }
            }
            if (selected[cur][next] < min_w) {
                min_w = selected[cur][next];
                min_from = cur;
                min_to = next;
            }
            cur = next;
        } while (cur != (int)i);
        selected[min_from][min_to] = 0;
    }

    free(visited);
    free(in_count);
}

/* Concatenate path strings (O(n²)) */
static char **
mgreedy_build_sequences(const FragmentArray *fa, int **selected, size_t n, size_t *out_count)
{
    char **sequences = malloc(n * sizeof(char *));
    size_t count = 0;

    for (size_t j = 0; j < n; j++) {
        bool is_start = true;
        for (size_t k = 0; k < n; k++) {
            if (selected[k][j] > 0) { is_start = false; break; }
        }
        if (!is_start) continue;

        size_t total_len = fa->items[j].len;
        {
            int from = (int)j;
            while (1) {
                int next = -1;
                for (size_t k = 0; k < n; k++) {
                    if (selected[from][k] > 0) { next = (int)k; break; }
                }
                if (next == -1) break;
                total_len += fa->items[next].len - (size_t)selected[from][next];
                from = next;
            }
        }

        char *seq = malloc(total_len + 1);
        memcpy(seq, fa->items[j].str, fa->items[j].len);
        size_t pos = fa->items[j].len;

        int from = (int)j;
        while (1) {
            int next = -1;
            for (size_t k = 0; k < n; k++) {
                if (selected[from][k] > 0) { next = (int)k; break; }
            }
            if (next == -1) break;
            size_t ov = (size_t)selected[from][next];
            memcpy(seq + pos, fa->items[next].str + ov, fa->items[next].len - ov);
            pos += fa->items[next].len - ov;
            from = next;
        }
        seq[pos] = '\0';
        sequences[count++] = seq;
    }

    *out_count = count;
    return sequences;
}

/* --- A.6: Branch & Bound primitives (search state, min-heap, lower bound) --- */

typedef struct {
    int      lb;       /* lower bound on final SCS length from this state */
    uint64_t mask;     /* bitmask of visited fragments                    */
    int      last;     /* index of last visited fragment                  */
    int      curlen;   /* total SCS length accumulated so far             */
    int     *path;     /* fragment-index sequence (owned by this state)   */
    int      path_len;
} BBState;

typedef struct {
    BBState *data;
    int      size;
    int      cap;
} MinHeap;

static void
mh_push(MinHeap *h, BBState s)
{
    if (h->size == h->cap) {
        h->cap  = h->cap ? h->cap * 2 : 64;
        h->data = realloc(h->data, (size_t)h->cap * sizeof(BBState));
    }
    int i = h->size++;
    h->data[i] = s;
    while (i > 0) {
        int p = (i - 1) / 2;
        if (h->data[p].lb <= h->data[i].lb) break;
        BBState tmp = h->data[p]; h->data[p] = h->data[i]; h->data[i] = tmp;
        i = p;
    }
}

static BBState
mh_pop(MinHeap *h)
{
    BBState ret = h->data[0];
    h->data[0]  = h->data[--h->size];
    for (int i = 0;;) {
        int l = 2*i+1, r = 2*i+2, m = i;
        if (l < h->size && h->data[l].lb < h->data[m].lb) m = l;
        if (r < h->size && h->data[r].lb < h->data[m].lb) m = r;
        if (m == i) break;
        BBState tmp = h->data[m]; h->data[m] = h->data[i]; h->data[i] = tmp;
        i = m;
    }
    return ret;
}

/* Compute a lower bound on the remaining characters needed (O(n²)) */
static int
remaining_lb(uint64_t mask, int n, int **ov, int *flen)
{
    int lb = 0;
    for (int u = 0; u < n; u++) {
        if (mask & ((uint64_t)1 << u)) continue;
        int best_ov = 0;
        for (int v = 0; v < n; v++) {
            if (v != u && ov[v][u] > best_ov) best_ov = ov[v][u];
        }
        lb += flen[u] - best_ov;
    }
    return lb;
}

/* ============================================================================
 * Section B: Best-solution registry
 * ============================================================================ */

typedef struct {
    char  *str;
    size_t len;
} BestSolution;

static BestSolution g_best = { NULL, 0 };

/* Log algo metrics; if candidate beats g_best, publish it to stdout, else free it */
static int
try_record_solution(char *candidate,
                    size_t cand_len,
                    const char *algo_label,
                    double elapsed_sec)
{
    if (candidate == NULL) {
        fprintf(stderr,
                "[algo=%s] len=%zu elapsed=%.3fs (no improvement over g_best)\n",
                algo_label,
                g_best.len,
                elapsed_sec);
        return 0;
    }

    fprintf(stderr,
            "[algo=%s] len=%zu elapsed=%.3fs result=%s\n",
            algo_label,
            cand_len,
            elapsed_sec,
            candidate);

    bool improved = (g_best.str == NULL) || (cand_len < g_best.len);

    if (!improved) {
        free(candidate);
        return 0;
    }

    free(g_best.str);

    g_best.str = candidate;
    g_best.len = cand_len;
    fputs(candidate, stdout);
    fputc('\n', stdout);
    fflush(stdout);

    return 1;
}

/* ============================================================================
 * Section C: Step 1 — preprocess
 *
 * KMP-fused preprocessing (Tarhio & Ukkonen 1988, Steps A1+A2 combined):
 *   For each fragment x_i, build the KMP automaton. Scan every other fragment
 *   x_j through it. During the scan:
 *     - If the automaton reaches state |x_i|, x_i is contained in x_j → mark
 *       x_i for removal.
 *     - After the full scan, the final state gives the overlap(x_j, x_i).
 *   This fuses Steps A1 (overlap computation) and A2 (containment removal)
 *   into a single O(mn) pass, exactly as the paper prescribes.
 *
 *   However, since GREEDY needs to re-merge dynamically (changing the fragment
 *   set), we still separate preprocessing from the overlap matrix used by
 *   MGREEDY/TGREEDY/HK/BB. The preprocessing here uses KMP for containment
 *   detection; the matrix build also uses KMP for overlap computation.
 * ============================================================================ */

/* Load fragments from a file (or stdin if "-") into a FragmentArray */
static FragmentArray
read_all_fragments_array(const char *file_name)
{
    FILE *in = stdin;
    if (strcmp(file_name, "-") != 0) {
        in = fopen(file_name, "r");
        if (in == NULL) {
            fprintf(stderr, "Error: cannot open input file: %s\n", file_name);
            exit(1);
        }
    }
    FragmentArray fa = fa_new();
    while (1) {
        char *buf = NULL;
        size_t bufsz = 0;
        ssize_t nread = getline(&buf, &bufsz, in);
        if (nread <= 0) { free(buf); break; }
        for (ssize_t i = 0; i < nread; i++) {
            if (buf[i] == '\n' || buf[i] == '\r') buf[i] = '\0';
        }
        if (buf[0] == '\0') { free(buf); continue; }
        fa_push(&fa, buf);
    }
    if (in != stdin) fclose(in);
    return fa;
}

/* Remove any fragment that is a substring of another fragment.
 * Uses KMP for O(|text| + |pattern|) substring detection per pair.
 *
 * Per Tarhio & Ukkonen Section 2.2: "This step is naturally implemented
 * by embedding [it] into step (A1): whenever step (A1) finds out that x_i
 * is a substring of x_j, remove x_i from R."
 *
 * We use the same swap-with-last removal and equal-length duplicate rule
 * as the original implementation for identical output behaviour.
 *
 * Total time: O(n · m) where m = total length of all fragments. */
static void
remove_substring_fragments_array(FragmentArray *fa)
{
    size_t i = 0;
    while (i < fa->count) {
        bool is_substring = false;

        /* Build KMP failure function for fragment i (the candidate to test) */
        int *fail_i = kmp_build_failure(fa->items[i].str, fa->items[i].len);

        for (size_t j = 0; j < fa->count; j++) {
            if (i == j) continue;

            /* Skip: if items[i] is same length as items[j] and j > i,
               keep items[i] (remove the later duplicate, not the earlier) */
            if (fa->items[i].len == fa->items[j].len && j > i) continue;

            /* Can items[i] be a substring of items[j]?
               Only possible if |items[j]| >= |items[i]|. */
            if (fa->items[j].len < fa->items[i].len) continue;

            /* KMP scan: use items[j] as text, items[i] as pattern */
            bool contained = false;
            size_t q = 0;
            for (size_t k = 0; k < fa->items[j].len; k++) {
                while (q > 0 && fa->items[i].str[q] != fa->items[j].str[k])
                    q = (size_t)fail_i[q - 1];
                if (fa->items[i].str[q] == fa->items[j].str[k])
                    q++;
                if (q == fa->items[i].len) {
                    contained = true;
                    break;
                }
            }

            if (contained) {
                is_substring = true;
                break;
            }
        }

        free(fail_i);

        if (is_substring) {
            free(fa->items[i].str);
            fa->items[i] = fa->items[fa->count - 1];
            fa->count--;
        } else {
            i++;
        }
    }
}

static FragmentArray
run_step1(const char *file_name)
{
    FragmentArray fa = read_all_fragments_array(file_name);
    remove_substring_fragments_array(&fa);
    return fa;
}

/* ============================================================================
 * Section D: Step 2 — correct / sub-optimal / quick algorithms
 * ============================================================================ */

/* GREEDY: iteratively merge the max-overlap pair until one string remains.
 * Now uses KMP overlap (Complexity: worst O(n³ · L), average O(n³ · L)). */
static int
run_greedy(const FragmentArray *fa)
{
    double t0 = mono_seconds();
    size_t n = fa->count;

    Fragment *frags = malloc(n * sizeof(Fragment));
    for (size_t i = 0; i < n; i++) {
        frags[i].len = fa->items[i].len;
        frags[i].str = malloc(frags[i].len + 1);
        memcpy(frags[i].str, fa->items[i].str, frags[i].len + 1);
    }

    size_t result_len;
    char *result = greedy_merge_pairs(frags, n, &result_len);
    free(frags);

    return try_record_solution(result, result_len, "GREEDY", mono_seconds() - t0);
}

/* MGREEDY: model fragments as nodes of a weighted directed graph (edge weight
 * = overlap), build a max-weight path cover via degree-one greedy edge selection,
 * break cycles at weakest edge, concatenate paths.
 * Now uses KMP overlap matrix (Complexity: worst O(n² · (L + log n)),
 * average O(n² · (L + log n))). */
static int
run_mgreedy(const FragmentArray *fa)
{
    double t0 = mono_seconds();
    size_t n = fa->count;

    int **overlap  = build_overlap_matrix(fa);
    int **selected = mgreedy_select_edges(overlap, n);
    mgreedy_break_cycles(selected, n);

    size_t seq_count = 0;
    char **seqs = mgreedy_build_sequences(fa, selected, n, &seq_count);

    size_t result_len = 0;
    for (size_t i = 0; i < seq_count; i++) result_len += strlen(seqs[i]);

    char *result = malloc(result_len + 1);
    size_t pos = 0;
    for (size_t i = 0; i < seq_count; i++) {
        size_t len = strlen(seqs[i]);
        memcpy(result + pos, seqs[i], len);
        pos += len;
        free(seqs[i]);
    }
    result[result_len] = '\0';
    free_overlap_matrix(overlap, n);
    free_overlap_matrix(selected, n);
    free(seqs);

    return try_record_solution(result, result_len, "MGREEDY", mono_seconds() - t0);
}

/* TGREEDY: Run MGREEDY first to compress the n fragments into a small set
 * of path strings, then run GREEDY on those paths to produce one supersequence.
 * Now uses KMP throughout. */
static int
run_tgreedy(const FragmentArray *fa)
{
    double t0 = mono_seconds();
    size_t n = fa->count;

    int **overlap  = build_overlap_matrix(fa);
    int **selected = mgreedy_select_edges(overlap, n);
    mgreedy_break_cycles(selected, n);

    size_t seq_count = 0;
    char **seqs = mgreedy_build_sequences(fa, selected, n, &seq_count);
    free_overlap_matrix(overlap, n);
    free_overlap_matrix(selected, n);

    Fragment *frags = malloc(seq_count * sizeof(Fragment));
    for (size_t i = 0; i < seq_count; i++) {
        frags[i].str = seqs[i];
        frags[i].len = strlen(seqs[i]);
    }
    free(seqs);

    size_t result_len;
    char *result = greedy_merge_pairs(frags, seq_count, &result_len);
    free(frags);

    return try_record_solution(result, result_len, "TGREEDY", mono_seconds() - t0);
}

static void
run_step2(const FragmentArray *fa, uint8_t mask)
{
    if (mask & ALGO_GREEDY)  run_greedy(fa);
    if (mask & ALGO_MGREEDY) run_mgreedy(fa);
    if (mask & ALGO_TGREEDY) run_tgreedy(fa);
}

/* ============================================================================
 * Section E: Step 3 — correct / optimal / slow algorithms
 * ============================================================================ */

#define HK_THRESHOLD 20

/* HELD_KARP: Build a 2D table indexed by (subset of fragments, last fragment),
 * where each cell holds the shortest supersequence length for ordering that
 * subset ending at that fragment. Fill the table from smaller subsets to larger;
 * the answer is the minimum entry in the final row (subset = all fragments).
 * Restricted to n ≤ 20 (Complexity: worst O(2ⁿ · n²), average O(2ⁿ · n²)). */
static int
run_held_karp(const FragmentArray *fa)
{
    double t0 = mono_seconds();
    size_t n  = fa->count;

    int **ov   = build_overlap_matrix(fa);
    int  *flen = malloc(n * sizeof(int));
    for (size_t i = 0; i < n; i++) flen[i] = (int)fa->items[i].len;

    int    full   = (1 << n) - 1;
    size_t states = (size_t)1 << n;

    int *dp  = malloc(states * n * sizeof(int));
    int *par = malloc(states * n * sizeof(int));
    for (size_t i = 0; i < states * n; i++) { dp[i] = INT_MAX / 2; par[i] = -1; }

    for (int v = 0; v < (int)n; v++)
        dp[(1 << v) * n + v] = flen[v];

    size_t upper_bound = (g_best.str == NULL) ? SIZE_MAX : g_best.len;

    /* For each filled cell (S, v) in increasing-S order, extend by every
       u ∉ S, updating table[S∪{u}, u] if shorter (O(2ⁿ · n²)) */
    for (int S = 1; S <= full; S++) {
        for (int v = 0; v < (int)n; v++) {
            if (!(S & (1 << v))) continue;
            int dpSv = dp[S * n + v];
            if (dpSv >= INT_MAX / 2) continue;
            for (int u = 0; u < (int)n; u++) {
                if (S & (1 << u)) continue;
                int S2       = S | (1 << u);
                int new_cost = dpSv + flen[u] - ov[v][u];
                if ((size_t)new_cost >= upper_bound) continue;
                if (new_cost < dp[S2 * n + u]) {
                    dp[S2 * n + u]  = new_cost;
                    par[S2 * n + u] = v;
                }
            }
        }
    }

    /* In the final row (S = all fragments), pick the column with minimum length */
    int best_total = INT_MAX / 2, best_last = -1;
    for (int v = 0; v < (int)n; v++) {
        if (dp[full * n + v] < best_total) {
            best_total = dp[full * n + v];
            best_last  = v;
        }
    }
    if (best_last == -1) {
        free(dp); free(par); free(flen);
        free_overlap_matrix(ov, n);
        return try_record_solution(NULL, 0, "HELD_KARP", mono_seconds() - t0);
    }
    /* Walk parents backward from that column to recover the ordering (O(n)) */
    int *path = malloc(n * sizeof(int));
    {
        int S = full, curr = best_last;
        for (int i = (int)n - 1; i >= 0; i--) {
            path[i] = curr;
            int prev = par[S * n + curr];
            S ^= (1 << curr);
            curr = prev;
        }
    }
    /* Emit the string by concatenating fragments along the ordering with
       their overlaps (O(n · L)) */
    char  *result = malloc((size_t)best_total + 1);
    size_t pos    = fa->items[path[0]].len;
    memcpy(result, fa->items[path[0]].str, pos);
    for (int i = 1; i < (int)n; i++) {
        int    u    = path[i-1], v2 = path[i];
        size_t skip = (size_t)ov[u][v2];
        memcpy(result + pos, fa->items[v2].str + skip, fa->items[v2].len - skip);
        pos += fa->items[v2].len - skip;
    }
    result[pos] = '\0';

    free(path);
    free(dp); free(par); free(flen);
    free_overlap_matrix(ov, n);

    return try_record_solution(result, (size_t)best_total, "HELD_KARP", mono_seconds() - t0);
}

/* BRANCH_AND_BOUND: Build every possible ordering of the fragments. Discard
 * any branch whose lower bound cannot beat the current best length. The first
 * ordering that uses every fragment is the optimal answer.
 * Restricted to n ≤ 63 (Complexity: worst O(n! · n³), average O(n! · n³)). */
static int
run_branch_and_bound(const FragmentArray *fa)
{
    double t0 = mono_seconds();
    size_t n  = fa->count;

    int **ov   = build_overlap_matrix(fa);
    int  *flen = malloc(n * sizeof(int));
    for (size_t i = 0; i < n; i++) flen[i] = (int)fa->items[i].len;

    size_t upper_bound = (g_best.str == NULL) ? SIZE_MAX : g_best.len;
    bool improved = false;

    MinHeap pq = { NULL, 0, 0 };

    /* Seed with each fragment as the start */
    for (int start = 0; start < (int)n; start++) {
        uint64_t mask = (uint64_t)1 << start;
        int lb = flen[start] + remaining_lb(mask, (int)n, ov, flen);
        if ((size_t)lb >= upper_bound) continue;

        BBState s;
        s.lb       = lb;
        s.mask     = mask;
        s.last     = start;
        s.curlen   = flen[start];
        s.path_len = 1;
        s.path     = malloc(sizeof(int));
        s.path[0]  = start;
        mh_push(&pq, s);
    }

    uint64_t full_mask = ((uint64_t)1 << n) - 1;

    while (pq.size > 0) {
        BBState cur = mh_pop(&pq);

        if ((size_t)cur.lb >= upper_bound) {
            free(cur.path);
            continue;
        }

        if (cur.mask == full_mask) {
            /* Reconstruct the superstring */
            char  *result = malloc((size_t)cur.curlen + 1);
            size_t rlen   = fa->items[cur.path[0]].len;
            memcpy(result, fa->items[cur.path[0]].str, rlen);
            for (int i = 1; i < cur.path_len; i++) {
                int u = cur.path[i-1], v = cur.path[i];
                size_t skip = (size_t)ov[u][v];
                memcpy(result + rlen, fa->items[v].str + skip,
                       fa->items[v].len - skip);
                rlen += fa->items[v].len - skip;
            }
            result[rlen] = '\0';

            if (try_record_solution(result, rlen, "BRANCH_AND_BOUND",
                                    mono_seconds() - t0)) {
                improved = true;
                upper_bound = g_best.len;
            }
            free(cur.path);
            continue;
        }

        for (int next = 0; next < (int)n; next++) {
            if (cur.mask & ((uint64_t)1 << next)) continue;

            uint64_t new_mask   = cur.mask | ((uint64_t)1 << next);
            int      new_curlen = cur.curlen + flen[next] - ov[cur.last][next];
            int      new_lb     = new_curlen + remaining_lb(new_mask, (int)n, ov, flen);

            if ((size_t)new_lb >= upper_bound) continue;

            BBState s;
            s.lb       = new_lb;
            s.mask     = new_mask;
            s.last     = next;
            s.curlen   = new_curlen;
            s.path_len = cur.path_len + 1;
            s.path     = malloc((size_t)s.path_len * sizeof(int));
            memcpy(s.path, cur.path, (size_t)cur.path_len * sizeof(int));
            s.path[cur.path_len] = next;
            mh_push(&pq, s);
        }
        free(cur.path);
    }

    for (int i = 0; i < pq.size; i++) free(pq.data[i].path);
    free(pq.data);
    free(flen);
    free_overlap_matrix(ov, n);

    if (!improved) {
        try_record_solution(NULL, 0, "BRANCH_AND_BOUND", mono_seconds() - t0);
    }
    return improved;
}

static void
run_step3(const FragmentArray *fa, uint8_t mask)
{
    bool want_hk = (mask & ALGO_HK) && fa->count <= HK_THRESHOLD;
    bool want_bb = (mask & ALGO_BB) && fa->count <= 63;

    if      (want_hk) run_held_karp(fa);
    else if (want_bb) run_branch_and_bound(fa);
}

/* ============================================================================
 * Section F: main()
 * ============================================================================ */

int
main(int argc, char *argv[])
{
    uint8_t     mask       = ALGO_ALL;
    const char *input_path = NULL;

    for (int i = 1; i < argc; i++) {
        if (strcmp(argv[i], "--algos") == 0) {
            if (i + 1 >= argc) {
                fprintf(stderr, "Error: --algos requires an argument\n");
                exit(1);
            }
            int parsed = parse_algos(argv[++i]);
            if (parsed < 0) {
                fprintf(stderr,
                        "Error: invalid --algos value '%s'. "
                        "Valid tokens: greedy,mgreedy,tgreedy,hk,bb\n",
                        argv[i]);
                exit(1);
            }
            mask = (uint8_t)parsed;
        } else if (input_path == NULL) {
            input_path = argv[i];
        } else {
            fprintf(stderr, "Error: unexpected argument '%s'\n", argv[i]);
            exit(1);
        }
    }

    if (input_path == NULL) {
        fprintf(stderr,
                "Usage: %s [--algos LIST] [ <input_file> | - ]\n"
                "  --algos LIST: subset of {greedy,mgreedy,tgreedy,hk,bb} (default: all).\n"
                "  Input: one fragment per line (no blank lines).\n"
                "  stdout: reconstructed text (non-increasing length).\n"
                "  stderr: per-algorithm and quality/test metrics.\n",
                argv[0]);
        exit(1);
    }

    FragmentArray fa = run_step1(input_path);

    if (fa.count == 0) {
        fputc('\n', stdout);
        fflush(stdout);
        fa_free(&fa);
        return 0;
    }

    if (fa.count == 1) {
        char *copy = malloc(fa.items[0].len + 1);
        memcpy(copy, fa.items[0].str, fa.items[0].len + 1);
        try_record_solution(copy, fa.items[0].len, "trivial", 0.0);
    } else {
        run_step2(&fa, mask);
        run_step3(&fa, mask);
    }

    fa_free(&fa);
    free(g_best.str);
    return 0;
}