/* ============================================================================
 * IFN664 Advanced Algorithms & Computational Complexity
 * Assignment 1 — Shortest Common Supersequence
 *
 * Compile:
 *   gcc -o text_reconstruction text_reconstruction.c
 *
 * Usage:
 *   ./text_reconstruction <input_file>          (or '-' for stdin)
 *
 * Output protocol:
 *   stdout : reconstructed text, one solution per line, in non-increasing length order.
 *   stderr : per-algorithm metrics:  [algo=NAME] len=N elapsed=T.TTTs
 *
 * Architecture:
 *   A. Common types & utilities
 *   B. SIGINT handling & best-solution registry
 *   C. Step 1 — Preprocess
 *   D. Step 2 — correct / sub-optimal / quick algorithms
 *   E. Step 3 — correct / optimal / slow algorithms
 *   F. main()
 *
 * Step 3 is exponential — Ctrl+C interrupts cleanly via g_interrupted;
 * best-so-far stays on stdout.
 * ============================================================================ */

/* ============================================================================
 * Section A: Common types & utilities
 * ============================================================================ */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include <assert.h>
#include <signal.h>
#include <time.h>
#include <stdint.h>

typedef struct {
    char  *str;     /* owned, null-terminated */
    size_t len;     /* cached strlen(str) */
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

/* Append `owned_str` to fa. The array takes ownership; do not free it elsewhere. */
static void
fa_push(FragmentArray *fa, char *owned_str)
{
    if (fa->count == fa->capacity) {
        size_t new_cap = fa->capacity == 0 ? 8 : fa->capacity * 2;
        fa->items = realloc(fa->items, new_cap * sizeof(Fragment));
        fa->capacity = new_cap;
    }
    fa->items[fa->count].str = owned_str;
    fa->items[fa->count].len = strlen(owned_str);
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

/* Max k > 0 such that suffix(a, k) == prefix(b, k). Returns 0 if no overlap.
 * Never returns k == la or k == lb (substring relations are removed in Step 1). */
static size_t
overlap_chars(const char *a, size_t la, const char *b, size_t lb)
{
    size_t max_k = (la < lb ? la : lb);
    if (max_k > 0) max_k--;     /* exclude full-prefix / full-suffix matches */
    for (size_t k = max_k; k > 0; k--) {
        if (memcmp(a + la - k, b, k) == 0) return k;
    }
    return 0;
}

/* Allocate a + b[ov..]. Result length = la + lb - ov. Caller frees. */
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

/* n x n matrix; M[i][j] = overlap_chars(items[i], items[j]); diagonal is 0. */
static int **
build_overlap_matrix(const FragmentArray *fa)
{
    size_t n = fa->count;
    int **m = malloc(n * sizeof(int *));
    for (size_t i = 0; i < n; i++) m[i] = calloc(n, sizeof(int));
    for (size_t i = 0; i < n; i++) {
        for (size_t j = 0; j < n; j++) {
            if (i == j) continue;
            m[i][j] = (int)overlap_chars(fa->items[i].str, fa->items[i].len,
                                         fa->items[j].str, fa->items[j].len);
        }
    }
    return m;
}

static void
free_overlap_matrix(int **m, size_t n)
{
    for (size_t i = 0; i < n; i++) free(m[i]);
    free(m);
}

/* ============================================================================
 * Section B: SIGINT handling & best-solution registry
 * ============================================================================ */

static volatile sig_atomic_t g_interrupted = 0;

typedef struct {
    char  *str;     /* owned; NULL means "no solution yet" */
    size_t len;
} BestSolution;

static BestSolution g_best = { NULL, 0 };

static void
sigint_handler(int signo)
{
    (void)signo;
    g_interrupted = 1;   /* async-signal-safe: only this assignment */
}

static void
install_sigint_handler(void)
{
    struct sigaction sa;
    sa.sa_handler = sigint_handler;
    sigemptyset(&sa.sa_mask);
    sa.sa_flags = 0;    /* no SA_RESTART: blocking syscalls return EINTR */
    sigaction(SIGINT, &sa, NULL);
}

static double
mono_seconds(void)
{
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return (double)ts.tv_sec + (double)ts.tv_nsec * 1e-9;
}

/* Submit a candidate solution. Takes ownership of `candidate`.
 *   - If strictly shorter than g_best (or g_best is empty):
 *       publish to stdout, log to stderr, store in g_best, return 1.
 *   - Otherwise: free `candidate` and return 0.
 * Centralising this is what guarantees the PDF's "non-increasing length" rule. */
static int
try_record_solution(char *candidate, size_t cand_len,
                    const char *algo_label, double elapsed_sec)
{
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
    fflush(stdout);    /* ensure visibility even if interrupted mid-Step 3 */
    fprintf(stderr, "[algo=%s] len=%zu elapsed=%.3fs\n",
            algo_label, cand_len, elapsed_sec);
    return 1;
}

/* ============================================================================
 * Section C: Step 1 — Preprocess (COMPLETE)
 *
 *   1. Load fragments line-by-line into a FragmentArray.
 *   2. Remove any fragment that is already a substring of another.
 *   3. Extract simple-path components from the overlap graph: any chain
 *      where every internal edge is "forced" (out-degree 1 ∧ in-degree 1)
 *      is collapsed into a single merged fragment. Cycles are left intact.
 * ============================================================================ */

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
        if (buf[0] == '\0') { free(buf); continue; }   /* skip blank lines */
        fa_push(&fa, buf);    /* fa takes ownership of buf */
    }
    if (in != stdin) fclose(in);
    return fa;
}

/* Remove any fragment that is a substring of another. Order is not preserved
 * (swap-with-last for O(1) deletion). Duplicate fragments keep exactly one. */
static void
remove_substring_fragments_array(FragmentArray *fa)
{
    size_t i = 0;
    while (i < fa->count) {
        bool is_substring = false;
        for (size_t j = 0; j < fa->count; j++) {
            if (i == j) continue;
            if (strstr(fa->items[j].str, fa->items[i].str) != NULL) {
                /* If equal-length duplicate, only the higher-index copy is
                 * removed; the j>i guard ensures we keep one representative. */
                if (fa->items[i].len == fa->items[j].len && j > i) continue;
                is_substring = true;
                break;
            }
        }
        if (is_substring) {
            free(fa->items[i].str);
            fa->items[i] = fa->items[fa->count - 1];
            fa->count--;
            /* re-check the swapped-in element at index i, do not advance */
        } else {
            i++;
        }
    }
}

/* Collapse forced-merge chains in the overlap graph.
 *
 * Definition: edge (i, j) is "forced" iff out_count[i] == 1 AND in_count[j] == 1.
 * Forced edges chained together form a simple directed path; collapsing such
 * paths is loss-free (any optimal SCS that uses any edge from i must use (i,j),
 * and symmetrically for j).
 *
 * Note: forced edges are NOT vertex-disjoint in general — they share endpoints
 * within a chain. The traversal below handles chains of arbitrary length.
 * Pure cycles (every vertex has next AND prev set, no head exists) are skipped
 * — they would lose information if collapsed; Step 2/3 will resolve them. */
static void
extract_simple_paths(FragmentArray *fa)
{
    if (fa->count < 2) return;

    size_t n = fa->count;
    int **m = build_overlap_matrix(fa);

    int *out_count = calloc(n, sizeof(int));
    int *in_count  = calloc(n, sizeof(int));
    int *out_only  = malloc(n * sizeof(int));
    int *in_only   = malloc(n * sizeof(int));
    for (size_t i = 0; i < n; i++) { out_only[i] = -1; in_only[i] = -1; }
    for (size_t i = 0; i < n; i++) {
        for (size_t j = 0; j < n; j++) {
            if (m[i][j] > 0) {
                out_count[i]++; out_only[i] = (int)j;
                in_count[j]++;  in_only[j]  = (int)i;
            }
        }
    }

    /* next_v[i] is i's forced successor (or -1); prev_v[j] is j's forced
     * predecessor (or -1). Either both endpoints participate, or neither. */
    int *next_v = malloc(n * sizeof(int));
    int *prev_v = malloc(n * sizeof(int));
    for (size_t i = 0; i < n; i++) { next_v[i] = -1; prev_v[i] = -1; }
    for (size_t i = 0; i < n; i++) {
        if (out_count[i] == 1) {
            int j = out_only[i];
            if (in_count[j] == 1) {
                next_v[i] = j;
                prev_v[j] = (int)i;
            }
        }
    }

    /* Find chain heads: next_v != -1 ∧ prev_v == -1.
     * Walk each chain to mark all its members; vertices not visited but with
     * pointers set belong to cycles, which we leave intact. */
    bool *in_chain = calloc(n, sizeof(bool));
    bool *is_head  = calloc(n, sizeof(bool));
    bool  any_chain = false;
    for (size_t i = 0; i < n; i++) {
        if (next_v[i] != -1 && prev_v[i] == -1) {
            is_head[i] = true;
            any_chain  = true;
            int cur = (int)i;
            while (cur != -1) {
                in_chain[cur] = true;
                cur = next_v[cur];
            }
        }
    }

    if (!any_chain) {
        free(out_count); free(in_count); free(out_only); free(in_only);
        free(next_v); free(prev_v); free(in_chain); free(is_head);
        free_overlap_matrix(m, n);
        return;
    }

    FragmentArray next_fa = fa_new();

    /* Emit one merged fragment per chain. */
    for (size_t i = 0; i < n; i++) {
        if (!is_head[i]) continue;
        char  *acc = malloc(fa->items[i].len + 1);
        memcpy(acc, fa->items[i].str, fa->items[i].len + 1);
        size_t acc_len = fa->items[i].len;
        int cur = (int)i;
        while (next_v[cur] != -1) {
            int nxt = next_v[cur];
            int ov  = m[cur][nxt];   /* tail of acc == tail of fragment `cur` */
            char *merged = merge_with_overlap(acc, acc_len,
                                              fa->items[nxt].str,
                                              fa->items[nxt].len,
                                              (size_t)ov);
            free(acc);
            acc = merged;
            acc_len = acc_len + fa->items[nxt].len - (size_t)ov;
            cur = nxt;
        }
        fa_push(&next_fa, acc);
    }

    /* Keep every vertex not absorbed into a chain (isolated vertices and
     * cycle members). Duplicate the strings so fa_free below is safe. */
    for (size_t k = 0; k < n; k++) {
        if (in_chain[k]) continue;
        char *copy = malloc(fa->items[k].len + 1);
        memcpy(copy, fa->items[k].str, fa->items[k].len + 1);
        fa_push(&next_fa, copy);
    }

    fa_free(fa);
    *fa = next_fa;

    free(out_count); free(in_count); free(out_only); free(in_only);
    free(next_v); free(prev_v); free(in_chain); free(is_head);
    free_overlap_matrix(m, n);
}

/* Step 1 orchestrator. */
static FragmentArray
preprocess(const char *file_name)
{
    FragmentArray fa = read_all_fragments_array(file_name);
    remove_substring_fragments_array(&fa);
    extract_simple_paths(&fa);
    return fa;
}

/* ============================================================================
 * Section D: Step 2 — correct / sub-optimal / quick algorithms (TEAM TODO)
 *
 * Each algorithm in this section is guaranteed to return a CORRECT SCS
 * (every input fragment appears in the output) and to run QUICKLY (polynomial
 * time), but the result may be SUB-OPTIMAL (longer than the true minimum).
 *
 * Each function should:
 *   1. Build any local data structures it needs (e.g., overlap matrix).
 *   2. Compute its candidate SCS string.
 *   3. Submit via try_record_solution(result, strlen(result), "NAME", elapsed).
 *      Do NOT call printf/fputs directly — that bypasses the non-increasing
 *      length guarantee.
 *   4. Return 0 on success.
 * ============================================================================ */

static int
run_greedy(const FragmentArray *fa)
{
    double t0 = mono_seconds();
    /* TODO: GREEDY — repeatedly merge the pair with the maximum overlap.
     *
     * Sketch:
     *   m = build_overlap_matrix(fa)
     *   working_set = copy of fa
     *   while working_set.count > 1:
     *       find (i, j) maximising m[i][j]
     *       replace items i and j with merge_with_overlap(items[i], items[j], m[i][j])
     *       rebuild affected rows/cols of m (or rebuild matrix)
     *   result = concatenation of remaining fragments
     *   try_record_solution(result, strlen(result), "GREEDY", mono_seconds()-t0)
     */
    (void)fa; (void)t0;
    return 0;
}

static int
run_mgreedy(const FragmentArray *fa)
{
    double t0 = mono_seconds();
    /* TODO: MGREEDY — modified greedy. */
    (void)fa; (void)t0;
    return 0;
}

static int
run_tgreedy(const FragmentArray *fa)
{
    double t0 = mono_seconds();
    /* TODO: TGREEDY — tour/path-based greedy. */
    (void)fa; (void)t0;
    return 0;
}

/* ============================================================================
 * Section E: Step 3 — correct / optimal / slow algorithms (TEAM TODO)
 *
 * Each algorithm in this section is guaranteed to return a CORRECT and
 * OPTIMAL SCS (minimum length amongst all correct reconstructions), but its
 * complexity is EXPONENTIAL — it may be SLOW or fail to terminate in
 * reasonable time on large inputs.
 *
 * Both algorithms must periodically check `g_interrupted` and return cleanly
 * when it is set. Any improvement already published via try_record_solution
 * is safe on stdout (fflush'd line by line).
 * ============================================================================ */

#define HK_THRESHOLD 20    /* 2^20 * 20^2 ≈ 4e8 ops, borderline feasible */

/* Helper: decide which optimal algorithm to dispatch.
 * Default heuristic — refine as the team gains insight. */
static bool
use_held_karp(const FragmentArray *fa)
{
    /* TODO: team to refine. Could also factor in fragment lengths, available
     * memory, alphabet size, or empirical timings. */
    return fa->count <= HK_THRESHOLD;
}

static int
run_held_karp(const FragmentArray *fa)
{
    double t0 = mono_seconds();
    /* TODO: Held-Karp bitmask DP.
     *
     * Sketch:
     *   n = fa->count
     *   m = build_overlap_matrix(fa)
     *   dp[mask][last] = min total length to cover set `mask` ending at `last`
     *   transition:
     *     dp[mask | (1<<j)][j] = min over i in mask of
     *         dp[mask][i] + len[j] - m[i][j]
     *
     *   Interrupt check: at the top of the outer subset-size loop:
     *     if (g_interrupted) { ... free local state ... return 0; }
     *
     *   On completion: reconstruct the order, build the string, submit:
     *     try_record_solution(result, result_len, "HELD_KARP", mono_seconds()-t0);
     */
    (void)fa; (void)t0;
    return 0;
}

static int
run_branch_and_bound(const FragmentArray *fa)
{
    double t0 = mono_seconds();
    /* TODO: Branch & Bound.
     *
     * Sketch:
     *   initial_bound = sum of fragment lengths (concatenation length)
     *   state = (chosen prefix order, remaining set, current length)
     *   recurse / use an explicit stack:
     *     prune when current_length + lower_bound(remaining) >= best_known
     *     on reaching a leaf with a complete order:
     *       try_record_solution(...)   ← publishes improvement incrementally
     *
     *   Interrupt check: at the top of each node expansion, before computing
     *   children's bounds:
     *     if (g_interrupted) return 0;
     */
    (void)fa; (void)t0;
    return 0;
}

static void
run_step3(const FragmentArray *fa)
{
    if (g_interrupted) return;
    if (use_held_karp(fa)) run_held_karp(fa);
    else                   run_branch_and_bound(fa);
}

/* ============================================================================
 * Section F: main()
 * ============================================================================ */

int
main(int argc, char *argv[])
{
    if (argc != 2) {
        fprintf(stderr,
                "Usage: %s [ <input_file> | - ]\n"
                "  Input: one fragment per line (no blank lines).\n"
                "  stdout: reconstructed text (non-increasing length).\n"
                "  stderr: per-algorithm metrics.\n",
                argv[0]);
        exit(1);
    }

    install_sigint_handler();

    FragmentArray fa = preprocess(argv[1]);

    if (fa.count == 0) {
        fputc('\n', stdout);
        fflush(stdout);
        fa_free(&fa);
        return 0;
    }

    if (fa.count == 1) {
        /* Preprocessing collapsed everything (or the input had one fragment):
         * the remaining fragment IS the optimal SCS. Skip Step 2/3. */
        char *copy = malloc(fa.items[0].len + 1);
        memcpy(copy, fa.items[0].str, fa.items[0].len + 1);
        try_record_solution(copy, fa.items[0].len, "trivial", 0.0);
    } else {
        run_greedy(&fa);
        run_mgreedy(&fa);
        run_tgreedy(&fa);
        run_step3(&fa);
    }

    if (g_interrupted) {
        fprintf(stderr, "interrupted; best len=%zu\n", g_best.len);
    }

    fa_free(&fa);
    free(g_best.str);
    return 0;
}
