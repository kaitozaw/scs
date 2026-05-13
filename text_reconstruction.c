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
 *   B. Best-solution registry
 *   C. Step 1 — Preprocess
 *   D. Step 2 — correct / sub-optimal / quick algorithms
 *   E. Step 3 — correct / optimal / slow algorithms
 *   F. main()
 * ============================================================================ */

/* ============================================================================
 * Section A: Common types & utilities
 * ============================================================================ */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include <assert.h>
#include <time.h>
#include <stdint.h>

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

/* Max k > 0 such that suffix(a, k) == prefix(b, k); 0 if no overlap. */
static size_t
overlap_chars(const char *a, size_t la, const char *b, size_t lb)
{
    size_t max_k = (la < lb ? la : lb);
    if (max_k > 0) max_k--;
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
            m[i][j] = (int)overlap_chars(fa->items[i].str, fa->items[i].len, fa->items[j].str, fa->items[j].len);
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
 * Section B: Best-solution registry
 * ============================================================================ */

typedef struct {
    char  *str;
    size_t len;
} BestSolution;

static BestSolution g_best = { NULL, 0 };

static double
mono_seconds(void)
{
    return (double)clock() / CLOCKS_PER_SEC;
}

/* Publish `candidate` if shorter than g_best; enforces non-increasing length on stdout. */
static int
try_record_solution(char *candidate, size_t cand_len, const char *algo_label, double elapsed_sec)
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
    fflush(stdout);
    fprintf(stderr, "[algo=%s] len=%zu elapsed=%.3fs\n", algo_label, cand_len, elapsed_sec);
    return 1;
}

/* ============================================================================
 * Section C: Step 1 — Preprocess
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
        if (buf[0] == '\0') { free(buf); continue; }
        fa_push(&fa, buf);
    }
    if (in != stdin) fclose(in);
    return fa;
}

/* Remove any fragment that is a substring of another fragment. */
static void
remove_substring_fragments_array(FragmentArray *fa)
{
    size_t i = 0;
    while (i < fa->count) {
        bool is_substring = false;
        for (size_t j = 0; j < fa->count; j++) {
            if (i == j) continue;
            if (strstr(fa->items[j].str, fa->items[i].str) != NULL) {
                if (fa->items[i].len == fa->items[j].len && j > i) continue;
                is_substring = true;
                break;
            }
        }
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
preprocess(const char *file_name)
{
    FragmentArray fa = read_all_fragments_array(file_name);
    remove_substring_fragments_array(&fa);
    return fa;
}

/* ============================================================================
 * Section D: Step 2 — correct / sub-optimal / quick algorithms
 * ============================================================================ */

static int
run_greedy(const FragmentArray *fa)
{
    double t0 = mono_seconds();
    /* TODO: GREEDY — repeatedly merge the pair with the maximum overlap. */
    /* End with: try_record_solution(result, result_len, "GREEDY", mono_seconds() - t0); */
    (void)fa; (void)t0;
    return 0;
}

static int
run_mgreedy(const FragmentArray *fa)
{
    double t0 = mono_seconds();
    /* TODO: MGREEDY — modified greedy. */
    /* End with: try_record_solution(result, result_len, "MGREEDY", mono_seconds() - t0); */
    (void)fa; (void)t0;
    return 0;
}

static int
run_tgreedy(const FragmentArray *fa)
{
    double t0 = mono_seconds();
    /* TODO: TGREEDY — tour/path-based greedy. */
    /* End with: try_record_solution(result, result_len, "TGREEDY", mono_seconds() - t0); */
    (void)fa; (void)t0;
    return 0;
}

static void
run_step2(const FragmentArray *fa)
{
    run_greedy(fa);
    run_mgreedy(fa);
    run_tgreedy(fa);
}

/* ============================================================================
 * Section E: Step 3 — correct / optimal / slow algorithms
 * ============================================================================ */

#define HK_THRESHOLD 20

static bool
use_held_karp(const FragmentArray *fa)
{
    return fa->count <= HK_THRESHOLD;
}

static int
run_held_karp(const FragmentArray *fa)
{
    double t0 = mono_seconds();
    /* TODO: Held-Karp. */
    /* End with: try_record_solution(result, result_len, "HELD_KARP", mono_seconds() - t0); */
    (void)fa; (void)t0;
    return 0;
}

static int
run_branch_and_bound(const FragmentArray *fa)
{
    double t0 = mono_seconds();
    /* TODO: Branch & Bound. */
    /* End with: try_record_solution(result, result_len, "BRANCH_AND_BOUND", mono_seconds() - t0); */
    (void)fa; (void)t0;
    return 0;
}

static void
run_step3(const FragmentArray *fa)
{
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

    FragmentArray fa = preprocess(argv[1]);

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
        run_step2(&fa);
        run_step3(&fa);
    }

    fa_free(&fa);
    free(g_best.str);
    return 0;
}
