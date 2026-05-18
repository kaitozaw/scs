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
 *   C. Step 1 — preprocess
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
#include <time.h>
#include <stdint.h>
#include <limits.h>

/* --- A.1: Fragments & string utilities --- */

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

/* Find the longest suffix-prefix overlap between a and b (O(L²)) */
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

/* --- A.2: Greedy pairwise merging --- */

/* Loop until one string remains (n − 1 iterations): scan all pairs to find the maximum-overlap pair; merge them into a single string (O(n² · L²) per iteration). */
static char *
greedy_merge_pairs(Fragment *frags, size_t n, size_t *out_len)
{
    while (n > 1) {
        size_t best_i = 0, best_j = 1, best_ov = 0;
        for (size_t i = 0; i < n; i++) {
            for (size_t j = 0; j < n; j++) {
                if (i == j) continue;
                size_t ov = overlap_chars(frags[i].str, frags[i].len, frags[j].str, frags[j].len);
                if (ov > best_ov) { best_ov = ov; best_i = i; best_j = j; }
            }
        }
        char *merged = merge_with_overlap(frags[best_i].str, frags[best_i].len, frags[best_j].str, frags[best_j].len, best_ov);
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

/* --- A.3: Overlap matrix & edge sort --- */

/* Build the overlap matrix (O(n² · L²)) */
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

/* --- A.4: Cycle-cover graph helpers --- */

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

/* --- A.5: Branch & Bound primitives (search state, min-heap, lower bound) --- */

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

static double
mono_seconds(void)
{
    return (double)clock() / CLOCKS_PER_SEC;
}

/* Publish `candidate` if it beats g_best; NULL only logs to stderr. */
static int
try_record_solution(char *candidate, size_t cand_len, const char *algo_label, double elapsed_sec)
{
    if (candidate == NULL) {
        fprintf(stderr, "[algo=%s] len=%zu elapsed=%.3fs (no improvement over g_best)\n", algo_label, g_best.len, elapsed_sec);
        return 0;
    }

    bool improved = (g_best.str == NULL) || (cand_len < g_best.len);
    fprintf(stderr, "[algo=%s] len=%zu elapsed=%.3fs result=%s\n", algo_label, cand_len, elapsed_sec, candidate);
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
run_step1(const char *file_name)
{
    FragmentArray fa = read_all_fragments_array(file_name);
    remove_substring_fragments_array(&fa);
    return fa;
}

/* ============================================================================
 * Section D: Step 2 — correct / sub-optimal / quick algorithms
 * ============================================================================ */

/* GREEDY: iteratively merge the max-overlap pair until one string remains (Complexity: worst O(n³ · L²), average O(n³ · L)). */
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

/* MGREEDY: model fragments as nodes of a weighted directed graph (edge weight = overlap), build a max-weight path cover via degree-one greedy edge selection, break cycles at weakest edge, concatenate paths (Complexity: worst O(n² · (L² + log n)), average O(n² · (L + log n))). */
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

/* TGREEDY: Run MGREEDY first to compress the n fragments into a small set of path strings, then run GREEDY on those paths to produce one supersequence (Complexity: worst O(n² · (L² + log n) + p³ · L′²), average O(n² · (L + log n) + p³ · L′)). */
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

/* HELD_KARP: Build a 2D table indexed by (subset of fragments, last fragment), where each cell holds the shortest supersequence length for ordering that subset ending at that fragment. Fill the table from smaller subsets to larger; the answer is the minimum entry in the final row (subset = all fragments). Restricted to n ≤ 20 (Complexity: worst O(2ⁿ · n²), average O(2ⁿ · n²)). */
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

    size_t upper_bound = g_best.len;

    /* For each filled cell (S, v) in increasing-S order, extend by every u ∉ S, updating table[S∪{u}, u] if shorter (O(2ⁿ · n²)) */
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

    /* In the final row (S = all fragments), pick the column with the minimum length (O(n)) */
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
    /* Emit the string by concatenating fragments along the ordering with their overlaps (O(n · L)) */
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

/* BRANCH_AND_BOUND: Build every possible ordering of the fragments. Discard any branch whose lower bound cannot beat the current best length. The first ordering that uses every fragment is the optimal answer. Restricted to n ≤ 63 (Complexity: worst O(n! · n³), average O(n! · n³)). */
static int
run_branch_and_bound(const FragmentArray *fa)
{
    double t0 = mono_seconds();
    size_t n  = fa->count;

    if (n > 63) return 0;

    int **ov   = build_overlap_matrix(fa);
    int  *flen = malloc(n * sizeof(int));
    for (size_t i = 0; i < n; i++) flen[i] = (int)fa->items[i].len;

    size_t   upper_bound = g_best.len;
    uint64_t full_mask   = ((uint64_t)1 << n) - 1;
    MinHeap  pq          = { NULL, 0, 0 };
    int      improved    = 0;

    /* For each fragment, push an initial state into the min-heap (state = lower bound, visited fragments, last fragment, current length) (O(n³)) */
    for (int start = 0; start < (int)n; start++) {
        uint64_t mask   = (uint64_t)1 << start;
        int      curlen = flen[start];
        int      lb     = curlen + remaining_lb(mask, (int)n, ov, flen);
        if ((size_t)lb >= upper_bound) continue;

        BBState s;
        s.lb       = lb;
        s.mask     = mask;
        s.last     = start;
        s.curlen   = curlen;
        s.path_len = 1;
        s.path     = malloc(sizeof(int));
        s.path[0]  = start;
        mh_push(&pq, s);
    }

    /* Loop until every possible ordering is explored (n! iterations): pop the lowest-bound state; if complete, return it as the answer; else extend by every unvisited fragment, pushing children that beat the current best (O(n³) per iteration) */
    while (pq.size > 0) {
        BBState cur = mh_pop(&pq);

        if ((size_t)cur.lb >= upper_bound) { free(cur.path); continue; }

        if (cur.mask == full_mask) {
            size_t rlen = (size_t)cur.curlen;
            if (rlen < upper_bound) {
                int   *path   = cur.path;
                char  *result = malloc(rlen + 1);
                size_t pos    = fa->items[path[0]].len;
                memcpy(result, fa->items[path[0]].str, pos);
                for (int i = 1; i < cur.path_len; i++) {
                    int    u    = path[i-1], v2 = path[i];
                    size_t skip = (size_t)ov[u][v2];
                    memcpy(result + pos, fa->items[v2].str + skip, fa->items[v2].len - skip);
                    pos += fa->items[v2].len - skip;
                }
                result[pos] = '\0';
                upper_bound  = rlen;
                improved    |= try_record_solution(result, rlen, "BRANCH_AND_BOUND", mono_seconds() - t0);
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

    FragmentArray fa = run_step1(argv[1]);

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
