# Hybrid SCS Solver

> **Algorithm:** `Hybrid_SCS_Solver`
> **Input:** A set *F* of text fragments
> **Output:** A shortest common superstring of all fragments in *F*, output progressively *(anytime behaviour)*

---

## Phase 0 — Input & Initialisation

```
FUNCTION main(filename):
    fragments = READ_FRAGMENTS(filename)
    // Read each line of file as a fragment string
    // Skip empty lines
    // Store in array fragments[0..m-1]

    IF |fragments| == 0:
        OUTPUT "No fragments provided"
        RETURN

    // Remove exact duplicate fragments immediately
    fragments = REMOVE_DUPLICATES(fragments)
    m = |fragments|

    // Begin solving — anytime output starts here
    CALL Hybrid_SCS(fragments, m)
```

---

## Phase 1 — Greedy Solution
> *Fast, sub-optimal — output immediately for anytime behaviour*

```
FUNCTION Hybrid_SCS(fragments[], m):

    // Step 1a: Remove subsumed fragments
    // A fragment f is subsumed if it appears as a substring of another
    // fragment g (g already covers f for free)
    active_frags = REMOVE_SUBSTRINGS(fragments, m)
    n = |active_frags|      // n <= m; potentially much smaller

    // Step 1b: Build overlap and cost matrices
    FOR i = 0 TO n-1:
        FOR j = 0 TO n-1:
            IF i != j:
                ov[i][j]   = OVERLAP(active_frags[i], active_frags[j])
                cost[i][j] = |active_frags[j]| - ov[i][j]
            ELSE:
                ov[i][i]   = 0
                cost[i][i] = INFINITY

    // Step 1c: Greedy overlap merging
    // Analogous to Kruskal's MST:
    // repeatedly merge the pair with greatest overlap
    greedy_result = GREEDY_MERGE(active_frags, n, ov, cost)

    OUTPUT greedy_result        // *** ANYTIME OUTPUT #1 ***
                                // User can stop here if good enough

    // Proceed to exact search to improve the solution
    CALL Exact_Search(active_frags, n, ov, cost, |greedy_result|)
```

### Helper: `REMOVE_SUBSTRINGS`

```
FUNCTION REMOVE_SUBSTRINGS(fragments[], m):
    active = [TRUE] * m
    FOR i = 0 TO m-1:
        FOR j = 0 TO m-1:
            IF i != j AND active[i] AND active[j]:
                IF fragments[i] is substring of fragments[j]:
                    active[i] = FALSE
                    BREAK
    RETURN { fragments[i] : active[i] == TRUE }
```

### Helper: `OVERLAP`

```
FUNCTION OVERLAP(a, b):
    // Length of longest suffix of a that is a prefix of b
    max_k = min(|a|, |b|)
    FOR k = max_k DOWNTO 1:
        IF a[|a|-k .. |a|-1] == b[0 .. k-1]:
            RETURN k
    RETURN 0
```

### Helper: `GREEDY_MERGE`

```
FUNCTION GREEDY_MERGE(fragments[], n, ov[][], cost[][]):

    // Collect all directed edges (i→j) sorted by overlap descending
    edges = []
    FOR i = 0 TO n-1:
        FOR j = 0 TO n-1:
            IF i != j:
                APPEND (ov[i][j], i, j) TO edges
    SORT edges BY overlap DESCENDING

    // Track path structure to avoid cycles and branching
    // (each node: at most 1 in-edge, 1 out-edge)
    in_degree  = [0]  * n
    out_degree = [0]  * n
    next_node  = [-1] * n       // successor in path
    prev_node  = [-1] * n       // predecessor in path
    component  = [i for i in 0..n-1]   // union-find for cycle detection

    chains_remaining = n        // number of separate path chains

    FOR EACH (overlap_val, i, j) IN edges:
        IF chains_remaining == 1:
            BREAK               // single chain formed, done

        // Reject edge if it would cause:
        //   (a) branching at i  (i already has out-edge)
        //   (b) branching at j  (j already has in-edge)
        //   (c) a cycle (i and j in same component) UNLESS last merge
        IF out_degree[i] > 0:                               CONTINUE
        IF in_degree[j]  > 0:                               CONTINUE
        IF FIND(component, i) == FIND(component, j)
           AND chains_remaining > 1:                        CONTINUE

        // Accept edge i → j
        next_node[i] = j
        prev_node[j] = i
        out_degree[i] = 1
        in_degree[j]  = 1
        UNION(component, i, j)
        chains_remaining = chains_remaining - 1

    // Reconstruct the single merged superstring
    // Find the start node (no incoming edge)
    start = 0
    FOR i = 0 TO n-1:
        IF prev_node[i] == -1:
            start = i
            BREAK

    result = fragments[start]
    curr   = start
    WHILE next_node[curr] != -1:
        nxt    = next_node[curr]
        result = result + fragments[nxt][ov[curr][nxt] .. end]
        curr   = nxt

    RETURN result
```

---

## Phase 1 Variants
> *Three greedy strategies are implemented — all produce a valid superstring; quality varies*

### Variant A — GREEDY (Simple Pairwise)
> *Scan all pairs each round and pick the maximum-overlap merge. O(n³), straightforward.*

```
FUNCTION SIMPLE_GREEDY_MERGE(fragments[], n):
    WHILE n > 1:
        best_ov = 0
        best_i  = 0
        best_j  = 1

        FOR i = 0 TO n-1:
            FOR j = 0 TO n-1:
                IF i != j:
                    ov_ij = OVERLAP(fragments[i], fragments[j])
                    IF ov_ij > best_ov:
                        best_ov = ov_ij
                        best_i  = i
                        best_j  = j

        // Merge: fragments[best_i] is the prefix, fragments[best_j] the suffix
        fragments[best_i] = fragments[best_i]
                          + fragments[best_j][best_ov .. end]
        REMOVE fragments[best_j]
        n = n - 1

    RETURN fragments[0]
```

---

### Variant B — MGREEDY (Edge-Cover)
> *Builds a max-overlap cycle cover then linearises — matches the `GREEDY_MERGE` pseudocode in Phase 1.*

See `GREEDY_MERGE` pseudocode in Phase 1 above.

---

### Variant C — TGREEDY (Two-Phase)
> *Runs MGREEDY first to collapse fragments into chains, then SIMPLE_GREEDY_MERGE on those chains.*

```
FUNCTION TGREEDY(fragments[], n, ov[][], cost[][]):
    // Phase 1: reduce n fragments into fewer merged chains via MGREEDY
    chains = GREEDY_MERGE(fragments, n, ov, cost)

    // Phase 2: merge the remaining chains with SIMPLE_GREEDY_MERGE
    RETURN SIMPLE_GREEDY_MERGE(chains, |chains|)
```

> **Rationale:** MGREEDY is fast but may leave several disjoint chains.
> The second pass merges them, often producing a shorter result than
> either strategy applied alone.

---

## Phase 2 — Exact Search
> *Held-Karp DP with upper-bound pruning*

```
FUNCTION Exact_Search(fragments[], n, ov[][], cost[][], upper_bound):
    // upper_bound = length of greedy solution
    // Used to prune states that cannot improve on the best known solution

    IF n > HELD_KARP_THRESHOLD:
        // Too large for exact DP; fall back to branch-and-bound
        CALL Branch_and_Bound(fragments, n, ov, cost, upper_bound)
        RETURN

    // Held-Karp DP
    dp[][]     = INFINITY   // dp[S][v]: min overlap-cost to visit subset S ending at v
    parent[][] = -1         // for path reconstruction

    // Base case: single-node subsets
    FOR v = 0 TO n-1:
        S        = (1 << v)
        dp[S][v] = 0

    // Fill DP table (subsets in increasing size order)
    FOR S = 1 TO (1 << n) - 1:
        FOR v = 0 TO n-1:
            IF bit v NOT set in S:      CONTINUE
            IF dp[S][v] == INFINITY:    CONTINUE

            FOR u = 0 TO n-1:
                IF bit u set in S:      CONTINUE    // u already visited

                S2       = S | (1 << u)
                new_cost = dp[S][v] + cost[v][u]

                // Pruning: abandon if partial cost already exceeds best known
                partial_len = BASE_LENGTH(S, fragments) + new_cost
                IF partial_len >= upper_bound:
                    CONTINUE

                IF new_cost < dp[S2][u]:
                    dp[S2][u]     = new_cost
                    parent[S2][u] = v

    // Extract best complete path
    FULL      = (1 << n) - 1
    best_cost = INFINITY
    best_last = -1

    FOR v = 0 TO n-1:
        IF dp[FULL][v] < best_cost:
            best_cost = dp[FULL][v]
            best_last = v

    // Traceback and merge
    path   = TRACEBACK(parent, n, FULL, best_last)
    result = MERGE_PATH(path, fragments, ov)

    IF |result| < upper_bound:
        upper_bound = |result|
        OUTPUT result       // *** ANYTIME OUTPUT #2 ***
                            // Better solution found by exact DP

    RETURN result
```

### Helper: `BASE_LENGTH`

```
FUNCTION BASE_LENGTH(S, fragments[]):
    // Minimum possible superstring length for subset S
    // = length of longest fragment in S (tightest lower bound)
    max_len = 0
    FOR i = 0 TO n-1:
        IF bit i set in S:
            max_len = max(max_len, |fragments[i]|)
    RETURN max_len
```

### Helper: `TRACEBACK`

```
FUNCTION TRACEBACK(parent[][], n, FULL, last):
    path = empty stack
    S    = FULL
    curr = last
    WHILE curr != -1:
        PUSH curr ONTO path
        prev = parent[S][curr]
        S    = S XOR (1 << curr)
        curr = prev
    RETURN path     // in correct forward order
```

### Helper: `MERGE_PATH`

```
FUNCTION MERGE_PATH(path[], fragments[], ov[][]):
    result = fragments[path[0]]
    FOR i = 1 TO |path| - 1:
        u      = path[i-1]
        v      = path[i]
        result = result + fragments[v][ov[u][v] .. end]
    RETURN result
```

---

## Phase 3 — Branch and Bound
> *Fallback for large n; still anytime*
>
> Best-first search over all fragment orderings, pruned by a lower bound on the
> remaining SCS length. Outputs each improvement immediately (anytime).

```
FUNCTION Branch_and_Bound(fragments[], n, ov[][], cost[][], upper_bound):
    // upper_bound = best length found so far (from greedy / Held-Karp)

    best = ""

    // Priority queue ordered by lower bound (min first)
    // Each state: (lower_bound, mask, last, cur_len, path[])
    pq = empty min-heap

    // Seed one search state per possible starting fragment
    FOR start = 0 TO n-1:
        mask   = (1 << start)
        curlen = |fragments[start]|
        lb     = curlen + REMAINING_LB(mask, fragments, n, ov)
        PUSH (lb, mask, start, curlen, [start]) ONTO pq

    WHILE pq NOT empty:
        (lb, mask, last, curlen, path) = POP_MIN(pq)

        // Prune: this branch cannot beat the current best
        IF lb >= upper_bound:   CONTINUE

        // Complete solution — all fragments visited
        IF mask == (1 << n) - 1:
            result = MERGE_PATH(path, fragments, ov)
            IF |result| < upper_bound:
                upper_bound = |result|
                best        = result
                OUTPUT result       // *** ANYTIME OUTPUT ***
            CONTINUE

        // Branch: extend path with each unvisited fragment
        FOR next = 0 TO n-1:
            IF bit next set in mask:    CONTINUE    // already visited

            new_mask   = mask | (1 << next)
            new_curlen = curlen + cost[last][next]
            new_lb     = new_curlen + REMAINING_LB(new_mask, fragments, n, ov)

            IF new_lb < upper_bound:
                PUSH (new_lb, new_mask, next, new_curlen,
                      path + [next]) ONTO pq

    RETURN best
```

### Helper: `REMAINING_LB`

```
FUNCTION REMAINING_LB(mask, fragments[], n, ov[][]):
    // Lower bound on additional characters still needed for unvisited fragments.
    // For each unvisited fragment u, it must contribute at least:
    //   |fragments[u]| - max_ov_into_u
    // where max_ov_into_u is the best overlap any fragment can feed into u.
    // Using the global maximum is optimistic (ignores ordering), so it is a
    // valid (never over-estimating) lower bound.

    lb = 0
    FOR u = 0 TO n-1:
        IF bit u set in mask:   CONTINUE    // already placed

        best_ov = 0
        FOR v = 0 TO n-1:
            IF v != u:
                best_ov = max(best_ov, ov[v][u])

        lb = lb + (|fragments[u]| - best_ov)

    RETURN lb
```
