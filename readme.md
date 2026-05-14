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

```
FUNCTION Branch_and_Bound(fragments[], n, ov[][], cost[][], upper_bound):
    // TODO
```
