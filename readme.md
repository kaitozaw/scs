# SCS Hybrid Solver

Anytime solver for Shortest Common Supersequence. Each algorithm publishes its best result immediately, in non-increasing length order.

```sh
gcc -o text_reconstruction text_reconstruction.c
./text_reconstruction <input_file>      # or '-' for stdin
```

- **stdin/input** — one fragment per line.
- **stdout** — reconstructed text, one solution per line.
- **stderr** — `[algo=NAME] len=N elapsed=T.TTTs result=...`.

## Pipeline ([text_reconstruction.c](text_reconstruction.c))

- **A** Common types & utilities (fragments, greedy merge, overlap, cycle-cover & B&B primitives)
- **B** Best-solution registry (`try_record_solution` enforces the stdout contract)
- **C** Step 1 — Preprocess: read fragments, drop substrings of others
- **D** Step 2 — sub-optimal / quick: GREEDY, MGREEDY, TGREEDY
- **E** Step 3 — optimal / slow: HELD_KARP (n ≤ 20) or BRANCH_AND_BOUND (n ≤ 63)
- **F** `main()`
