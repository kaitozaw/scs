#!/usr/bin/env python3
"""
quality.py — Solution-quality measurement for the SCS solver.

Runs text_reconstruction on a single input and reports, per algorithm, how
good each reconstruction is: correctness, length versus direct concatenation,
and length versus a theoretical lower bound.

This logic previously lived in Section B.1 of text_reconstruction.c. It is pure
post-hoc analysis derived from the input fragments and the solver's [algo=...]
output, so it lives here instead of bloating the C source. The [test=...] lines
printed below reproduce the format the C solver used to emit directly.

Standard library only; reuses helpers from bench.py. No external dependencies.

Usage:
  python3 quality.py <input_file>            # or '-' for stdin
  python3 quality.py --algos greedy,hk <input_file>
  python3 quality.py --timeout 60 <input_file>

Output:
  stdout : reconstructed text (forwarded from the solver).
  stderr : [algo=NAME] ... lines (forwarded from the solver), followed by
           [test=NAME] correct=yes|no direct=N len=N reduction=N
                       improvement=P% lb=N gap=N est_ratio=R elapsed=T.TTTs
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

from bench import PER_RUN_TIMEOUT_S, _overlap, ensure_built, parse_stderr, run_solver

# ============================================================================
# Section 1: Input loading & preprocessing
# ============================================================================


def read_fragments(path: str) -> list[str]:
    """Load fragments from a file (or stdin if '-'), one per line, skipping
    blank lines. Mirrors read_all_fragments_array() in text_reconstruction.c."""
    text = sys.stdin.read() if path == "-" else Path(path).read_text()
    return [line for line in text.splitlines() if line]


def remove_substring_fragments(frags: list[str]) -> list[str]:
    """Drop any fragment that is contained in another fragment.

    Faithful port of remove_substring_fragments_array() in text_reconstruction.c,
    including the swap-with-last removal and the equal-length duplicate rule
    (keep the earlier index). The solver computes its quality metrics over this
    preprocessed set, so quality.py must preprocess identically.
    """
    items = list(frags)
    i = 0
    while i < len(items):
        is_substring = False
        for j in range(len(items)):
            if i == j:
                continue
            if items[i] in items[j]:
                if len(items[i]) == len(items[j]) and j > i:
                    continue
                is_substring = True
                break
        if is_substring:
            items[i] = items[-1]
            items.pop()
        else:
            i += 1
    return items


# ============================================================================
# Section 2: Quality metrics  (ports of the test_* helpers from Section B.1)
# ============================================================================


def direct_concat_len(frags: list[str]) -> int:
    """Baseline length: all fragments joined with no overlap."""
    return sum(len(f) for f in frags)


def longest_fragment_lb(frags: list[str]) -> int:
    """Lower bound 1: the result is at least as long as the longest fragment."""
    return max((len(f) for f in frags), default=0)


def overlap_based_lb(frags: list[str]) -> int:
    """Lower bound 2: total length minus the best incoming overlap per fragment.

    A theoretical estimate, not the true optimum. Clamped to stay >= 1.
    """
    total = direct_concat_len(frags)
    max_possible_saving = 0
    for j, fj in enumerate(frags):
        best_incoming = 0
        for i, fi in enumerate(frags):
            if i == j:
                continue
            ov = _overlap(fi, fj)
            if ov > best_incoming:
                best_incoming = ov
        max_possible_saving += best_incoming
    if max_possible_saving >= total:
        return 1
    return total - max_possible_saving


def final_lower_bound(frags: list[str]) -> int:
    """Stronger of the two lower bounds."""
    return max(longest_fragment_lb(frags), overlap_based_lb(frags))


def is_correct(frags: list[str], candidate: str | None) -> bool:
    """A reconstruction is correct iff every fragment occurs inside it."""
    if candidate is None:
        return False
    return all(f in candidate for f in frags)


# ============================================================================
# Section 3: Report
# ============================================================================


def report(frags: list[str], stderr_text: str) -> None:
    """Print a [test=...] line for every algorithm that produced a result."""
    direct = direct_concat_len(frags)
    lb = final_lower_bound(frags)

    for rec in parse_stderr(stderr_text):
        if rec.result is None:
            continue  # NULL candidate — the C solver skipped these too
        cand_len = rec.length
        correct = is_correct(frags, rec.result)
        reduction = direct - cand_len
        gap = cand_len - lb
        improvement = (reduction / direct * 100.0) if direct > 0 else 0.0
        est_ratio = (cand_len / lb) if lb > 0 else 0.0
        print(
            f"[test={rec.name}] correct={'yes' if correct else 'no'} "
            f"direct={direct} len={cand_len} reduction={reduction} "
            f"improvement={improvement:.2f}% lb={lb} gap={gap} "
            f"est_ratio={est_ratio:.4f} elapsed={rec.elapsed_s:.3f}s",
            file=sys.stderr,
        )


# ============================================================================
# Section 4: main()
# ============================================================================


def main() -> None:
    p = argparse.ArgumentParser(
        description="Solution-quality measurement for the SCS solver")
    p.add_argument("input",
                   help="input file (one fragment per line), or '-' for stdin")
    p.add_argument("--algos", default=None,
                   help="subset of {greedy,mgreedy,tgreedy,hk,bb} (default: all)")
    p.add_argument("--timeout", type=float, default=PER_RUN_TIMEOUT_S,
                   help=f"solver timeout in seconds (default: {PER_RUN_TIMEOUT_S})")
    args = p.parse_args()

    ensure_built(False)

    fragments = read_fragments(args.input)
    run = run_solver(fragments, algos=args.algos, timeout=args.timeout)

    # Forward the solver's own output, then append quality metrics.
    sys.stdout.write(run.stdout)
    sys.stdout.flush()
    sys.stderr.write(run.stderr)
    sys.stderr.flush()
    if run.timed_out:
        print(f"[quality] solver timed out after {args.timeout}s; "
              f"metrics below cover partial output only", file=sys.stderr)

    report(remove_substring_fragments(fragments), run.stderr)


if __name__ == "__main__":
    main()
