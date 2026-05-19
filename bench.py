#!/usr/bin/env python3
"""
bench.py — Empirical evaluation harness for the SCS solver.

Runs three experiments:
  Exp 1: Step 2 (GREEDY/MGREEDY/TGREEDY) optimality + time across charset x n.
  Exp 2: Step 3 (HK/BB) timing without UB pruning.
  Exp 3: Step 3 (HK/BB) timing with UB pruning (from Step 2).

Outputs formatted ASCII tables to stdout and per-experiment CSVs to results/.
Standard library only; no external dependencies.

Defaults:
  --exp        all                 (run Exp 1, 2, 3 in order)
  --ns         5,10,15,...,60      (12 values, step 5)
  --charsets   binary,dna,alphabet,ascii  (Exp 1)
               binary,ascii                (Exp 2 and 3)
  --timeout    30  seconds per solver invocation
  --rebuild    off (skips gcc unless binary missing)

Usage:
  # Full sweep: all 3 experiments with default grid
  python3 bench.py

  # Run a single experiment
  python3 bench.py --exp 1
  python3 bench.py --exp 2
  python3 bench.py --exp 3

  # Subset of charsets (override defaults)
  python3 bench.py --charsets binary
  python3 bench.py --charsets binary,dna

  # Subset of n values (any comma list works)
  python3 bench.py --ns 5,10,15
  python3 bench.py --ns 20,40,60

  # Tune per-run timeout (seconds). Increase for larger n on slow charsets.
  python3 bench.py --timeout 60
  python3 bench.py --exp 2 --charsets binary --timeout 300

  # Force gcc rebuild before benchmarking (after C source edits)
  python3 bench.py --rebuild

  # Combined: one experiment, small grid, custom timeout
  python3 bench.py --exp 3 --charsets binary,ascii --ns 5,10,20,40 --timeout 60

Output:
  - Formatted ASCII tables printed to stdout (one per experiment).
  - CSV files written to results/ (exp1.csv, exp2.csv, exp3.csv) for plotting.
  - DNF rows appear when a run exceeds --timeout; subsequent larger-n runs
    for the same charset are auto-skipped within an experiment.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import random
import re
import string
import subprocess
import sys
import time
from dataclasses import dataclass, field
from pathlib import Path

# ============================================================================
# Section 1: Config
# ============================================================================

CHARSETS: dict[str, str] = {
    "binary":   "01",
    "dna":      "ACGT",
    "alphabet": string.ascii_lowercase,
    # printable non-space ASCII (33..126); excludes space so empty-line parsing
    # in text_reconstruction's read_all_fragments_array stays safe.
    "ascii":    "".join(chr(c) for c in range(33, 127)),
}

DEFAULT_NS = list(range(5, 65, 5))   # 5, 10, ..., 60
L = 10
PER_RUN_TIMEOUT_S = 30
HK_THRESHOLD = 20
BB_LIMIT = 63

SCRIPT_DIR = Path(__file__).resolve().parent
BIN_PATH = SCRIPT_DIR / "text_reconstruction"
SRC_PATH = SCRIPT_DIR / "text_reconstruction.c"
RESULTS_DIR = SCRIPT_DIR / "results"


# ============================================================================
# Section 2: Generator
# ============================================================================

def stable_seed(charset: str, n: int) -> int:
    h = hashlib.sha256(f"{charset}|{n}".encode()).digest()
    return int.from_bytes(h[:8], "big")


def generate_fragments(charset: str, n: int, frag_len: int = L) -> list[str]:
    """Generate n fragments of fixed length by sampling random windows from a base string."""
    rng = random.Random(stable_seed(charset, n))
    chars = CHARSETS[charset]
    L_base = max(n * frag_len // 2, frag_len)
    base = "".join(rng.choice(chars) for _ in range(L_base))
    fragments = []
    for _ in range(n):
        start = rng.randrange(0, L_base - frag_len + 1)
        fragments.append(base[start:start + frag_len])
    return fragments


def _overlap(a: str, b: str) -> int:
    """Longest suffix-of-a matching prefix-of-b (excludes full-string match)."""
    max_k = min(len(a), len(b)) - 1
    for k in range(max_k, 0, -1):
        if a[-k:] == b[:k]:
            return k
    return 0


def compute_input_stats(fragments: list[str]) -> tuple[int, int]:
    """
    Return (overlap_edges, cycle_cover_size) — same convention as MGREEDY in
    text_reconstruction.c: greedy max-weight selection with in/out-degree <= 1,
    then count cycles in the resulting permutation graph.
    """
    n = len(fragments)
    # Build overlap matrix; count positive overlap pairs.
    edges = []
    overlap_edges = 0
    for i in range(n):
        for j in range(n):
            if i == j:
                continue
            w = _overlap(fragments[i], fragments[j])
            if w > 0:
                overlap_edges += 1
                edges.append((w, i, j))
    # Greedy max-weight selection with in/out-degree <= 1.
    edges.sort(reverse=True)
    succ = [-1] * n
    out_used = [False] * n
    in_used = [False] * n
    for _w, i, j in edges:
        if out_used[i] or in_used[j]:
            continue
        succ[i] = j
        out_used[i] = True
        in_used[j] = True
    # Count cycles. First walk all path-starting nodes (no predecessor).
    pred = [-1] * n
    for i in range(n):
        if succ[i] != -1:
            pred[succ[i]] = i
    visited = [False] * n
    for start in range(n):
        if pred[start] == -1:
            v = start
            while v != -1 and not visited[v]:
                visited[v] = True
                v = succ[v]
    # Remaining unvisited nodes belong to cycles.
    cycles = 0
    for start in range(n):
        if not visited[start]:
            v = start
            while not visited[v]:
                visited[v] = True
                v = succ[v]
            cycles += 1
    return overlap_edges, cycles


# ============================================================================
# Section 3: Runner
# ============================================================================

@dataclass
class SolverRun:
    stdout: str
    stderr: str
    returncode: int
    wall_s: float
    timed_out: bool


def run_solver(fragments: list[str], algos: str | None = None,
               timeout: float = PER_RUN_TIMEOUT_S) -> SolverRun:
    cmd = [str(BIN_PATH)]
    if algos:
        cmd.extend(["--algos", algos])
    cmd.append("-")
    input_text = "\n".join(fragments) + "\n"
    t0 = time.perf_counter()
    try:
        proc = subprocess.run(
            cmd,
            input=input_text,
            capture_output=True,
            text=True,
            timeout=timeout,
        )
        wall = time.perf_counter() - t0
        return SolverRun(proc.stdout, proc.stderr, proc.returncode, wall, False)
    except subprocess.TimeoutExpired as e:
        wall = time.perf_counter() - t0
        out = (e.stdout or b"").decode("utf-8", errors="replace") if isinstance(e.stdout, bytes) else (e.stdout or "")
        err = (e.stderr or b"").decode("utf-8", errors="replace") if isinstance(e.stderr, bytes) else (e.stderr or "")
        return SolverRun(out, err, -1, wall, True)


# ============================================================================
# Section 4: Parser
# ============================================================================

ALGO_LINE_RE = re.compile(
    r"^\[algo=(?P<name>[A-Z_]+)\]\s+len=(?P<len>\d+)\s+elapsed=(?P<t>[\d.]+)s"
    r"(?:\s+result=(?P<r>.+))?\s*$"
)
# Matches both:
#   [algo=HK] len=N elapsed=T.TTTs result=ABC
#   [algo=HK] len=N elapsed=T.TTTs (no improvement over g_best)
NO_IMPROVEMENT_RE = re.compile(
    r"^\[algo=(?P<name>[A-Z_]+)\]\s+len=(?P<len>\d+)\s+elapsed=(?P<t>[\d.]+)s\s+\(no improvement"
)


@dataclass
class AlgoRecord:
    name: str
    length: int
    elapsed_s: float
    result: str | None
    no_improvement: bool


def parse_stderr(text: str) -> list[AlgoRecord]:
    records = []
    for line in text.splitlines():
        stripped = line.strip()
        m = ALGO_LINE_RE.match(stripped)
        if m:
            records.append(AlgoRecord(
                name=m.group("name"),
                length=int(m.group("len")),
                elapsed_s=float(m.group("t")),
                result=m.group("r"),
                no_improvement=m.group("r") is None,
            ))
            continue
        m2 = NO_IMPROVEMENT_RE.match(stripped)
        if m2:
            records.append(AlgoRecord(
                name=m2.group("name"),
                length=int(m2.group("len")),
                elapsed_s=float(m2.group("t")),
                result=None,
                no_improvement=True,
            ))
    return records


# ============================================================================
# Section 5: OPT extraction
# ============================================================================

def find_record(records: list[AlgoRecord], name: str) -> AlgoRecord | None:
    for r in records:
        if r.name == name:
            return r
    return None


def opt_length(records: list[AlgoRecord]) -> tuple[int | None, str]:
    """Return (opt_len, source). Source: 'HK', 'BB', 'step2-tie', or 'none'."""
    hk = find_record(records, "HELD_KARP")
    if hk and not hk.no_improvement:
        return hk.length, "HK"
    bb = find_record(records, "BRANCH_AND_BOUND")
    if bb and not bb.no_improvement:
        return bb.length, "BB"
    # Fallback: tie to the best Step 2 result if any.
    step2 = [find_record(records, n) for n in ("GREEDY", "MGREEDY", "TGREEDY")]
    step2_lens = [r.length for r in step2 if r is not None]
    if step2_lens:
        return min(step2_lens), "step2-tie"
    return None, "none"


# ============================================================================
# Section 6: Table formatter (stdlib only)
# ============================================================================

def render_table(headers: list[str], rows: list[list], aligns: list[str] | None = None) -> str:
    """Render an ASCII table. aligns: list of 'l' or 'r' per column (default 'r')."""
    cols = len(headers)
    if aligns is None:
        aligns = ["r"] * cols
    str_rows = [[str(c) for c in r] for r in rows]
    widths = [
        max(len(headers[i]), *(len(r[i]) for r in str_rows)) if str_rows else len(headers[i])
        for i in range(cols)
    ]

    def fmt_cell(s: str, w: int, a: str) -> str:
        return s.ljust(w) if a == "l" else s.rjust(w)

    sep = "+" + "+".join("-" * (w + 2) for w in widths) + "+"
    out = [sep]
    out.append("| " + " | ".join(fmt_cell(headers[i], widths[i], aligns[i]) for i in range(cols)) + " |")
    out.append(sep)
    for r in str_rows:
        out.append("| " + " | ".join(fmt_cell(r[i], widths[i], aligns[i]) for i in range(cols)) + " |")
    out.append(sep)
    return "\n".join(out)


# ============================================================================
# Section 7: CSV writer
# ============================================================================

class CsvWriter:
    def __init__(self, path: Path, fieldnames: list[str]):
        path.parent.mkdir(parents=True, exist_ok=True)
        self.path = path
        self.fh = path.open("w", newline="")
        self.writer = csv.DictWriter(self.fh, fieldnames=fieldnames)
        self.writer.writeheader()
        self.fh.flush()

    def write(self, row: dict) -> None:
        self.writer.writerow(row)
        self.fh.flush()

    def close(self) -> None:
        self.fh.close()


# ============================================================================
# Section 8: Experiments
# ============================================================================

def fmt_time(s: float | None, dnf: bool) -> str:
    if dnf:
        return "DNF"
    if s is None:
        return "-"
    return f"{s:.3f}"


def fmt_int(v: int | None, dnf: bool = False) -> str:
    if dnf:
        return "---"
    if v is None:
        return "-"
    return str(v)


def fmt_ratio(v: float | None) -> str:
    if v is None:
        return "-"
    return f"{v:.3f}"


def experiment_1(charsets: list[str], ns: list[int], timeout: float) -> None:
    print("\n=== Experiment 1: Step 2 optimality & time ===")
    csv_path = RESULTS_DIR / "exp1.csv"
    writer = CsvWriter(csv_path, [
        "charset", "n", "cycle_cover", "overlap", "algo", "len", "elapsed_s",
        "opt_len", "opt_source", "optimality_ratio", "wall_s", "dnf",
    ])

    for charset in charsets:
        print(f"\n[charset = {charset}]")
        headers = ["n", "cycle_cover", "overlap", "OPT",
                   "g_len", "g_ratio", "g_t(s)",
                   "m_len", "m_ratio", "m_t(s)",
                   "t_len", "t_ratio", "t_t(s)"]
        table_rows = []

        for n in ns:
            frags = generate_fragments(charset, n)
            overlap_edges, cycle_cover = compute_input_stats(frags)
            run = run_solver(frags, algos=None, timeout=timeout)
            records = parse_stderr(run.stderr)
            opt_len, opt_src = opt_length(records)

            row_cells = [n, cycle_cover, overlap_edges, fmt_int(opt_len, run.timed_out)]
            for algo_name, prefix in (("GREEDY", "g"), ("MGREEDY", "m"), ("TGREEDY", "t")):
                r = find_record(records, algo_name)
                if r is None or r.no_improvement:
                    ratio = None
                    row_cells.extend([fmt_int(None), fmt_ratio(None), fmt_time(None, run.timed_out)])
                    writer.write({
                        "charset": charset, "n": n,
                        "cycle_cover": cycle_cover, "overlap": overlap_edges,
                        "algo": algo_name,
                        "len": "", "elapsed_s": "",
                        "opt_len": opt_len if opt_len is not None else "",
                        "opt_source": opt_src,
                        "optimality_ratio": "",
                        "wall_s": f"{run.wall_s:.3f}",
                        "dnf": str(run.timed_out).lower(),
                    })
                else:
                    ratio = r.length / opt_len if opt_len else None
                    row_cells.extend([
                        fmt_int(r.length), fmt_ratio(ratio),
                        fmt_time(r.elapsed_s, False),
                    ])
                    writer.write({
                        "charset": charset, "n": n,
                        "cycle_cover": cycle_cover, "overlap": overlap_edges,
                        "algo": algo_name,
                        "len": r.length,
                        "elapsed_s": f"{r.elapsed_s:.6f}",
                        "opt_len": opt_len if opt_len is not None else "",
                        "opt_source": opt_src,
                        "optimality_ratio": f"{ratio:.6f}" if ratio is not None else "",
                        "wall_s": f"{run.wall_s:.3f}",
                        "dnf": str(run.timed_out).lower(),
                    })
            table_rows.append(row_cells)

        print(render_table(headers, table_rows))

    writer.close()
    print(f"\nCSV: {csv_path}")


def _step2_lens(records: list[AlgoRecord]) -> dict:
    """Extract step-2 lengths for CSV (warm-start runs only)."""
    out = {}
    for name_s, key in (("GREEDY", "greedy_len"),
                        ("MGREEDY", "mgreedy_len"),
                        ("TGREEDY", "tgreedy_len")):
        rs = find_record(records, name_s)
        out[key] = rs.length if rs and not rs.no_improvement else ""
    return out


def _run_step3_algo(frags: list[str], algo_token: str, expected_algo: str,
                    warm_start: bool, timeout: float):
    """Single Step-3 measurement. Returns (record, wall_s, timed_out, step2_lens)."""
    algos = (f"greedy,mgreedy,tgreedy,{algo_token}" if warm_start else algo_token)
    run = run_solver(frags, algos=algos, timeout=timeout)
    records = parse_stderr(run.stderr)
    rec = find_record(records, expected_algo)
    s2 = _step2_lens(records) if warm_start else {}
    return rec, run.wall_s, run.timed_out, s2


def _experiment_step3(name: str, charsets: list[str], ns: list[int],
                      timeout: float, warm_start: bool) -> None:
    print(f"\n=== {name} ===")
    fname = "exp3.csv" if warm_start else "exp2.csv"
    csv_path = RESULTS_DIR / fname

    fields = ["charset", "n", "algo", "len", "elapsed_s", "wall_s", "dnf"]
    if warm_start:
        fields += ["greedy_len", "mgreedy_len", "tgreedy_len"]
    writer = CsvWriter(csv_path, fields)

    headers = ["charset", "n", "HK_len", "HK_t(s)", "BB_len", "BB_t(s)"]
    aligns = ["l", "r", "r", "r", "r", "r"]
    table_rows = []

    for charset in charsets:
        skip_bb = False
        for n in ns:
            hk_applicable = n <= HK_THRESHOLD
            bb_applicable = n <= BB_LIMIT and not skip_bb
            if not hk_applicable and not bb_applicable:
                continue  # nothing to measure for this (charset, n)

            frags = generate_fragments(charset, n)
            row = [charset, n]

            # --- HK measurement (independent process, never warm-started by BB) ---
            if hk_applicable:
                rec, wall, dnf, s2 = _run_step3_algo(
                    frags, "hk", "HELD_KARP", warm_start, timeout)
                csv_row = {
                    "charset": charset, "n": n, "algo": "HELD_KARP",
                    "len": rec.length if rec else "",
                    "elapsed_s": f"{rec.elapsed_s:.6f}" if rec else "",
                    "wall_s": f"{wall:.3f}",
                    "dnf": str(dnf).lower(),
                }
                if warm_start:
                    csv_row.update(s2)
                writer.write(csv_row)
                if dnf:
                    row.extend(["---", "DNF"])
                elif rec:
                    row.extend([rec.length, f"{rec.elapsed_s:.3f}"])
                else:
                    row.extend(["-", "-"])
            else:
                row.extend(["-", "-"])

            # --- BB measurement (independent process, no HK influence) ---
            if bb_applicable:
                rec, wall, dnf, s2 = _run_step3_algo(
                    frags, "bb", "BRANCH_AND_BOUND", warm_start, timeout)
                csv_row = {
                    "charset": charset, "n": n, "algo": "BRANCH_AND_BOUND",
                    "len": rec.length if rec else "",
                    "elapsed_s": f"{rec.elapsed_s:.6f}" if rec else "",
                    "wall_s": f"{wall:.3f}",
                    "dnf": str(dnf).lower(),
                }
                if warm_start:
                    csv_row.update(s2)
                writer.write(csv_row)
                if dnf:
                    row.extend(["---", "DNF"])
                    skip_bb = True  # don't waste time on larger n with same charset
                elif rec:
                    row.extend([rec.length, f"{rec.elapsed_s:.3f}"])
                else:
                    row.extend(["-", "-"])
            else:
                row.extend(["-", "-"])

            table_rows.append(row)

    print(render_table(headers, table_rows, aligns=aligns))
    writer.close()
    print(f"\nCSV: {csv_path}")


def experiment_2(charsets: list[str], ns: list[int], timeout: float) -> None:
    _experiment_step3("Experiment 2: Step 3 without UB pruning",
                      charsets, ns, timeout, warm_start=False)


def experiment_3(charsets: list[str], ns: list[int], timeout: float) -> None:
    _experiment_step3("Experiment 3: Step 3 with UB pruning",
                      charsets, ns, timeout, warm_start=True)


# ============================================================================
# Section 9: main()
# ============================================================================

def ensure_built(rebuild: bool) -> None:
    if rebuild or not BIN_PATH.exists():
        print(f"Compiling: gcc -O2 -Wall -Wextra -o {BIN_PATH.name} {SRC_PATH.name}",
              file=sys.stderr)
        subprocess.run(
            ["gcc", "-O2", "-Wall", "-Wextra", "-o", str(BIN_PATH), str(SRC_PATH)],
            check=True, cwd=SCRIPT_DIR,
        )


def parse_list(s: str, valid: set[str] | None = None) -> list[str]:
    items = [x.strip() for x in s.split(",") if x.strip()]
    if valid:
        bad = [x for x in items if x not in valid]
        if bad:
            raise argparse.ArgumentTypeError(f"unknown values: {bad}; valid: {sorted(valid)}")
    return items


def parse_int_list(s: str) -> list[int]:
    return [int(x.strip()) for x in s.split(",") if x.strip()]


def main() -> None:
    p = argparse.ArgumentParser(description="SCS empirical benchmark harness")
    p.add_argument("--exp", default="all", choices=["1", "2", "3", "all"],
                   help="Which experiment to run (default: all)")
    p.add_argument("--charsets", type=str, default=None,
                   help="Override charsets (comma-separated). Default per experiment.")
    p.add_argument("--ns", type=str, default=None,
                   help="Override n values (comma-separated). Default: 5..60 step 5.")
    p.add_argument("--timeout", type=float, default=PER_RUN_TIMEOUT_S,
                   help=f"Per-run timeout seconds (default: {PER_RUN_TIMEOUT_S})")
    p.add_argument("--rebuild", action="store_true",
                   help="Re-run gcc before benchmarking")
    args = p.parse_args()

    ensure_built(args.rebuild)

    ns = parse_int_list(args.ns) if args.ns else DEFAULT_NS

    exp1_charsets = parse_list(args.charsets, set(CHARSETS)) if args.charsets else list(CHARSETS.keys())
    exp23_charsets = parse_list(args.charsets, set(CHARSETS)) if args.charsets else ["binary", "ascii"]

    if args.exp in ("1", "all"):
        experiment_1(exp1_charsets, ns, args.timeout)
    if args.exp in ("2", "all"):
        experiment_2(exp23_charsets, ns, args.timeout)
    if args.exp in ("3", "all"):
        experiment_3(exp23_charsets, ns, args.timeout)


if __name__ == "__main__":
    main()
