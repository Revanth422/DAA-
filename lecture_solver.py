#!/usr/bin/env python3
"""
Lecture Supervision -- Efficient Open Domination on a Circular-Arc Graph
========================================================================

Problem
-------
n lectures repeat every day (period L = 1440 minutes). A lecture is a
circular arc [start, end) and may wrap past midnight. Choose a subset S such
that EVERY lecture overlaps EXACTLY ONE lecture of S (a lecture never counts
itself). Equivalently: S is an *efficient open dominating set* (total perfect
code) of the circular-arc "overlap" graph.

Overlap convention
------------------
Lectures are half-open arcs [start, end): a lecture ending at 13:00 and one
starting at 13:00 do NOT overlap. start == end means a full 24-hour lecture.

Usage
-----
python3 lecture_solver.py                 # type the lectures in interactively
python3 lecture_solver.py input.txt       # solve a file ("-" = read stdin)
python3 lecture_solver.py --demo          # the 8 built-in test cases
python3 lecture_solver.py --stress 500    # compare with brute force
python3 lecture_solver.py --bench         # running-time measurements
python3 lecture_solver.py --latex-report  # LaTeX fragment (used by report.tex)
python3 lecture_solver.py --latex-bench   # LaTeX timing table (used by report.tex)

Input format: first line n, then n lines "HH:MM HH:MM" (start end); "#" starts
a comment. Raw integers (minutes after midnight) are accepted as well.

Indexing: lectures are numbered 1..n in all output (0..n-1 inside the code).
"""
from __future__ import annotations

import argparse
import random
import sys
import time
from bisect import bisect_left

DAY = 1440  # minutes per day
INF = float("inf")


# --------------------------------------------------------------------------
# 1. Parsing and basic geometry
# --------------------------------------------------------------------------
def parse_time(token: str, period: int = DAY) -> int:
    """'HH:MM' (00:00..24:00) or a plain integer -> integer time in [0, period)."""
    token = token.strip()
    if ":" in token:
        parts = token.split(":")
        try:
            if len(parts) != 2:
                raise ValueError
            hh, mm = int(parts[0]), int(parts[1])
        except ValueError:
            raise ValueError(f"cannot read time {token!r}; use HH:MM") from None
        if not (0 <= mm <= 59 and 0 <= hh <= 24 and (hh < 24 or mm == 0)):
            raise ValueError(f"time out of range: {token!r}")
        return (hh * 60 + mm) % period
    try:
        return int(token) % period
    except ValueError:
        raise ValueError(f"cannot read time {token!r}; use HH:MM") from None


def fmt_time(t: int) -> str:
    return f"{t // 60:02d}:{t % 60:02d}"


def to_arcs(pairs, period=DAY):
    """[(start, end)] -> [(start, duration)], duration in 1..period."""
    arcs = []
    for s, e in pairs:
        s %= period
        d = (e - s) % period
        arcs.append((s, d if d else period))
    return arcs


def arc_mask(s: int, d: int, period: int) -> int:
    """Bit-mask of the unit cells covered by the arc (cyclic)."""
    full = (1 << period) - 1
    m = ((1 << d) - 1) << s
    return (m | (m >> period)) & full


def build_graph(arcs, period):
    """Return (cell masks, adjacency bitsets over lecture indices).

    This is the graph model used by the polynomial algorithm.
    """
    n = len(arcs)
    masks = [arc_mask(s, d, period) for s, d in arcs]
    adj = [0] * n
    for i in range(n):
        for j in range(i + 1, n):
            if masks[i] & masks[j]:
                adj[i] |= 1 << j
                adj[j] |= 1 << i
    return masks, adj


# --------------------------------------------------------------------------
# 2. Independent checkers (verification only; they do NOT use build_graph)
# --------------------------------------------------------------------------
def arcs_overlap(a, b, period) -> bool:
    """Do arcs a=(s1,d1) and b=(s2,d2) share a unit cell?

    Two cyclic arcs intersect iff one of them contains the start cell of
    the other. Pure modular arithmetic -- no bit-masks.
    """
    (s1, d1), (s2, d2) = a, b
    return (s2 - s1) % period < d1 or (s1 - s2) % period < d2


def adjacency_arithmetic(arcs, period):
    """Adjacency bitsets computed with arcs_overlap only."""
    n = len(arcs)
    adj = [0] * n
    for i in range(n):
        for j in range(i + 1, n):
            if arcs_overlap(arcs[i], arcs[j], period):
                adj[i] |= 1 << j
                adj[j] |= 1 << i
    return adj


def verify(pairs, chosen, period=DAY) -> bool:
    """True iff every lecture overlaps exactly one chosen lecture.

    Uses its own arithmetic overlap test, so it is independent of the
    bit-mask graph that the solver works on.
    """
    arcs = to_arcs(pairs, period)
    n = len(arcs)
    if not chosen or any(not (0 <= c < n) for c in chosen):
        return False
    adj = adjacency_arithmetic(arcs, period)
    cmask = 0
    for c in chosen:
        cmask |= 1 << c
    return all(bin(adj[v] & cmask).count("1") == 1 for v in range(n))


def brute_force(pairs, period=DAY):
    """Exponential reference solver (only for small n); independent of solve()."""
    n = len(pairs)
    adj = adjacency_arithmetic(to_arcs(pairs, period), period)
    for subset in range(1, 1 << n):
        if all(bin(adj[v] & subset).count("1") == 1 for v in range(n)):
            return [i for i in range(n) if subset >> i & 1]
    return None


# --------------------------------------------------------------------------
# 3. The polynomial algorithm
# --------------------------------------------------------------------------
def _run_start_len(umask: int, period: int):
    """Start cell and length of the (single) run of 1-cells in umask."""
    full = (1 << period) - 1
    if umask == full:
        return None
    rot = ((umask << 1) | (umask >> (period - 1))) & full  # cell c-1 -> bit c
    start_bit = umask & ~rot
    return start_bit.bit_length() - 1, bin(umask).count("1")


def candidate_blocks(arcs, masks, adj, period):
    """
    A *block* is an overlapping pair (a, b) with N(a) & N(b) = {} .
    Every solution is a disjoint union of blocks, so only these pairs matter.
    Returns list of (a, b, union_mask, run) ; run = None if union is the circle.
    """
    n = len(arcs)
    blocks = []
    for a in range(n):
        for b in range(a + 1, n):
            if adj[a] >> b & 1 and (adj[a] & adj[b]) == 0:
                u = masks[a] | masks[b]
                blocks.append((a, b, u, _run_start_len(u, period)))
    return blocks


def solve(pairs, period=DAY):
    """
    Return a sorted list of 0-based lecture indices forming a valid subset,
    or None if none exists. Time O(n^3) (plus O(n^2 L / w) for the masks).
    """
    n = len(pairs)
    if n < 2:
        return None
    L = period
    arcs = to_arcs(pairs, L)
    masks, adj = build_graph(arcs, L)
    blocks = candidate_blocks(arcs, masks, adj, L)

    # Special case: one block already wraps the whole circle.
    for a, b, _, run in blocks:
        if run is None:
            return [a, b]

    # Lecture 0 is touched by exactly one block in any solution -> anchor it.
    anchors = [k for k, blk in enumerate(blocks) if blk[2] & masks[0]]

    for x0 in anchors:
        chain = _try_anchor(x0, blocks, arcs, L)
        if chain is not None:
            return sorted(v for k in chain for v in blocks[k][:2])
    return None


def _try_anchor(x0, blocks, arcs, L):
    """Linear DP over blocks, with block x0 placed at coordinate 0."""
    s0, len0 = blocks[x0][3]

    # lecture arcs in coordinates where x0 starts at 0
    sig = [(s - s0) % L for s, _ in arcs]
    end = [sg + d for sg, (_, d) in zip(sig, arcs)]
    order = sorted(range(len(arcs)), key=lambda i: sig[i])
    sorted_sig = [sig[i] for i in order]
    pref_max = [-INF]                      # pref_max[k] = max end of first k
    for i in order:
        pref_max.append(max(pref_max[-1], end[i]))
    suf_min = [INF] * (len(arcs) + 1)      # suf_min[k] = min end from k on
    for k in range(len(arcs) - 1, -1, -1):
        suf_min[k] = min(suf_min[k + 1], end[order[k]])
    wrap_max = max(e - L for e in end)     # copies shifted by -L

    def Q(p):  # furthest right end of a lecture that starts before p
        return max(pref_max[bisect_left(sorted_sig, p)], wrap_max)

    def T(p):  # earliest end of a lecture that starts at or after p
        return suf_min[bisect_left(sorted_sig, p)]

    # nodes: (start, end, block-id); node 0 is the anchor itself
    nodes = [(0, len0, x0)]
    others = []
    for k, blk in enumerate(blocks):
        if k == x0:
            continue
        s_k, l_k = blk[3]
        ys = (s_k - s0) % L
        if ys >= len0 and ys + l_k <= L:   # disjoint from the anchor
            others.append((ys, ys + l_k, k))
    nodes += sorted(others)

    m = len(nodes)
    reach = [False] * m
    parent = [-1] * m
    reach[0] = True
    Qe = [Q(nd[1]) for nd in nodes]
    Te = [T(nd[1]) for nd in nodes]

    def chain_to(j):
        out = []
        while j != -1:
            out.append(nodes[j][2])
            j = parent[j]
        return out[::-1]

    # k = 1 : the anchor alone covers everything
    if L < Te[0]:
        return [x0]

    for j in range(1, m):
        ys = nodes[j][0]
        for i in range(j):
            if reach[i] and nodes[i][1] <= ys and Qe[i] <= ys < Te[i]:
                reach[j], parent[j] = True, i
                break
        # closing edge j -> anchor shifted by L (needs k >= 2)
        if reach[j] and Qe[j] <= L < Te[j]:
            return chain_to(j)
    return None


# --------------------------------------------------------------------------
# 4. Input helpers
# --------------------------------------------------------------------------
def parse_tokens(tok, period=DAY):
    """Token list [n, s1, e1, s2, e2, ...] -> list of (start, end)."""
    if not tok:
        raise ValueError("empty input")
    try:
        n = int(tok[0])
    except ValueError:
        raise ValueError(f"first value must be the number of lectures, got {tok[0]!r}") from None
    if n < 0 or len(tok) != 1 + 2 * n:
        raise ValueError(f"expected n and then {n} start/end pairs "
                         f"({1 + 2 * n} values), got {len(tok)} values")
    return [(parse_time(tok[1 + 2 * i], period), parse_time(tok[2 + 2 * i], period))
            for i in range(n)]


def read_file(path, period=DAY):
    """Read a file, or standard input when path is '-'. '#' starts a comment."""
    fh = sys.stdin if path == "-" else open(path)
    with fh:
        text = "\n".join(line.split("#")[0] for line in fh)
    return parse_tokens(text.split(), period)


def interactive(period=DAY):
    """Ask the user to type the lectures one by one."""
    n = int(input("Number of lectures n: "))
    pairs = []
    for i in range(1, n + 1):
        while True:
            try:
                s, e = input(f"Lecture {i} - start end (HH:MM HH:MM): ").split()
                pairs.append((parse_time(s, period), parse_time(e, period)))
                break
            except ValueError:
                print("  please enter two times such as 09:30 11:00")
    return pairs


def hm(s, e):
    return (parse_time(s), parse_time(e))


# --------------------------------------------------------------------------
# 5. Built-in test cases (name, lectures, expected: solvable?)
# --------------------------------------------------------------------------
TESTS = [
    ("T1 two overlapping lectures",
     [hm("09:00", "10:00"), hm("09:30", "10:30")], True),
    ("T2 night lecture wraps midnight",
     [hm("23:00", "01:00"), hm("00:30", "02:00")], True),
    ("T3 two chains around the day",
     [hm("01:00", "03:00"), hm("02:00", "04:00"),
      hm("13:00", "15:00"), hm("14:00", "16:00")], True),
    ("T4 four blocks, private lectures",
     [hm("00:00", "02:00"), hm("01:00", "03:00"), hm("00:30", "01:00"),
      hm("08:00", "10:00"), hm("09:00", "11:00"), hm("08:30", "09:00"),
      hm("16:00", "18:00"), hm("17:00", "19:00")], True),
    ("T5 lecture touched by two blocks",
     [hm("01:00", "03:00"), hm("02:00", "04:00"),
      hm("05:00", "07:00"), hm("06:00", "08:00"),
      hm("03:30", "05:30")], False),
    ("T6 full-day lecture + private partner",
     [hm("00:00", "00:00"), hm("10:00", "11:00")], True),
    ("T7 isolated lecture",
     [hm("09:00", "10:00"), hm("09:30", "10:30"), hm("15:00", "16:00")], False),
    ("T8 touching endpoints do not overlap",
     [hm("09:00", "10:00"), hm("10:00", "11:00")], False),
]


# --------------------------------------------------------------------------
# 6. Random, planted and mutated instances (for testing)
# --------------------------------------------------------------------------
def random_instance(rng, n, period):
    return [(rng.randrange(period), rng.randrange(period)) for _ in range(n)]


def planted_instance(rng, period, k, extras):
    """Instance that is solvable BY CONSTRUCTION.

    The circle is cut into k windows (each >= 3 cells). Each window holds one
    block {a, b} (crossing or nested) and no lecture leaves its region, so
    every lecture overlaps exactly one region. Up to `extras` private lectures
    are added inside a \\ b or b \\ a. Returns a shuffled, rotated list.
    """
    if 3 * k > period:
        raise ValueError("too many blocks for this period")
    spare = period - 3 * k
    cuts = sorted(rng.randrange(spare + 1) for _ in range(k - 1))
    pieces = [b - a for a, b in zip([0] + cuts, cuts + [spare])]
    lectures, privates = [], []      # arcs as (start, length)
    pos = 0
    for piece in pieces:
        width = 3 + piece
        ell = rng.randint(3, width)
        rs = pos + rng.randint(0, width - ell)
        if rng.random() < 0.6:       # crossing pair
            q = rng.randint(1, ell - 2)
            p = rng.randint(q + 1, ell - 1)
            lectures.append((rs, p))
            lectures.append((rs + q, ell - q))
            privates.append((rs, q))             # inside a \ b
            privates.append((rs + p, ell - p))   # inside b \ a
        else:                        # nested pair b inside a
            v = rng.randint(1, ell)
            u = rng.randint(0, ell - v)
            lectures.append((rs, ell))
            lectures.append((rs + u, v))
            if u > 0:
                privates.append((rs, u))
            if u + v < ell:
                privates.append((rs + u + v, ell - u - v))
        pos += width
    for _ in range(extras):
        if not privates:
            break
        s, d = rng.choice(privates)
        x = rng.randint(0, d - 1)
        y = rng.randint(x + 1, d)
        lectures.append((s + x, y - x))
    shift = rng.randrange(period)
    pairs = [((s + shift) % period, (s + d + shift) % period) for s, d in lectures]
    rng.shuffle(pairs)
    return pairs


def mutate(rng, pairs, period):
    """Random small change: move an endpoint, add or drop a lecture."""
    pairs = list(pairs)
    kind = rng.randrange(3)
    if kind == 0 and pairs:
        i = rng.randrange(len(pairs))
        s, e = pairs[i]
        pairs[i] = (rng.randrange(period), e) if rng.random() < .5 else (s, rng.randrange(period))
    elif kind == 1:
        pairs.append((rng.randrange(period), rng.randrange(period)))
    elif len(pairs) > 2:
        del pairs[rng.randrange(len(pairs))]
    return pairs


def _compare(inst, period):
    """Return (solvable, ok) comparing solve() with the brute force."""
    got = solve(inst, period)
    ref = brute_force(inst, period)
    ok = (got is None) == (ref is None)
    if got is not None:
        ok = ok and verify(inst, got, period)
    # the solver's bit-mask graph must agree with the arithmetic graph
    arcs = to_arcs(inst, period)
    ok = ok and build_graph(arcs, period)[1] == adjacency_arithmetic(arcs, period)
    return ref is not None, ok


def stress_category(kind, trials, seed, period, max_n):
    """Run one category; return (instances, solvable, disagreements)."""
    rng = random.Random(seed)
    solvable = bad = 0
    for _ in range(trials):
        if kind == "random":
            inst = random_instance(rng, rng.randint(2, max_n), period)
        else:
            k = rng.randint(1, min(max_n // 2, period // 3))
            inst = planted_instance(rng, period, k, rng.randint(0, max_n - 2 * k))
            if kind == "planted":
                if solve(inst, period) is None:    # must be solvable
                    bad += 1
            else:                                  # "mutated"
                inst = mutate(rng, inst, period)
        sol, ok = _compare(inst, period)
        solvable += sol
        bad += (not ok)
    return trials, solvable, bad


STRESS_PLAN = [  # (label, kind, period, max_n)
    ("random, period 24", "random", 24, 9),
    ("planted, period 24", "planted", 24, 10),
    ("mutated, period 24", "mutated", 24, 10),
    ("planted, period 1440", "planted", DAY, 10),
    ("mutated, period 1440", "mutated", DAY, 10),
]


def run_stress(trials_each, seed=1):
    return [(label,) + stress_category(kind, trials_each, seed + i, period, max_n)
            for i, (label, kind, period, max_n) in enumerate(STRESS_PLAN)]


# --------------------------------------------------------------------------
# 7. Benchmark
# --------------------------------------------------------------------------
def bench_instance(n, worst, seed=7):
    """Benchmark instance with about n lectures.

    worst=False: planted (solvable) instance.
    worst=True : adversarial for the number of anchors. Lecture 1 is a
        half-day arc that contains n/2 tiny private lectures, so every pair
        {lecture 1, tiny lecture} is a block meeting lecture 1 (n/2 anchors);
        the other half of the day holds ~n/4 crossing blocks, and one extra
        isolated lecture can never be supervised, so no solution exists and
        the DP has to run for every anchor.
    """
    if not worst:
        rng = random.Random(seed)
        k = max(1, n // 4)
        return planted_instance(rng, DAY, k, n - 2 * k)
    half = DAY // 2
    tiny = n // 2
    inst = [(0, half)]
    inst += [(2 * i, 2 * i + 1) for i in range(tiny)]
    kb = max(2, (n - tiny - 2) // 2)
    for i in range(kb):
        w = half + 4 * i
        inst += [(w, w + 2), (w + 1, w + 3)]
    inst.append((DAY - 3, DAY - 2))          # isolated lecture
    return inst


def _best_time(inst, repeat=3):
    best = INF
    for _ in range(repeat):
        t0 = time.perf_counter()
        solve(inst)
        best = min(best, time.perf_counter() - t0)
    return best


def run_bench(sizes=(100, 200, 400, 800)):
    rows = []
    for n in sizes:
        easy = bench_instance(n, False)
        hard = bench_instance(n, True)
        rows.append((len(easy), _best_time(easy), len(hard), _best_time(hard, 1)))
    return rows


# --------------------------------------------------------------------------
# 8. Reporting
# --------------------------------------------------------------------------
def report(pairs, name="input"):
    """Print the answer for one instance; return the chosen list or None."""
    sol = solve(pairs)
    print(f"--- {name} (n = {len(pairs)}) ---")
    if sol is None:
        print("No valid subset exists.")
        return None
    print("Selected lectures (1-indexed):", [c + 1 for c in sol])
    cmask = sum(1 << c for c in sol)
    arcs = to_arcs(pairs)
    adj = adjacency_arithmetic(arcs, DAY)
    for v, (s, d) in enumerate(arcs):
        sup = [u + 1 for u in range(len(arcs)) if (adj[v] >> u & 1) and (cmask >> u & 1)]
        print(f"lecture {v + 1:>3} {fmt_time(s)}-{fmt_time((s + d) % DAY)}"
              f" supervised by {sup}")
    print("Independent verification:", "PASSED" if verify(pairs, sol) else "FAILED")
    return sol


def _tex(s):
    return s.replace("&", r"\&").replace("_", r"\_").replace("#", r"\#")


def latex_report():
    out = [r"\begin{longtable}{p{4.6cm}p{5.6cm}p{2.6cm}c}",
           r"\toprule Test & Lectures (start--end) & Algorithm output & OK \\ \midrule \endhead"]
    all_ok = True
    for name, pairs, expect in TESTS:
        got = solve(pairs)
        ref = brute_force(pairs)
        ok = ((got is None) == (ref is None) == (not expect)) and \
             (got is None or verify(pairs, got))
        all_ok &= ok
        lec = ", ".join(f"L{i + 1}\\,{fmt_time(s)}--{fmt_time(e)}"
                        for i, (s, e) in enumerate(pairs))
        res = "none" if got is None else "$S=\\{" + ",".join(str(c + 1) for c in got) + "\\}$"
        mark = "\\checkmark" if ok else "\\texttimes"
        out.append(f"{_tex(name)} & {{\\footnotesize {lec}}} & {res} & {mark} \\\\")
    out.append(r"\bottomrule \end{longtable}")

    rows = run_stress(200)
    out.append(r"\begin{center}\begin{tabular}{lrrr}")
    out.append(r"\toprule Instance family & Instances & Solvable & Disagreements \\ \midrule")
    total_bad = 0
    for label, trials, sol, bad in rows:
        total_bad += bad
        out.append(f"{_tex(label)} & {trials} & {sol} & {bad} \\\\")
    out.append(r"\bottomrule \end{tabular}\end{center}")
    good = all_ok and total_bad == 0
    out.append("\\par\\noindent\\textbf{" +
               ("All built-in tests and all cross-checks passed." if good else
                "WARNING: a test or cross-check failed.") + "}")
    return "\n".join(out)


def latex_bench():
    out = [r"\begin{center}\begin{tabular}{rrrrr}",
           r"\toprule $n$ & solvable: time (s) & $n$ & worst case: time (s) & growth \\ \midrule"]
    prev = None
    for n1, t1, n2, t2 in run_bench():
        growth = "--" if prev is None or prev == 0 else f"$\\times{t2 / prev:.1f}$"
        out.append(f"{n1} & {t1:.3f} & {n2} & {t2:.3f} & {growth} \\\\")
        prev = t2
    out.append(r"\bottomrule \end{tabular}\end{center}")
    return "\n".join(out)


# --------------------------------------------------------------------------
# 9. Command line
# --------------------------------------------------------------------------
def _run_demo():
    all_ok = True
    for name, pairs, expect in TESTS:
        sol = report(pairs, name)
        all_ok &= (sol is not None) == expect and (sol is None or verify(pairs, sol))
    print("All demo cases passed." if all_ok else "SOME DEMO CASES FAILED.")
    return 0 if all_ok else 1


def _run_stress(trials):
    if trials < 1:
        print("--stress must be a positive integer.", file=sys.stderr)
        return 2
    each = max(1, trials // len(STRESS_PLAN))
    bad_total = 0
    for label, n, sol, bad in run_stress(each):
        print(f"{label:<24} {n:>5} instances, {sol:>4} solvable, {bad} disagreements")
        bad_total += bad
    print("Stress test PASSED." if bad_total == 0 else "Stress test FAILED.")
    return 0 if bad_total == 0 else 1


def _run_bench():
    print(f"{'n':>5} {'solvable (s)':>13} | {'n':>5} {'worst case (s)':>15} {'growth':>8}")
    prev = None
    for n1, t1, n2, t2 in run_bench():
        growth = "" if not prev else f"x{t2 / prev:.1f}"
        print(f"{n1:>5} {t1:>13.3f} | {n2:>5} {t2:>15.3f} {growth:>8}")
        prev = t2
    return 0


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("file", nargs="?", help='input file ("-" = standard input)')
    ap.add_argument("--demo", action="store_true", help="run the 8 built-in tests")
    ap.add_argument("--stress", type=int, metavar="N",
                    help="compare with brute force on about N random/planted/mutated cases")
    ap.add_argument("--bench", action="store_true", help="print running times")
    ap.add_argument("--latex-report", action="store_true")
    ap.add_argument("--latex-bench", action="store_true")
    a = ap.parse_args()

    if a.latex_report:
        print(latex_report())
        return 0
    if a.latex_bench:
        print(latex_bench())
        return 0
    if a.bench:
        return _run_bench()
    if a.stress is not None:
        return _run_stress(a.stress)
    if a.demo:
        return _run_demo()
    try:
        pairs = read_file(a.file) if a.file else interactive()
    except (ValueError, OSError) as err:
        print(f"Input error: {err}", file=sys.stderr)
        return 2
    report(pairs, a.file or "typed input")
    return 0


if __name__ == "__main__":
    sys.exit(main())
