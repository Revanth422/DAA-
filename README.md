# DAA-
# Lecture Supervision on a Circular Timetable

Design and Analysis of Algorithms -- course project.

`n` lectures repeat every day on a circular timetable of `L = 1440` minutes (a lecture may
wrap past midnight). The program selects a subset `S` of the lectures such that **every
lecture overlaps exactly one lecture of `S`** (a lecture never counts itself), or reports
that no such subset exists.

In graph terms this is an *efficient open dominating set* (total perfect code) of the
circular-arc overlap graph. The algorithm runs in **O(n^3)** time (for fixed `L`) and
O(n^2) bits of space; correctness proofs, the complexity analysis and experiments are in
[`report.pdf`](report.pdf) (LaTeX source: [`report.tex`](report.tex)).

## Files

| File | Purpose |
|------|---------|
| `lecture_solver.py` | solver, independent verifier, brute-force checker, test generators, benchmark, command-line tool |
| `sample.txt` | example input (8 lectures, 3 blocks) |
| `mydata.txt` | a small input file to edit and try your own data |
| `report.tex`, `report.pdf` | the project report (source and compiled PDF) |
| `generated_tests.tex`, `generated_bench.tex`, `demo_output.txt`, `sample_output.txt` | files the report builds from the program output; kept in the repository so the PDF also builds where `-shell-escape` is unavailable |
| `.github/workflows/test_solver.yml` | continuous integration: runs the tests on every push |

## Requirements

Python 3.8 or newer (standard library only, nothing to install).

## Input format

First line: `n`. Then `n` lines, one per lecture: `start end` as `HH:MM` (24-hour clock,
`24:00` means midnight). Text after `#` is ignored. Plain integers are read as minutes after
midnight.

```
4
23:00 01:00     # lecture 1 wraps past midnight
00:30 02:00
13:00 15:00
14:00 16:00
```

Conventions:

* Intervals are half-open: a lecture ending at 10:00 and one starting at 10:00 do **not** overlap.
* `start end` equal (e.g. `00:00 00:00`) means a full 24-hour lecture.
* Lectures are numbered `1..n` in the order given.

## Usage

```
python3 lecture_solver.py                  # type the lectures in interactively
python3 lecture_solver.py mydata.txt       # solve a file
cat mydata.txt | python3 lecture_solver.py -   # read from standard input
python3 lecture_solver.py --demo           # the 8 built-in test cases (T1..T8)
python3 lecture_solver.py --stress 1000    # compare with brute force (random, planted, mutated)
python3 lecture_solver.py --bench          # running-time measurements
```

Exit status: `0` on success (also when the answer is "no valid subset"), `1` if a test or
stress run fails, `2` on bad input.

## Example

```
$ python3 lecture_solver.py mydata.txt
--- mydata.txt (n = 4) ---
Selected lectures (1-indexed): [1, 2, 3, 4]
lecture   1 23:00-01:00 supervised by [2]
lecture   2 00:30-02:00 supervised by [1]
lecture   3 13:00-15:00 supervised by [4]
lecture   4 14:00-16:00 supervised by [3]
Independent verification: PASSED
```

If no subset exists the program prints `No valid subset exists.`

## How it works (short)

1. A **block** is a pair of overlapping lectures `{a, b}` with no third lecture overlapping
   both. Every solution is a disjoint union of blocks such that every lecture hits exactly
   one block region, and there are at most `2n` candidate blocks.
2. Lecture 1 hits exactly one block in any solution; each candidate block hitting it is tried
   as an **anchor**. The circle is rotated so the anchor starts at 0 and unrolled.
3. A dynamic programme chains the remaining blocks from left to right; two consecutive
   blocks are compatible if no lecture meets both and no lecture lies strictly in the gap
   (tests with a prefix-maximum `Q` and a suffix-minimum `T`).

See the report for the full proofs.

## Testing

* `--demo` runs the 8 hand-written cases (wrap-around, touching endpoints, full-day lecture, ...).
* `--stress N` compares the polynomial algorithm with exhaustive search on random,
  *planted* (solvable by construction) and *mutated* instances, for periods 24 and 1440.
* The answers are checked by `verify`, and the exhaustive search uses its own arithmetic
  overlap test, so neither shares code with the solver's bit-mask graph.
* GitHub Actions (`.github/workflows/test_solver.yml`) runs all of this on every push.

## Rebuilding the report

```
pdflatex -shell-escape report.tex      # run it twice for the table of contents
```

`-shell-escape` lets LaTeX run `lecture_solver.py` so that the tables and raw outputs in the
report are the real program output. Without it (for example on Overleaf) the committed
`generated_*.tex` and `*_output.txt` files are used. Before building, fill in your name and
roll number at the top of `report.tex`.
