# Lab Assignment 7 - Multiobjective Optimization (C++17)

All 37 algorithms, the shared toolkit (`moo_utils`), problems P1-P14 and the experiment driver are written
from scratch in **C++17 using only the standard library** (no Eigen/Armadillo/NLopt, no optimisation library).

## 1. What you need to install

Only a C++17 compiler and `make`. Nothing else is required to build and run.

| OS | Command |
|---|---|
| Ubuntu / Debian / WSL | `sudo apt update && sudo apt install -y g++ make` |
| macOS | `xcode-select --install` (gives clang++ and make) |
| Windows | Easiest: install WSL (Ubuntu) and use the Ubuntu line. Or MSYS2/MinGW-w64: `pacman -S mingw-w64-ucrt-x86_64-gcc make` |

Check: `g++ --version` (GCC 9+ or clang 10+ works).

**Optional** (only for drawing figures): `python3` + `pip install matplotlib numpy`, then `python3 plot_results.py`.
Plotting is not part of any algorithm; you may use any other tool (Excel, gnuplot) on the CSV files instead.

## 2. Build

```
cd moo_lab7
make            # builds ./moo   (or: g++ -std=c++17 -O2 main.cpp moo_utils.cpp -o moo)
```
On Windows/MinGW the binary is `moo.exe`; replace `./moo` by `moo` below.

## 3. Run

All results are written to `results/` (CSV + text). Seeds are fixed in the code.

| Command | What it does | Approx. time |
|---|---|---|
| `./moo unit` | Unit tests: toolkit (hand-checked cases), Deb's 3 rules on P6, problem sanity checks from the sheet, one test per algorithm | 5 s |
| `./moo log` | **Task A3** implementation log (`results/implementation_log.csv`) + fronts in `results/fronts/` (budget 10,000) | 20 s |
| `./moo trace` | **Task C hand traces**: TOPSIS (alg 8), IBEA (alg 32), MOEA/D-ACO (alg 20), every number printed | 1 s |
| `./moo sens` | **Task C sensitivity** (kappa for IBEA, rho for MOEA/D-ACO, weights for TOPSIS), 10 seeds | 1-2 min |
| `./moo compare` | **Task D** head-to-head, 10 seeds: P3, P5, P10 (M=5), P8; ranked tables in `results/taskD_tables.csv` | 5-15 min |
| `./moo mcdm` | TOPSIS and VIKOR on P12 and on the NSGA-II front of P4, weight-shift test | 1 s |
| `./moo p11` | Goal Programming (unit / normalised weights) and Lexicographic on P11 | 1 s |
| `./moo mtl` | MGDA, GradNorm vs fixed weights on P13(a) (`results/P13_trajectories.csv`) | 1 s |
| `./moo investigate` | Numbers answering the "Investigate" questions (P1, P2, P3, P5, P6, P10, P14, P7) | 10 s |
| `./moo all` | everything above | ~15 min |

Options: `--budget N` (default 10000 evaluations), `--seeds K` (default 10), e.g. `./moo compare --seeds 3 --budget 5000` for a quick test.

## 4. Files

| File | Content |
|---|---|
| `moo_utils.hpp/.cpp` | RNG, dominance, Deb's constrained dominance, fast non-dominated sort, crowding distance, HV (exact 2-D, Monte-Carlo M>=3), IGD, Spacing, Das-Dennis points, linear solver |
| `problems.hpp` | P1-P7, P10, P11, P13, P14 (objectives + constraint violation, normalisation bounds, true fronts) |
| `classical.hpp` | Algs 1-9 (own scalar solver: (1+1)-ES + Nelder-Mead), TOPSIS, VIKOR |
| `evo.hpp` | Algs 10-19, 21, 30-35 |
| `swarm.hpp` | Algs 22-24, 26-29 |
| `aco.hpp` | P8 (TSP, seed 7), algs 20 and 25, permutation NSGA-II and MOMA |
| `grad.hpp` | Algs 36 MGDA and 37 GradNorm |
| `main.cpp` | Registry of the catalogue, metrics, unit tests, all experiments |

## 5. Conventions and honest limitations 

* All objectives are minimised internally (maximisation is negated). HV uses objectives normalised by fixed per-problem bounds and reference point 1.1; Spacing is Schott's; IGD uses the analytic front where available (P1, P2, P3, P6, P7, P10, P13, P14), otherwise the pooled non-dominated union of all runs.
* Constraints: algorithm 14 (plain NSGA-II) ignores constraints on purpose (to show what happens without them); algorithm 21 and all others use Deb's rules (or an equivalent penalty for scalarising methods).
* Algorithms 1-4 sweep 20 settings; 5-7 report the solution for each preference input; 8-9 are ranking methods (not optimisers).
* Not done / simplified: **P9** (feature selection) needs the instructor's CSV and a hand-coded 5-NN + 5-fold CV, so MOGWO and MOMA are tested on P5 instead; **P13(b)** (neural network) is not implemented; P12 and P13(a) are handled by `mcdm` and `mtl`. The lecture "car" example is not in the sheet, so the toolkit tests use small hand-checkable cases instead. Add your lecture example to `unitTests()`.
* MOEA/D-ACO, MOBA, MOSA, MO-TLBO are compact versions of the published methods (see the comments on each function); say so in the report.
* Task B tables, the study cards, the written hand-trace explanations, plots and the PDF report are yours to write from the numbers in `results/`. The assignment grades your own explanations in the viva.

## 6. Deep-dive algorithms (Task C)

Assumption: roll number 241070045 -> R = 45 (last two digits). Algorithm 1 = ((45-1) mod 37)+1 = **8 TOPSIS**, Algorithm 2 = ((45+11) mod 37)+1 = **20 MOEA/D-ACO**, Algorithm 3 = ((45+23) mod 37)+1 = **32 IBEA**. If your *registration* number ends in different digits, recompute and tell me; the code for all 37 is already present, only the traces/sensitivity choices change.

## 7. Academic integrity

The sheet allows AI help for debugging but requires that you can explain every line in the viva and that explanations come from your own experiments. Read each module, run it, change a line and predict the effect before submitting.
