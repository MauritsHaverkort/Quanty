# Stage-A prototype (Python) — reference implementation of the checkpoints

Scratch implementation of the equations of `QuantyDFT_StageA_DiracAtom.md`, used to produce the prototype numbers in §12. Requires numpy, scipy and sympy. Keep all files in one directory (`dkb.py` is imported by the others).

| File | Checkpoint |
|---|---|
| `dkb.py` | DKB B-spline Dirac solver (point nucleus); `python3 dkb.py` runs check 1a |
| `check1b_2_nucleus.py` | 1b (point-nucleus intruder, Z = 92, κ = +1) and 2 (Gaussian nucleus, U) |
| `check3_confinement.py` | 3 (Lorentz-scalar vs vector confinement) |
| `check4_poisson.py` | 4 (Galerkin Y^k); argument = number of knot intervals, e.g. `python3 check4_poisson.py 160` |
| `check5_aoc_vs_cowan.py` | 5 (configuration average = Cowan, symbolic) |
| `check6_fock_gradient.py` | 6 (Fock operator vs finite-difference gradient) |
| `check6_fock_gradient.cpp` | 6 for the C code (`DFT/AtomHF.cpp`); a C++ harness that includes Quanty, build command in its header |

Note: `dkb.py` uses the order k (degree k−1); the Quanty C core uses the degree.
