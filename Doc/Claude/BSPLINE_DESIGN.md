# B-spline consolidation: design note

Status: proposal, 2026-11-09. Companion to `QuantyDFT_Step0_BSplineConsolidation.md`
(steps 0.1–0.5 and 0.7). The CODATA update (0.6) is done separately and first.

## 1. Decisions taken

- **C core, class kept.** New plain-array C functions do the work. The class
  `BSpline::TBSpline : Interpolation::InterpolatingFunction` stays and calls them.
  Replacing the class hierarchy with the 2018 `InterpolatingFunctionType` struct
  (commented out in `Core/StructDef.cpp:1838`) is a later, separate project.
- **Location.** Everything reachable from a main Quanty function (`NewOperator`,
  `Interpolate`, `ReadFPLOBasisFunctions`, `GetSlaterIntegrals`, …) moves to
  `Core/BasicMath/Interpolation/`. Code reachable only as `Martin.<Name>` stays in
  `Martin/` and calls the core.
- **Old code.** `EvaluateOld` and `EvalDerivativeOld` stay in `BSpline.cpp` inside
  `#if 0 … #endif`: kept for the record, not compiled.
- **Outside the knot range** a basis function, and a B-spline function built from
  them, is 0, value and derivatives alike.
- **Tests** use analytic functions only; an FPLO test comes later, with the tutorials.
- `Martin.HF` need not reproduce its old output.

## 2. Files after the change

| File in `Core/BasicMath/Interpolation/` | Contents |
|---|---|
| `GaussLegendre.{h,cpp}` | The 34-digit table and `GaussLegendreQuadrature(n)`, `…Node`, `…Weight`, moved from `BasicMath_StandardFunctions`. New `GaussLegendreNodesWeights(N,a,b,x,w)`: table for N ≤ 127, Golub–Welsch (moved from `Martin/Math/Quad.cpp`) above. The old `Quad::` entry points become thin wrappers. |
| `BSpline.{h,cpp}` | The C core (§3) and the C-named successors of the namespace functions (§4). |
| `BSplineMatrix.{h,cpp}` | The matrix builders now in `Martin/Hartree-Fock/Hartree-Fock.cpp` that `AtomicOperators` reaches: `CreateLaplacian`, `CreatePotentialMatrix`, `PutBSplinesOnQuadGrid`, `CalcOffDiagSlaterInts`, `CalcDiagSlaterInts`, `AssembleSlaterInts2`. They are rebuilt on `BSplineMatrixBanded`. |
| `TBSpline.{h,cpp}` | The class; its internals call the core. |

`Core` then no longer includes anything from `Martin/`. `AtomicOperators/GaussQuadrature.cpp`
keeps `GetMultipoleIntegral`, which is physics; its `GetKnots`/`GetWeights` become
wrappers.

## 3. The C core

Degree convention throughout (order = degree + 1), matching Lua `"Degree"`. Plain
arrays, return 0 on success and 1 on error, no printing inside the core. The caller
allocates what it can size; where the core allocates (`double **`), the caller frees.

```c
uint32_t BSplineFindInterval(double t, const double *knot, uint32_t nKnot, uint32_t degree);
int BSplineValues(double t, uint32_t interval, const double *knot, uint32_t degree,
                  uint32_t nDeriv, double *val);     /* val[(nDeriv+1)*(degree+1)] */
int BSplineKnotsExponential(uint32_t degree, uint32_t nInterval, double r1, double rMax,
                  double **knot, uint32_t *nKnot);
int BSplineQuadratureSetup(const double *knot, uint32_t nKnot, uint32_t degree, uint32_t nQ,
                  double **rq, double **wq, uint32_t *nPoint);
int BSplineMatrixBanded(const double *knot, uint32_t nKnot, uint32_t degree,
                  const double *rq, const double *wq, uint32_t nPoint, const double *fq,
                  uint32_t a, uint32_t b, double *band);
```

- **`BSplineFindInterval`**: binary search for μ with `knot[μ] ≤ t < knot[μ+1]`. At
  `t == knot[last]` it returns the last non-empty interval, which is the convention of
  the present `Evaluate`.
- **`BSplineValues`**: de Boor's `bsplvd`. The triangular recurrence gives the
  degree+1 nonzero splines `B[μ-degree … μ]` at `t`. Derivatives come from the same
  table by differencing the coefficients. Cost O(degree²) per point, against
  O(2^degree) per spline for the present recursion. Only convex combinations enter the
  values, so the routine is stable for repeated knots too.
- **`BSplineQuadratureSetup`**: on every non-empty knot interval, nQ Gauss–Legendre
  points from the table. Exact for polynomials of degree 2nQ−1 on each interval, so
  nQ = degree + 1 integrates a product of two splines exactly.
- **`BSplineMatrixBanded`** is what the old note called `BSplineBandedAssemble`. It
  computes

  M_ij = ∫ B_i^(a)(r) f(r) B_j^(b)(r) dr ≈ Σ_q w_q f(r_q) B_i^(a)(r_q) B_j^(b)(r_q),

  for derivative orders a, b ∈ {0,1,2}. Only |i−j| ≤ degree is nonzero, so it is
  stored in LAPACK band format. f enters as its values `fq` on the quadrature points,
  with no callback, which keeps it OpenMP-safe. One routine then replaces Martin's
  separate builders:

  | f, a, b | Matrix |
  |---|---|
  | 1, 0, 0 | overlap S |
  | 1/2, 1, 1 | kinetic energy ½∫B_i′B_j′ |
  | −Z/r + l(l+1)/(2r²), 0, 0 | potential and centrifugal |
  | r^k or r^(−k−1), 0, 0 | Slater-integral pieces |

## 4. Namespace functions → C names

| Old (`namespace BSpline`) | New | Change |
|---|---|---|
| `Evaluate(i,n,x,t)` | `BSplineEvaluate` | FindInterval + Values, 0 outside the support |
| `EvalDerivative` | `BSplineEvaluateDerivative` | same |
| `Derivative(idx,…)` | `BSplineDerivative` | sum over the degree+1 nonzero splines only |
| `deBoor` | `BSplineDeBoor` | unchanged |
| `GenerateGrid`, `GetKnotMatrix` | `BSplineKnotsAveraged`, `BSplineCollocationMatrix` | collocation matrix in band storage (0.7) |
| `GenerateFisherGrid` | `BSplineKnotsFischer` | declared in the header |
| `Overlap`, `OverlapMatrix` | via `BSplineMatrixBanded` | |
| `Components` | `BSplineProjection` | debug prints removed |
| `CalculateCellIntegrals`, commented-out `OverlapMatrix` | removed | stubs, no callers |

The arguments become plain arrays. There are 17 live call sites in 3 files: 2 in
`AtomicOperators/AtomicAux.cpp`, 5 in `Martin/functions_martin.cpp` and 10 in
`Martin/Hartree-Fock/Hartree-Fock.cpp`. The last 10 are rewritten anyway when the
builders move. The class name `BSpline::TBSpline` (33 uses) stays until the class
project.

## 5. Behaviour to keep, and two bugs

- `TBSpline::operator()(t)` outside the knot range prints an error and returns NaN
  today. **Decided (Maurits, 2026-11-09): it returns 0**, without a message, for the
  value and for every derivative. Basis functions already return 0 there.
- The vector `operator()` resets its interval index to `1` instead of `degree` when
  the points are not sorted. For degree > 1, `deBoor` then returns NaN. Using
  `BSplineFindInterval` fixes this; it changes results only where they were NaN.

## 6. Verification

- **Before any change:** `TestScripts/BSpline/` with `InterpolatingFunction.BSpline`
  on analytic functions (values and derivatives 0–3, both end points, repeated knots,
  degree 1–9). Record the reference output with the unmodified binary. Register the
  scripts in `.gitlab-ci.yml` in a new job `test-jobs-bspline`, built like the
  `test-jobs-math` job.
- **At an interior knot** the old code computes the value (`deBoor`) on the interval to
  the left. The derivatives (`Derivative` → `EvalDerivative`) come from the half-open
  degree-0 splines [x_i, x_i+1), which is the interval to the right. Only the
  derivative of order = degree jumps at a knot, so only it depends on the convention.
  The new code uses the interval to the right everywhere except at the last knot, which
  reproduces the old derivatives. The baseline test evaluates at every knot.
- **Integrate is not exact for B-splines.** `Interpolation::Integrate` uses 2 Gauss
  points per knot interval, exact up to cubics. For a B-spline of degree 9 the integral
  of a function accurate to 2.6e-12 is off by 1.6e-5. Fix: ⌈(degree+1)/2⌉ points for
  a B-spline. This is a separate commit after the quadrature move, which re-records
  only the integrals.
- **Old vs new kernels:** a C harness in the scratchpad, not in the repository. It
  compiles `EvaluateOld`/`EvalDerivativeOld` against the new code and requires a
  maximum relative difference < 1e-13 for degree 1–9, random and repeated knots, and
  both ends.
- `TestScripts/Math/GaussLegendreQuadrature.Quanty` unchanged.

## 7. Commit order on `Claude`

1. CODATA 2022 and `cAU` (step 0.6)
2. Baseline tests `TestScripts/BSpline/`
3. Gauss–Legendre and `Quad` into `Interpolation/`; Core no longer includes `Martin/`
3b. `Integrate` exact for B-splines
4. C core plus the renamed namespace functions
5. Dead code out
6. Matrix builders out of `Martin/Hartree-Fock/` into `BSplineMatrix`
7. Banded collocation in `TBSpline::Create` (0.7)

Later, in separate commits as preparation for the next sessions: 0.7b (`Z` as a
double, the HF `eps` option, the `float` check) and 0.8 (`DFT/`, the `DKBSpline`
skeleton).
