# Quanty DFT/HF add-on — Stage A design: relativistic atom on B-splines

Status: design, ready for implementation. Step 0 is done (branch `Claude`, up to commit 9a14a753; interfaces in `Doc/Claude/BSPLINE_DESIGN.md`). Its open remainder is the first task of stage A (§0.2). Date: 2026-11-11, revision 3.
- Revision 2 added: nucleus models including Lua functions and non-spherical nuclei; $V_\beta$ notation; on-site orthonormalization; Lua interface aligned with `ReadFPLOBasisFunctions`, `NewOperator("Dirac")` and `NewOperator("AU")`.
- Revision 3: §2, §5, §7 and §11 refer to the B-spline core that now exists (`BSplineFindInterval`, `BSplineValues`, `BSplineKnotsExponential`, `BSplineQuadratureSetup`, `BSplineMatrixBanded`, `GaussLegendreNodesWeights`). It adds §0.2 (status of step 0 and stage-A task A0) and adopts the step-0 conventions: 0 outside the knot range, and `cAU` read at call time.
Scope of stage A: atomic eigenstates and radial 4-spinor functions on a B-spline basis, with a Lorentz-scalar confining potential; Dirac–Hartree–Fock (configuration average, Cowan-compatible) and Dirac–Kohn–Sham (LDA, PBE). The output is a set of contracted radial spinors that become the local basis for stages B–E.

## 0. Fixed decisions

| Item | Decision |
|---|---|
| Units | Hartree atomic units. The speed of light in atomic units is the entry **`cAU`** of the `PhysicalConstants` table (`Core/StructDef.cpp`), $=\alpha^{-1}=137.035999177$ (CODATA 2022). It is a pure number like `alpha`. To change it, the user sets `PhysicalConstants.cAU.value` in Lua; for example, `1e8` gives the nonrelativistic limit. Lua bindings read `PhysicalConstants.cAU.value` at call time and pass $c$ as an explicit argument to the C kernels. C code never uses a global or literal $c$. Energies measured from $mc^2$. |
| Precision | Double precision everywhere (no `float`). Every tolerance is a user option. Designed for high-precision atomic physics; see §0.1. |
| Angular basis | 4-spinors $(n,\kappa,\mu)$, jj order: per $l$-shell first $j=l-\tfrac12$ then $j=l+\tfrac12$, $\mu=-j\ldots j$ ascending. A single fixed unitary maps to tesseral $S_{lm}\otimes\chi_\sigma$ (Quanty phase convention) at the Operator boundary. |
| Relativity | Full 4-component Dirac; noncollinear (noncollinearity is irrelevant for the spherical atom of stage A; enters stage D). |
| Radial functions | Numerical only: every basis function is a contraction over a B-spline dual-kinetic-balance (DKB) set. |
| Confinement | Lorentz scalar $V_\beta(r)$, entering as $\beta V_\beta$. The symbol $S$ is reserved for overlap matrices. Code name `Vbeta`. |
| Nucleus | `"Point"`, `"Gaussian"`, `"Fermi"`, `"Function"` (Lua function $\rho_N(r)$); `"Function3D"` ($\rho_N(x,y,z)$, multipole-expanded; see §3). |
| On-site orthonormality | All local basis orbitals of one atom are orthonormal, also when generated with orbital-dependent confinement or ionicity (§6.1). Inter-site overlap is handled later by $S(\mathbf R)$ in TightBinding. |
| B-spline convention | The Lua and C interfaces use the polynomial **degree** $d$, matching the existing Quanty option `"Degree"` and Martin's `BSpline::` code where `n` = degree. This document uses the order $k=d+1$. |
| Lua integration | A new function `CalculateAtomicBasisFunctions` returns `G, F, x` exactly as `ReadFPLOBasisFunctions` does, plus a fourth return value, so that `NewOperator("Dirac",…)` and `NewOperator("AU",…)` consume the result unchanged (§11). |
| Prerequisite | Step 0 (B-spline, quadrature and constants consolidation) is done; see §0.2. The remainder (0.7b, 0.8 and the three conventions) is task A0, the first step of stage A. |
| Atomic HF | Configuration average over the **nonrelativistic** configuration (= trace over all determinants of $nl^q$; equals GRASP "EAL" over all CSFs generated from the configuration). Reduces exactly to Cowan's $E_{av}$ for $c\to\infty$ (§8.6, verified). |
| XC | Native: Slater exchange, PW92 and PBE. Additionally VWN5 and the MacDonald–Vosko relativistic exchange factor, needed only for validation against NIST (§12). No libxc. |
| TightBinding $S(\mathbf R)$ | Separate session between stage A and stages B–E. |

## 0.1 Precision requirements

One future application is high-precision physics. The code must therefore never be the limiting factor, and its accuracy must be controllable and reported.

1. **Scalar type.** `double` throughout, and `uint32_t` for indices.
   - No `float`, no single-precision math functions (`sqrtf`, `expf`, …) and no truncated literals such as `137.0`.
   - All new code uses a `typedef double DFTReal;` so that a later `long double` or `__float128` build of the LAPACK-free kernels (quadrature, Poisson solver, Slater integrals) stays possible.
   - Physical inputs are `double`, including $Z$ and $A$: fractional $Z$ is used for VCA-type tests. The existing `NewOperator("Dirac")` still takes `int Z` (`DiracH(…, const int Z)`); that is fixed in task A0.
2. **No hard-coded tolerances.** Every threshold is an option with a documented default:
   - SCF energy and orbital gradient (§8.5), potential mixing residual (§9);
   - orthonormality check and redundancy threshold (§6.1);
   - normalization check of Lua nuclear densities (§3);
   - quadrature points per interval and number of knot intervals (§2).

   A preset `{"Precision","High"}` tightens all of them together. Defaults: $\Delta E<10^{-10}$ Ha, gradient $<10^{-8}$. High: $\Delta E<10^{-13}$ Ha relative to $|E|$, gradient $<10^{-11}$, $n_q=k+8$, and a knot-convergence report.
3. **Summation accuracy.** Accumulate quadrature sums for energies, Slater integrals, normalizations and multipole moments with compensated (Neumaier) summation. Matrix assembly may use BLAS.
4. **Reported accuracy.** The `Atom` table returns:
   - eigen-residuals $\|Hc-\varepsilon Sc\|$;
   - orthonormality error;
   - the Brillouin residual (§8.5);
   - the difference between the two total-energy expressions (§9);
   - the virial ratio in the nonrelativistic limit;
   - the change of $E$ and $\varepsilon$ under one knot refinement (optional, `{"ConvergenceCheck",true}`).
5. **Constants.** Use one consistent CODATA set (2022) from the core table (`Core/StructDef.cpp`, updated in step 0), with conversion 1 bohr = 52917.7210544 fm. Echo the constants used into the `Atom` table, so every result records its $c$, $a_0$ and nuclear parameters.
6. **Quadrature tables.** The Gauss–Legendre table (`Core/BasicMath/Interpolation/GaussLegendre`) has 27–34 significant digits for N ≤ 127. That is sufficient for double and long double; for a future `__float128` build, regenerate it with ≥36 digits. Above 127 points the Golub–Welsch fallback is accurate only to about 4e-10 in the weights. Stage A must either stay at N ≤ 127, which is always enough per knot interval, or refine the fallback first.
7. **Extensibility for high-precision physics** (not in stage A, but not blocked by the design):
   - user-supplied additional radial potentials, `{"ExtraPotentials",{…}}`, sampled in the binding as in §3. This covers QED model potentials such as Uehling vacuum polarization or effective self-energy potentials;
   - finite nuclear mass;
   - the Breit interaction: Martin's `BreitInteraction.cpp` exists and should be reviewed then.

## 0.2 Status of step 0 and tasks A0, A1, A2

### What stage A can build on (branch `Claude`, commits f50b204f … 9a14a753)

| Item | Where | Notes for stage A |
|---|---|---|
| CODATA 2022 constants, `cAU` | `Core/StructDef.cpp`; test `TestScripts/PhysicalConstants.Quanty` | `alpha` is defined as `1/cAU`. `NewOperator("Dirac")` and `Martin.NewDirac` read `PhysicalConstants.cAU.value` at every call and raise an error if it is not a number. New bindings do the same. |
| One Gauss–Legendre source | `Core/BasicMath/Interpolation/GaussLegendre.{h,cpp}` | `GaussLegendreNodesWeights(N,a,b,x,w)` gives the 34-digit table for N ≤ 127. Above that, the Golub–Welsch fallback has weights accurate only to about 4e-10, so stage A must not request more than 127 points per interval (raise an error), or that fallback must first be refined (Newton polishing of the nodes plus the closed-form weights). |
| B-spline C core | `Core/BasicMath/Interpolation/BSpline.{h,cpp}` | `BSplineFindInterval`, `BSplineValues` (de Boor, derivatives as in Piegl and Tiller A2.3, O(d²) per point), `BSplineKnotsExponential` (added in A1), `BSplineDerivative` (coefficient differencing). Degree convention. Basis functions and splines are **0 outside the knot range**, value and derivatives alike. Unclamped knot vectors are handled. |
| Quadrature and banded matrices | `Core/BasicMath/Interpolation/BSplineMatrix.{h,cpp}` | `BSplineQuadratureSetup(knot,nKnot,nQ,…)` returns nQ points per non-empty interval. `BSplineMatrixBanded(…, fq, a, b, first, n, band, ldBand)` returns $\int B_i^{(a)}fB_j^{(b)}dr$ for $a,b\in\{0,1,2\}$, with $f$ given at the quadrature points (no callback), for basis functions `first…first+n-1`, in LAPACK band storage for $i\ge j$. Also Martin's Slater cell method (`BSplineCellValues`, `BSplineSlaterCellTables`, `BSplineSlaterIntegral`). |
| Exact integration | `InterpolatingFunction.Integrate`; `GetSlaterIntegrals`, `GetMultipoleIntegral` | `Integrate` is exact for B-splines of any degree. The atomic operators now use 34-digit rather than 6-digit quadrature. |
| `TBSpline` | banded collocation (`dgbsv`); arithmetic and `==` for B-splines | `ReadFPLOBasisFunctions(…,{{"Method","BSpline"}})` is 5× faster. |
| Tests | `TestScripts/InterpolatingFunction/` (one CI job) | CI passes for all step-0 pushes (checked by Maurits, 2026-11-11). |

### Tasks A0 and A1: done (2026-11-11, branch `Claude`, commits 5f74f9ae … a8983969)

**A0**
- $Z$ is a `double` along `LuaOperatorNewAtomicDirac` → `CreateOperatorDirac` (both overloads) → `CalculateDiracMatrixElements` → `DiracH`, and in `Martin.NewDirac` (`luaL_checknumber`). Test `TestScripts/DiracAtom/NewOperatorDiracZ`: $h(Z)$ linear to 2e-16; the old code rounded 1.5 to an integer.
- **HF `eps`: nothing to do.** `Martin.AtomicHartreeFock` has had it as the option `{"Precision", eps}` (default 1e-6) since 2019. Note the name clash with the planned preset `{"Precision","High"}` of §0.1 (a different function, so harmless).
- $O(h^2)$ note: as comments at `LuaOperatorNewAtomicDirac` and `LuaOperatorNewAtomicUGeneralFunction`. Correction to the step-0 plan: only `"Dirac"` uses `CreateWeights`; the general `"AU"` path uses a trapezoid-type sum (relativistic, `Slater_Integrate`) or 3 Gauss points per interval with the kinked kernel $r_<^k/r_>^{k+1}$ (non-relativistic); both are $O(h^2)$. A user-level note on quanty.org is still to be written (no page for `NewOperator("Dirac")` exists in the wiki drafts).
- `TestScripts/NoSinglePrecision.sh` (own CI job, no build): `float`, `sqrtf`, `expf`, `powf`, `logf`, `fabsf`, `erff`, `sinf`, `cosf`, `tanf`, `atanf` in `Core/BasicMath`, `AtomicOperators`, `DFT`, `Martin`. Passes.
- `DFT/` (`DFT.h` with `DFTReal`, `DFT.cpp`), included from `header.h` and `main.cpp`. `Interpolation::DKBSpline` with `lua_check` overload. Explicit errors for `"DKBSpline"` in `+ - * / ^ .. == unary-` and in `NewOperator("AU")` (any `G` or `F`).
- The three conventions: answered in §11.2 — `G` = $P$, `F` = $Q$, FPLO sign = this document's, $\mu$ ascending; $l=4$ works, $l=5,6$ labels are misparsed (plan there).
- Not done: the bookkeeping item (`grep "Quad::" Core/`), it was not in this session's scope.

**A1**
- `BSplineKnotsExponential` (§2.1). It was listed here and in `BSPLINE_DESIGN.md` as existing; it did not exist, and was written to the §2.1 formula.
- `BSplinePoissonSetup`/`BSplinePoissonSolve` in `BSplineMatrix.cpp` (§7, general $k$, Robin term from the spline values at $R$). Check 4 with a C harness: identical to `check4_poisson.py` for $k=0,1,2,4$; hydrogen 1s $Y^0$ to 4.0e-9 (80 intervals) and 4.5e-12 (160). $k>0$ has no CI test until A2.
- `DFT/DiracRadial.cpp`: `RadialBasisCreate`; nucleus models Point, Gaussian, Fermi, Function (Poisson, $k=0$), compensated normalization; `DiracDKBMatrices` from the explicit block list (in the source, checked against `dkb.py` to 1e-16 relative for $\kappa=-4\ldots3$); `DiracRadialEigen` (`dsygvd`); root selection; `DiracComponents`; sign convention.
- `DKBSpline` evaluation of any derivative order, exact; `Integrate` per knot interval.
- Temporary Lua entry `E, G, F, Info = DiracRadialSolve(Z, kappa, {options})` (`DFT/LuaDiracRadial.cpp`; options in its header comment).
- Tests `TestScripts/DiracAtom/` (CI job `test-jobs-diracatom`): `PointNucleus` (1a, 1b), `GaussianNucleus` (2), `FunctionNucleus` (13, spherical), `DKBSplineEvaluate`, `NewOperatorDiracZ`. All pass on macOS/Accelerate; the CI result is pending.

Results (degree 8, 13 points per interval, $c$ = 137.035999177):

| Check | Result |
|---|---|
| 1a | all roots ≤ 3.7e-9 Ha (most ≤ 1e-9); U 1s₁/₂ = −4861.1979032416 (−2.4e-8, discretization of $r^\gamma$); $Z=1$, $n=5$: 2.6e-6 (box); 126 of 252 states below $-c^2$ for every case |
| 1b | intruder at −16542.6751 Ha, $\int P^2=0.507$, $\int Q^2=0.493$, $\langle r\rangle=0.021$; then 2p₁/₂–4p₁/₂ to ≤ 1.6e-9 |
| 2 | −4854.1589200338, −1257.2395087575, −1089.6114149875 Ha, the same for 80, 120, 160 intervals to ≤ 7e-11; the prototype's 9 decimals to 4e-9. Fermi ($t$ = 2.30 fm): 1s₁/₂ = −4854.0905053681, 0.0684 Ha above Gaussian |
| 13 | `"Function"` with the Gaussian $\rho$: charge to 1.4e-14, $V_N$ to 1.8e-7 / 2.2e-9 / 8.0e-11 relative for 80 / 120 / 160 intervals (the test uses 160), orbital energies to 1.4e-14 relative already at 80 |

**Deviations from the design, decided or to be decided**
1. **Rayleigh quotient** (decided by Maurits, 2026-11-11: option, default on). The eigenvalues of a dense generalized solver err by about $\varepsilon\,\max|E|$; with a point nucleus and $r_1=10^{-6}/Z$ the DKB spectrum reaches $|E|\approx10^{12}$ Ha, so `dsygvd` alone gives ±1e-8 … 1.5e-7 Ha (both in C and in scipy; the prototype's quoted errors were one draw of this noise). `DiracRadialEigen` replaces each eigenvalue by $c^THc/c^TSc$; option `{"RayleighQuotient", false}` gives the plain values, which `Info.Edsygvd` always reports. The 1e-13 relative goal of `{"Precision","High"}` is reachable only with it.
2. **Singular integrands.** §2.3 said none is singular. For the first working spline ($B_1\propto r$) the DKB products contain $B_1^2/r^3$, $B_1'B_1/r^2$, $B_1'^2/r$ (and $B_1^2V/r^2$ for a point nucleus), which are logarithmically divergent for $|\kappa|\ge2$. The Gauss sum is finite, and the eigenvalues do not change (≤ 2e-10 Ha) when $B_1$ is left out for $|\kappa|\ge2$, or when the quadrature goes from $k+4$ to $k+16$ points (Gaussian nucleus). Only the point-nucleus $\kappa=+1$ intruder moves (by 1.5e4 Ha): it is a quadrature artefact. §2.3 is corrected.
3. **DKB functions at the ends.** $P(0)=h\,b_1B_1'(0)(1-\kappa)\neq0$, and $P(R_{\max})$, $Q(R_{\max})$ contain $S'(R_{\max})\neq0$: the boundary conditions hold only approximately (1e-12 of $\max|P|$ for hydrogen). The boundary term of the integration by parts in §5, $-c\,P_\alpha(0)Q_\beta(0)$, is therefore not 0 for the $B_1$ functions when $|\kappa|\ge2$ ($\propto1-\kappa^2$); H stays symmetric by construction, and by item 2 the eigenvalues do not see it.
4. **Spline count.** The knots of §2.1 give $N=M+k-1$ splines, not $M+k-2$ (corrected).
5. **`DKBSpline` data** (§10 said: pointer to the shared knot set). Each object owns a copy of its knots, $\kappa$, $c$ (the balance $1/2c$ is part of the function) and coefficients, so that Lua's garbage collector needs no reference counting. A knot set is ~140 doubles.
6. **Sign convention** (not in the design): $P>0$ on its innermost lobe, as FPLO's files and Grant. The innermost lobe is the first local extremum of $|P|$ above $\sqrt{\varepsilon}\max|P|$ — a rounding-noise floor, below which $P\propto r^{l+1}$ is noise for large $l$; not a physics tolerance.
7. **`"Function"` option form:** `{"Nucleus","Function"},{"Density",rho}` (`LuaCheckOptions` cannot take a string or a function under one name); §11.1 shows `{"Nucleus", rho}`. To settle for `CalculateAtomicBasisFunctions`. The density is sampled on the atomic quadrature points (the first ~40 knot intervals lie inside the nucleus with $r_1=R_{\rm rms}/50$), not on a separate grid; it is always rescaled to $Z$, with a warning above `DensityTolerance` (1e-10).
8. **Fermi nucleus for light atoms:** $c_F^2=\tfrac53(R_{\rm rms}^2-\tfrac75\pi^2a^2)$ is negative for small $A$ (e.g. $A=1$); `DiracRadialSolve` then raises an error. Visscher–Dyall use Gaussian nuclei there; a default rule is needed in A5.
9. **Defaults** of the temporary entry, not in the design: $r_{\min}=R_{\rm rms}/50$ (finite) or $10^{-6}/Z$ (point), $R_{\max}=40$; a finite nucleus needs `A` or `Rrms` (no default mass table).

**Found on the way, not fixed (no piecemeal guards)**
- `NewOperator("H0")` (`LuaOperatorNewAtomicH0`) still reads $Z$ as `int`.
- `lua_check(L, arg, BSpline::TBSpline*&)` (`luaInOut.cpp`) ignores `arg` and reads index −1, and casts without a type check; `NewOperator("H0")` given a `DKBSpline` (or a `CubicSpline`) would be misread. `NewOperator("Dirac")` given `DKBSpline`s runs through the $O(h^2)$ grid path on the quadrature points (the exact path is A5).
- `LuaOperatorNewDirac`/`LuaOperatorNewAtomicDirac` read $\kappa$ with `luaL_checkunsigned`, which works for negative $\kappa$ only through wrap-around.

**Remaining for A3 and later**
- Geometric DKB blocks are recomputed per call; cache them per knot set for the SCF (A3/A4), only the $V$, $V_\beta$ blocks change.
- Eigen-residuals are reported (`Info.Residual`); their size (≈1e-8 … 1e-4 absolute) reflects the same $\varepsilon\max|E|$ scale and is not yet a useful accuracy measure.

### Task A2: done (2026-11-11, commits 5017d0ec, 141893ac)

- **Confinement** (§6): `ConfinementPotential` $(r/r_c)^N$ in `DiracRadial.cpp`; `DiracRadialSolve` options `ConfinementRadius`, `ConfinementPower` (default 8), `ConfinementType` (`"Scalar"`, or `"Vector"` for comparison), `ConfinementVmax`; `Info.BetaVbeta` $=\langle\beta V_\beta\rangle$, `Info.NormQ` $=\int Q^2$. The split $r_c=P_{nl}P_{\rm at}\rho$ is bookkeeping for A5.
- **$R_{\max}$ follows the confinement** (decided by Maurits): $R_{\max}=r_c\,V_{\max}^{1/N}$ with option $V_{\max}$ (default $10^5$ Ha; $4.2\,r_c$ for $N=8$), an explicit `rMax` wins. Reason, measured: the eigenvalue rounding grows with $V_\beta(R_{\max})$. With $R_{\max}=40$, $r_c=0.5$: $\max|E|=2\times10^{15}$ Ha, the shift is not converged (1.271509 vs 1.271465 Ha for 140 and 200 intervals) and even the Rayleigh quotient scatters by 1.6e-8 Ha. With $4r_c$: $\max|E|=2\times10^{10}$, 140 and 200 intervals agree to 3e-8. $3r_c$ is already too small (−2.5e-4 Ha).
- **Check 3 targets replaced** by converged values (Rayleigh quotient, automatic $R_{\max}$, 140 intervals; 200 agree to 1e-9): §12. The old ones were the prototype's `dsygvd` values at $R_{\max}=40$; the two small ones were off in the third digit (2.2239e-7 against 2.2412e-7 = $\langle\beta V_\beta\rangle$). The C code agrees with the new values to ≤ 1.2e-10 relative (2.2e-7 Ha to its rounding floor $\varepsilon|E|$). With the automatic $R_{\max}$ the vector variant has an intruder for every $r_0$ (ΔE ≈ −1.5e4 Ha, $\int Q^2\approx1$), since $V_{\rm conf}$ exceeds $2c^2$ inside the box.
- **$Y^k$** (§7): temporary `RadialPoisson(rho, k, {opts})`, $Y^k$ as an exact `"BSpline"`. Check 4 in CI: hydrogen 1s $Y^0$ to 6.8e-9 / 7.1e-12 (80 / 160 intervals); $\rho=r^6e^{-2r}(1+0.3r^3)$, $k=0,1,2,4$, against incomplete gamma functions to 2e-9 of $\max|Y|$ (discretization at $r\approx8$, where the knots are 1 bohr apart).
- **On-site orthonormalization** (§6.1): `DiracOrthonormalize` (C, `DiracRadial.cpp`) and the temporary Lua `DiracOrthonormalize(G, F, {opts})` (S rebuilt from the `DKBSpline`s' knots, $\kappa$, $c$, nQ: exactly the solver's S). **Two passes, always** (decided by Maurits): one Cholesky or Löwdin pass leaves $\varepsilon\,{\rm cond}(O)$, 3e-10 for two 3d₅/₂ of $Z=25$, 26 (smallest eigenvalue of $O$ 2.8e-7); the second gives 2e-16 and keeps the Gram–Schmidt order and Löwdin's symmetry. §6.1 corrected.
- **Check 14 with stand-ins** (not decided explicitly; my default): different $Z$ for $Q^{\rm ion}$ and different $r_c$ for $P$, since $Q^{\rm ion}$ needs the SCF atom. Results: $\|C^TSC-1\|_{\max}$ = 2.2e-16 (Gram–Schmidt), 1.3e-15 (Löwdin); first orbital unchanged to 0; Löwdin independent of the input order to 1.2e-11 ($\varepsilon\,{\rm cond}(O)$). The real $Q^{\rm ion}$ test belongs to A3/A4.
- `luaL_error` knows no `%g`/`%e`; several A1 error messages were garbled. All DFT Lua errors now go through `DFTLuaError` (vsnprintf).
- Prototype bug, not fixed there: `dkb.py`'s `knots()` can give a last interior knot one ulp above `rmax` and so a decreasing knot vector (scipy refuses, e.g. `rmax = 40, nint = 200`). The C code sets $t_M=R_{\max}$ exactly.

### Order of work in stage A (one commit series per task, each closed by its checkpoints in §12)

| Task | Content | Checks |
|---|---|---|
| A0 | Remainder of step 0 (above) — done | step-0 tests plus the float check |
| A1 | Nucleus models (§3, without `Function3D`), DKB assembly from `BSplineMatrixBanded` (§5), eigenproblem, root selection; `DKBSpline` evaluation — done | 1a, 1b, 2, 13 (spherical part) |
| A2 | Poisson $Y^k$ in `BSplineMatrix` (§7), Lorentz-scalar confinement (§6), on-site orthonormalization (§6.1) — done | 3, 4, 14 |
| A3 | XC functionals and Dirac–Kohn–Sham atom (§9) | 9, 10 |
| A4 | Configuration-average DHF (§8) | 5, 6, 7, 8 |
| A5 | `CalculateAtomicBasisFunctions` binding, exact `"Dirac"`/`"AU"` paths (§11), `Function3D` (§3.1), Fourier–Bessel (§10) | 11, 12, 13 |

### Known limits that stay (not part of stage A)
- Martin's continuum code (`ContinuumStates`, `Eigendifferential`, reachable only as `Martin.*`) uses the compiled `cAU`; a change in a script does not reach it.
- Old code with `Rydberg = 13.6057` hard-coded (`Model/`, `Yi/`, `Robert/`, `MartinDresden/`) is untouched.
- There is no `BSplineProjection` yet. Stage A does not need one: orbitals come from eigenproblems, not from projection.
- quanty.org still shows CODATA 2014 and nothing yet on interpolating functions being 0 outside their grid.

## 1. Notation

Spinor:
$$\psi_{n\kappa\mu}(\mathbf r)=\frac1r\begin{pmatrix}P_{n\kappa}(r)\,\Omega_{\kappa\mu}(\hat r)\\ i\,Q_{n\kappa}(r)\,\Omega_{-\kappa\mu}(\hat r)\end{pmatrix},\qquad \int_0^\infty (P^2+Q^2)\,dr = 1 .$$
$l(\kappa)=\kappa$ for $\kappa>0$, $-\kappa-1$ for $\kappa<0$; $j=|\kappa|-\tfrac12$; $g_\kappa=2j+1=2|\kappa|$. Large component carries $l=l(\kappa)$, small component $\bar l=l(-\kappa)$.

## 2. B-spline module

### 2.1 Knots and basis
Order $k$ (polynomial degree $d=k-1$), $k=8$–$10$. Knot vector from `BSplineKnotsExponential(degree, nInterval, r1, rMax, …)`: $k$-fold knot at $0$, interior knots $t_j = r_1 (R_{\max}/r_1)^{(j-1)/(M-1)}$, $j=1..M$, $(k-1)$ further knots at $R_{\max}$ ($t_M=R_{\max}$ is the $k$-th). Number of splines $N = M + k - 1$ ($M$ knot intervals including $[0,r_1]$; corrected in A1, it said $M+k-2$). `BSplineKnotsExponential` implements exactly this (A1).
- $r_1$: well inside the nucleus, $r_1 \approx R_{\rm rms}/50$ for a finite nucleus.
- $R_{\max}$: with confinement $(r/r_c)^N$, where it reaches $V_{\max}$ (option, default $10^5$ Ha): $R_{\max}=r_cV_{\max}^{1/N}$, $4.2\,r_c$ for $N=8$ (A2). Much larger $R_{\max}$ makes $V_\beta(R_{\max})$, the largest eigenvalue and with it the rounding of all eigenvalues huge (§0.2, A2).
- Boundary conditions $P=Q=0$ at $r=0$ and $r=R_{\max}$: remove $B_1$ (the only spline nonzero at $0$) and $B_N$ (the only spline nonzero at $R_{\max}$). Working set $\{B_i\}_{i=2}^{N-1}$, size $n=N-2$.

### 2.2 Evaluation
`BSplineFindInterval` plus `BSplineValues(t, interval, knot, degree, nDeriv=2, val)` give $B_i$, $B_i'$, $B_i''$ of the $d+1$ nonzero splines at each point. ($B''$ is needed only in the DKB small components, §5.) Precompute these once per knot set at all quadrature points: $(3)\times(d+1)$ values per point plus the index of the first nonzero spline, as in `BSplineCellValues`. Contracted functions are evaluated through these tables, or for arbitrary $t$ through `BSplineValues`. They are 0 outside $[0,R_{\max}]$, following the step-0 convention.

### 2.3 Quadrature
`BSplineQuadratureSetup(knot, nKnot, nQ, …)` with $n_Q = k+4$ points per non-empty knot interval (High precision: $k+8$; never more than 127, §0.2). This is exact for $\langle B_i|r^m|B_j\rangle$ polynomial integrands and converges spectrally for $1/r^n$, $V$, $V_\beta$ and Coulomb terms. The nodes are those of the 34-digit table in `Core/BasicMath/Interpolation/GaussLegendre`. The working splines vanish at $r=0$ as $r^{\ge1}$, but the DKB products of the first one contain terms like $B_1^2/r^3$ that are logarithmically divergent for $|\kappa|\ge2$; the Gauss sum stays finite and the eigenvalues are insensitive to it (§0.2, A1 deviation 2).

## 3. Nucleus

Use a finite nucleus for production. A point nucleus with $\kappa>0$ and large $Z$ produces a spurious intruder state in the DKB/B-spline spectrum (§12, check 1b). The point nucleus remains available for tests against analytic hydrogen-like results.

Option `Nucleus` (default `"Gaussian"`):

| Value | $\rho_N$ | $V_N$ |
|---|---|---|
| `"Point"` | $Z\delta(\mathbf r)$ | $-Z/r$ |
| `"Gaussian"` | $Z(\zeta/\pi)^{3/2}e^{-\zeta r^2}$, $\zeta=3/(2R_{\rm rms}^2)$ | $-Z\,\mathrm{erf}(\sqrt\zeta r)/r$ (analytic) |
| `"Fermi"` | $\rho_0/(1+e^{(r-c_F)/a})$, $a=t/(4\ln3)$, $t=2.30$ fm; $c_F$ from $R_{\rm rms}^2=\tfrac35c_F^2+\tfrac75\pi^2a^2$; $\rho_0$ from $\int\rho_N=Z$ | radial Poisson (§7) |
| `"Function"` | Lua function `rho(r)` | radial Poisson (§7) |
| `"Function3D"` | Lua function `rho(x,y,z)` | multipole expansion, see below |

- **Nuclear radius:** default $R_{\rm rms}=(0.836\,A^{1/3}+0.570)$ fm (Visscher & Dyall). Options `A` (mass number) or `Rrms` (explicit, in bohr). $1\ \text{bohr}=52917.7210544$ fm (CODATA 2022).
- **Lua callbacks:** C never calls back into Lua. The Lua binding samples the user function before entering the C code: `rho(r)` on a dedicated fine radial grid (B-spline quadrature points of the nuclear region plus the atomic grid). The binding checks the normalization $4\pi\int r^2\rho_N=Z$, warns on a deviation larger than $10^{-10}$, and rescales. $V_N$ then follows from §7 with $k=0$: $V_N(r)=-4\pi\,Y^0[r^2\rho_N](r)/r$.

### 3.1 Non-spherical nuclei: `"Function3D"`

This is possible in the present design. Project the density onto real spherical harmonics, for even $L\le L_{\max}$ (default $L_{\max}=4$):
$$\rho_{LM}(r)=\int d\Omega\,\rho_N(r,\Omega)\,S_{LM}(\Omega).$$
The Lua binding samples $\rho_N(x,y,z)$ on a product grid: Gauss–Legendre in $\cos\theta$, uniform in $\phi$, with $n_\theta\ge L_{\max}+8$ and $n_\phi\ge 2L_{\max}+16$, times the radial nuclear grid. The electron potential energy of each multipole follows from §7:
$$V_{LM}(r)=-\frac{4\pi}{2L+1}\,\frac{Y^L[r^2\rho_{LM}](r)}{r},\qquad V_N(\mathbf r)=\sum_{LM}V_{LM}(r)\,S_{LM}(\hat r).$$
For $L=0$ this reduces to the spherical case.

How it is used:
- **$L=0$** enters the atomic SCF.
- **$L\ge2$** is exported as a one-particle Quanty Operator, e.g. the electric quadrupole hyperfine interaction including finite-size and penetration effects. The radial integrals are $\int V_{LM}(PP'+QQ')\,dr$. The angular factor is the standard jj matrix element of the rank-$L$ tensor, identical for large and small components; the parity condition $l+l'+L$ even holds for both.
- The energies involved are many orders of magnitude below the SCF tolerance, so first-order treatment is exact for all practical purposes.
- A self-consistent non-spherical treatment would need $LM$-resolved electron densities. That machinery arrives in stage D and can then take $V_{LM}$ directly.

Magnetic hyperfine interaction needs the nuclear magnetization distribution, not the charge distribution, and enters as $-\boldsymbol\alpha\cdot\mathbf A$. It is left for later, but the `Nucleus` object should reserve fields for the magnetic moment and its distribution (point, uniform sphere, or Lua function).

## 4. Radial Dirac equation with vector $V$ and Lorentz-scalar $V_\beta$

Hamiltonian $H = c\,\boldsymbol\alpha\cdot\mathbf p + \beta\,(c^2+V_\beta) + V - c^2$. Radial form:
$$\begin{pmatrix} V+V_\beta & c\big(-\frac{d}{dr}+\frac{\kappa}{r}\big)\\ c\big(\frac{d}{dr}+\frac{\kappa}{r}\big) & V-V_\beta-2c^2\end{pmatrix}\begin{pmatrix}P\\Q\end{pmatrix}=E\begin{pmatrix}P\\Q\end{pmatrix}.$$

$V_\beta$ acts as a position-dependent mass $c^2+V_\beta$. The local decay constant $\sqrt{(c^2+V_\beta)^2-(E+c^2-V)^2}/c$ stays real and grows with $V_\beta$, so both components are confined and there is no Klein tunnelling.

A vector confinement ($V\to V+V_{\rm conf}$) instead pushes negative-continuum states into the bound-state region once $V_{\rm conf}\gtrsim 2c^2$ inside the box. This was observed in the prototype (§12, check 3).

## 5. Dual-kinetic-balance basis and matrix elements

DKB (Shabaev et al. 2004). For each working spline $B_i$:
$$u_i^{(1)}=\begin{pmatrix}B_i\\ \frac{1}{2c}\big(B_i'+\frac{\kappa}{r}B_i\big)\end{pmatrix},\qquad u_i^{(2)}=\begin{pmatrix}\frac{1}{2c}\big(B_i'-\frac{\kappa}{r}B_i\big)\\ B_i\end{pmatrix}.$$
This gives $2n$ functions per $\kappa$, collected as $u_\alpha=(P_\alpha,Q_\alpha)$. The balance uses the constant mass $c$, not $c^2+V_\beta$. That is harmless because $V_\beta\approx0$ wherever the wavefunction is non-negligible, and completeness is unaffected.

Matrix elements: integrate the off-diagonal derivative terms by parts (boundary terms vanish). This makes $H$ manifestly symmetric and needs only first derivatives of $P$ and $Q$:
$$H_{\alpha\beta}=\int_0^{R_{\max}}\!\Big[P_\alpha(V+V_\beta)P_\beta+Q_\alpha(V-V_\beta-2c^2)Q_\beta + c\,(P_\alpha'Q_\beta+Q_\alpha P_\beta') + \frac{c\kappa}{r}(P_\alpha Q_\beta+Q_\alpha P_\beta)\Big]dr,$$
$$S_{\alpha\beta}=\int_0^{R_{\max}}(P_\alpha P_\beta+Q_\alpha Q_\beta)\,dr .$$
The required derivatives are $Q^{(1)\prime}=\frac1{2c}(B''+\kappa B'/r-\kappa B/r^2)$ and $P^{(2)\prime}=\frac1{2c}(B''-\kappa B'/r+\kappa B/r^2)$.

Generalized symmetric eigenproblem $Hc=ESc$ (LAPACK `dsygvd`), dimension $2n$; the returned eigenvalues are the Rayleigh quotients $c^THc/c^TSc$ of the eigenvectors (default; §0.2, A1 deviation 1).
- **Spectrum:** $n$ eigenvalues below $-2c^2$ (negative continuum) and $n$ above. Electronic states are those with $E>-c^2$.
- **Root selection** for subshell $(n,\kappa)$: index $n-l(\kappa)-1$ within the electronic branch. In self-consistent open-shell iterations, use maximum overlap with the previous iterate instead.

Implementation with the step-0 core:
- **Building blocks.** All blocks of $H$ and $S$ are linear combinations of
$$M^{(a,b)}[f]_{ij}=\int B_i^{(a)}\,f\,B_j^{(b)}\,dr,\qquad a,b\in\{0,1,2\},\quad f\in\{r^{-m},\ V r^{-m},\ V_\beta r^{-m}\},\ m=0\ldots3,$$
  each from one `BSplineMatrixBanded` call with `first = 1` and `n` = $N-2$, where `fq` holds $f$ at the quadrature points.
  - The powers of $1/r$ come from expanding the DKB functions (§5) in $B$, $B'$, $B''$ and $\kappa B/r$. For example, $\frac{c\kappa}{r}P^{(2)}Q^{(1)}$ contains $\kappa^3B_iB_j/r^3$, and $Q^{(1)}VQ^{(1)}$ contains $\kappa^2B_iVB_j/r^2$. Write out the explicit list of $(a,b,f)$ combinations in A1, and check it against the prototype `dkb.py`, which evaluates $P,Q$ directly.
  - The $\kappa$-independent geometric blocks ($f=r^{-m}$) are computed once per knot set. The $V$ and $V_\beta$ blocks are recomputed whenever the potential changes (each SCF iteration). Per $\kappa$, only the coefficients ($\kappa$, $\kappa^2$, $\kappa^3$) of the combination change.
- **Non-symmetric blocks.** `BSplineMatrixBanded` stores only $i\ge j$, which is complete for $a=b$. For $a\ne b$ the full matrix needs both triangles. They follow from two calls, using $M^{(a,b)}_{ij}=M^{(b,a)}_{ji}$: the lower triangle of $(a,b)$ and the lower triangle of $(b,a)$, transposed, which also supplies the upper triangle of $(a,b)$.
- **Eigenproblem.** Assemble the $2n\times2n$ matrices densely for `dsygvd`; with $n\approx150$ this is cheap. A banded-block solver is an optimization for later.

## 6. Basis-generating potential and parameters

For basis generation (FPLO-like, but with a Lorentz-scalar confinement):
$$V_{nl}(r) = V^{\rm at}[Q^{\rm ion}_{nl}](r),\qquad V_{\beta,nl}(r)=\Big(\frac{r}{P_{nl}\,P_{\rm at}\,\rho}\Big)^{N}.$$
- **$V^{\rm at}[Q^{\rm ion}]$:** self-consistent atomic potential (§8 or §9), optionally recomputed with charge $Q^{\rm ion}$ removed. This generates more compact polarization orbitals. (The FPLO symbol $Q$ is renamed $Q^{\rm ion}$ to avoid a clash with the small component.)
- **$\rho$:** compression radius. In stage A it is a free parameter; in stage B it is made environment-dependent.
- **$P_{\rm at}$, $P_{nl}$:** atom-wide and orbital-specific scale factors.
- **$N$:** power, of order 8–14.
- Both $j$ sub-shells of an $nl$ shell use the same $(Q^{\rm ion},P)$.
- **All orbitals of one atom share one knot set**, whatever their confinement. This makes the on-site overlaps of §6.1 exact matrix operations in coefficient space.

Energy bookkeeping:
- The confinement contributes $\langle\beta V_\beta\rangle=\int V_\beta\,(P^2-Q^2)\,dr$. This is the scalar density, not $P^2+Q^2$.
- Report the confined-atom total energy and $\langle\beta V_\beta\rangle$ separately.
- The solid/molecule Hamiltonian of stages D–E contains no $V_\beta$.

### 6.1 On-site orthonormalization

Orbitals generated in different potentials ($Q^{\rm ion}_{nl}$, $P_{nl}$ differ per orbital) are not mutually orthogonal. All local orbitals of an atom are made orthonormal as follows.

- **Different $\kappa$:** orthogonal by the spin-angular parts, for large and small components alike. Nothing to do.
- **Equal $\kappa$:** all orbitals live in the same DKB space (same knot set, same $\kappa$). Collect their coefficient vectors as the columns of $C_\kappa$ ($2n\times m_\kappa$) and form the overlap $O_\kappa=C_\kappa^{T}S_\kappa C_\kappa$, where $S_\kappa$ is the DKB overlap matrix of §5.
- **Default: ordered Gram–Schmidt**, implemented as a Cholesky factorization $O_\kappa=LL^T$, then $C_\kappa\leftarrow C_\kappa L^{-T}$.
  - Order: core, then semicore, then the chemical valence orbitals (first member of each multi-orbital), then the remaining multi-orbital members and polarization orbitals in increasing $n$. This is FPLO's ordering.
  - Gram–Schmidt leaves the first orbital of each $\kappa$ unchanged and deforms the others minimally with respect to those before them. The atomic ground-state orbitals therefore keep their shape.
  - The order can be overridden per atom.
- **Option: Löwdin** ($C\leftarrow CO^{-1/2}$) within user-defined groups of equal rank.
- **Checks:**
  - Stop with an error if an eigenvalue of $O_\kappa$ falls below $10^{-8}$; the basis is then redundant and must be pruned.
  - After orthonormalization verify $\|C^TS_\kappa C-1\|_{\max}<10^{-13}$. The operation is exact in coefficient space (unlike FPLO, no orthogonalization on a grid), but in floating point one pass leaves $\varepsilon\,{\rm cond}(O)$ (3e-10 for similar orbitals); it is therefore done **twice**, which gives 1e-16 whenever ${\rm cond}(O)<1/\sqrt\varepsilon$ — guaranteed by the redundancy check (A2).
- **Consequence:** the orthonormalized orbitals are no longer eigenfunctions of their generating Hamiltonians, so $h_D$ is a full matrix within each $\kappa$ block. `NewOperator("Dirac")` already builds the full $\kappa$ block (`CalculateDiracMatrixElements`).
- **Inter-site overlaps** (stages B–E) remain non-zero and enter through $S(\mathbf R)$.

## 7. Radial Poisson equation (Hartree and exchange potentials)

For a radial pair density $\rho_{cd}(r)=P_cP_d+Q_cQ_d$:
$$Y^k_{cd}(r) = r\int_0^\infty \rho_{cd}(r')\,\frac{r_<^k}{r_>^{k+1}}\,dr',\qquad R^k(ab;cd)=\int_0^\infty \rho_{ac}(r)\,\frac{Y^k_{bd}(r)}{r}\,dr .$$
$Y^k$ satisfies $Y''-\frac{k(k+1)}{r^2}Y=-\frac{(2k+1)}{r}\rho$ with $Y(0)=0$ and $Y\propto r^{-k}$ outside the density.

Galerkin solution in the same B-spline set: keep $B_N$ and remove only $B_1$. With $Y=\sum_j a_jB_j$, and the sign chosen so that the matrix is symmetric positive definite:
$$\sum_j\Big[\int B_i'B_j' + k(k+1)\!\int\frac{B_iB_j}{r^2} + \frac{k}{R_{\max}}B_i(R_{\max})B_j(R_{\max})\Big]a_j = (2k+1)\int\frac{B_i\,\rho}{r}.$$
The last matrix term is the Robin condition $Y'(R_{\max})=-kY(R_{\max})/R_{\max}$, exact when $\rho=0$ beyond $R_{\max}$; for $k=0$ it gives $Y'=0$.
- **Assembly:** two `BSplineMatrixBanded` calls, $(a,b)=(1,1)$ with $f=1$ and $(0,0)$ with $f=1/r^2$, with `first = 1` and `n` = $N-1$, i.e. keeping $B_N$. Add the Robin term to the last diagonal element, since only $B_N(R_{\max})=1$ is nonzero there.
- The matrix is banded (bandwidth $d$) and SPD. Factorize once per $k$ (`dpbtrf`) and reuse for all pair densities (`dpbtrs`).
- The right-hand side $\int B_i\rho/r$ is a quadrature sum with the values of $\rho$ at the quadrature points.
- Evaluate $Y^k/r$ at quadrature points; that is all the potentials need.
- Accuracy is verified in §12 (check 4).

## 8. Configuration-average Dirac–Hartree–Fock

### 8.1 Bookkeeping
Nonrelativistic shells $s=(n,l)$ with occupation $q_s$ and degeneracy $G_s=4l+2$. Relativistic subshells $a\in s$, $\kappa_a\in\{l, -l-1\}$ (one subshell for $l=0$), degeneracy $g_a$.

The average over all determinants of the configuration, which is basis-invariant and therefore equal to Cowan's average, gives:
- **mean occupation** $\bar q_a = q_s\,g_a/G_s$;
- **number of electron pairs** (unordered):
$$N_{aa}=\frac{q_s(q_s-1)}{2}\,\frac{g_a(g_a-1)}{G_s(G_s-1)},\quad N_{ab}\big|_{a\ne b\in s}=q_s(q_s-1)\frac{g_ag_b}{G_s(G_s-1)},\quad N_{ab}\big|_{a\in s,b\in t\ne s}=q_sq_t\frac{g_ag_b}{G_sG_t}.$$

**Do not** use fixed fractional subshell occupations $q_a=\bar q_a$ in the pair counts. Pair counts $\bar q_a(\bar q_a-1)/2$ do not reproduce Cowan. For $d^3$ the error is $\tfrac{44}{1875}F^2-\tfrac{8}{375}F^4$ (§12, check 5).

### 8.2 Angular coefficients
$$\Gamma^k_{ab}=\begin{pmatrix} j_a & k & j_b\\ \tfrac12 & 0 & -\tfrac12\end{pmatrix}^2\ \text{if } l_a+l_b+k \text{ even, else } 0 .$$
The parity condition is the same for large and small components, since $\bar l_a+\bar l_b+k$ has the same parity.

### 8.3 Energy functional
$$E_{av}=\sum_a \bar q_a\, I_a + \sum_a N_{aa}\Big[F^0(aa)-\frac{g_a}{g_a-1}\sum_{k>0}\Gamma^k_{aa}F^k(aa)\Big] + \sum_{a<b}N_{ab}\Big[F^0(ab)-\sum_k \Gamma^k_{ab}\,G^k(ab)\Big],$$
with $I_a=\langle a|h_D|a\rangle$ ($h_D$ = §4 with $V=V_N$, plus $\beta V_\beta$ if confined), $F^k(ab)=R^k(ab;ab)$ and $G^k(ab)=R^k(ab;ba)$. The $aa$ term is absent for $g_a=1$.

Derivation of the $aa$ bracket: average over pairs of distinct spinors $m\neq m'$ in one subshell. The direct sum gives $g_a^2F^0$; the exchange sum gives $g_a^2\sum_k\Gamma^k_{aa}F^k$ with $\Gamma^0_{aa}=1/g_a$; the $m=m'$ terms cancel. Dividing by $g_a(g_a-1)$ gives the bracket.

### 8.4 Fock operators (orbital-specific)
Varying $E_{av}$ with respect to $(P_a,Q_a)$ and dividing by $2\bar q_a$:
$$f_a\psi_a = h_D\psi_a + \frac{1}{r}\Big\{\frac{2N_{aa}}{\bar q_a}\Big[Y^0_{aa}-\frac{g_a}{g_a-1}\sum_{k>0}\Gamma^k_{aa}Y^k_{aa}\Big] + \sum_{b\ne a}\frac{N_{ab}}{\bar q_a}Y^0_{bb}\Big\}\psi_a \;-\; \sum_{b\ne a}\frac{N_{ab}}{\bar q_a}\sum_k\Gamma^k_{ab}\,\frac{Y^k_{ab}}{r}\,\psi_b .$$
The potentials multiply both components ($\psi=(P,Q)$).

Check: a closed subshell alone ($\bar q_a=g_a$, $N_{aa}=g_a(g_a-1)/2$) gives $(g_a-1)Y^0_{aa}-g_a\sum_{k>0}\Gamma^k Y^k_{aa}$, i.e. the standard self-interaction-free closed-shell form.

In the DKB basis:
$$(F_a)_{\alpha\beta}=H_{\alpha\beta}+\int (P_\alpha P_\beta+Q_\alpha Q_\beta)\,v_a\,dr,$$
where $v_a$ is the local part (the curly bracket divided by $r$). The exchange part enters as the vector $x_{a,\alpha}=-\sum_b\frac{N_{ab}}{\bar q_a}\sum_k\Gamma^k_{ab}\int (P_\alpha P_b+Q_\alpha Q_b)\,Y^k_{ab}/r\,dr$, so it does not have to be stored as a matrix.

The energy gradient with respect to the coefficient vector $c_a$ of the normalized orbital is
$$\partial E_{av}/\partial c_a = 2\bar q_a\,(F_ac_a + x_a - \varepsilon_a S c_a).$$
This was verified against finite differences (§12, check 6).

### 8.5 SCF procedure
1. **Initial guess:** a self-consistent LDA atom (§9) with the same configuration. This is robust and also gives the guess for $V^{\rm at}$.
2. **Each iteration:** compute all $Y^k_{bb}$ and $Y^k_{ab}$ (§7). For each subshell, build $F_a$ and the exchange vector $x_a$ from the current orbitals. Then either (i) diagonalize $F_a$ with the exchange approximated by a local-equivalent potential (the exchange term divided by the previous orbital, per component, regularized at nodes, as in RCN), or (ii) solve the inhomogeneous linear system $(F_a-\varepsilon_aS)c_a=-x_a$ with $\varepsilon_a$ updated from normalization; see open point 3. Pick the root by maximum overlap with the previous iterate, as in Froese Fischer's method and Cowan's RCN. Orthogonalize orbitals of equal $\kappa$ (Löwdin, or Schmidt starting from the most strongly bound).
3. **Off-diagonal Lagrange multipliers / Brillouin condition:** stationarity under rotation of two occupied orbitals $a,b$ of equal $\kappa$ requires
$$\bar q_a\langle b|f_a|a\rangle=\bar q_b\langle a|f_b|b\rangle .$$
   It is automatically satisfied (redundant rotation) when both subshells are closed. Otherwise enforce it with a rotation step per pair, or use DIIS on the projected orbital gradient of §8.4, which contains it.
4. **Convergence:** max |orbital gradient| $<10^{-8}$ and $\Delta E<10^{-10}$ Ha. Acceleration: DIIS over the stacked $(c_a)$ with error vectors equal to the projected gradients.

### 8.6 Nonrelativistic limit
With $j$-independent radial functions ($R^k$ equal for both $j$ subshells), $E_{av}$ of §8.3 equals Cowan's
$$E_{av}=\sum_s q_sI_s+\sum_s\tfrac{q_s(q_s-1)}{2}\Big[F^0-\tfrac{2l+1}{4l+1}\sum_{k>0}\begin{pmatrix}l&k&l\\0&0&0\end{pmatrix}^2F^k\Big]+\sum_{s<t}q_sq_t\Big[F^0-\tfrac12\sum_k\begin{pmatrix}l_s&k&l_t\\0&0&0\end{pmatrix}^2G^k\Big]$$
identically. This was verified symbolically for $p,d,f$ intra-shell and $sp,pd,df,sd$ inter-shell terms (§12, check 5).

## 9. Dirac–Kohn–Sham atom (LDA, PBE)

Spherical density from mean occupations (only $\bar q_a$ enters; the pair counts of §8.1 are HF-specific):
$$n(r)=\frac{1}{4\pi r^2}\sum_a\bar q_a\,(P_a^2+Q_a^2).$$
- **Hartree:** $V_H(r)=Y^0[\,4\pi r^2 n\,](r)/r$ from §7.
- **LDA:** $v_{xc}=\partial(n\,\epsilon_{xc})/\partial n$, with Slater exchange plus PW92 correlation (production) or VWN5 (NIST validation). The MacDonald–Vosko factor applies to exchange only and is a validation option.
- **PBE, in weak form** to avoid second derivatives of $n$ in the potential. With $\sigma=|\nabla n|^2=n'^2$ and $t_{\alpha\beta}=(P_\alpha P_\beta+Q_\alpha Q_\beta)/(4\pi r^2)$:
$$\langle u_\alpha|v_{xc}|u_\beta\rangle=\int_0^{R_{\max}}4\pi r^2\Big[\frac{\partial e_{xc}}{\partial n}\,t_{\alpha\beta}+2\frac{\partial e_{xc}}{\partial\sigma}\,n'\,t_{\alpha\beta}'\Big]dr.$$
  Here $n'$ and $t'$ are analytic from the B-spline derivatives; $e_{xc}$ is the energy density per volume.
- **Total energy:** $E=\sum_a\bar q_a\langle a|h_D|a\rangle+\tfrac12\int 4\pi r^2 nV_H+E_{xc}[n]$. Equivalently, the eigenvalue sum minus $\tfrac12\int nV_H-\int nv_{xc}+E_{xc}$, which serves as a consistency check.
- **Mixing:** Pulay/Anderson mixing on $V_H+v_{xc}$ at the quadrature points.

## 10. Output of stage A

1. **Contracted radial spinors:** `{kappa, n, coefficient vector c in R^{2n} over u^{(1)},u^{(2)}, knot set}`, orthonormalized per §6.1.
   - $P$ and $Q$ are not pure B-spline expansions: they contain $B'$ and $B/r$. The object must therefore always be evaluated through the DKB definition; that is exact and needs no refitting.
   - To Lua they are handed out as two `InterpolatingFunction_Type` userdata per orbital (large part `G` = $P$, small part `F` = $Q$) of a new subclass with `type = "DKBSpline"`. Each holds a pointer to the shared knot set, $\kappa$, its component flag, and the full coefficient vector $c$. Its `operator()(t,k)` evaluates exactly, including derivatives.
   - Existing consumers that sample on `GetX()` therefore keep working. `GetX()` returns the quadrature points.
   - The `DKBSpline` class follows the `TBSpline` conventions settled in step 0:
     - it is 0 outside the knot range, value and all derivatives;
     - `Integrate` is exact (per-interval Gauss–Legendre);
     - arithmetic with other `InterpolatingFunction` types raises an explicit error until it is needed, since the sum of two DKB functions with different $\kappa$ is not a DKB function.
2. **Atomic quantities:**
   - orbital energies $\varepsilon_{n\kappa}$, total energy, $\langle\beta V_\beta\rangle$, $\langle r^k\rangle$;
   - $R^k(ab;cd)$, $F^k$ and $G^k$ jj-resolved;
   - effective spin-orbit $\zeta_{nl}=2(\varepsilon_{n,l+1/2}-\varepsilon_{n,l-1/2})/(2l+1)$;
   - $j$-averaged radial quantities for LS-type Quanty Hamiltonians (weights $g_a/G_s$).
3. **Fourier–Bessel transforms**, prepared for stage C (normalization preserving $\int k^2|\tilde\phi|^2dk=\int P^2dr$):
$$\tilde P(k)=\sqrt{\tfrac{2}{\pi}}\int_0^\infty r\,P(r)\,j_{l}(kr)\,dr,\qquad \tilde Q(k)=\sqrt{\tfrac{2}{\pi}}\int_0^\infty r\,Q(r)\,j_{\bar l}(kr)\,dr,$$
   computed by quadrature per knot interval, with $n_q\ge k+k_{\max}\Delta r/\pi$ points per interval, on a uniform $k$ grid. $j_l$ is evaluated by downward recurrence for $kr<l$.

## 11. Code organisation and Lua interface

Based on branch `Claude` as of 2026-11-11, after step 0 (`LuaInterface/main.cpp` includes everything; `AtomicOperators/` holds the atomic Lua operators; `Core/BasicMath/Interpolation/` holds the B-spline, quadrature and interpolating-function code). Quanty conventions apply: single compilation unit, unprefixed and non-`static` helpers, return 0/1, `malloc`/`free`, C++11, Lua as the interface only. Core routines do not print.

| File | Content |
|---|---|
| `Core/BasicMath/Interpolation/BSpline.cpp`, `BSplineMatrix.cpp`, `GaussLegendre.cpp` (exist, step 0) | knots, values and derivatives, quadrature, banded matrices; used as is |
| `Core/BasicMath/Interpolation/BSplineMatrix.cpp` (extend) | banded Galerkin Poisson $Y^k$ (§7): assembly from `BSplineMatrixBanded`, Robin term, `dpbtrf`/`dpbtrs`. Suggested name: `BSplinePoissonSetup`/`BSplinePoissonSolve`. It belongs next to the other matrix builders because it is generic. |
| `Core/BasicMath/Interpolation/DKBSpline.cpp` (skeleton in A0, filled in A1) | `InterpolatingFunction` subclass `"DKBSpline"` (§10) |
| `DFT/DiracRadial.cpp` (new directory, included from `main.cpp` like `Tight_Binding/`) | nucleus models (§3), DKB assembly, generalized eigenproblem, root selection, on-site orthonormalization (§6.1), Fourier–Bessel transform |
| `DFT/XCFunctional.cpp` | Slater, PW92, VWN5, PBE ($e$, $\partial e/\partial n$, $\partial e/\partial\sigma$), MacDonald–Vosko factor |
| `DFT/AtomSCF.cpp` | configuration parsing, $\bar q_a$ and $N_{ab}$, LDA/PBE SCF, configuration-average DHF SCF with DIIS, Slater integrals |
| `DFT/LuaAtom.cpp` | Lua binding `CalculateAtomicBasisFunctions`; sampling of Lua nucleus functions (§3) |
| `AtomicOperators/AtomicDirac.cpp`, `AtomicU.cpp` | exact fast paths for `"DKBSpline"` (below) |

### 11.1 Lua interface: drop-in replacement for `ReadFPLOBasisFunctions`

The existing workflow, from `Debug examples/AtomicOperators/AtomicOperators.Quanty`:
```lua
orbitalLabels = {"1s1/2","2s1/2","2p1/2","2p3/2"}
indices, NF = CreateAtomicIndicesList(orbitalLabels,{{"Kappas",true}})
G, F = ReadFPLOBasisFunctions(orbitalLabels,"+fval.001.1","+fval.001.2",{{"Dirac",true}})
H0 = NewOperator("Dirac",NF,Z,indices,G,F)
U  = NewOperator("AU",NF,indices,G,F)
```
The proposed workflow changes only the third line:
```lua
orbitalLabels = {"1s1/2","2s1/2","2p1/2","2p3/2","3s1/2","3p1/2","3p3/2","3d3/2","3d5/2","4s1/2"}
indices, NF = CreateAtomicIndicesList(orbitalLabels,{{"Kappas",true}})
G, F, x, Atom = CalculateAtomicBasisFunctions(orbitalLabels, {
    {"Z", 28}, {"A", 58}, {"Nucleus", "Gaussian"},          -- or "Point","Fermi", a function rho(r), {"Function3D", rho(x,y,z), Lmax}
    {"Configuration", "1s2 2s2 2p6 3s2 3p6 3d8 4s2"},
    {"Method", "DHF"},                                       -- "DHF" (configuration average), "LDA", "PBE", "LDA-VWN", "RLDA-MV"
    {"Degree", 8}, {"NIntervals", 120}, {"rMin", 1e-6}, {"rMax", 40},
    {"Confinement", {rho = 6.0, P = 1.0, N = 10,
                     Orbitals = {["4s"] = {P = 0.9, Qion = 0}, ["4p"] = {P = 0.85, Qion = 2}}}},
    {"Orthonormalize", "GramSchmidt"},                       -- or "Lowdin"; order = orbitalLabels unless {"OrthoOrder",{...}}
    {"Precision", "High"} })                                 -- or individual tolerances, §0.1
H0 = NewOperator("Dirac",NF,Z,indices,G,F,{{"Nucleus",Atom.Nucleus}})
U  = NewOperator("AU",NF,indices,G,F)
```
- Nonrelativistic limit: set `PhysicalConstants.cAU.value = 1e8` before the call. This affects the atom solver and `NewOperator("Dirac")` identically, because both read the same entry.
- Options use the existing Quanty `{{"Name",value},…}` list form (`LuaTestOptionsList` / `LuaCheckOptions`).
- **Orbital labels:** labels in `orbitalLabels` that are not in the configuration are generated as empty or polarization orbitals in the SCF potential, with their own confinement. Labels absent from `orbitalLabels` are still part of the SCF but are not returned.
- **Fourth return value `Atom`:** a Lua table with total energy, $\langle\beta V_\beta\rangle$, orbital energies per label, $\bar q_a$, `Nucleus` (model and parameters, for reuse in `NewOperator("Dirac")`), Slater integrals (jj-resolved and $j$-averaged), $\zeta_{nl}$, and the potentials as `InterpolatingFunction`s. For `"Function3D"` it also contains the $L\ge2$ nuclear multipole potentials, ready for a one-particle Operator.

### 11.2 Required changes to existing operators

- **`NewOperator("Dirac",…)`** (`LuaOperatorNewAtomicDirac` → `CreateOperatorDirac` → Martin's `DiracH`):
  - Done in step 0: $c$ is read from `PhysicalConstants.cAU.value` at every call.
  - Still a point nucleus ($-Z/r$) with `int Z` (the `double` conversion is task A0). Integrals use the $O(h^2)$ `CreateWeights` grid quadrature with dual functions on `G[0]->GetX()`.
  - Add the option `{"Nucleus", …}` with the same models as §3.
  - Add an exact path when all functions are `"DKBSpline"`: $\langle a|h_D|b\rangle=c_a^TH^{\rm DKB}_\kappa c_b$ with $H^{\rm DKB}$ from §5 (without $V_\beta$).
  - With the solver's nucleus and $c$, this must reproduce the solver's one-electron integrals to machine precision (checkpoint 12).
- **`NewOperator("AU",…)`** (`LuaOperatorNewAtomicU`): it dispatches on `InterpolatingFunction::type`.
  - `"BSpline"` goes to `LuaOperatorNewAtomicUBSpline`, which since step 0 uses `BSplineCellValues`, `BSplineSlaterCellTables` and `BSplineSlaterIntegral` (Martin's cell method, in Core).
  - Everything else goes to `LuaOperatorNewAtomicUGeneralFunction`, which now has 34-digit nodes but still converges as $O(h^2)$ because of `CreateWeights`.
  - Add a third branch for `"DKBSpline"` that computes $R^k(ab;cd)$ with the Galerkin Poisson solver of §7 on the shared knot set, then calls the existing `CreateOperatorAtomicURelativisticOneIntegral…` assembly unchanged. The cell method cannot be reused directly: DKB functions are not pure B-spline expansions ($B'$ and $B/r$ terms). In the nonrelativistic limit, however, it serves as an independent check of the Poisson route.
  - Until that branch exists (A0 gives an explicit error; see §0.2), `"DKBSpline"` objects must not silently fall through to the general-function path.
- **Conventions, verified in task A0** (2026-11-11):
  - (i) **`G` = $P=rg$, `F` = $Q=rf$**, normalized as $\int(P^2+Q^2)\,dr=1$.
    - Code: `DiracH` (`Martin/DiracHamiltonian.cpp`) integrates with the weights of `CreateWeights`, which approximate $\int\ldots dr$ (no $r^2$), and applies exactly the radial operator of §4 with $P$ = `G`, $Q$ = `F`: `H11G` $=-Z/r\,P$, `H12F` $=\kappa Q/r-Q'$, `H21G` $=\kappa P/r+P'$, `H22F` $=(-2c^2-Z/r)\,Q$. The relativistic Slater integrals (`Slater_Integrate`) use $P_aP_c+Q_aQ_c$ with $dr$ weights as well.
    - Data: the FPLO pair `Debug examples/AtomicOperators/+fval.001.{1,2}` (Be, 88 orbitals up to 10f) has $\int(G^2+F^2)\,dr=1.000000$ for every orbital ($\int r^2(\ldots)dr$ = 0.24 … 115), $G\propto r^{l+1}$ and $F\propto r^{\bar l+1}$ at the origin.
  - (ii) **The FPLO small component has this document's sign**: `F` $=Q$ with $\psi=\frac1r(P\,\Omega_{\kappa\mu},\,iQ\,\Omega_{-\kappa\mu})$.
    - Data: at the maximum of $|G|$, $F/[(G'+\kappa G/r)/2c]$ = 0.9998 (1s), 1.0000 (2s, 2p₁/₂, 2p₃/₂, 3d₃/₂, 3d₅/₂, 4f₅/₂, 4f₇/₂); with the opposite sign it would be −1. $F/G<0$ for $\kappa<0$ and $>0$ for $\kappa>0$ there.
    - Code: `DiracH` has $+\kappa/r$ and $+d/dr$ in the lower-left element, the same operator as §4, so FPLO functions and this document's $Q$ enter it identically.
  - (iii) **$\mu$ order: ascending, $\mu=-j\ldots j$**, per orbital; orbitals in the order given to `CreateAtomicIndicesList` (each list holds $2j+1$ consecutive indices).
    - The order is fixed by the consumers, not by `CreateAtomicIndicesList`: `CreateOperatorCoulombFourRelativisticShells` (`Model/StandardOperators.cpp`) reads `Jindex[(2j+2μ)/2]`, i.e. position $j+\mu$; `CreateOperatorDirac` is diagonal in $\mu$.
    - The transformation in use is **`YToJJzMatrix`** (`Tight_Binding_lib.cpp`, Lua `YtojjzMatrix`), not `SphericalHarmonicKTojjz`: the latter has **no callers** and is a fixed table for $l\le3$. `YToJJzMatrix` is general in $l$: rows $j=l-\tfrac12$ ($\mu$ ascending) then $j=l+\tfrac12$ ($\mu$ ascending); columns $(m_l,\downarrow),(m_l,\uparrow)$ for $m_l=-l\ldots l$; Condon–Shortley Clebsch–Gordan coefficients, $|l-\tfrac12,\mu\rangle=-\sqrt{\tfrac{l-\mu+1/2}{2l+1}}|\mu-\tfrac12,\uparrow\rangle+\sqrt{\tfrac{l+\mu+1/2}{2l+1}}|\mu+\tfrac12,\downarrow\rangle$. This is the standard $\Omega_{\kappa\mu}$ (Grant), consistent with (ii). `YToKMatrix` (tesseral) has tables up to $l=4$ and a general fallback.
    - **$l=4$ works**: `CreateAtomicIndicesList({"5g7/2","5g9/2"},{{"kappas",true}})` gives $\kappa=4,-5$ with 8 and 10 indices.
    - **$l=5,6$ do not**: the label parser `relAngularMomentumZNumbers` reads one digit before `/2`, so `6h11/2` gives $\kappa=-1$ with 2 indices, `7i11/2` and `7i13/2` give $\kappa=-1,-2$ — silently wrong. Plan for $l\le6$: parse the whole integer $2j$ before `/2` in `relAngularMomentumZNumbers` (and `kappa()`); `NamesToKappa` in `AtomicU.cpp` knows only s–g and needs the same; `YToJJzMatrix` needs no change; check `YToKMatrix`'s general fallback against the Quanty tesseral phase convention for $l=5,6$ before use.

### 11.3 C data structures (sketch)
```c
typedef struct { uint32_t degree, nKnot, nSpline, nQuadPerInterval, nQuad;
                 double *knot, *rq, *wq;            /* quadrature points/weights            */
                 double *B, *dB, *d2B;              /* (degree+1) nonzero splines per point  */
                 uint32_t *firstSpline;             /* index of first nonzero spline per pt  */ } BSplineBasisType;
typedef struct { int model; double Z, A, Rrms, cFermi, aFermi;
                 uint32_t Lmax; double *VLM;        /* tabulated V_LM(r) for Function/Function3D */ } NucleusType;
typedef struct { uint32_t n; int kappa; BSplineBasisType *basis; double *c; /* 2*(nSpline-2) */
                 double energy, occupation, Qion, Pconf; } RadialSpinorType;
typedef struct { NucleusType nucleus; BSplineBasisType basis; uint32_t nSub; RadialSpinorType *orb;
                 double *qbar, *Npair; double speedOfLight, confRho, confP, confN; } DiracAtomType;
```

## 12. Checkpoints

Prototype results (Python, from a scratch implementation of exactly the equations above, double precision, $c$ = 137.035999177) are given where available; the C implementation must reproduce them. The prototype scripts are in `Doc/Claude/QuantyDFT_StageA/QuantyDFT_StageA_proto/` (one script per check, see its README).

| # | Test | Target / prototype result |
|---|---|---|
| 1a | Point nucleus, hydrogen-like; Z = 1, 80, 92; κ = −1, ±2, −3, +1 (+1 except Z = 92); order 9, 120 intervals | Exact $E_{n\kappa}=c^2[(1+(Z/c)^2/(n-|\kappa|+\gamma)^2)^{-1/2}-1]$: lowest three roots per κ: errors ≤ 3e-8 Ha (Z = 92, 1s), ≤ 5e-9 Ha otherwise. Z = 92, 1s₁/₂ = −4861.1979032 Ha with $c$ = 137.035999177. Exception: Z = 1, n = 5 has 2.6e-6, which is the finite box ($R_{\max}$ = 80 bohr), not the basis. No spurious states. |
| 1b | Point nucleus, Z = 92, κ = +1 | One intruder state at −16542.7 Ha (⟨r⟩ = 0.02 bohr, half large/half small component). All other roots remain exact. This is the reason a finite nucleus is required (§3). |
| 2 | Gaussian nucleus, Z = 92, A = 238 ($R_{\rm rms}$ = 1.0867e-4 bohr) | With $c$ = 137.035999177 and 1 bohr = 52917.7210544 fm: 1s₁/₂ = −4854.15892003, 2p₁/₂ = −1257.23950876, 2p₃/₂ = −1089.61141499 Ha. Converged to ~1e-9 Ha (relative 2e-13) for 80–160 intervals; no intruder. In `{"Precision","High"}` mode the C code should reach the relative 1e-13 level and report it. |
| 3 | Lorentz-scalar confinement, Fe²⁵⁺, $V_\beta=(r/r_0)^8$ | $E-E_{\rm free}$ on the same knots ($A=56$ Gaussian nucleus, 140 intervals, $R_{\max}=r_0\,(10^5)^{1/8}$, Rayleigh quotients), $r_0$ = 0.5, 1, 2, 4: 3d₅/₂ 1.2714638668, 1.3936838314e-2, 5.7358135983e-5, 2.2412199030e-7 Ha; 4s₁/₂ 24.515661797, 1.5071830268, 1.3505189202e-2, 5.4031644439e-5 Ha (A2; the prototype at $R_{\max}=40$ with `dsygvd` gave 1.27, 1.4e-2, 5.7e-5, 2.2e-7). Small component decays with large component. $\langle\beta V_\beta\rangle\to\Delta E$ for weak confinement. The vector variant shows negative-continuum intruders (ΔE ~ 10⁴ Ha; with $R_{\max}=40$ at $r_0$ = 1–2, with the automatic $R_{\max}$ at every $r_0$). |
| 4 | Galerkin $Y^k$, k = 0–4 vs. direct quadrature; hydrogen 1s $Y^0=1-e^{-2r}(1+r)$ | 80 intervals: 7e-9; 160 intervals: 8e-12 (order 9) |
| 5 | Configuration average vs. Cowan (symbolic) | identical for $p^q,d^q,f^q$ (q = 2, 3, 5) and $sp,pd,df,sd$ pairs; the fixed-occupation variant fails ($d^3$: $\tfrac{44}{1875}F^2-\tfrac{8}{375}F^4$) |
| 6 | Fock gradient (§8.4) vs. finite differences, Ni-like $3d_{3/2},3d_{5/2},4s$, $3d^84s^2$ | relative agreement 2e-11 – 8e-11 |
| 7 | DHF closed shells, point nucleus: He, Ne, Ar | Compare with GRASP/literature DHF totals (He −2.8618133 Ha, Ne −128.6919 Ha — confirm to µHa with GRASP2018 at the same $c$). |
| 8 | Nonrelativistic limit ($c\to 10^4$) of DHF-AOC vs. Cowan RCN (HF, no relativistic corrections), e.g. Ni 3d⁸4s², C 2p² | $E_{av}$, ε, $F^k(3d,3d)$, $G^k$: agreement to RCN's numerical accuracy |
| 9 | Dirac–Kohn–Sham vs. NIST atomic reference data | NIST uses VWN correlation; its RLDA adds the MacDonald–Vosko factor to exchange. Compare LDA with VWN5 (nonrelativistic, $c\to\infty$) and RLDA with VWN5 + MV to µHa. First pin down NIST's nucleus model and $j$-subshell occupations from the documentation. |
| 10 | PBE atom | Compare to an independent all-electron code (e.g. FHI-aims or PySCF with a large basis); consistency of the two total-energy expressions (§9) to 1e-9 Ha |
| 11 | Fourier–Bessel | Plancherel $\int k^2(\tilde P^2+\tilde Q^2)dk=1$; inverse transform reproduces $P,Q$ |
| 12 | `NewOperator("Dirac")` with `"DKBSpline"` input, same nucleus and $c$ | one-electron matrix equals the solver's to $10^{-12}$ relative; diagonal equals $\varepsilon$ for unorthogonalized eigen-orbitals of the bare nucleus |
| 13 | Nucleus `"Function"` given the Gaussian $\rho_N(r)$; `"Function3D"` given a spherical $\rho_N$ | $V_N$ equal to the analytic erf form, and orbital energies equal to the `"Gaussian"` run, to $10^{-10}$ relative; `"Function3D"` spherical: all $V_{L>0}<10^{-12}$ |
| 14 | On-site orthonormalization (§6.1) with deliberately different $Q^{\rm ion},P$ per orbital | $\|C^TS_\kappa C-1\|_{\max}<10^{-13}$; first orbital per $\kappa$ unchanged to $10^{-14}$ (A2: with $Z$ and $r_c$ as stand-ins, 2.2e-16 and 0; the $Q^{\rm ion}$ version in A3/A4) |

## 13. Open points

1. **Slater-integral export convention for Quanty:** jj-resolved $R^k$ (four radial variants for $d$–$d$), the $g$-weighted $j$-average, or both. Proposal: both, with LS-Quanty scripts using the average plus $\zeta$.
2. **NIST validation details** (check 9): nucleus model, $c$ value and subshell occupations of the RLDA reference.
3. **Exchange in the DHF eigenproblem:** treat $x_a$ through the previous orbital (simple, linear convergence) or solve the inhomogeneous linear system per iteration (faster). Decide at implementation, after check 6.

## References
- V. M. Shabaev et al., Phys. Rev. Lett. 93, 130405 (2004) — dual kinetic balance.
- L. Visscher, K. G. Dyall, At. Data Nucl. Data Tables 67, 207 (1997) — nuclear models.
- I. P. Grant, *Relativistic Quantum Theory of Atoms and Molecules* (Springer 2007) — DHF, angular coefficients.
- P. Jönsson et al., Comput. Phys. Commun. 237, 184 (2019) — GRASP2018 (EAL).
- R. D. Cowan, *The Theory of Atomic Structure and Spectra* (1981), ch. 6–7 — average-of-configuration HF.
- C. Froese Fischer et al., B-spline HF methods; C. Froese Fischer, *The Hartree–Fock Method for Atoms* (1977).
- A. H. MacDonald, S. H. Vosko, J. Phys. C 12, 2977 (1979); NIST Atomic Reference Data for Electronic Structure Calculations (Kotochigova et al.), https://math.nist.gov/DFTdata/atomdata/
- K. Koepernik, H. Eschrig, Phys. Rev. B 59, 1743 (1999) — FPLO basis concept (confined, compressed atomic orbitals).
