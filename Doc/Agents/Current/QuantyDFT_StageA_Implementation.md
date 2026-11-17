# Atomic radial functions in Quanty: the stage-A implementation

This document describes the code of stage A of the DFT/HF add-on as it is. Stage A gives Quanty its own routine for atomic radial functions: relativistic (Dirac) or nonrelativistic (Schrödinger) atoms, self-consistent in configuration-average Hartree–Fock or in Kohn–Sham density functional theory, on a B-spline basis. The results feed Quanty's existing atomic operators directly.

Companion documents:

- `QuantyDFT_StageA_DiracAtom.md` — the design record, finished with task A5e on 2026-11-11: equations (sections 1–11), checkpoints (section 12), and in section 0.2 the history of every task (A0–A5e) with the decisions and the deviations from the design. This document does not repeat that history; it describes the result.
- the folder `QuantyDFT_StageA_proto`, next to the design record — the Python prototype of the design checks, and the C harness `check6_fock_gradient.cpp`.

Contents:

1. What it does, in one script
2. Lua interface
3. Conventions
4. Numerical methods
5. Code map
6. Verification
7. Known limits and open points

---

## 1. What it does, in one script

```lua
-- the minimal call: the ground state of the neutral Ni atom (Dirac-Fock, Gaussian nucleus of 58Ni)
G, F, Atom = CalculateAtomicBasisFunctions(28)

-- a basis: orbitals up to 5s, 4p and 4d, all confined, the empty 4p and 4d contracted
G, F, Atom = CalculateAtomicBasisFunctions(28, {
    {"Configuration", "[Ar] 3d8 4s2"},
    {"NMax", {s = 5, p = 4, d = 4}},
    {"Confinement", {Rho = 3.5, N = 8}},
    {"Q", {["4p"] = 2, ["4d"] = 3}} })

indices, NF = CreateAtomicIndicesList(Atom.Labels, {{"Kappas", true}})
H0 = NewOperator("Dirac", NF, 28, indices, G, F, Atom.Nucleus)   -- exact one-electron part
U  = NewOperator("AU", NF, indices, G, F)                        -- exact Coulomb part
print(Atom.Etot, Atom.Zeta["3d"], Atom.SlaterIntegralsNL["3d 3d 3d 3d"][2])
```

`CalculateAtomicBasisFunctions` replaces `ReadFPLOBasisFunctions`: it returns the large and small components `G` and `F` of every orbital and a table `Atom` with energies, Slater integrals, spin–orbit constants and diagnostics. `NewOperator("Dirac")` and `NewOperator("AU")` recognise the returned functions (type `"DKBSpline"`) and compute their matrix elements exactly on the same basis, instead of on a grid.

---

## 2. Lua interface

### 2.1 `CalculateAtomicBasisFunctions(Z [, {options}])`

Returns `G, F, Atom`. Every returned function (G, F and the potentials) gives the quadrature points of the radial basis as `f.x`, so `x = G[1].x` where the points are needed. This is the only Lua entry to the atom solvers: the temporary `DiracKohnSham` and `DiracHartreeFock` of stages A3 and A4 were merged into it on 2026-11-15 (their extra fields moved to `Atom`).

**Arguments:** `Z`, the nuclear charge (a real number), and the options. The option `Configuration` (an argument before 2026-11-16) is the configuration, e.g. `"[Ar] 3d8 4s2"`, per relativistic subshell `"[Kr] 4d3/2(4) 4d5/2(5)"`, or both kinds mixed. An nl entry is averaged over its 4l+2 spinors (Cowan), an n κ entry over its 2j+1 spinors; fractional occupations are allowed, and an occupation 0 (`"4p0"`, `"5p1/2(0)"`) adds an empty orbital. Without a configuration the ground state of the neutral atom is used (`AtomGroundStateConfiguration`, from the NIST Atomic Spectra Database, Z = 1–103). No option is required.

**Which orbitals:** every orbital of the configuration and all `n` up to `NMax` are calculated and returned, ordered by n, l and j. `Atom.Labels` lists them in that order, ready for `CreateAtomicIndicesList`; a selection is made in Lua.

Options:

| Option | Default | Meaning |
|---|---|---|
| `Method` | `"DHF"` | `"DHF"` (configuration-average Hartree–Fock), `"LDA"` (Slater + PW92), `"LDA-VWN"`, `"RLDA-MV"`, `"PBE"`, `"PBE-VWN"` |
| `Relativistic` | `true` | `false`: Schrödinger equation; `G` are then B-splines and `F` is empty |
| `NMax` | none | `{s = 7, p = 6, d = 6, f = 5}`; per κ: `["p1/2"] = 7` |
| `Q` | none | `{["4p"] = 1, ["5s1/2"] = 2}`: extra nuclear charge of the equation of an orbital; an nl key sets both j (n κ wins); only calculated orbitals |
| `P` | 1 | `{["4d5/2"] = 0.8}`: confinement scale of an orbital, likewise; needs `Confinement` |
| `Confinement` | none | `{Rho = 6, N = 10}`: confinement radius r_c = P(n κ) · Rho for every orbital |
| `ConfinementVmax` | 1e5 Ha | with a confinement and no `rMax`: rMax = max r_c · Vmax^(1/N) |
| `Orthonormalize` | `"GramSchmidt"` | per κ from the lowest n upwards; or `"Lowdin"` |
| `RedundancyTolerance` | 1e-8 | error if an overlap eigenvalue is smaller (linearly dependent functions) |
| `FourierBessel` | none | `{KMax = 30, NK = 1001}`: momentum-space transforms in `Atom.FourierBessel` |
| `Precision` | `"Normal"` | `"High"`: EnergyTolerance 1e-13 |E|, OrbitalTolerance 1e-9, NQuadrature Degree+9, ConvergenceCheck on |
| `ConvergenceCheck` | `false` | repeat with twice the knot intervals and report the changes |
| `Nucleus` | `"Gaussian"` | `"Point"`, `"Gaussian"`, `"Uniform"`, `"Fermi"`, `"Function"`, `"Function3D"` |
| `A`, `Rrms` | most abundant isotope | mass number, or the rms radius in bohr (Rrms = (0.836 A^(1/3) + 0.570) fm from A); default `AtomDefaultMassNumber` (NIST isotopic compositions) |
| `FermiThickness` | 2.30 fm | skin thickness of the Fermi nucleus |
| `Density` | | Lua function ρ(r) (`"Function"`) or ρ(x, y, z) (`"Function3D"`), charge per volume, rescaled to Z |
| `DensityTolerance` | 1e-10 | relative charge error above which a warning is printed |
| `NucleusLmax`, `NucleusNTheta`, `NucleusNPhi` | 4, Lmax+8, 2Lmax+16 | multipole expansion and angular grid of `"Function3D"` |
| `Degree`, `NIntervals` | 8, 120 | B-spline degree and number of knot intervals |
| `rMin`, `rMax` | Rrms/50 or 1e-6/Z; 40 | first interior knot and the end of the knots (bohr) |
| `NQuadrature` | Degree+5 | Gauss points per knot interval (at most 127) |
| `MaxIterations`, `History` | 200, 8 | SCF iterations and Pulay history |
| `Mixing`, `OrbitalMixing` | 0.3, 0.5 | Kohn–Sham potential mixing; orbital (DIIS) mixing |
| `EnergyTolerance` | 1e-10 Ha | change of Etot between iterations |
| `ResidualTolerance` | 1e-8 | Kohn–Sham: density-weighted potential residual |
| `OrbitalTolerance` | 1e-8 | Hartree–Fock and generating SCFs: the orbital step |
| `MixDensity` | GGA: true | Kohn–Sham: mix the density instead of the potential |
| `DensityCutoff` | 1e-30 | density below which e_xc = 0 |
| `RayleighQuotient`, `OriginBoundary` | `true` | see section 4.2 |

The returned `Atom` table:

- **Ground state:** `Etot`, `Ekin`, `Eenuc`, and `Edirect`, `Eexchange` (DHF) or `Ecoul`, `Exc` (Kohn–Sham); `EtotEigenvalueSum`, Etot through the eigenvalues (Kohn–Sham: Σ q̄ε − ∫(f v + f′ g) + Ecoul + Exc, equal at self-consistency; DHF: ½ Σ q̄ (ε + I), equal for any orbitals); `Converged`, `Iterations`, `Residual`, `dE`; Kohn–Sham: `PotentialResidual` (|EtotEigenvalueSum − Etot| ≲ √N times it) and `MixDensity`; DHF: `Gradient`, `PairGradient`; Schrödinger DHF: `EavJJ` (E_av with the j-subshell angular factors, = Etot, design check 5); Schrödinger: `VirialRatio` = −(Etot − Ekin)/Ekin.
- **Echo of the input:** `Method`, `Relativistic`, `c`, `Z`, `Configuration`, `A`, `NElectron`, `rMin`, `rMax`, and `Nucleus` as an options list that `NewOperator("Dirac", ..., Atom.Nucleus)` accepts.
- **Per orbital, `Orbitals[i]`:** `Label`, `n`, `l`, `kappa`, `j`, `Occupation`, `Q`, `P`, `rc`, `Energy` (eigenvalue of its generating equation), `EnergyGround` (⟨a|F|a⟩ of the returned function with the ground-state operator), `Vbeta` (⟨βV_β⟩), `RExpectation` (⟨r⁻³⟩ for l > 0, ⟨r⁻¹⟩, ⟨r⟩, ⟨r²⟩), `EtotGenerating`, `Generated`, `Converged`.
- **Slater integrals:** `SlaterIntegrals` in the format of `GetSlaterIntegrals` — `t["3d5/2 4s1/2 3d5/2 4s1/2"] = {[k] = F^k}` and `t["3d5/2 4s1/2 4s1/2 3d5/2"] = {[k] = G^k}` — jj-resolved in Dirac mode; `SlaterIntegralsNL`, the same for nl shells with the radial integrals averaged over j with weights (2j+1)/(4l+2); `Zeta[nl]` = 2(⟨F⟩_{l+1/2} − ⟨F⟩_{l−1/2})/(2l+1).
- **Diagnostics:** `MinOverlapEigenvalue`, `OrthonormalityError`, `Potential = {rVN, rVH, rVxc}` ($rV$ as BSplines through the quadrature points; `rVxc` Kohn-Sham only, the local potential also for a GGA: $v - g'$, $g'$ the derivative of a BSpline of $g$), `ConvergenceCheck = {NIntervals, dEtot, dEnergy}`, `FourierBessel = {k, P, Q}`, and for `"Function3D"` `NuclearMultipoleMoments[i] = {L, M, Value}`.

### 2.2 Operators on the returned functions

These dispatch to the exact DKB paths when `G` or `F` contains a `DKBSpline` (then all functions must be DKB spinors of one knot set, made with the present c):

- `NewOperator("Dirac", NF, Z, indices, G, F, {nucleus options})` — ⟨a|h_D|b⟩ = c_aᵀ H_κ c_b with the DKB matrix of the nucleus given (default point nucleus). For grid functions (FPLO) the old path is unchanged; a `Nucleus` option is an error there.
- `NewOperator("AU", NF, indices, G, F [, kappas1, kappas2] [, {options}])` — the Coulomb interaction with R^k from the Galerkin Poisson solution, assembled exactly as `Create_Relativistic_Coulomb_Operator` does with grid integrals. The forms and options of the grid path (`AtomicU.cpp`), with the same tuple generators and assembly routines:
  - `{"conserve", true}`: only the terms that conserve every orbital's occupation, (ab|ab) and (ab|ba) (`CreateSlaterTupleTwoFilled`);
  - `kappas1, kappas2`: the first #kappas1 orbitals are a core, the rest valence; only the core-valence terms (`CreateSlaterTupleOneFilled`). With `{"1P-Operator", true}` (the default) their mean field with the core filled (`MeanFieldOperatorOld`): the Hartree–Fock potential of the core on the valence. The DKB path fills the fermion indices of the core orbitals; the grid path fills the first fermion indices, which is the same only when the core comes first in the indices.
- `NewOperator("NuclearMultipole", NF, Z, indices, G, F, {nucleus options})` — the L ≥ 2 part of a `"Function3D"` nucleus as a one-particle operator (e.g. the electric quadrupole hyperfine interaction).

### 2.3 Other Lua functions

The temporary test entries of stages A1–A4 (`DiracRadialSolve`, `RadialPoisson`, `DiracOrthonormalize`, `DiracKohnSham`, `DiracHartreeFock`) were merged into section 2.1 on 2026-11-15: a configuration of empty orbitals only (`"1s0 2s0 2p0"`) is the bare nucleus, `Atom.Potential.rVH` the Poisson solution, the option `Orthonormalize` the on-site orthonormalization.

| Function | Returns |
|---|---|
| `XCFunctional(method, n, sigma)` | `e, de/dn, de/dsigma` (kept: the only Lua access to the functionals) |
| `math.SphericalBesselJ(l, x, {{"Method", ...}})` | j_l(x); `"Recurrence"` (default) or `"Series"` |

---

## 3. Conventions

- **Units:** Hartree atomic units. The speed of light is `PhysicalConstants.cAU.value`, read at every call and passed down; no C routine uses a fixed c. Energies are measured from mc².
- **Radial functions:** `G` = P = r g and `F` = Q = r f, with ψ = (P Ω_κμ, i Q Ω_−κμ)/r, normalized as ∫(P² + Q²) dr = 1. This is FPLO's convention, and the one `NewOperator("Dirac")` (Martin's `DiracH`) and the relativistic Slater integrals already used.
- **κ and labels:** κ = −(l+1) for j = l+½ and κ = l for j = l−½; labels `"4d5/2"` (Dirac) or `"4d"` (Schrödinger), l up to 6 (`spdfghi`). μ runs from −j to j in the index lists of `CreateAtomicIndicesList`.
- **Sign:** P > 0 on its innermost lobe (`InnermostLobeSign`): the first local extremum of |P| above ε^(1/3) max|P|.
- **Occupations:** an nl shell is shared between its j subshells in proportion to 2j+1 (the j-average of NIST and of Cowan); an n κ entry fixes the subshell occupation.
- **Real spherical harmonics** (Function3D): S_LM = (Y_L,−M + (−1)^M Y_LM)/√2, S_L,−M = i(Y_L,−M − (−1)^M Y_LM)/√2, without the Condon–Shortley phase inside; the same S_LM define the multipole densities and the operator.
- **Code:** C written in `.cpp` files, one compilation unit (`DFT/DFT.cpp` is included from `LuaInterface/main.cpp`), `malloc`/`free`, routines return 0 or 1 and do not print, `double` throughout (`typedef double DFTReal`), no hard-coded tolerances: every threshold is an option with a documented default.

---

## 4. Numerical methods

### 4.1 Radial basis

- **Knots:** exponential knot sequence (`BSplineKnotsExponential`): rMin, the first interior knot, to rMax with NIntervals intervals, and degree+1 coinciding knots at 0 and at rMax.
- **Quadrature:** Gauss–Legendre with NQuadrature points per interval from a 34-digit table (at most 127 points). All integrals are sums over these points, with Neumaier (compensated) summation where it matters.
- **Dual kinetic balance (DKB):** with the working B-splines B_i (all but the first and the last) and h = 1/(2c), the 2n basis functions are
  u_i⁽¹⁾ = (B_i, h(B_i' + κB_i/r)) and u_i⁽²⁾ = (h(B_i' − κB_i/r), B_i).
  A spinor is a coefficient vector of length 2n; P and Q are evaluated exactly through these formulas (`DKBSplineEvaluate`, class `Interpolation::DKBSpline`), never refitted. In Lua arithmetic and `==` (`LuaInterpolateAsBSpline`, `LuaInterface/Interpolate.cpp`) a DKBSpline acts as the BSpline of its degree through its values on its quadrature points; the result is a BSpline, and integrals of products such as $\int(P_aP_b+Q_aQ_b)\,dr$ are exact to rounding (Hg: 3e-15).
- **Schrödinger mode** uses the B_i alone (n functions, Q = 0).

### 4.2 One-electron problem

- **Matrices:** H and S of the DKB basis are built per κ by quadrature (`DiracDKBMatrices`); a potential enters as M[V]_ab = ∫ V (P_a P_b + Q_a Q_b) dr, a Lorentz-scalar confinement as ∫ V_β (P_a P_b − Q_a Q_b) dr.
- **Eigenproblem:** the generalized symmetric problem is solved with LAPACK `dsygv` (QR iteration), not `dsygvd`: divide and conquer gives the small eigenvalues of this strongly graded spectrum only to ε·max|E|. With `RayleighQuotient` the energies are the Rayleigh quotients of the eigenvectors.
- **Origin boundary:** with `OriginBoundary` the DKB functions on B_1 that are nonzero at r = 0 are left out for κ ≠ −1, which removes the spurious state of a point nucleus for κ = +1.
- **Root selection:** state (n, κ) is the electronic root n − l − 1, counted above −c².
- **Nucleus models** (`NucleusPotential`): point; Gaussian (analytic erf); uniform sphere (analytic, R = √(5/3) Rrms); Fermi, a Lua function ρ(r), and the spherical part of `"Function3D"` through the Poisson equation.
- **Confinement:** V_β = (r/r_c)^N as a Lorentz scalar (β V_β, a position-dependent mass).

### 4.3 Poisson equation and Slater integrals

- **Method:** Y^k for any radial pair density is the Galerkin solution of Y'' − k(k+1)Y/r² = −(2k+1)ρ/r on the B-splines (with a Robin condition at rMax), factorized once per k (`BSplinePoissonSetup`, `BSplinePoissonSolve`).
- **Symmetry:** the discrete R^k(ab;cd) = (2k+1) b_acᵀ A⁻¹ b_bd is exactly symmetric, so the discrete Hartree–Fock energy is a polynomial in the coefficients and the Fock operator is its exact derivative.
- **Slater integrals:** the same route gives the Slater integrals in `Atom.SlaterIntegrals` and the `"AU"` operator.

### 4.4 Kohn–Sham atom

- **Functionals:** `XCFunctional.cpp` holds Slater exchange, PW92, VWN5, the MacDonald–Vosko relativistic exchange factor, and PBE (on PW92 or on VWN5), spin unpolarized, with e, ∂e/∂n and ∂e/∂σ analytic.
- **GGA:** enters in the weak form, v = V_H + ∂e/∂n − 4(∂e/∂σ)n'/r, g = 2(∂e/∂σ)n'.
- **SCF** (`AtomSCF`): Thomas–Fermi start, Pulay (DIIS) mixing of the potential (LDA) or of the density (GGA). The residual is weighted with the density, because the tail density is rounding noise.
- **Energies:** Etot = Ekin + Eenuc + Ecoul + Exc; the eigenvalue-sum expression is reported for the first-order consistency check.

### 4.5 Configuration-average Hartree–Fock

- **Energy functional** (`AtomHF.cpp`): E_av averaged over all determinants of the configuration, with the pair numbers N_ab of design section 8.1 and the angular factors Γ^k = (j_a k j_b; ½ 0 −½)² (Dirac) or ½(l_a k l_b; 0 0 0)² (Cowan).
- **Fock operator:** F_a = h + Σ_b C_ab (J_b − Σ_k Γ^k_ab K^k[b]), with exact exchange matrices K^k[b] = (2k+1) Bᵀ A⁻¹ B. The self term is written in this same J − K form, so every closed subshell of a κ has the closed-shell (Roothaan) operator.
- **SCF** (`HFSCF`): starts from the LDA-VWN atom. Each iteration:
  1. solves one eigenproblem per class of subshells with the same operator, with the other occupied orbitals of the κ projected out exactly (a level shift on their span);
  2. takes Newton steps for the rotations between occupied subshells of one κ (Brillouin), only where the energy is not invariant under them;
  3. applies DIIS on the stacked coefficients, followed by Löwdin orthonormalization.

  It converges on the orbital step ‖c_out − c_in‖_S, which contains both stationarity conditions; the gradient in its dual norm has a rounding floor and is reported only.
- **Identity:** E = ½ Σ q̄ (ε + I) holds for any orbitals and checks the Fock matrices against the energy.
- **Empty subshells** (q = 0) have the Fock operator of V^N, the q → 0 limit.
- **Kohn–Sham mode:** the same orbital SCF has one (h + v_KS[ρ]), used for the generating calculations below.

### 4.6 Basis functions (`AtomBasis.cpp`)

1. **Ground state:** Kohn–Sham, or Hartree–Fock from the LDA start. Empty orbitals come from the ground-state operator (V^N Fock operator, or the Kohn–Sham potential), root n − l − 1.
2. **Generating calculations:** for every orbital with Q ≠ 0 or a confinement, its own generating calculation, in which only its own equation gets ΔH = M[(Q/Z)V_N] + M_β[V_β]; all other electrons see Z and no confinement.
   - An occupied target: the orbital SCF with the target as its own class.
   - An empty target: one eigensolve.
3. **Orthonormalization:** the functions of different generating calculations are not orthogonal; they are orthonormalized per κ, Gram–Schmidt from the lowest n upwards (Cholesky of the overlap in coefficient space, two passes).

**Zero mode.** If the target and another occupied orbital of its κ are both closed (HF) or equally occupied (KS), E is invariant under their rotation. Only the extra term would decide it, and a contracted target would turn into the 1s-like combination. Such rotations are therefore not taken: the target stays the eigenvector of its own operator, orthogonal to the others.

### 4.7 Further outputs

- **Function3D** (`NucleusMultipole.cpp`): ρ(x, y, z) is sampled on a Gauss–Legendre × uniform angular grid at every radial point and projected on S_LM (even L ≤ NucleusLmax). The potentials are V_LM = −4π/(2L+1) Y^L[r²ρ_LM]/r, and ⟨κμ|S_LM|κ'μ'⟩ comes from Wigner–Eckart in Grant's convention (the same for large and small components).
- **Fourier–Bessel** (`FourierBessel.cpp`): P̃(k) = √(2/π) ∫ r P j_l(kr) dr and Q̃ with l̄, by Gauss–Legendre per knot interval with (degree+1) + ⌈k_max Δr/π⌉ points.
- **Spherical Bessel functions:** `SphericalBesselJ` in Core uses the recurrence method `SphericalBesselJRecurrence`: upward for x > l, Miller downward with the sum rule for x ≤ l, a series for tiny x. It supersedes the former series method, which is kept as `SphericalBesselJSeries`.

---

## 5. Code map

Everything is in `DFT/` unless noted. Header files declare the routines of their `.cpp`.

| File | Contents |
|---|---|
| `DFT.h`, `DFT.cpp` | include order of the add-on; `DFTReal` |
| `DiracRadial.cpp` | radial basis (`RadialBasisCreate`, `RadialBasisFromKnots`); nucleus models (`NucleusPotential`, `NucleusDensity`, `NucleusCharge`, `NucleusFermiParameters`); DKB matrices (`DiracDKBMatrices`, `DiracActiveFunctions`); eigenproblem (`DiracRadialEigen`, `DiracGeneralizedEigen`, `DiracElectronicFirst`, `DiracRootIndex`); spinor values (`DiracComponents`, `DiracExpectation`); sign (`InnermostLobeSign`, `DiracSignConvention`); confinement (`ConfinementPotential`, `ConfinementRMax`); orthonormalization (`DiracOrthonormalize`, `DiracOrthonormalityError`); `NeumaierAdd` |
| `XCFunctional.cpp` | `XCEvaluate`, `XCPW92`, `XCVWN`, `XCMacDonaldVosko`, `XCMethodFromName` |
| `AtomSCF.cpp` | configuration parser (`AtomParseConfiguration`, `AtomShellType`); NIST tables (`AtomGroundStateConfiguration`, `AtomDefaultMassNumber`); the atom (`DiracAtomType`, `AtomSetup`, `AtomFree`); orbitals on the points (`AtomOrbitalOnPoints`, `AtomLocalFunctions`, `AtomPotentialMatrix`); Kohn–Sham density and potential (`AtomDensity`, `AtomPotential`); the Kohn–Sham SCF (`AtomSCF`); `AtomSignConvention` |
| `AtomHF.cpp` | `DiracHFType`; angular factors and pair numbers (`HFGamma`, `HFSetup`, `HFCoupling`); energy (`HFEnergy`); Fock and exchange matrices (`HFFockMatrices`, `HFExchangeMatrix`); the orbital SCF (`HFGradient`, `HFClassOrbitals`, `HFOrthonormalize`, `HFSCF`); `RadialYkOnPoints`, `HFIntegral` |
| `AtomBasis.cpp` | generating calculations (`AtomTargetMatrices`, `AtomRootSolve`, `AtomBasisFunctions`, `AtomBasisType`); `AtomOrthonormalizeBasis` |
| `NucleusMultipole.cpp` | real spherical harmonics, multipole projection, potentials, moments, `RealHarmonicMatrixElement` |
| `FourierBessel.cpp` | `FourierBesselTransform` |
| `LuaDiracRadial.cpp` | nucleus and basis options of every Lua entry (`DFTLuaNucleusOptions`, `DFTLuaNucleus`, `DFTLuaNuclearPotential`, `DFTLuaNucleusBasis`) |
| `LuaAtom.cpp` | `XCFunctional`; `DFTLuaPushRadialFunctions`, `DFTLuaHasOption` |
| `LuaAtomicBasis.cpp` | `CalculateAtomicBasisFunctions`; labels, Slater-integral tables |
| `LuaAtomicOperators.cpp` | the DKB paths of `NewOperator("Dirac")`, `("AU")`, `("NuclearMultipole")` |
| `Core/BasicMath/Interpolation/` | `BSplineKnotsExponential` (`BSpline.cpp`); `BSplinePoissonSetup`/`Solve` (`BSplineMatrix.cpp`); the `DKBSpline` class (`DKBSpline.cpp`); the Gauss–Legendre table |
| `Core/BasicMath/BasicMath_StandardFunctions.cpp` | `SphericalBesselJ`, `SphericalBesselJSeries`, `SphericalBesselJRecurrence` |
| `AtomicOperators/AtomicDirac.cpp`, `AtomicU.cpp` | dispatch to the DKB paths; `OperatorAddKappaBlock` |
| `LuaInterface/ssp.cpp` | registration of the Lua functions and of `NewOperator("NuclearMultipole")` |

**Data flow of `CalculateAtomicBasisFunctions`:**
1. parse Z, the configuration (default `AtomGroundStateConfiguration`) and the options;
2. collect the orbitals (configuration and `NMax`), add empty shells for the `NMax` orbitals not in the configuration, and attach Q and P;
3. nucleus and basis (`DFTLuaNucleusBasis`);
4. `AtomSetup`;
5. `AtomBasisFunctions`: the ground state through `AtomSCF` and `HFSetup`/`HFSCF`, then the generating calculations;
6. `AtomOrthonormalizeBasis`;
7. output: `DFTLuaPushRadialFunctions`, the `Atom` table, the Slater tables, and optionally the Fourier–Bessel transforms and the convergence check (a second call with twice the intervals).

---

## 6. Verification

All tests are in `TestScripts/DiracAtom/` (CI job `test-jobs-diracatom`), plus `TestScripts/SphericalBesselJ.Quanty` and the source check `TestScripts/NoSinglePrecision.sh`. Thresholds come from measurements on two platforms or from bounds, never from a single run.

| What | Reference | Agreement | Test |
|---|---|---|---|
| Hydrogen-like ions, point nucleus | analytic Dirac energies | ≤ 5e-9 Ha (U 1s: 3e-8 Ha) | `PointNucleus` |
| U, Gaussian nucleus | converged prototype values | ≤ 1e-8 Ha (three grids agree to 1e-9 Ha) | `GaussianNucleus` |
| Nucleus `"Function"` | the `"Gaussian"` energies | ≤ 1e-10 relative | `FunctionNucleus` |
| Lorentz-scalar confinement | converged prototype shifts (Fe²⁵⁺) | ≤ 1e-8 relative | `Confinement` |
| Galerkin Y^k | analytic hydrogen Y⁰ (`rVH`); hydrogen Slater integrals k = 0, 1, 2 | 1e-8 (80 intervals), 1e-11 (160); 1.8e-11 | `RadialPoisson` |
| On-site orthonormalization | C^T S C = 1 | < 1e-13 (two passes) | `Orthonormalize` |
| LDA, RLDA atoms He–U | NIST SRD 141 | Etot ≤ 5.3e-7 Ha, eigenvalues ≤ 1.2e-6 Ha | `NISTAtoms` |
| PBE-VWN and VWN atoms | Lehtola, PRA 101, 012516 | ≤ 6.1e-7 Ha | `PBEAtoms` |
| HF limits He–Xe (nonrelativistic) | Lehtola 2019/2020 | ≤ 5.1e-10 Ha (table rounding) | `HartreeFockLimits` |
| Dirac–Fock He | design value | −2.8618133 Ha, 7 digits | `HartreeFockLimits` |
| j-subshell energy = Cowan's | identity | ≤ 1.8e-12 | `HartreeFockAverage` |
| Fock operator = energy gradient | finite differences | 2.5e-13 – 2.8e-10 | C harness (check 6) |
| Exact `"Dirac"` operator | DKB eigenvalues | 3.7e-15 relative | `NewOperatorDKB` |
| `"AU"` operator | mean of ⟨det\|H\|det⟩ over all determinants = E_av | 2.4e-15 relative | `NewOperatorDKB` |
| `"AU"` options | conserve = full within every occupation sector; core form = full − E_core on core + one valence electron | 5e-15 relative | `NewOperatorAUCore` |
| Generating calculations | Hellmann–Feynman in Q and in r_c | ≤ 3.5e-9, ≤ 3.6e-8 | `AtomicBasisFunctions` |
| V^N virtual orbitals | Koopmans with Quanty operators | 3.4e-14 | `AtomicBasisFunctions` |
| Galerkin Slater integrals | Martin's cell method | 9.5e-13 | `AtomicBasisFunctions` |
| Function3D | Gaussian nucleus; analytic moments; NewOperator("CF") rotated to jj | 2.1e-11; 2.7e-13; 5.7e-15 | `Function3DNucleus` |
| Fourier–Bessel | analytic hydrogen; Plancherel | 8.1e-10; 2.3e-10 | `FourierBessel` |
| Report fields | hydrogen ⟨rᵏ⟩, virial; Ne HF virial | 12 digits; 3.8e-10 | `AtomicBasisReport` |
| Re, PW92, Dirac | FPLO 21.00 eigenvalues and G, F | ≤ 8.7e-8 Ha; ≤ 2.7e-6 | `FPLOAtom` |
| Spherical Bessel functions | mpmath (30 digits) | ≤ 1.3e-13, l ≤ 60 | `SphericalBesselJ` |

The FPLO comparison uses FPLO's own occupations (5d3/2 full, one electron in 5d5/2) and one fitted uniform-sphere radius (6.90765721 fm, set by the 1s). It confirms the PW92 parameters and the Dirac–Kohn–Sham chain against an independent code.

---

## 7. Known limits and open points

- **SCF floors** (measured):
  - the Hartree–Fock orbital step bottoms out at 1.6e-10 for Hg with a finite nucleus and 2.3e-9 with a point nucleus;
  - the Kohn–Sham DIIS residual stalls at 1.7e-9 (Re), 7.5e-9 (U, finite nucleus) and 8.7e-8 (U, point nucleus).

  Point nuclei are for comparisons; use a finite nucleus for precision work.
- **Spin:** spin-unpolarized Kohn–Sham only; no spin-polarized atom in stage A.
- **Performance:** a heavy Dirac atom costs about 1 s per Hartree–Fock iteration; a basis in which every orbital has its own generating calculation costs one SCF per occupied orbital (Ni, 15 orbitals: 43 s).
- **Operators:**
  - `NewOperator("Dirac")` has no finite nucleus for grid (FPLO) functions;
  - the multipole operator covers even L only and no magnetic hyperfine fields.
- **Interfaces:** the temporary entries of section 2.3 have no `Precision` preset. The interfaces for the later stages (environment-dependent compression, inter-site overlaps S(R)) do not exist yet; the stage-A design lists them.
- **Labels:** l is limited to 6.
