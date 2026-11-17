# DMFT and Green's-function building blocks

What a dynamical mean-field (DMFT) loop in Quanty is built from: which Lua function
does each step, what it returns and in which representation, and the physics and
numerics that decide which route to take. A DMFT loop is a Lua script; Quanty
supplies the steps, not the loop. The worked loops are tutorials in the separate
tutorials repository (`06_Dynamical_Mean_Field_Theory`); their settings and run
histories stay there.

The response-function representations (`l`, `t`, `a`, `n`, `d`) and how `A0`, `B0`
compose are in `QUANTY_NOTES_ARCHITECTURE.md`. Conjugation of the head, character
versus residue, `BlockListOfPoles` storing residues and what `ToMatrix` returns are in
`CONVENTIONS.md`, §4. Known defects in the conversions are in `TODO.md` and
`TODO_ConversionAccuracy.md`. Code is cited by file and function name. Notation as in
`CONVENTIONS.md`, §1: `^*` conjugate, `^T` transpose, `^dag` both.


## 1. The loop, one function per step

| step | Lua call | returns |
|---|---|---|
| 1. lattice | `G, Aug = CalculateG(TB, Sigma, {..., {"HybridizationFunction",true}})` | `G`: block list of poles, the local Green's function. `Aug`: one augmented propagator per atom (§3.3) |
| 2. hybridisation | `Delta = ResponseFunction.Downfold(Aug[a], R)` | `Delta(w)`, the hybridisation without the on-site energy |
| 3. Weiss field | `Gbath = ResponseFunction.FromEffectiveHamiltonian(Delta + eps0)` | `Gbath` with `Gbath^-1 = w - eps0 - Delta` |
| 4. bath geometry | `GNat = ResponseFunction.ChangeType(Gbath, "Nat", {...})` | block natural-impurity-orbital form: a valence and a conduction chain |
| 5. impurity model | `H0, BitField, Tails = ResponseFunction.ToOperator(GNat, {{"AddSpin",true},{"Nval",Nv},{"Ncon",Nc}})` | the one-particle impurity operator, the bit masks and `T`, the cut-off chain tails |
| 6. interaction, double counting | `MeanFieldOperator`, `Eigensystem` | the cluster ground state(s) |
| 7. impurity spectrum | `CreateSpectra(H, {a, q}, psi, Hyb, {{"Tensor",true}, ...})` with `q = [a, H1]` | the augmented propagator of `(a, q)`, a block chain |
| 8. self-energy | `Sigma = ResponseFunction.Downfold(Gtilde, R)` | `Sigma(w) - Sigma^HF` |
| 9. feed back | `ChangeType(Sigma, "Tri", {{"NTriMax",N}})`, mix, back to 1 | |

Steps 2–3 have two older alternatives, `CalculateGbathfromGandSigma` and
`CalculateHybridizationFunction`; §3 says why the downfold is the one to use. Step 8
is the symmetric improved estimator of Zacinskis et al.: the downfold of the
augmented propagator onto the `q` half. `Sigma^HF` is not in it; it is read off the
spectral weight of `Gtilde` itself, `B0^* B0^T`, which `ChangeType(Gtilde,"Dense")`
returns (header of `LuaResponseFunctionDownfold`). `Downfold` returns a second value,
`S`, the spectral weight of the result: the sum rule to check it against.

`R` selects the subspace: `R = ( 0 | 1 )`, `n` rows and `2n` columns, takes the lower
right block of a two by two block structure; for one orbital `Matrix.New({{0,1}})`.
Any `R` of full row rank works, its rows need not be orthonormal.

A constant response function, for a static shift, is
`ResponseFunction.New({{A0}, {}, mu=0, type="Tri"})`; block-ness has to match what
it is added to.


## 2. The lattice Green's function: `CalculateG`

`LuaCalculateG` in `LuaInterface/TightBinding.cpp`, also registered as the global
`CalculateG` in `LuaInterface/ssp.cpp`.

- **Options and defaults.** `Nk` (`{40,40,40}`), `Emin`/`Emax`/`NE` (`-10`, `10`,
  `2000`) or `EnergyGrid` (strictly increasing, not combinable with the three),
  `Type` (`"ListOfPoles"`; anything else goes through `ChangeType`),
  `HybridizationFunction` (`false`), `Method` (left out, or `"Mesh"`), `Lambda`
  (`"Mesh"` only, default 2, at least 1).
- **Sigma.** Any of the five types, block or not (a single-valued one is copied to
  block size 1), or a table with one entry per atom (`LuaReadSelfEnergyPerAtom`). Its
  block size must equal `NFermion` of the tight-binding model.
- **The default route** diagonalises `H(k)` plus the self-energy embedded as auxiliary
  states. Its cost grows steeply with the number of poles of `Sigma`. Measured, `30^3`
  k-points, `NE = 1000`: 5 poles 0.3 s, 20 poles 2.4 s, 60 poles 24 s, 120 poles
  213 s. So cut `Sigma` back before every call, `ChangeType(Sigma,"Tri",{{"NTriMax",N}})`.
  Mixing two pole lists concatenates them (`ResponseFunctionAdd` →
  `AddBlockListOfPoles`), so without the cut the count grows every iteration.
- **`{"Method","Mesh"}`** (since b36a85aa), `TightBindingHamiltonianPlusSelfEnergyOnMesh`
  in `Core/ResponseFunctionBathExcitations.cpp`. It evaluates
  `G(z) = mean_k (z - H(k) - Sigma(z))^-1` off the real axis and takes the residue at
  each grid point without a broadening (§6.2). Its cost does not grow with the poles
  of `Sigma`, about `Ny * NGrid * (Nk n^3 + NSigma n^2)`. The grid is extended until it
  holds every pole of `G`, with spacings growing by `Lambda` (f24a1c90); the result
  has poles on grid points only. Agrees with the default route to about 1e-6
  relative (`TestScripts/TightBindingDef/CalculateGMesh.Quanty`).
- **`HybridizationFunction`** (cb8a3a02) adds the second result `Aug`, filled in the
  same k-loop that builds `G` (345c0a58, `AugmentedPerAtomType` with Init/Add/Finish/
  Free in `Core/BasicMathBlock/BasicMath_BlockTypeConversion.cpp`, driven by
  `TightBindingHamiltonianPlusEmbeddedSelfEnergyToBlockListOfPolesGrid`). The embedded
  self-energy comes from the `Block*InversedToCompactMatrix` of its type, so all five
  types share one routine. Asking for it routes `G` through the grid form; `G` then
  agrees with the `Emin`/`Emax` route to 1e-13, not bit for bit. Over ten DMFT
  iterations it cost 31 s on top of 2:58 for `G` alone; as a separate k-loop it had
  cost 3.5 minutes.

**The on-site energy is in `Aug`.** `eps0_a` is the total weight (sum of residues) of
the upper right `n × n` block of `Aug[a]`, exactly, since `sum_n u (H_k u)^dag = H_k`.
It is the bare lattice on-site energy, including `ShiftOnsiteEnergy`, excluding the
`a0` of `Sigma`. The downfold leaves it out; without it `FromEffectiveHamiltonian`
puts the impurity level at zero (an error of 1.9 in a measured case, 1e-4 with it).

**Trap for anyone editing the k-loop.** `ComplexEigenFunctionToCharacter` writes the
character into `H`, so after the first eigenvector `H` no longer holds the Hamiltonian.
`H u` needs it, hence the `MatrixCopy` of `H(k)` kept for the eigenvector loop.
Without it the hybridisation came out wrong by 0.33 instead of 4e-5.


## 3. The Weiss field: three routes

### 3.1 Subtraction — `CalculateGbathfromGandSigma(G, Sigma)`

`Gbath = (G^-1 + Sigma)^-1`, `CalculateGBathFromGandSigma` in
`Core/BasicMath/BasicMath_GreenFunctionArithmetic.cpp`. `Method` is `"Subtract"`
(default) or `"Eigensystem"`; they gave bit-identical results in a large-U run. The
`Format` option is parsed but no conversion is done (a TODO in the routine): the
result has whatever type the method produced.

It is the route that fails. With a Mott pole in `Sigma`, `zeta = w - Sigma(w)` sweeps
through infinity and `G^-1 + Sigma` is a difference of two huge numbers. Against the
exact law `Delta -> <(eps - <eps>)^2> / zeta`, the subtraction is lost at
`|zeta| = 1e8` (5.4e-7 against 1.9e-9 expected); the Schur complement of §3.3 tracks
the law to `1e10`. In a controlled test with an exactly symmetric `Sigma`, where the
Weiss occupation must be 0.5, the subtraction gives 0.142: wrong by a factor 3.5
**for a perfect `Sigma`**.

### 3.2 Composition — `CalculateHybridizationFunction(Glat0, Sigma)`

Despite its name it returns the Weiss field, `1 / ( Glat0(w - Sigma)^-1 + Sigma )`,
as an Anderson matrix (`CalculateGBath`: "adds a self energy to all sites except for
the impurity site"). `Glat0` is the **non-interacting** local `G`, computed once.
Both arguments block or both single valued, same block size.

Why it is stable: with a local `Sigma` the lattice `G` depends on `w` only through
`zeta`, so `Delta(w) = R(zeta)` with `R(z) = z - 1/Glat0(z)`. For a chain
`Glat0 = b0^2/(z - a0 - b1^2/(z - a1 - ...))`, `R(z) = a0 + b1^2/(z - a1 - ...)`: the
chain with its head removed. The cancellation is done analytically, once.
`R` is a Jacobi continued fraction and `Delta = R(w - Sigma)` a composition of
Herglotz functions, so negative spectral weight is structurally impossible. Measured:
Weiss occupation 0.500000 where subtraction gives 0.142. Cut the result to a chain
(`ChangeType(...,"Tri",{{"NTriMax",N}})`); it has thousands of poles otherwise.

**Exact only for scalar `zeta`.** For a matrix `G` the substitution replaces
`<d zeta^-1 d>` by `B1^dag zeta^-1 B1`, which holds only when `zeta` commutes with
every `eps_k`. Measured on random `2 × 2` `eps_k`: scalar `zeta` 4e-13,
`zeta = diag(4i, 4.2i)` 1.2e-2, `diag(4i, 6i)` 8.4e-2. It fails already for a
constant anisotropic `zeta`. Single band: exact. Multi-orbital: use §3.3.

### 3.3 Downfold of the augmented propagator — the route to use

Take the operator pair `c` and `Hc`. For a one-body `H`, `[c_i, H] = sum_j H_ij c_j`:
the second operator is the commutator with the Hamiltonian, as `q = [a, H1]` is with
the interaction. With `Q = (w - P)^-1`, `P` the frequency-independent extended
Hamiltonian (`Sigma` as auxiliary states),

    Aug = [ <c|Q|c>     <c|Q H|c>   ]  =  [ G   Kbar    ]
          [ <c|H Q|c>   <c|H Q H|c> ]     [ K   M - <H> ]

and the Schur complement of the `(1,1)` block is `zeta - G^-1 - <H> = Delta`. That is
exactly `ResponseFunction.Downfold(Aug, ( 0 | 1 ))`, the operation already used for
the self-energy; the constant it leaves out is `<H>` here and `Sigma^HF` there.
Then `G0^-1 = w - <H> - Delta`: an Anderson matrix with head `<H>` and bath the poles
of `Delta`, a relabelling, no inversion.

- Exact for matrices, including off-diagonal `zeta` (1e-15 on random `2 × 2` cases).
- One diagonalisation of `P_k` gives all four blocks as pole lists on the same poles;
  the extra cost is two matrix products per k-point.
- `Delta` has its poles at the zeros of `G`, so some inversion is unavoidable; but
  `1 + K = zeta G` tends to 1 near the Mott pole, so the inversion inside the Schur
  complement is trivially conditioned exactly where the subtraction cancels.
- Test: `TestScripts/TightBindingDef/HybridizationFunction.Quanty`.

The recipe, with `eps0` from §2:

    G, Aug = CalculateG(TB, Sigma, {..., {"HybridizationFunction",true}})
    Delta  = ResponseFunction.Downfold(Aug[1], Matrix.New({{0,1}}))
    Gbath  = ResponseFunction.FromEffectiveHamiltonian(Delta + Eps0)

`FromEffectiveHamiltonian(X)` returns `Y` with `Y^-1 = w - X`. A `t` stays a `t`;
`l`, `a`, `n`, `d` come back as `a` (`ResponseFunctionFromEffectiveHamiltonian`).
Converting a `Delta` of about 1000 grid poles straight to `Nat` cost 5 s per loop;
cutting `Aug` or `Delta` to an `NTriMax` chain first keeps that off the bill.


## 4. The impurity: natural orbitals and `ToOperator`

`ChangeType(Gbath, "Nat", {...})` splits the bath into a valence and a conduction
chain. A pole counts as occupied by a Fermi function at `EffectivekBT`, whose default
for `Nat` and `DoubleTri` is `sqrt(DBL_MIN)`, a step. `NTriMax`, `NTriMaxVal` and
`NTriMaxCon` cap the chains (header of `LuaChangeTypeResponseFunction`).

`ToOperator` on a `Nat` returns three things:

- the operator;
- a bit-field table: `Core`, `Impurity`, `ImpurityExt`, `Bath`, `BathVal`, `BathCon`,
  `BathExt`, each with an `N...` count and, with `AddSpin`, `...Up`/`...Dn`; `NF`,
  `NB`; and **`T`**, the normalised transformation
  (`GetNormalizedTFromValandConInBlockNaturalImpurityOrbitalMatrixType`) that makes
  every orbital full or empty at `U = 0`. `Rotate(H, Matrix.Expand(BitField.T, NF, 1))`
  puts the occupied impurity combination on the `Impurity` bits and its empty partner
  on `BathExt`; scripts need not build that rotation by hand from `b0`;
- when `Nval`/`Ncon` cut a chain, a table of the discarded tails, `GvalBath` and
  `GconBath` (type `t`, `A0` set to zero), with `GvalBathIndex`/`GconBathIndex`, the
  fermion indices of the last kept sites they couple to. These are the hybridisation
  functions handed to `CreateSpectra`.

`Nval = 0` or `Ncon = 0` is refused for `Nat`: the first site of each chain is part of
the impurity. `Nval`/`Ncon` also apply to a `DoubleTri`, where an empty branch is
legal (`QUANTY_NOTES_ARCHITECTURE.md`); an insulator can keep only the valence chain.

**The split is ill posed when a pole sits on `mu`.** In a Mott insulator one pole of
the Weiss field sits essentially at `mu` (§8). Measured by feeding the analytic Weiss
field to `ChangeType(...,"Nat")` + `ToOperator` with `Nval = Ncon = 2`: a chain
empties (`NBathVal = 0` or `NBathCon = 0`) exactly when the pole's distance `d` from
`mu` exceeds `EffectivekBT`. Choose `EffectivekBT` separately from the physical `kBT`
used for fillings. It is not monotonic: 0.3 was worse than 0.001, 1.0 much better but
shifted the Hubbard peaks (±2.87 against ±3.04). Anneal it rather than fix a hot value.


## 5. The impurity spectrum and the self-energy

**`CreateSpectra` with a hybridisation.** The fourth argument is a list
`{{Delta_1, {ops_1}}, {Delta_2, {ops_2}}, ...}`, each a response function and as many
operators as its block size (`LuaCreateSpectra` in `LuaInterface/ssp.cpp`). With a
hybridisation, only `{"Tensor",true}` on a sparse route is implemented: a dense route
or `Tensor = false` stops with "not yet implemented". A small cluster falls below
`DenseBorder` and `"Auto"` picks dense, so ask for `{"Method","Lanczos"}` or
`{"Method","Mesh"}`.

- `NTriFull` (blocks of the returned chain) defaults to `NTri` times the number of
  hybridisation functions. `GridPoints` sets the grid the tail poles are binned on in
  `CalculateGforHwithSigmaAsOperator` (`Core/ResponseFunction.cpp`).
- **`GridPoints` is on the scale of `H`, not of the spectrum.** A bath state is binned
  at `e_tail + E_n`, `E_n` an eigenvalue of `H` in the extended space; `<psi|H|psi>`
  (or the `E0` option) is subtracted only afterwards. A grid built on the spectrum
  axis is off by `E0`: with `E0 = -5757` everything fell into one bin.
- Without `GridPoints`, poles are kept exactly up to 20000 bath states
  (`SigmaAsOperatorMaxBathExact`) and otherwise binned on `8 * NTriFull` points spread
  evenly over the full pole range, cut back by a memory limit. Far tail poles
  (thousands of eV) make that spacing coarse; supply a grid fine over the window and
  coarse outside. A supplied grid is used as given, plus one point at each end where
  the poles reach beyond it; this route has no `Lambda`.

**The augmented propagator** is the particle part plus the mirrored hole part,
`Gtilde = Gp + InvertEnergy(Gm)`. Both are chains, and their sum is a chain (§7.2), so
building `Gtilde` and downfolding it cost nothing measurable. Measured on a five-site
impurity: all remaining time is `CreateSpectra`, 75 s of 82 s at `NTri = 160`.

**Chain length.** For the bare cluster `Sigma` is converged once the Krylov space is
exhausted (20 blocks for a five-site impurity; quasiparticle weight `Z` identical at
20, 40, 80, 160). The bath size mattered far more: `Z = 0.126` at `Nval = Ncon = 2`,
`0.225` at 3. Once the hybridisation tails are in, that no longer holds: a chain of
60 against 240 differed by 2.9 on a scale of 4, and the exhausted Lanczos needed about
5000 blocks. The mesh route (§6.2) avoids the chain length altogether.

**Cutting `Sigma`** (`NTriMax`, moment-preserving) replaces an isolated pole at
`w = 0` by a pair straddling it, so `Sigma` vanishes at 0 instead of diverging.
`ResponseFunction.MergePolesInInterval(X, Emin, Emax)` merges the poles in a window
into one at their trace-weighted centroid (`MergeBlockListOfPolesInInterval`), which
undoes such a pair faithfully; it does not snap to zero.


## 6. Bath excitations and the mesh methods

### 6.1 `NBathExcitations`

`CreateSpectra(H, T, psi, Hyb, {{"NBathExcitations",p}})` allows `p` particles plus
holes at once in the baths the hybridisation functions stand for: the `p` of Lu, Cao,
Hansmann and Haverkort, PRB 100, 115134 (2019); their `L` is `Nval`/`Ncon` of
`ToOperator`. Since 30a1dfdb.

- `p = 1` without `NonCrossing` goes through `CalculateGforHwithSigmaAsOperator`, as
  always. Otherwise `CalculateGforHwithSigmaAsOperatorBathExcitations` in
  `Core/ResponseFunctionBathExcitations.cpp`: an explicit auxiliary Hamiltonian and a
  block Lanczos, crossed terms included.
- `{"NonCrossing",true}`: the nested approximation, excitations distinguishable. At
  `p = 1` it changes only the route, which is how the tests compare old and new.
- The coupling operators may be any operators; bath statistics follow from their
  fermion parity (`OperatorFermionParity`). The ground state stays at `p = 0`.
- `TestScripts/CreateSpectraWithBathExcitations.Quanty` is exact to 1e-14 against
  explicit baths. On a tutorial-sized impurity `p = 2` changed a peak of 1.05 by 0.72,
  of which the crossed terms were 0.41.

### 6.2 `Method "Mesh"` — residues without a broadening

Three routines put a response function on an energy grid without a Lorentzian:
`CreateSpectra` (`CalculateGforHwithSigmaAsOperatorBathExcitationsMesh`, e5bba26f),
`ResponseFunction.Downfold` (`ResponseFunctionDownfoldOnMesh`, 195f2582) and
`CalculateG` (§2). Each pole of the exact function is split over its two neighbouring
grid points (tent), conserving weight and first moment. The residue uses only `Im G`:
with `A_y(x) = -(1/pi) Im G(x + iy)`, harmonic in `(x, y)`,

    R_k = -(1/dw) int_0^inf y [ A_y(w_{k-1}) - 2 A_y(w_k) + A_y(w_{k+1}) ] dy

on 8-point Gauss–Legendre panels in `ln y`, plus an analytic tail from the moments
`M0`, `M1`, `M2` beyond. Truncating the `y` integral instead loses weight as `1/ymax`.

- The three routines carry their own copies of this quadrature; merging them into one
  helper is planned. Only `TightBindingHamiltonianPlusSelfEnergyOnMesh` has the extra
  node on `[0, ymin]` that keeps the weight of a pole lying exactly on a grid point
  (it lost `2 ymin / (pi h)` without it).
- `CreateSpectra` `"Mesh"` takes `Emin`/`Emax`/`NE` only: uniform, no `Lambda`, no
  window. Its mesh runs at that spacing out to a bound on the spectrum, which can reach
  ±10–13 for a [-4, 4] window and widens as `Sigma` widens; that, not the window, sets
  its cost. With `p > 1` it needs `NonCrossing`.
- `Downfold` `"Mesh"` takes `Emin`/`Emax`/`NE` or a non-uniform `EnergyGrid`, and
  `Lambda`; `G` must have no constant part. Mesh `G` plus mesh downfold converges to
  the exact `Sigma` as `dw^2`. Splitting `G` into `Gp` and `Gm` before downfolding adds
  an `O(dw^2)` error.
- `Lambda` exists because a uniform extension spread a few tiny far poles (2e-4 of the
  weight) over hundreds of grid points, and a cut chain then spent its nodes on them:
  a comb in the Hubbard bands.
- Validation: with every tail energy on a grid point the binning is exact, so mesh
  and Lanczos (both non-crossing) agree to 3e-8. For a discrete toy spectrum binning
  moves poles by up to `dw`, and comparing residues is meaningless.


## 7. Arithmetic the loop relies on

### 7.1 `CalculateGfromG0andSigma` is a substitution, not Dyson

`G(w) = G0(w - Sigma(w))`: every site of the lattice gets `Sigma`. It agrees with the
Dyson equation only when `G0` is a single pole. A constant `Sigma` shifts every pole by
`+a0` and leaves the weights alone; the `Sigma->N == 0` branch subtracted it until
3cfe66b4 (`GListListListRepresentation`; test
`TestScripts/Responsefunction/CalculateGfromG0AndSigma.Quanty`).

### 7.2 Chain plus chain is a chain

Since 528411fc, `ResponseFunctionAdd` adds two `t` as a chain:
`BlockTriDiagonalMatrixAdd` in `Core/BasicMathBlock/BasicMath_BlockTridiagonalMatrix.cpp`
does a Householder band reduction with bulge chasing of the direct sum
(`HermitianBlockBandReduce`, via `BlockChainsSeenThroughV`), then makes every coupling
Hermitian positive semi-definite; on a rank drop it continues at the smaller block size
(`BlockSizeList`). Measured, block size 10: chain length 2000 in 3.0 s against 7.8 s
through poles, 4000 in 14 s against 59 s.

- **Not minimal.** For `G + G`, rounding splits each degenerate pair by ~1e-16 with
  O(1) weight shares, and the coupling where the Krylov space "ends" is 1e-3–1e-2 for
  chains longer than about 10. No threshold can fix that; the chain comes out up to
  `N_A + N_B` long but still evaluates to rounding. Any deflation deep in the chain is
  lost the same way, since whether it occurs at block `k` is itself conditioned like
  `eps (|H|/beta)^k`. Adding decoupled sectors separately would deflate exactly; it is
  deferred.
- **Deep coefficients are ill conditioned.** A symmetry broken at 1e-16 in the head
  reaches ~1e-9 at block 14 (the pole route breaks it at O(1)). `G` itself stays right
  to 1e-15; test chain coefficients at 1e-8. In `TestScripts/Responsefunction/AddTri.Quanty`
  the deflation case passes with a singular value of ~1e-13 against a 2e-12 threshold,
  a margin of ten.
- **Lua `__eq` fires only for table == table or userdata == userdata.** A result
  pushed as a table (capital type letter) compares false with a userdata, silently.
  Arithmetic results therefore use the lower-case type, which `LuaPush_ResponseFunction`
  pushes as a userdata.

### 7.3 What `Downfold` costs

The route is chosen by the input type (`ResponseFunctionDownfold`):

- **List of poles, Anderson, `Nat`, `DoubleTri`** → a pole list, `Sigma` comes back
  an Anderson matrix. Block size 20, 200 poles: 105 s, of which `EffectiveHamiltonian`
  is 98 s; the projection adds about 7 s. `ChangeType(Sigma,"ListOfPoles")` afterwards
  cost another 420 s, which is why the result is left an Anderson matrix.
- **Tri** → `BlockTriDiagonalMatrixDownfoldChain`, by the same band reduction as §7.2
  since 3cb7f0b0 (the block Lanczos survives as `BlockTriDiagonalMatrixDownfoldChainOld`).
  `NTriMax` cuts the chain before the reduction. Block size 20, chain 100, subspace 10:
  exact 2.7 s; `NTriMax` 50 0.67 s (1e-16 from exact), 20 0.10 s (2e-9), 10 0.02 s
  (1e-5).

The chain route is not faster at equal accuracy; it gives a ten times more compact
answer and is the only one that can be truncated. `NTriMax` is the speed knob. Old and
new chain routes agree to 1e-12 except deep in a full 199-block chain, where they part
by up to 0.3 while both reproduce `(R G^-1 R^dag)^-1` to 5e-15: the gauge sensitivity
of deep coefficients again.

### 7.4 `InvertEnergy` is a mirror only for self-adjoint operator pairs

`InvertEnergy` maps `G(mu + w) -> -G(mu - w)^T` (`ResponseFunctionInvertEnergy`). That
is the particle-hole mirror of a propagator whose operators map into their own
adjoints. Particle-hole symmetry does that to `a`, not to `q`: with
`H1 = U n_up n_dn` (no Hartree–Fock part subtracted), `q = [a, H1] = U n_-s a` and

    a -> a^dag        q -> U a^dag - q^dag

so the mirror of the augmented propagator is `M InvertEnergy(Gtilde) M^T` with
`M = {{1,0},{U,-1}}` on the `(a, q)` blocks: a shear, not a sign. Verified at 142
frequencies, 1.2e-10 against 2.1 for the plain `InvertEnergy`, wrong in `aq`, `qa` and
`qq`. Averaging `Gtilde` with its plain `InvertEnergy` changed `Z` from 0.126 to 0.481
on a solution that was already symmetric to 1e-10. **Symmetrise `Sigma` instead**,
after the downfold: built from `a` alone, with `A0 = 0` for the dynamic part,
`InvertEnergy` is the mirror there, and the average is the no-op it should be.


## 8. Large U: the Mott insulator

**The Weiss field.** For a single band with `Sigma(w) = W/w` near zero (`W = U^2/4`)
and a symmetric band of second moment `s2`:

    zeta  = (w^2 - W)/w
    G0^-1 = G^-1 + Sigma = w (w^2 - W - s2)/(w^2 - W)

`G0` has a pole exactly at `w = 0` with residue `W/(W + s2)`, plus two satellites at
`±sqrt(W + s2)`; `Delta = s2 w/(w^2 - W)` vanishes linearly at `mu`: the impurity
decouples. Checked on a simple cubic lattice at `U = 6` (`s2 = 1/24`, `W = 9`):
predicted residue 0.99539, measured 0.9953. Its occupation, 0.5007, is
`f(-2.7e-6/kBT)`: the orbital sits on `mu`, half filled. One-sided weight and an
occupation of 0.5 are both right.

**The particle-hole flat direction.** Once the solution is gapped, nothing pins the
centre of the gap. Measured in the converged insulator: `dn/d(dE) = -3.9e-4`, so a
filling accuracy of 1e-4 pins the double counting `dE` only to ±0.26. Moving `dE`
inside the gap translates the spectrum rigidly and does not dope: it is a gauge. The
gauge is not harmless, because it decides whether the big Weiss pole at `mu` is
counted valence or conduction (§4). With the subtraction route, `Sigma`'s Mott pole
drifted 1e-14 → 1e-12 → 1e-9 → 5e-6 → 7e-3, a factor of about 1e3 per iteration; the
chain emptied when the drift passed `EffectivekBT`, and the spectrum never converged
(10–90 % asymmetric). Ruled out by test: the double counting, `kBT`, the mesh, the
`NSigma` cut, the bath size, and the `Method` of `CalculateGbathfromGandSigma`.

**Cures, measured on the same `U = 6` run, ten iterations:**

| route | asymmetry | in-gap weight | symmetry assumed |
|---|---|---|---|
| subtraction (§3.1) | 1.27 | 0.053 | none |
| composition (§3.2) | 0.0042 | 0.00000 | none, but scalar `zeta` only |
| downfold (§3.3) | 0.0092 | 0.00000 | none |
| symmetrising `Sigma` (§7.4) | 0.0000 | — | particle-hole |

With the downfold the Weiss occupation converged to 0.4999, `dE = -3.000000` every
iteration, no negative-weight warnings. A lattice that is not particle-hole symmetric
has the same flat direction in the insulator; there `dE` needs a condition other than
the filling. That is open.


## 9. Converging a loop: checks, settings, failure modes

From one- and two-band loops; they hold for anyone using these functions.

- **Restrictions per group, the same for ground state and spectra.** After the
  rotation by `T` the cluster has four groups: `Impurity` (full at `U = 0`), `BathExt`
  (its empty partner), `BathVal`, `BathCon`. Restrict holes in the full groups and
  electrons in the empty ones per group. Too tight lost up to 10 % of the spectral
  weight at large U. A larger space for the excited states than for the ground state
  puts weight on the wrong side of `mu` (`E(N-1)` fell 0.3 below `E(N)`).
- **Print three checks every iteration**; they are cheap and change nothing:
  `E(N-1) - E(N)` and `E(N+1) - E(N)` of the cluster at `mu = 0` (both positive);
  the particle and hole weight of the impurity spectrum against `psi`, and the part
  below zero; the density matrices of lattice, `psi` and `Gtilde` (equal at
  self-consistency). A charge gap going to zero warned one to three iterations before
  failure. Start the `N ± 1` searches from all determinants of that sector; a single
  configuration can land in the wrong symmetry sector.
- **Energy window**: the band plus `U` on each side, not `U/2`. `Sigma` has weight
  beyond the Hubbard peaks of `G`; with too small a window most of it (63 % in one
  case) is lumped onto a few wing poles, which put sharp satellites into `G` and widen
  the range every iteration.
- **Broadening on a mesh** must exceed the spacing (about 1.6 times) or the mesh shows
  as a comb. A coarse spacing breaks the filling near `mu`; a large U wants a grid fine
  near `mu` and coarse far out, which `CalculateG` and `Downfold` take (`EnergyGrid`)
  and `CreateSpectra` `"Mesh"` does not (its Lanczos route with `GridPoints` does).
- **Chemical potential in a gap**: bisection on the filling lands on a gap edge; take
  the middle of the plateau.
- **Lattice and impurity filling** are fixed by different things. `dE` pins the
  impurity; nothing pins the lattice, which drifted 2.00 → 2.21 in a two-band case.
  Adding a constant to `Sigma` so that the lattice keeps its filling held it.
- **Large U is not solved** (in a two-band chain, `U` at twice the band width and
  above): the impurity filling jumps with `dE` (a level crossing of the fixed-`N`
  cluster), the cluster prefers `N - 1` with the hole in the impurity's bath partner,
  and a Mott start loses its gap within a few iterations. Tried without success:
  pinning the filling, centring the charge gap, an `E0` correction, mixing 0.25–0.5,
  `EffectivekBT` 0.1, a logarithmic grid.


## 10. Measuring and testing

- **Particle-hole asymmetry**: use the relative L1 measure
  `int |A(w) - A(-w)| dw / int A dw`. A pointwise maximum is useless next to sharp
  Hubbard bands: a shift of 0.04 already gives `max|dG| ~ 3`.
- **Self-energy**: compare the quasiparticle weight `Z = 1/(1 - dRe Sigma/dw|_0)`, not
  `max|dSigma|` on a mesh; a near-divergent pole makes the latter report only how far
  that pole moved.
- **Binned functions agree only to the binning error.** The downfold and
  `zeta - G^-1 - eps` are the same function binned differently: 2.4e-2 apart at
  `Gamma = 0.05`, 4e-5 at `Gamma = 1`. Compare at a broadening well above the spacing.
- **Trust `CalculateRho`** for occupations over a hand pairing of `ToTable` entries.
- Tests that pin the pieces: `TestScripts/TightBindingDef/` `CalculateG`,
  `CalculateGMesh`, `CalculateGSigmaTypes`, `CalculateGSigmaPerAtom`,
  `HybridizationFunction`, `CalculateGAgainstCreateSpectra`;
  `TestScripts/Responsefunction/` `Downfold`, `DownfoldMesh`, `AddTri`, `InvertEnergy`,
  `CalculateHybridizationFunction`, `CalculateGBathFromGAndSigma`,
  `CalculateGfromG0AndSigma`, `MergePolesInInterval`, `PlaceOnMesh`;
  `TestScripts/CreateSpectraWithSelfEnergy.Quanty` and
  `TestScripts/CreateSpectraWithBathExcitations.Quanty`. `CalculateG.Quanty` (about
  10 minutes) and `HybridizationFunction.Quanty` (about 3) are slow on a laptop; the
  CI runs them.


## 11. Open points

- A window option for `CreateSpectra` `"Mesh"` (compute only `[Emin, Emax]`, keep `M0`
  and `M1` exact by collapsing the outside weight into a few poles) is agreed and not
  built.
- One helper for the three copies of the mesh quadrature (§6.2).
- Removing the `O(dw^2)` error of splitting `G` before a mesh downfold needs `Gp` and
  `Gm` at the same `z`.
- A compression step for chains whose sum is not minimal (§7.2).
- Positivity of the weights of the Schur complement in §3.3 is expected (`Aug` is a
  Gram matrix of resolvent elements) but not checked; nor is the case where `Sigma`
  lives on a correlated sub-block only.
- The double counting of a correlated insulator without particle-hole symmetry (§8),
  and the large-U failure modes of §9.
- Near-zero pole weights make `CompactMatrixSqrt` print "Negative value in
  CompactMatrixSqrt" for rounding-level negative eigenvalues, for any augmented
  propagator. A weight threshold would cure it; `MinimalPoleWeigth` is refused by every
  `ChangeType` conversion today.
