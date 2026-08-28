# Plan: response-function type `d`/`D` — double tri-diagonal

Goal: `ResponseFunction.ToOperator{Nv=…, Nc=…}` on a response function held as a *double*
tri-diagonal matrix — one impurity site with a valence chain and a conduction chain hanging
off it — so `Nv` and `Nc` truncate the two chains independently and the two tails come back
as `GvalBath` / `GconBath`.

Scope decided: **the single-valued (non-block) type first.** The block type
(`BlockDoubleTriDiagonalMatrixType`) comes after, as a Block-style reimplementation of the
`Blok…` function list — see §8.

---

## 0. Status as of `b57cd101` (2026-08-28)

The type exists, evaluates, prints, converts to a dense matrix, produces spectra
and can now be **built** from any of the other four non-block formats, all
verified against each other. What is missing is the conversions *out* of `d`,
the block conversions, and the operator.

### The `-> d` column of `ConvertResponsefunction` is done (non-block)

`ResponseFunction.ChangeType(G,"DoubleTri")` works from `a`, `t`, `l` and `n`,
and reproduces the original response function to 1e-13 relative.
`TestScripts/Responsefunction/ChangeTypeDoubleTri.Quanty` is the test: a 999
level bath at finite `kBT`, at `kBT = 0`, with `mu` off centre, with `mu` outside
the band on either side so one chain comes out empty, a six level bath, a single
bath level, a bath with threefold degenerate levels, and a table round trip.

`a -> d` is the real one; `t -> d`, `l -> d` and `n -> d` route through it.
Note that `l -> d` **cannot** split the poles at `mu` and tri-diagonalize each
half, which is what §4 item 1 used to say: that is the `n` recipe and it works
there only because `n` adds its two chains where `d` nests them. A list of poles
is the function itself, so it first has to become the chain whose spectrum it is,
i.e. the Anderson form.

### The chain construction is shared, and no longer O(N^3)

`PolesToTriDiagonalCoefficients` in
[BasicMath_MatrixTypeConversion.cpp](../../Core/BasicMath/BasicMath_MatrixTypeConversion.cpp)
builds the Jacobi matrix of a discrete measure by the method of Rutishauser,
Kahan, Pal and Walker (Gragg and Harrod, Numer. Math. 44 (1984) 317). It is used
by `AndersonMatrixToTriDiagonalMatrix`, `ListOfPolesToTridiagonalMatrix` and both
chains of `AndersonMatrixToDoubleTriDiagonalMatrix`.

It replaces Lanczos with full re-orthogonalization, which needed the whole
`N x N` basis. Measured in Quanty on `l -> t`:

| poles | before | after |
|---|---|---|
| 500 | 0.054 s | 0.002 s |
| 1000 | 0.530 s | 0.007 s |
| 2000 | 5.46 s | 0.028 s |

at O(N) instead of O(N^2) memory. The coefficients agree with the old routine to
1e-14 relative; the resulting `G` is accurate to ~1e-13 where the old one reached
~1e-14, which is the honest cost of not re-orthogonalizing. Plain Lanczos without
re-orthogonalization, for contrast, is wrong by 1e-2 at 1000 poles.

Two guards in that routine are load bearing rather than cosmetic: poles at the
same energy are **merged**, and the chain is capped at the number of distinct
poles carrying weight. RKPW is not robust against a measure whose support is
smaller than the number of poles handed to it -- on a measure with every node
exactly doubled it is wrong by 0.37 relative where re-orthogonalized Lanczos is
right to 1e-14 -- and degenerate bath levels are common.

### Conventions §3.0 and §3.1 are settled on this path

`AndersonMatrixToDoubleTriDiagonalMatrix` writes `val.a[0] = con.a[0] =
M.a[1]/2`, so their **sum** is the impurity energy, and it splits with a true
**Fermi function**, storing `kBT` in the struct. The old code wrote the full
`M.a[1]` in both, which is a relative error of 190 on `G` for an impurity energy
of 1.7 -- invisible only because every existing test has a zero impurity energy.
The dummy bath site the old code invented for an empty branch is gone: an empty
branch is `val.N = 0` with `val.b[0] = 0`, which is exact.

The other consumers listed in §3.0 still read `(val.a[0] + con.a[0]) / 2` and are
now inconsistent with what the constructors produce. That sweep is the next thing
to do, and it is no longer optional.

### Working and verified

| | |
|---|---|
| struct, both forms | `DoubleTriDiagonalMatrixType` (with `kBT`), `BlockDoubleTriDiagonalMatrixType` |
| type letter | `d`/`D`; Lua `"DoubleTri"` / `"DTri"` |
| Lua table ⇄ struct | read and write, `kBT` round-trips, `ToTable`, `New`, `ChangeType(…,"Table"/"Userdata")`, `print` |
| `Evaluate` | non-block and block, checked against an independent continued fraction to ~1e-15 |
| dense matrix | all four `…ToCompactMatrix[WithDeflation]`, checked against `Evaluate` for real, complex and mixed, `BlockSize` 1 and 2, and empty chains |
| spectra | `ToSpectra` with broadening, non-block and block, agrees with `Evaluate` to the last digit |
| accessors | `GetMu`, `GetQComplex`, `GetBlockSize`, `GetA0`, `GetB0`, block and non-block |
| `Free` | `FreeResponsefunction`, `FreeBlockMatrixFromLua` |

### The conversions that are still stubs

Non-block (`BasicMath_MatrixTypeConversion.cpp`), all of them *out* of `d`:
`DoubleTriDiagonalMatrixToListOfPoles`,
`DoubleTriDiagonalMatrixToAndersonMatrixMatrix`,
`DoubleTriDiagonalMatrixToBADoubleTriDiagonalMatrix`.

`ChangeType` now refuses these with a catchable lua error instead of pushing a
response function that was never converted: `ConvertResponsefunction` sets
`Mout->data` to `NULL` up front and returns non-zero if no cell filled it in.
Before that, `ChangeType(Gd,"Tri")` pushed an uninitiallized pointer.

Block (`BasicMath_BlockTypeConversion.cpp`):
`BlockListOfPolesToBlockDoubleTriDiagonalMatrix`,
`BlockDoubleTriDiagonalMatrixToBlockListOfPoles`,
`BlockAndersonMatrixToBlockDoubleTriDiagonalMatrix`,
`BlockDoubleTriDiagonalMatrixToBlockAndersonMatrix`,
`BlockTriDiagonalMatrixToBlockDoubleTriDiagonalMatrix`,
`BlockDoubleTriDiagonalMatrixToBlockTriDiagonalMatrix`,
`BlockDoubleTriDiagonalMatrixToBlockNaturalImpurityMatrix`,
`BlockNaturalImpurityOrbitalMatrixToBlockDoubleTridiagonalMatrix`, and
`BlockDoubleTriDiagonalMatrixToCompactMatrix` **by value** — see the trap below.

Also stubs, silently returning 0: `BlockDoubleTriDiagonalMatrixChop` and
`…ChopReduceDimension`. A chop that does nothing is harmless in a way a
conversion that does nothing is not, so those were left; they want a comment
saying the no-op is deliberate.

### A trap worth removing

There are **two overloads** of `BlockDoubleTriDiagonalMatrixToCompactMatrix`:

- `(BlockDoubleTriDiagonalMatrixType const M, CompactMatrixType * const H)` —
  by value, a TODO stub at `BasicMath_BlockTypeConversion.cpp:2416`, declared in
  that header;
- `(BlockDoubleTriDiagonalMatrixType const * A, CompactMatrixType * const B)` —
  by pointer, the working implementation at
  `BasicMath_CompactMatrixVector.cpp:5312`.

The live call site passes a pointer and so gets the working one, but a caller who
dereferences silently gets the stub. Delete the by-value declaration and stub.

### Dispatch sites still without a `d` case

`ResponseFunctionToDensityMatrix`, `ResponseFunctionEq`, the three
`ResponseFunctionAddSelfEnergy…`, `ResponseFunctionInvertEnergy`,
`ResponseFunctionTranspose`, `ResponseFunctionEffectiveHamiltonian`,
`ResponseFunctionFromEffectiveHamiltonian`, `ResponsefunctionChop`,
`ResponsefunctionMul`, `FreeBlockMatrix`, `BlockResponseFunctionRotate`,
`LuaResponseFunctionEigensystem`, `ResponseFunctionToLuaBitMask`,
`LuaResponseFunctionToTightBinding`, `LuaResponseFunctionToOperator`,
`LuaResponseFunctionLen`, `LuaResponseFunctionIndexMu`, `…IndexType`,
`…IndexName`, `…NewIndexMu`, `…NewIndexName`, `lua_isBlockGreensfunction`.

Two of those bite today: `FreeBlockMatrix` leaks and complains for a block `d`
(it is called from `FreeSelfEnergyArrayFromLua` and `LuaImpurityGreensFunction`),
and `lua_isBlockGreensfunction` rejects `"DoubleTri"`, which silently excludes
`d` from the self-energy reader. The six index metamethods mean `G.type`, `G.mu`
and `G.name` do not work on a `d` userdata.

### Still open from §2, §3.0 and §3.1

- §2(b) `DoubleTriDiagonalMatrixInverse` builds the conduction block from
  `M.val`, and uses `MatrixType`.
- §2(c) `DoubleTriDiagonalMatrixAdd` loop bound, and the `uint32_t` underflow at
  `val.N == 0`.
- §3.0, the averaged head. `DoubleTriDiagonalMatrixToSpectrum` is fixed, but the
  four matrix–vector products and `…MatrixInverse`
  (`BasicMath_DoubleTriDiagonalMatrixVector.cpp:137, 266, 298, 339, 367`) and the
  eight sites in the two `…ToAndersonImpurityOperator` builders
  (`MatrixOperator.cpp:635…840`) still read `(val.a[0] + con.a[0]) / 2`.
- §3.1, the split rule. `kBT` replaced `Epsilon` in `DoubleTriDiagonalMatrixType`
  and `BlockDoubleTriDiagonalMatrixType`, but
  `AndersonMatrixToDoubleTriDiagonalMatrixPlusT:338` still uses the linear ramp,
  and the other five structs still carry `Epsilon`.

### Not testable yet

The deflated block path. Nothing can build a deflated block `d`: the Lua reader
derives `BlockSize` from the matrices and takes no per-block list, and both
`BlockListOfPolesToBlockDoubleTriDiagonalMatrix` and
`BlockDoubleTriDiagonalMatrixChopReduceDimension` are stubs. `WithDeflation` is
verified only to agree with the plain version when there is nothing to deflate.

### Tooling notes, learned the hard way

- **Quanty Lua options are a list of `{"name", value}` pairs.** A table with
  named keys (`{Emin=-4}`) has an empty array part, passes `LuaTestOptionsList`
  and is then silently ignored. Every option in this document must be written
  `{{"Emin",-4},{"GammaL",0.3}}`.
- Piped stdout is block buffered and the buffer is lost when the process dies, so
  a crash looks silent. Re-run under `script -q /dev/null` to see the last
  `printf`.
- `Matrix.Sub` is 1-based. `Complex.Abs` and `Complex.Conjugate` do not exist.
  `Spectra.ToTable(sp)` gives `sp[ispectrum][ipoint] = {energy, value}`.
- `ToSpectra` broadening is `GammaL`, not `Gamma`.

---

## 1. The two structs combine differently — this is the whole point

`DoubleTriDiagonalMatrixType` and `BADoubleTriDiagonalMatrixType` have near-identical
members (`mu`, `kBT` resp. `Epsilon`, `a0`, `b0`, `TriDiagonalMatrixType val`, `con`) and
completely different meanings. Verified from the code, not from the comments:

### `n` — `BADoubleTriDiagonalMatrixType`, **additive**

`EvaluateBADoubleTriDiagonalMatrix`
([BasicMath_GreenFunctionArithmetic.cpp:4071](../../Core/BasicMath/BasicMath_GreenFunctionArithmetic.cpp#L4071))
and `BADoubleTriDiagonalMatrixToSpectrum`
([Spectral.cpp:1609](../../Core/Spectral.cpp#L1609)):

```
G(w) = a0 + b0² · [ g_val(w) + g_con(w) ]

g_val(w) = val.a[0] + val.b[0]² / (w − val.a[1] − val.b[1]²/(w − val.a[2] − …))
```

Two *independent* chains; `val` and `con` are each a standalone `TriDiagonalMatrixType` with
the standard convention (`a[0]` = additive offset, `b[0]` = weight). The observed orbital is
a linear combination of the two chain heads — that is what `T` in the struct comment
([StructDef.cpp:586](../../Core/StructDef.cpp#L586)) does, and what the cached
`acvB`/`acvA`/`acvBA`/`bvB0`/`bcB0`/`bvA0`/`bcA0` fields hold.

### `d` — `DoubleTriDiagonalMatrixType`, **nested**

From `RealDoubleTriDiagonalMatrixVector`
([BasicMath_DoubleTriDiagonalMatrixVector.cpp:263](../../Core/BasicMath/BasicMath_DoubleTriDiagonalMatrixVector.cpp#L263)),
`AndersonMatrixToDoubleTriDiagonalMatrixPlusT`
([BasicMath_MatrixTypeConversion.cpp:277](../../Core/BasicMath/BasicMath_MatrixTypeConversion.cpp#L277))
and `DoubleTriDiagonalMatrixToAndersonImpurityOperator`
([MatrixOperator.cpp:741](../../Core/MatrixOperator.cpp#L741)):

```
G(w) = a0 + b0² / ( w − G.val(w) − G.con(w) )

G.val(w) = val.a[0] + val.b[0]² / (w − val.a[1] − val.b[1]²/(w − val.a[2] − …))
G.con(w) = con.a[0] + con.b[0]² / (w − con.a[1] − con.b[1]²/(w − con.a[2] − …))
```

`val` and `con` each hold *one bath branch as seen from the impurity*, in the ordinary
tri-diagonal convention: `a[0]` a site energy, `b[0]` the hopping into that branch,
`a[1…N]` / `b[1…N-1]` the chain. `G.val` and `G.con` are therefore each an ordinary
type-`t` response function, evaluated by the existing `TriDiagonalMatrixToSpectrum` /
`EvaluateTriDiagonalMatrix` with no special casing, and `d` composes them by putting **both
in the same denominator**.

The impurity site energy is `val.a[0] + con.a[0]` — the plain **sum**, no half. Matrix
dimension `1 + val.N + con.N`.

Site layout of the corresponding Hamiltonian:

| site | energy | hopping to previous |
|---|---|---|
| 0 (impurity) | `val.a[0] + con.a[0]` | — |
| 1 … val.N | `val.a[1…val.N]` | `val.b[0]` then `val.b[1…val.N-1]` |
| val.N+1 … val.N+con.N | `con.a[1…con.N]` | `con.b[0]` then `con.b[1…con.N-1]` |

Contrast `n`, which adds the two *Green's functions* rather than the two hybridisations.

### On `a0` / `b0`

`a0` is a site energy and `b0` a hopping — at every level, in all of these types. That is
what they mean whenever a response function holds an effective Hamiltonian, a self-energy or
a hybridisation function rather than `G` itself: `G = 1/(w − Heff)`, `Heff = a0 + b0²·g`, so
`a0` is the static on-site part and `b0²·g` the hybridisation. The ASCII art in each struct
comment is the auxiliary matrix whose resolvent is the Green's function; `a0`/`b0` are its
site-0 energy and its first hopping.

Read that way there is nothing irregular about `d`: the head of `val` is a site energy and a
hopping like everywhere else. What is specific to `d` is only that the *same* site is the
head of both branches, so its energy is split across the two head slots.

Maurits is writing these formulas into the struct comments in
[StructDef.cpp](../../Core/StructDef.cpp); that file is the definitive statement of intent,
and this document follows it.

### Decision: the shared site is `val.a[0] + con.a[0]`, and there is no separate member

Two candidates were considered:

1. `G = a0 + b0²/(w − G.val − G.con)`, impurity energy `= val.a[0] + con.a[0]`. **Chosen.**
2. `G = a0 + b0²/(w − a1 − G.val − G.con)` with a new member `a1`, requiring
   `val.a[0] == con.a[0] == 0`.

Option 1 wins because it introduces no *reserved-zero* slot, hence no invariant for future
code to violate. Every slot is live, and every existing tri-diagonal routine that reads
`a[0]` as an additive offset contributes the right term once the result lands in the
denominator.

Option 2 is the design the block NIO type already has — an outer `a0` plus two head slots
that are supposed to be zero — and it is measurably fragile there: the invariant is enforced
in only one place, and three groups of functions disagree about it. See §10 for the evidence.

Consequences, all in our favour:

- `DoubleTriDiagonalMatrixToSpectrum` becomes correct by *deleting* a term (§2a).
- `G.val` / `G.con` remain self-contained response functions, so `ToSpectra`,
  `ToListOfPoles`, `Part` etc. apply to them unchanged.
- The `/2` that currently appears in eight places moves into the constructors, where it is
  written once and visibly (§3.0).

Cost, stated plainly: "the impurity site energy" is a sum of two fields, and how it is split
between them is a gauge with no canonical value. Normalise to equal halves in the
constructors — as a *convention*, not an invariant, so nothing breaks silently if a caller
shifts weight from one to the other.

### Consequence

Every conversion and every action must be written fresh for `d`. Nothing can be copied from
the `n`/`BADouble…` or `BlockNaturalImpurityOrbital…` twins beyond the function signature
and the file it lives in. (My earlier note claiming `n` is `d` plus a head rotation was
wrong: `n` sums two continued fractions, `d` nests them in one denominator. Different
functions, not a change of basis on the same one.)

Where it *doesn't* matter: the two structs lay out the chains identically —
`a[1…N]` energies, `b[1…N-1]` hoppings, `a[0]`/`b[0]` at the head — so truncation and tail
extraction are shared code. `DoubleTriDiagonalMatrixPart(M, A, iEndval, iEndcon)` truncates,
and `TriDiagonalMatrixPart` / `TriDiagonalMatrixPartAddSpin` extract a tail, unchanged; see
§3.

---

## 2. Three defects in the existing `d` functions — (a) fixed, (b) and (c) open

The `DoubleTriDiagonalMatrixType` family is dormant — `DoubleTriDiagonalMatrixAdd`,
`…ScalarMultiply`, `…Copy`, `…Chop`, `ComplexDoubleTriDiagonalMatrixVector` and
`MatrixToDoubleTriDiagonalMatrixNoOverlap` have no callers at all, and the only consumer
outside `Test/` is [DMFT/SingleSiteXAS.cpp](../../DMFT/SingleSiteXAS.cpp), which still
`#include`s `../Core/StructDef.c` and therefore no longer builds. So these are cheap to fix
— but they cannot be reused as-is:

**(a) `DoubleTriDiagonalMatrixToSpectrum` — FIXED in `23018386`.** It double-counted the
shared impurity site: `GV.Re` from `TriDiagonalMatrixToSpectrum` already contains `val.a[0]`
and `GC.Re` contains `con.a[0]`, so `GV + GC` is the whole `w`-denominator and the extra
`avc0` term had to go, from both the `denom` and the `G->Re` line. Verified against
`Evaluate` to the last digit. Kept here because it is the reference case for the same
mistake elsewhere: **any** consumer that adds `val.a[0] + con.a[0]` on top of an evaluated
`val`/`con` is counting the impurity twice.

**(b) `DoubleTriDiagonalMatrixInverse`**
([BasicMath_DoubleTriDiagonalMatrixVector.cpp:125](../../Core/BasicMath/BasicMath_DoubleTriDiagonalMatrixVector.cpp#L125))
builds the conduction block from `M.val` instead of `M.con` — `MI->VRe[i][0] = M.val.b[0]`,
`MI->VRe[i][i] = M.val.a[1]`, and `M.val.a/b` throughout the con loop. A copy-paste slip.
It also uses `MatrixType`, which should be `CompactMatrixType`.

**(c) `DoubleTriDiagonalMatrixAdd`**
([…cpp:159](../../Core/BasicMath/BasicMath_DoubleTriDiagonalMatrixVector.cpp#L159))
loops `i < M->val.N - 1`, so `val.a[val.N-1]` and `val.b[val.N-1]` are never added; with
`val.N == 0` the `uint32_t` underflows into a runaway loop. It also ignores `a0`/`b0`.

Related, lower stakes: `DoubleTriDiagonalMatrixScalarMultiply` scales every `a` and `b` and
then *resets* `a0 = 0`, `b0 = 1` — scaling the matrix rather than the Green's function, and
discarding the offset. Which of the two is intended needs deciding before it is used from
`ResponsefunctionMul`.

Also worth knowing: `ResponseFunctionToSpectra`
([Core/ResponseFunction.cpp:383](../../Core/ResponseFunction.cpp#L383)) already casts the
non-block `case 'n'` payload to `DoubleTriDiagonalMatrixType*` when it should be
`BADoubleTriDiagonalMatrixType*`. It reads only `->Name`, the first member of both structs,
so it works by luck today — but it will read as intentional once `d` exists, so it is worth
correcting while touching that switch.

---

## 3. Struct and core algebra (non-block)

The struct itself needs no change, and the §1 formulas are now in its comment in
[StructDef.cpp](../../Core/StructDef.cpp#L529). Of the lists below, items 5 (`Evaluate`) and
the `…ToCompactMatrix` pair are done; the two sweeps §3.0 and §3.1 are not. See §0.

### 3.0 Convention sweep: the eight `/2` sites

The struct itself needs no new member under option 1, but the current code encodes the *old*
average reading — `(val.a[0] + con.a[0]) / 2` — in eight places. Each becomes the plain sum,
i.e. the `/ 2` is deleted:

| # | function | file:line |
|---|---|---|
| 1 | `DoubleTriDiagonalMatrixToSpectrum` | [Spectral.cpp:1583–1584](../../Core/Spectral.cpp#L1583) — delete the whole term, see §2a |
| 2 | `DoubleTriDiagonalMatrixInverse` | [BasicMath_DoubleTriDiagonalMatrixVector.cpp:137](../../Core/BasicMath/BasicMath_DoubleTriDiagonalMatrixVector.cpp#L137) |
| 3 | `RealDoubleTriDiagonalMatrixVector` | [:266](../../Core/BasicMath/BasicMath_DoubleTriDiagonalMatrixVector.cpp#L266) |
| 4 | `ComplexDoubleTriDiagonalMatrixVector` | [:298–299](../../Core/BasicMath/BasicMath_DoubleTriDiagonalMatrixVector.cpp#L298) |
| 5 | `RealVectorDoubleTriDiagonalMatrixVector` | [:339](../../Core/BasicMath/BasicMath_DoubleTriDiagonalMatrixVector.cpp#L339) |
| 6 | `ComplexVectorDoubleTriDiagonalMatrixVector` | [:367](../../Core/BasicMath/BasicMath_DoubleTriDiagonalMatrixVector.cpp#L367) |
| 7 | `DoubleTriDiagonalMatrixToAndersonImpurityOperator` | [MatrixOperator.cpp:765, 790, 815, 840](../../Core/MatrixOperator.cpp#L765) |
| 8 | `SecondDoubleTriDiagonalMatrixToAndersonImpurityOperator` | [MatrixOperator.cpp:635, 660, 685, 710](../../Core/MatrixOperator.cpp#L635) — operates on the wrapped `Dt` |

`DoubleTriDiagonalMatrixToTriDiagonalMatrix`
([BasicMath_MatrixTypeConversion.cpp:1668](../../Core/BasicMath/BasicMath_MatrixTypeConversion.cpp#L1668))
does not read the head directly — it goes through `RealDoubleTriDiagonalMatrixVector`
(line 1697) — so item 3 fixes it too.

And the *producing* side must write halves. `AndersonMatrixToDoubleTriDiagonalMatrixPlusT`
([:277](../../Core/BasicMath/BasicMath_MatrixTypeConversion.cpp#L277)) currently writes
`H->val.a[0] = MVTri.a[1]` and `H->con.a[0] = MCTri.a[1]`, both equal to the full impurity
energy; under option 1 each takes half of it. Net effect on every consumer is unchanged —
this is a re-partition of where the factor lives, not a change of physics — which is why it
can land as one commit with a round-trip test either side.

`Test/LanczosDyson.cpp:290` also hard-codes the average and will need the same edit; it is
the only place outside `Core/` that does.

Note the two `…ToAndersonImpurityOperator` functions are DMFT-specific two-spin + `U`
builders, not the general `…ToImpurityHamiltonian` of §6.1 — but they are the working
reference for the site layout, so keeping them consistent is worth the four-line change.

### 3.1 Convention sweep: `Epsilon` → `kBT`, linear ramp → Fermi function

Decided: the valence/conduction split uses a **true Fermi–Dirac function**, and the
temperature it was built with is **stored in the struct** — the representation genuinely
depends on it, so it is real provenance, not decoration. Scope for now is
`DoubleTriDiagonalMatrixType` only; the other five structs carrying the field follow later
(§11).

**Rename the member.** `double Epsilon` → `double kBT` in
`DoubleTriDiagonalMatrixType`, and the corresponding `epsilon` parameters of the
constructing functions → `kbT`. The member is currently write-only (§11), so nothing reads
it that could break.

**Replace the ramp with the Fermi function.** The live splitting code in
`AndersonMatrixToDoubleTriDiagonalMatrixPlusT`
([BasicMath_MatrixTypeConversion.cpp:277](../../Core/BasicMath/BasicMath_MatrixTypeConversion.cpp#L277))
is a linear ramp over `[mu − epsilon, mu + epsilon]`:

```c
MV.b[MV.N-1] = M.b[i-1] * sqrt(0.5 - (M.a[i] - mu) / (DBL_MIN + 2*epsilon));
MC.b[MC.N-1] = M.b[i-1] * sqrt(0.5 + (M.a[i] - mu) / (DBL_MIN + 2*epsilon));
```

becomes

```c
f = FermiFunction(M.a[i], mu, kbT);
MV.b[MV.N-1] = M.b[i-1] * sqrt(f);
MC.b[MC.N-1] = M.b[i-1] * sqrt(1.0 - f);
```

`b` is a hopping here, so the weight is `b²` and the `sqrt` is required — same as
`AndersonMatrixToBADoubleTriDiagonalMatrixTemperature`
([:1002](../../Core/BasicMath/BasicMath_MatrixTypeConversion.cpp#L1002)), which is the
existing reference implementation of exactly this.

Two properties of `FermiFunction`
([BasicMath_StandardFunctions.cpp:134](../../Core/BasicMath/BasicMath_StandardFunctions.cpp#L134))
matter and are worth relying on deliberately:

- it returns **exact** `1.0` for `E − mu <= −40·kbT` and **exact** `0.0` for
  `E − mu >= 40·kbT`. So the "hard zero outside a finite window" property that motivated the
  linear approximation is retained; the window is just `±40 kbT` instead of `±epsilon`. The
  existing three-way branch (val only / con only / both) therefore stays, with the cutoff
  test becoming `f == 1.0` / `f == 0.0`;
- `kbT == 0` is safe — it returns `0.5` at `E == mu` (checked first) and then a clean step.
  So the `T = sqrt(DBL_MIN)` default in `LuaChangeTypeResponseFunction`
  ([ResponseFunction.cpp:33](../../LuaInterface/ResponseFunction.cpp#L33)) can simply
  become `0`.

**Blast radius is nil.** Outside
[BasicMath_MatrixTypeConversion.cpp](../../Core/BasicMath/BasicMath_MatrixTypeConversion.cpp)
itself, the `d` constructors are called only from `DMFT/SingleSiteXAS.cpp` (which no longer
builds — it `#include`s `../Core/StructDef.c`) and from `Test/TestBasicMath.cpp`,
`Test/TestGreenFunctions.cpp`, `Test/Old_Functions.cpp`. There is no live production caller,
so both the rename and the rule change are free.

**Everything new follows the same rule.** `ListOfPolesToDoubleTriDiagonalMatrix` (§4 item 1)
uses `FermiFunction` on the pole residues directly — no `sqrt`, since a `ListOfPoles` `b` is
already a weight — exactly as `CopyBlockListOfPolesVal` / `…Con` do. `Copy` and `Part` carry
`kBT` across unchanged. `d → n` and `n → d` (§4 items 5–6) take `kBT` from the source
struct rather than as an argument, which is the point of storing it.

### Reuse as-is

`InitDoubleTriDiagonalMatrix`, `InitIdentityDoubleTriDiagonalMatrix`,
`FreeDoubleTriDiagonalMatrix`, `DoubleTriDiagonalMatrixCopy`,
`DoubleTriDiagonalMatrixPart`, `DoubleTriDiagonalMatrixChop` — untouched by the convention
change; they copy or test elements without interpreting the head. The four matrix–vector
products are correct apart from the one-line `/2` edit of §3.0.

`…Part(M, A, iEndval, iEndcon)` is the `Nv`/`Nc` truncation primitive.

**Tail extraction needs no new code either.** `TriDiagonalMatrixPart(M, A, iStart, iEnd)`
([BasicMath_TriDiagonalMatrixVector.cpp:214](../../Core/BasicMath/BasicMath_TriDiagonalMatrixVector.cpp#L214))
sets `M->a[0] = A.a[iStart-1]`, `M->b[0] = A.b[iStart-1]` and then copies
`A.a[iStart…iEnd]` — i.e. it emits a chain whose head is *the last site before the cut and
its coupling into the tail*, which is exactly the `a0`-is-a-site-energy /
`b0`-is-a-hopping convention of §1. Since `d` and `n` lay the chain out identically
(`a[1…N]`, `b[1…N-1]`), the call the `n` path already uses —
`TriDiagonalMatrixPart(tail, val, Nv+1, MNv)`, and `TriDiagonalMatrixPartAddSpin` for
`AddSpin` — is correct for `d` unchanged. Only the *edge cases* differ and need checking:
`Nv == 0` (no valence bath, but the impurity still carries `val.a[0]`) and `Nv >= val.N`
(nothing to extract — the existing `MNv > Options.Nv` test already handles this).

### Rewrite

| # | function | why |
|---|---|---|
| 1 | `DoubleTriDiagonalMatrixAdd` | §2(c) |
| 2 | `DoubleTriDiagonalMatrixScalarMultiply` | semantics of `a0`/`b0` |
| 3 | `DoubleTriDiagonalMatrixToSpectrum` | §2(a) |
| 4 | `DoubleTriDiagonalMatrixInverse` → `DoubleTriDiagonalMatrixTypeToCompactMatrix` | §2(b), plus `MatrixType` → `CompactMatrixType` and the §3.0 edit |

### Write new

| # | function | model / note |
|---|---|---|
| 5 | `EvaluateDoubleTriDiagonalMatrix` (typed + `void const*` overload) | the **nested** form — contrast `EvaluateBADoubleTriDiagonalMatrix` |
| 6 | `DoubleTriDiagonalMatrixTypeToEigenSystemCompact` | cf. `BADoubleTriDiagonalMatrixTypeToEigenSystemCompact` |
| 7 | `DoubleTriDiagonalMatrixToEigenvalues` | cf. `TriDiagonalMatrixToEigenvalues` |
| 8 | `SymmetrizeDoubleTriDiagonalMatrix` | cf. `SymmetrizeBADoubleTriDiagonalMatrix`; must keep `val.a[0] == con.a[0]` |
| 9 | `ShiftDoubleTriDiagonalMatrix` | shifts `val.a[*]`, `con.a[*]`, `mu` |
| 10 | `DoubleTriDiagonalMatrixTranspose` | real symmetric — a no-op, but the dispatch needs it |
| 11 | `DoubleTriDiagonalMatrixEffectiveHamiltonian` | strip `a0`/`b0`, keep the nested denominator |
| 12 | `DoubleTriDiagonalMatrixFromEffectiveHamiltonian` | inverse of 11 |
| 13 | `SetDoubleTriDiagonalMatrixConVal` | re-impose `val.a[0] == con.a[0]` after edits; cf. `SetBADoubleTriDiagonalMatrixConVal` |

Items 3, 5 and 11 are the ones that have to encode `Heff = Heff_val + Heff_con − eimp`
correctly, and item 3 is where it is currently coded wrong (§2a) — so those three deserve
the round-trip test of §7 rather than an eyeball check.

---

## 4. Conversions

`ConvertResponsefunction`
([BasicMath_GreenFunctionArithmetic.cpp:4130](../../Core/BasicMath/BasicMath_GreenFunctionArithmetic.cpp#L4130))
is a complete N×N table. A fifth type takes 4×4 = 16 cells to 5×5 = 25 — **9 new cells**.

Done:

- `AndersonMatrixToDoubleTriDiagonalMatrix` — splits the Anderson bath at `mu` with a Fermi
  function and builds the two chains with `PolesToTriDiagonalCoefficients`. Rewritten; §3.0
  and §3.1 are implemented here.
- `TriDiagonalMatrixToDoubleTriDiagonalMatrix` — routes via Anderson
- `ListOfPolesToDoubleTridiagonalMatrix` — routes via Anderson, see §0
- `BADoubleTriDiagonalMatrixToDoubleTriDiagonalMatrix` — routes via the poles, see §0
- `DoubleTriDiagonalMatrixToTriDiagonalMatrix`
- **`DoubleTriDiagonalMatrixToCompactMatrix[WithDeflation]`** and
  **`BlockDoubleTriDiagonalMatrixToCompactMatrix[WithDeflation]`** — done, see §0
- the `d → d` copy cell in `ConvertResponsefunction` and `ConvertBlockMatrix`, which is what
  `ChangeType(…,"Table"/"Userdata")` uses

To write:

| # | function | note |
|---|---|---|
| 1 | `ListOfPolesToDoubleTriDiagonalMatrixOccupation` | only if the occupation-constrained variant is wanted, cf. `ListOfPolesToBADoubleTridiagonalMatrixOccupation` — note that one stores `Abs(mu − muOccupation)` in the field, a third meaning again (§11) |
| 2 | `DoubleTriDiagonalMatrixToListOfPoles` | diagonalise; the one real algorithm still missing. `TriDiagonalMatrixToEigensystemCompact` on the dense form is the obvious route, and `DoubleTriDiagonalMatrixToCompactMatrix` already builds it |
| 3 | `DoubleTriDiagonalMatrixToAndersonMatrix` | route via `l` |
| 4 | `DoubleTriDiagonalMatrixToBADoubleTriDiagonalMatrix` | `d → n`: route via `l`, as `n → d` does. A rotation of the shared site into the natural pair is not available — the two forms combine their chains differently (§1) |
| 5 | `DoubleTriDiagonalMatrixCopy` | exists |

Item 2 unlocks 3 and 4 and closes the `d →` row, which is all that is left of the non-block
conversion table.

`AndersonMatrixToDoubleTriDiagonalMatrixNonOrthogonal` is a `/* To DO */` stub that returns
1. Nothing calls it; leave it.

### The three `…PlusT` variants are kept, and are now the odd ones out

`AndersonMatrixToDoubleTriDiagonalMatrixPlusT`,
`TriDiagonalMatrixToDoubleTriDiagonalMatrixPlusT` and
`TriDiagonalMatrixToSecondDoubleTriDiagonalMatrixPlusT` still hold the original
implementation, to be rewritten rather than replaced. Nothing outside `Test/` calls them, and
`Test/` is not in the build.

They no longer agree with the routine that `ChangeType` uses: they write the whole impurity
energy in both chains instead of half in each, split with the linear ramp instead of a Fermi
function, invent a zero-coupling bath site for an empty branch, and tri-diagonalize with
re-orthogonalized Lanczos. A comment on each says so, so that a caller who binds to one by
accident finds out. Rewriting them on top of `PolesToTriDiagonalCoefficients` means deciding
how `T` is produced — the recurrence gives it in O(NPoles*N) but only to 3e-8 orthogonality,
so machine precision still costs one Gram–Schmidt sweep and the O(NPoles*N) basis. `T` is not
exposed by any `ChangeType` on the lua level, so that is a project of its own rather than a
loose end here.

Cells 4 and 6 can legitimately route through `l` (`ConvertResponsefunction` already recurses
that way elsewhere) — only 1, 3 and 5 are new mathematics.

---

## 5. Response-function dispatch sites

Each needs a `case 'd': case 'D':`, or a new `else if((32|A->type)=='d')` branch. Almost all
mechanical.

**The live list is in §0** — it is regenerated from the source rather than maintained here,
because it shrinks every session. Regenerate it with a scan for functions that switch on the
type character and contain no `'d'` case.

Two notes that do not change:

- Several of these `switch`es end in a `default:` that `printf`s and continues rather than
  returning an error. A half-finished `d` therefore produces wrong numbers rather than
  stopping. Either finish a function completely or add `d` to its *guard* first.
- The order that matters: `FreeBlockMatrix` and `lua_isBlockGreensfunction` are the two that
  misbehave today rather than merely being absent (§0), so they are worth doing regardless of
  what else is picked up.

---

## 6. The actual goal: operator + `Nv`/`Nc`

### 6.1 Operator

New in [Core/MatrixOperator.cpp](../../Core/MatrixOperator.cpp), next to
`BADoubleTriDiagonalMatrixToImpurityHamiltonian` (line 167):

```c
int DoubleTriDiagonalMatrixToImpurityHamiltonian(OperatorType *O, DoubleTriDiagonalMatrixType *M,
                                                bool AddSpin, uint32_t FermionShift, uint32_t NF);
```

Take the **signature and the `AddSpin`/`FermionShift`/`NF` handling** from
`BADoubleTriDiagonalMatrixToImpurityHamiltonian`, and the **site layout and coefficients**
from `DoubleTriDiagonalMatrixToAndersonImpurityOperator`
([MatrixOperator.cpp:741](../../Core/MatrixOperator.cpp#L741)) — that function already
emits exactly the right Hamiltonian for this struct (impurity energy as
`val.a[0]/2 + con.a[0]/2`, `val.b[0]` and `con.b[0]` as the two impurity→chain hoppings),
just wrapped in a two-spin + `U` DMFT-specific interface.

Fermion ordering, `AddSpin == false`: `[FermionShift core] [impurity ×1] [val ×val.N] [con ×con.N]`.

### 6.2 `LuaResponseFunctionToOperator`

In [LuaInterface/ResponseFunction.cpp:2672](../../LuaInterface/ResponseFunction.cpp#L2672):

- **lines 2710–2711** — the `Nv`/`Nc` guards currently reject anything but `N`/`n`. Widen to
  include `D`/`d`.
- **new `case 'D': case 'd':`** in the switch at line 2714, structured like the `n` case:
  shorten `val.N`/`con.N`, build the operator, push the bitmask table, restore `val.N`/`con.N`,
  push `{GvalBath, GconBath}`. Two inner branches only (`AddSpin` true/false) rather than
  four, since only the non-block type exists at this stage. Uses `TriDiagonalMatrixPart` /
  `TriDiagonalMatrixPartAddSpin` unchanged (§3).
- the `n` case reaches into `M->data` to overwrite `val.N`/`con.N` and restores it afterwards.
  Same trick works for `d`; worth a comment, since a `luaL_error` between shorten and restore
  would leave the caller's object truncated.

### 6.3 `ResponseFunctionToLuaBitMask`

[LuaInterface/ResponseFunction.cpp:959](../../LuaInterface/ResponseFunction.cpp#L959),
switch at line 976. A genuinely new variant, not a copy of `n`:

| | `n` (non-block) | `d` |
|---|---|---|
| `type` string | `"Nat"` | new, e.g. `"DoubleTri"` |
| `NImpurity` | 2 (impurity + the extra natural orbital) | 1 |
| `BathExt` | present | **absent** |
| `NBathVal` | `val.N - 1` | `val.N` |
| `NBathCon` | `con.N - 1` | `con.N` |

Both `AddSpin` values; `Deflate` has no meaning without `BlockSizeList`, so the non-block
path ignores it as it already does for `n`.

---

## 7. Lua interface, IO, tests

### Type-name string

Pick one and use it at all five sites:

| site | file:line |
|---|---|
| `LuaChangeTypeResponseFunction` (+ its error message) | [ResponseFunction.cpp:109](../../LuaInterface/ResponseFunction.cpp#L109) |
| `luaL_checkResponsefunction` (+ its error message) | [TableToArray.cpp:2343](../../LuaInterface/lua_extensions/TableToArray.cpp#L2343) |
| `lua_isBlockGreensfunction` | [TableToArray.cpp:2419](../../LuaInterface/lua_extensions/TableToArray.cpp#L2419) |
| `LuaResponseFunctionIndexType` | [ResponseFunction.cpp:3309](../../LuaInterface/ResponseFunction.cpp#L3309) |
| `ResponseFunctionToLuaBitMask` | [ResponseFunction.cpp:976](../../LuaInterface/ResponseFunction.cpp#L976) |

### Table ⇄ struct — DONE

`luaL_checkDoubleTriDiagonalMatrixType`, `luaL_checkBlockDoubleTriDiagonalMatrixType`,
`LuaPushTable_DoubleTriDiagonal` and `LuaPushTable_BlockDoubleTriDiagonal` all exist, with
their `void*` overloads, and `LuaPush_ResponseFunction` and `LuaResponseFunctionToTable`
dispatch on `d`. The block push takes `B0` from `M->b0` and not as an identity matrix, unlike
its `n` model, because the `d` struct has a real `b0`.

Still missing: the cross-type branches of the two readers, i.e. reading a `Tri` / `And` /
`Nat` / `ListOfPoles` table *as* a `DoubleTri`. Those need the §4 conversions.

Use `kBT` as the Lua table field name for `d` from the start. The `n` type reads and writes
this as `epsilon`
([TableToArray.cpp:1452, 1694, 3270, 3310](../../LuaInterface/lua_extensions/TableToArray.cpp#L1452)),
which is user-facing and therefore stuck until `n` is migrated (§11) — but `d` has no legacy
to carry, so it should not inherit the old name.

The table layout can otherwise match `n`'s — the chains are laid out identically (§1). What
must not be interchangeable is the `type` field: `n` and `d` combine the same two chains differently, so a
table written as `Nat` and read back as `DoubleTri` gives silently wrong numbers. Reject a
missing or mismatched `type` rather than defaulting.

### Lua dispatch sites

`LuaResponseFunctionToTable` (173), `LuaResponseFunctionEigensystem` (394),
`LuaResponseFunctionToTightBinding` (2605), `LuaResponseFunctionToOperator` (2714),
`LuaResponseFunctionLen` (3201), `LuaResponseFunctionIndexMu` (3242),
`LuaResponseFunctionIndexType` (3309), `LuaResponseFunctionIndexName` (3350),
`LuaResponseFunctionNewIndexMu` (3519), `LuaResponseFunctionNewIndexName` (3587), all in
[LuaInterface/ResponseFunction.cpp](../../LuaInterface/ResponseFunction.cpp); plus
`LuaPush_ResponseFunction` ([TableToArray.cpp:3401](../../LuaInterface/lua_extensions/TableToArray.cpp#L3401)).

Everything else on the `ResponseFunction` metatable — `ReducePoles`,
`RemovePolesInInterval`, `MergePolesInInterval`, `Rotate`, `Chop`, `Elements`, `ToMatrix`,
`ToSpectra`, `CalculateRho`, `CalculateSelfEnergy`, `CalculateHybridizationFunction`,
`CalculateGfromG0andSigma`, `CalculateGbathfromGandSigma`, `InvertEnergy`, `Transpose`, and
the arithmetic metamethods — needs **no** change of its own once §4 and §5 are complete.
That is the payoff for doing §4/§5 properly.

`CalculateGbathfromGandSigma`
([ResponseFunction.cpp:559](../../LuaInterface/ResponseFunction.cpp#L559)) already takes
`Nv`, `Nc` and a `GBathFormat` string; check that its format parsing accepts the new name.

### IO

`PrintDoubleTriDiagonalMatrix` and `StreamPrintDoubleTriDiagonalMatrix` exist — check they
print the index-shifted convention meaningfully (label the shared site) rather than dumping
`val` and `con` as if they were standalone chains. No `FilePrint` variant is needed (only
`n` and `l` have one).

### Tests

In [Test/TestGreenFunctions.cpp](../../Test/TestGreenFunctions.cpp) (which already declares a
`DoubleTriDiagonalMatrixType` at line 138):

1. **`l → d → l` round trip** — the spectrum must be reproduced. Catches the §2(a) class of
   error and any convention slip in §4.
2. **`d` vs `n` on the same input** — convert one list of poles to both and check the spectra
   agree. This is the test that pins the additive-vs-nested distinction.
3. **`Nv`/`Nc` consistency** — `ToOperator{Nv,Nc}` plus the returned `GvalBath`/`GconBath`
   must reproduce the untruncated result. This is what guards the `Nv == 0` and
   `Nv >= val.N` edge cases of §3.
4. **`eimp ≠ 0`** — every one of the above with a non-zero impurity energy, since that is
   exactly what §2(a) gets wrong.

Note `Test/TestSpectral.cpp:103` calls `DoubleTriDiagonalMatrixToSpectrum` with an 8-argument
signature that no longer exists, so that file is stale and not built — worth checking what in
`Test/` still compiles before relying on it.

### Documentation

The type name and the widened `Nv`/`Nc` applicability are user-facing → quanty.org wiki page
for `ResponseFunction`. See [WIKI_EDITING_WORKFLOW.md](WIKI_EDITING_WORKFLOW.md).

---

## 8. The block type — mostly done

`BlockDoubleTriDiagonalMatrixType` exists with `Init`, `Free`, `Copy`, `Evaluate`,
`…ToCompactMatrix[WithDeflation]` and `…ToSpectra[WithBroadening]`. What is left of this
section is `Chop` / `ChopReduceDimension` (stubs), the conversions in §0, and the retirement
of the `Blok` family. The function-list mapping below is kept as the reference for what a
complete `Block…` family looks like.

### Original plan

`BlockDoubleTriDiagonalMatrixType` — same functionality as
`BlokDoubleTriDiagonalMatrixType`, in the `Block…` storage scheme: `BlockSize`, `QComplex`,
and flat `double*` / `Complex16*` arrays, instead of an array of separately-allocated
`MatrixType` blocks.

The `Blok` ASCII art
([BasicMath_BlokDoubleTriDiagonalMatrixVector.cpp:1](../../Core/BasicMath/BasicMath_BlokDoubleTriDiagonalMatrixVector.cpp#L1))
confirms the same **nested** structure as the scalar type — `acv0 = (a0val + a0con)/2`,
`bv0`/`bc0` as impurity↔chain hoppings, `A[1;;3][1;;3]` and `A[4;;6][4;;6]` as ordinary
block-tri-diagonal chains — so §1 carries over unchanged, with `a0`/`b0` becoming
`CompactMatrixType` and `val`/`con` becoming `BlockTriDiagonalMatrixType`.

Function list, from
[BasicMath_BlokDoubleTriDiagonalMatrixVector.h](../../Core/BasicMath/BasicMath_BlokDoubleTriDiagonalMatrixVector.h):

| Blok function | Block equivalent |
|---|---|
| `InitBlokDoubleTriDiagonalMatrix` | `InitBlockDoubleTridiagonalMatrix` |
| `FreeBlokDoubleTriDiagonalMatrix` | `FreeBlockDoubleTridiagonalMatrix` |
| `InitFrameBlokDoubleTriDiagonalMatrix` | **drop** |
| `FreeFrameBlokDoubleTriDiagonalMatrix` | **drop** |
| `InitIdentityBlokDoubleTriDiagonalMatrix` | `InitIdentityBlockDoubleTridiagonalMatrix` |
| `BlokDoubleTriDiagonalMatrixChop` | `BlockDoubleTridiagonalMatrixChop` (+ `…ChopReduceDimension`, new — `BlockSizeList` deflation has no `Blok` analogue) |
| `BlokDoubleTriDiagonalMatrixAdd` | `BlockDoubleTridiagonalMatrixAdd` |
| `BlokDoubleTriDiagonalMatrixScalarMultiply` | `…ScalarMultiply` |
| `BlokDoubleTriDiagonalMatrixPart` | `…Part` — the `Nv`/`Nc` primitive |
| `BlokDoubleTriDiagonalMatrixCopy` | `CopyBlockDoubleTridiagonalMatrix` |
| `BlokDoubleTriDiagonalMatrixTranspose` | `BlockDoubleTriDiagonalMatrixTypeTranspose` |
| `BlokDoubleTriDiagonalMatrixComplexConjugate` | `…ComplexConjugate` |
| `BlokDoubleTriDiagonalMatrixConjugateTranspose` | `…ConjugateTranspose` |
| `Real/ComplexBlokDoubleTriDiagonalMatrixVector` | `Real/ComplexBlockDoubleTridiagonalMatrixVector` |
| `Real/ComplexVectorBlokDoubleTriDiagonalMatrixVector` | `Real/ComplexVectorBlockDoubleTridiagonalMatrixVector` |
| `MatrixToBlokDoubleTriDiagonalMatrixNoOverlap` | `CompactMatrixToBlockDoubleTridiagonalMatrixNoOverlap` |

`Init/FreeFrame` have no analogue by construction: the whole point of the `Block` layout is
one flat allocation instead of a descriptor array plus per-block allocations, so the
frame/data split disappears. That is 2 fewer functions than the `Blok` list, and 1 more
(`ChopReduceDimension`) for the deflation the `Block` family supports.

Plus, mirroring §3–§7 for the block case: `ComplexBlockDoubleTridiagonalMatrixToReal` /
`Real…ToComplex`, `Shift…`, `Evaluate…` (×2), `…EigenValues`, `…EigenSystem`,
`…ToCompactMatrix` (+`WithDeflation`), `…ToSpectra` (+`WithBroadening`), `…ToTightBinding`,
`…EffectiveHamiltonian` / `…FromEffectiveHamiltonian`, `…GetTail` /
`…GetTailAndAddSpin`; the 9 new cells of `ConvertBlockMatrix`
([BasicMathBlock.cpp:153](../../Core/BasicMathBlock/BasicMathBlock.cpp#L153));
`BlockDoubleTridiagonalMatrixToImpurityHamiltonian`; the block halves of every dispatch site
in §5; and `PrintBlockDoubleTriDiagonalMatrix` /
`StreamPrintBlockDoubleTriDiagonalMatrix`.

Under option 1 (§1) the block type needs **no** shared-site member either: the impurity block
is `val.a[0] + con.a[0]`, and

```
G(w) = A0 + B0*·( w − G.val(w) − G.con(w) )⁻¹·B0^T
```

with `G.val` evaluated by the existing `EvaluateBlockTriDiagonalMatrix`, which already
includes `val.a[0]` and `val.b[0]` per the `BlockTriDiagonalMatrixType` convention. So the
struct is `Name`, `BlockSize`, `QComplex`, `mu`, `Epsilon`, `CompactMatrixType a0`,
`CompactMatrixType b0`, `val`, `con` — and nothing is reserved-zero.

Note the block tail extraction needs *no* new code: `BlockTriDiagonalMatrixGetTail` and
`BlockTriDiagonalMatrixGetTailAndAddSpin` already act on a `BlockTriDiagonalMatrixType`, so
they apply to `.val` and `.con` directly, as `TriDiagonalMatrixPart` does in the non-block
case (§3).

### Build registration

The build is a unity build (`.cpp` files are `#include`d), so new files are invisible until
listed:

1. `#include "BasicMath_BlockDoubleTridiagonalMatrix.h"` in
   [BasicMathBlock.h](../../Core/BasicMathBlock/BasicMathBlock.h), **before**
   `BasicMath_BlockTypeConversion.h`;
2. `#include "BasicMath_BlockDoubleTridiagonalMatrix.cpp"` in
   [BasicMathBlock.cpp](../../Core/BasicMathBlock/BasicMathBlock.cpp);
3. struct into [StructDef.cpp](../../Core/StructDef.cpp);
4. check whether [CMakeLists.txt](../../CMakeLists.txt) and
   [Model.xcodeproj/project.pbxproj](../../Model.xcodeproj/project.pbxproj) list headers
   explicitly.

### Naming

`BlockDoubleTriDiagonalMatrixType` will coexist with `BlokDoubleTriDiagonalMatrixType`,
one letter apart, for as long as the `Blok` family survives — same as
`BlockTriDiagonalMatrixType` vs `BlokTriDiagonalMatrixType` today. Say so in the struct
comment. The `Blok` double-tri-diagonal family has exactly three entry points outside its
own file (`MatrixToBlokDoubleTriDiagonalMatrixNoOverlap`, `PrintBlok…`,
`StreamPrintBlok…`), so retiring it once the `Block` version lands is cheap.

---

## 9. Staging

Stages 0 to 2 of the original plan are done. What is left, in the order I would take it:

| # | content | why this order |
|---|---|---|
| ~~0~~ | ~~the `→ d` column, non-block~~ | **done**, see §0 |
| 1 | the `§3.0` head sweep | now urgent rather than tidy: the constructors write halves, the five matrix–vector products and the eight sites in the two `…ToAndersonImpurityOperator` builders still read the average, so they disagree with what `ChangeType` produces |
| 2 | `DoubleTriDiagonalMatrixToListOfPoles` | closes the `d →` row: it unblocks the zero-broadening spectrum path, and by routing through `l` it gives `d → t`, `d → a` and `d → n` without writing each one |
| 3 | the two that misbehave today: `FreeBlockMatrix`, `lua_isBlockGreensfunction` | small, and they are wrong rather than absent |
| 4 | §2(b) and §2(c) | dormant functions, but they are reachable now that the conversions exist |
| 5 | the six Lua index metamethods | `G.type`, `G.mu`, `G.name` on a `d` userdata |
| 6 | §6, the operator and the bitmask | the original goal |
| 7 | the rest of the §0 dispatch list, the block conversions, and the deflated path once something can build a deflated `d` | |

This is no longer a purely additive change. `AndersonMatrixToTriDiagonalMatrix` and
`ListOfPolesToTridiagonalMatrix` were rewritten onto the shared kernel, so `a → t`, `l → t`,
`l → a` and the non-block `→ n` path all move by about 1e-14; the 510 assertions in
`TestScripts/Responsefunction/All.Quanty` still pass, and the chain coefficients were compared
against a binary built from the previous commit. The three `…PlusT` variants keep their
original implementation and now disagree with it (§4). `ListOfPolesToTridiagonalMatrix` now
propagates `mu` from the list of poles, where it used to leave it at 0. And the `default:` branches of `ConvertResponsefunction` no longer print and
continue: an unimplemented cell leaves `Mout->data` at `NULL` and the function returns
non-zero (§0).

---

## 10. Adjacent findings: reserved-zero slots in the block NIO type

Recorded because it is the evidence behind the option-1 decision in §1, and because it is a
latent defect set in live code — `n` is a working, Lua-exposed type.

`BlockNaturalImpurityOrbitalMatrixType` has an outer `CompactMatrixType a0` and no `b0`; the
scaling lives in `val.b[0]` / `con.b[0]`. Its convention is additive,

```
G(w) = a0 + G.val(w) + G.con(w)
```

so the total constant offset is `a0 + val.a[0] + con.a[0]` — the plain sum. The head slots
`val.a[0]` / `con.a[0]` are *intended* to be zero, but nothing enforces that except the
constructors, and the code does not agree on what they mean:

| path | `val.a[0]` non-zero → |
|---|---|
| `EvaluateBlockNaturalImpurityOrbital` | **honoured** — delegates to `EvaluateBlockTriDiagonalMatrix`, which adds its own `a[0]` |
| `…MatrixToSpectra` | **honoured** — `Gval[i] + Gcon[i] + A->a0` |
| `…MatrixToBlockListOfPoles` | **honoured** — and hence `…ToBlockAndersonMatrix` / `…ToBlockTridiagonalMatrix`, which route through it |
| `BlockResponseFunctionGetA0` | **honoured** — returns `a0 + val.a[0] + con.a[0]`, carrying a `//TODO: test if this is correct` |
| `…MatrixChop` | touched — loops all `N+1` a-blocks, block 0 included |
| `…MatrixToCompactMatrix`, `…WithDeflation` | **dropped** — built from `TMT` + `a[k≥2]` + `b[k≥2]` |
| `…MatrixToImpurityHamiltonian` | **dropped** — same |
| `GetTFromValandCon…`, `SetTMTFromValandCon…` | **dropped** — read only `b[0]` and `a[1]`, `a[2]`, … |
| `ResponsefunctionMul` (block `n`) | **wrong** — scales `a0` by `x` and `val.b[0]`/`con.b[0]` by `√x`, but not `val.a[0]`/`con.a[0]`, so `x·G ≠ Mul(G, x)` |
| `ResponseFunctionInvertEnergy` (block `n`) | skips block 0, i.e. treats it as an offset |

The "dropped" rows are not wrong in themselves — an additive constant on `G` has no home in
a matrix or a Hamiltonian. The problem is that `Evaluate` and `ToMatrix` can disagree about
the same object.

Two places do enforce zero:

- `BlockListOfPolesToBlockNaturalImpurityOrbitalMatrix` `memset`s `Aval.a0` / `Acon.a0` to
  zero after copying the poles' offset into `B->a0`;
- `luaL_checkBlockNatrualImpurityOrbitalMatrixType`
  ([TableToArray.cpp](../../LuaInterface/lua_extensions/TableToArray.cpp)) does the full
  normalisation, under its own comment
  `// set A0 of val and con to zero  /  A0 = A0 + B^* A0con B^T + B^* A0val B^T`: it
  transforms both head blocks with the table's `B0`, adds them into `a0`, zeroes them in
  place, then folds `B0` into `val.b[0]` / `con.b[0]`.

So a `n` built from Lua or from a list of poles is fine. One built or edited through the C
API is not guaranteed to be.

Worth doing at some point, independent of `d`:

1. Resolve the `BlockResponseFunctionGetA0` TODO — it *is* correct, in that it agrees with
   `Evaluate` and `ToSpectra`. Delete the TODO and say so.
2. Fix block-`n` `ResponsefunctionMul` to scale `val.a[0]` / `con.a[0]` by `x` as well —
   currently a no-op in practice, but it is the one operation that is outright wrong rather
   than merely inconsistent.
3. State the reserved-zero invariant in the struct comment, since two constructors already
   depend on it.

---

## 11. Adjacent findings: the `Epsilon` split parameter — partly resolved

`Epsilon` became `kBT` in `DoubleTriDiagonalMatrixType` and
`BlockDoubleTriDiagonalMatrixType`, and the convention is settled: a true Fermi function, the
temperature stored with the object. The **implementation** is not done —
`AndersonMatrixToDoubleTriDiagonalMatrixPlusT` still uses the linear ramp — and the other five
structs still carry `Epsilon`. The mixed real/complex half of this has moved to its own entry
in [TODO.md](TODO.md). The survey below is kept as the record of why.

### Survey

Survey behind the §3.1 decision, and the follow-up work it implies for the other types.

### `Epsilon` is write-only

Every use of the member outside `Test/`, across all six structs that carry it:

| kind | where |
|---|---|
| set to `0` | `Init…` for `DoubleTriDiagonal`, `BADouble`, `BlockNIO`, `Blok` |
| propagated | `Copy` / `Part` (and `*B = A` struct assignment in the NIO copy) |
| set from an argument | only `AndersonMatrixToDoubleTriDiagonalMatrixPlusT`, `ListOfPolesToBADoubleTridiagonalMatrix`, and `ListOfPolesToBADoubleTridiagonalMatrixOccupation` — the last storing `Abs(mu − muOccupation)`, a different quantity |
| Lua | read from / written to the table field `epsilon` |
| print | `StreamPrint…` |

It never appears on the right-hand side of a computation. Every splitting site takes the
parameter as a function argument instead. `BlockNaturalImpurityOrbitalMatrixType.Epsilon` is
the extreme case: set to `0` in `Init` and never set again, because the block constructors
take `kbT` and do not record it — so it is *always* 0 while the comment above it purports to
describe the split.

### The documented rule was never the implemented one

`StructDef.cpp` says *"values in a range of [Ef-1/2 Epsilon , Ef+1/2 Epsilon] contribute half
to val and half to con"*. What `ListOfPolesCopyVal` / `…Con`
([BasicMath_ListOfPoles.cpp](../../Core/BasicMath/BasicMath_ListOfPoles.cpp)) actually do:

```c
if      (M->a[i] >= Ef + Epsilon)  M->b[i] = 0;
else if (M->a[i] >  Ef - Epsilon)  M->b[i] *= ((Ef - M->a[i]) + Epsilon) / (2 * Epsilon);
```

That is a **linear ramp**, 1 at `Ef − Epsilon` down to 0 at `Ef + Epsilon`, with the con
mirror summing to 1. So the comment is wrong twice over: half/half holds only exactly at
`Ef`, and the window is `±Epsilon` (full width `2·Epsilon`), not `±Epsilon/2`.

### Three rules coexist today

| rule | implemented in | used by |
|---|---|---|
| linear ramp on the residue | `ListOfPolesCopyVal` / `…Con` | non-block `l/t/a → n` |
| linear ramp on the hopping (`sqrt`) | `AndersonMatrixToDoubleTriDiagonalMatrixPlusT` | non-block `→ d` |
| Fermi function | `CopyBlockListOfPolesVal` / `…Con`; `AndersonMatrixToBADoubleTriDiagonalMatrixTemperature` | all block `→ n`; the DMFT / `Model/AndersonImpurity` / `Yi` paths |

The two ramps agree with each other on the weight. The Fermi routine
`AndersonMatrixToBADoubleTriDiagonalMatrixTemperature` exists but is **not reachable from
`ConvertResponsefunction`** — only from `DMFT/SingleSiteFunctions.cpp`,
`Model/AndersonImpurity.cpp` and `Yi/*`.

### `EffectivekBT` currently means two different things

`ResponseFunction.ChangeType{EffectivekBT = …}` — documented as *"used to determine if a
state belongs to the valence or conduction bands (occupied or empty using Fermi Dirac
statistics)"* — becomes the `T` argument of `ConvertResponsefunction`, and then:

| path | receives it as | rule |
|---|---|---|
| block `→ n` | `kbT` | Fermi function ✓ matches the name |
| non-block `→ n` | `Epsilon` | linear ramp half-width ✗ |

Same number, different physics, no warning. The default `sqrt(DBL_MIN)` degenerates to a step
at `mu` in both, so it does not bite today; `EffectivekBT = 0.025` gives a Fermi tail of
10–90% width ≈ 4.4 kbT ≈ 0.11 eV in one path and a ramp confined to exactly ±0.025 eV in the
other.

### Follow-up beyond `d`

Once `d` is consistent (§3.1), the same treatment applies to the remaining five structs —
`BADoubleTriDiagonalMatrixType`, `SecondDoubleTriDiagonalMatrixType`,
`SecondAndersonMatrixType`, `BlockNaturalImpurityOrbitalMatrixType`,
`BlokDoubleTriDiagonalMatrixType`:

1. rename `Epsilon` → `kBT`;
2. replace `ListOfPolesCopyVal` / `…Con` with the Fermi weighting, which makes the non-block
   `→ n` path agree with the block one and with `EffectivekBT`'s documented meaning;
3. have the block NIO constructors actually *store* the `kbT` they are given, instead of
   leaving the field at 0;
4. decide the Lua field migration: `epsilon` is user-facing for `n`, so accept both names on
   read and emit `kBT`, keeping `epsilon` as a deprecated alias rather than breaking existing
   scripts;
5. `ListOfPolesToBADoubleTridiagonalMatrixOccupation` stores `Abs(mu − muOccupation)` in the
   field — decide whether that is a `kBT` at all, or belongs somewhere else.
