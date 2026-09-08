# Plan: response-function type `d`/`D` — double tri-diagonal

Goal: `ResponseFunction.ToOperator{Nv=…, Nc=…}` on a response function held as a *double*
tri-diagonal matrix — one impurity site with a valence chain and a conduction chain hanging
off it — so `Nv` and `Nc` truncate the two chains independently and the two tails come back
as `GvalBath` / `GconBath`.

Scope decided: **the single-valued (non-block) type first.** The block type
(`BlockDoubleTriDiagonalMatrixType`) comes after, as a Block-style reimplementation of the
`Blok…` function list — see §8.

---

## 0. Status as of `eb91c974` (2026-09-07)

**Regenerated from the source, not from memory.** The lists below were produced by
scanning the tree; redo the scans rather than trusting this section after a few
more commits.

### Starting a new session from here

**`BlockDoubleTriDiagonalMatrixToBlockAndersonMatrix` is written**, as §8
decided it: one Anderson block per merged pole, `A_{k+1} = e_k·I` and
`B_k = sqrt(W_k)` kept whole and dense. It reproduces `G` to 2.4e-13, every bath
block is bit-exactly `e·I`, and a degenerate bath comes back as one block per
distinct energy. Its signature gained `MinimalPoleDistance` and `relativeDistance`
and takes `M` by pointer; the two call sites are `ConvertBlockMatrix` and
`luaL_checkBlockAndersonMatrixType`, which passes `-1.0, false` for the standard
default.

**`d → t` is written too, through `a`** — Maurits chose that route: `d → a → t`,
with `SingularValue` for the deflation of the second step and
`MinimalPoleDistance` for the merge of the first. Same signature shape as
`BlockAndersonMatrixToBlockTridiagonalMatrix`, and the `NTriDiag` argument the
stub carried is gone, since `a → t` has no maximum chain length to give it to.

**`d → l` is written, and it does not go through `a`** — the whole matrix is
band-reduced and diagonalized directly, `MinimalPoleDistance` and
`relativeDistance` only.

**`d → n` is written as `d → l → n`**, the route the scalar version takes, with
`kBT` for the Fermi split, `SingularValue` for the deflation of the two block
Lanczos runs and `MinimalPoleDistance` for the merge.
`BlockListOfPolesToBlockNaturalImpurityOrbitalMatrix` gained `SingularValue` so
that it could be threaded; it used to hard-code `1000·DBL_EPSILON`, which is what
its other callers now pass, so nothing changed for them.

**The conversion table and §6 are both done.** `ToOperator` works for a `d` in
both forms, with `Nv`/`Nc` truncating the two chains and the tails coming back as
`GvalBath`/`GconBath`, and `CreateBitMasks` counts the right orbitals. See §6 and
the `d` sections of `ToOperator.Quanty` and `CreateBitMasks.Quanty`.

**The next work is the list of broken entry points** in "The rest of the lua
surface" just below: two of them segfault, four are silent no-ops that return the
input unchanged, and the rest fail in the middle without raising. That subsection
carries the order to take them in and which of them share a root cause.

Read §1 for what `d` means and how it differs from `n`, and
[TODO_ConversionAccuracy.md](TODO_ConversionAccuracy.md) for what the accuracy
numbers do and do not say.

Three test files, all registered in `All.Quanty` **and** in `.gitlab-ci.yml` — a
new test needs both:

| file | covers |
|---|---|
| `ChangeTypeDoubleTri.Quanty` | into `d`, non-block |
| `ChangeTypeFromDoubleTri.Quanty` | out of `d`, non-block |
| `ChangeTypeDoubleTriBlock.Quanty` | the block row, both directions |
| `ToOperator.Quanty` | `d` to an impurity Hamiltonian, both forms, and `FermionShift` with `AddSpin` for all four block types |
| `CreateBitMasks.Quanty` | the orbital count and the masks of a `d`, both forms |
| `ChangeType.Quanty` | that `SingularValue` reaches the six cells that take it |
| `LengthZero.Quanty` | a block list of poles with no poles at all, through every conversion |
| `__add.Quanty` | a real block list of poles plus a complex one |

**Build out of the source tree** — `g++ … -o <somewhere else>/Quanty` — and run
from the repo root as CI does, `./Quanty "TestScripts/Responsefunction/X.Quanty"`.
Nothing depends on where the binary lands, and a stray `Quanty` in the root is not
in `.gitignore`.

### The state of the block row

| | `→ l` | `→ t` | `→ a` | `→ n` | `→ d` |
|---|---|---|---|---|---|
| from `l` | — | ✓ | ✓ | ✓ | ✓ |
| from `t` | ✓ | — | ✓ | ✓ | ✓ |
| from `a` | ✓ | ✓ | — | ✓ | ✓ |
| from `n` | ✓ | ✓ | ✓ | — | ✓ |
| from `d` | ✓ | ✓ | ✓ | ✓ | ✓ copy |

**The block row is complete**, both directions, and every cell is tested against
the function it came from. Each outbound cell is its own reading of the same
matrix, which is the point: `a` from the two chain spectra side by side, `t`
through `a`, `l` by band-reducing the whole matrix, `n` through `l`.
`Chop` works both ways with and without deflation, the six index metamethods work,
and both lua readers convert on all four type branches.

`MinimalPoleDistance` is now an option on `ChangeType`, threaded through
`ConvertResponsefunction`/`ConvertBlockMatrix` to every conversion that merges.
Unset arrives negative, and **`BlockListOfPolesMergeDegeneratePoles` is the single
place that turns that into the standard value** `10·DBL_EPSILON·(BandWidth + DBL_EPSILON)`.
Do not add a second default anywhere; forward what you were given.
`relativeDistance` and `relativeWeigth` are refused at the lua layer — neither is
implemented. Note `MinimalPoleDistance == 0` merges **nothing**, not even
bit-identical poles: the test in the merge is a strict `<`.

### Done: the whole non-block conversion table

Every cell of `ConvertResponsefunction` involving `d` works, both ways —
`l`,`t`,`a`,`n` → `d` and `d` → `l`,`t`,`a`,`n`, plus the copy. Each direction
uses a different reading of the same matrix, and none was a template for the next:

| | how | cost |
|---|---|---|
| `a → d` | Fermi split of the bath, RKPW on each half | O(N·k), the shared `PolesToTriDiagonalCoefficients` |
| `t → d`, `l → d`, `n → d` | route through `a` | — |
| `d → a` | the minor is the two chains, **disconnected** | two tri-diagonal eigenproblems |
| `d → t` | interleaved `imp, val_1, con_1, val_2, con_2, …` is a band matrix, `kd = 2`, observed site first | one `dsbtrd` |
| `d → l` | `val_N … val_1, imp, con_1 … con_N` is plain tri-diagonal, `kd = 1`, weights from the middle row | one `dstevd` |
| `d → n` | through `d → l`, the only thing the two forms agree on | — |

Also done and verified: the struct in both forms, `Evaluate`, the four
`…ToCompactMatrix[WithDeflation]`, `ToSpectra`, the accessors, `Free`, the lua
table round trip, §3.0 (the shared site is the **sum** of the two head slots,
everywhere) and §3.1 on the `a → d` path (a true Fermi function, `kBT` stored).

### The rest of the lua surface, audited at `eb91c974`

Re-derived from the source and **verified by running**, not by reading: every entry of the
`ResponseFunction` metatable was called on a single-valued and on a block `d`, with a `t`
and an `n` of the same shape as controls, each probe in its own process so a crash only
took that probe with it. Controls matter here — several things that look like `d` gaps are
how every type behaves.

**Works.** `ChangeType` (every cell, both directions), `Evaluate` / `__call`, `ToTable`,
`ToUserdata`, `ToMatrix`, `ToSpectra`, `ToOperator`, `CreateBitMasks`, `Elements`, `Chop`,
`__tostring`, reading and writing `mu` / `type` / `name`, `G+G`, `G-G`, `G==G`,
`CalculateSelfEnergy`, `ReducePoles`, `RemovePolesInInterval`, `MergePolesInInterval`,
`PlaceOnMesh`, and both free paths (`FreeResponsefunction`, `FreeBlockMatrixFromLua`).

**Broken, worst first.** `Rotate` is **done** — the `d` case goes through the list of poles
like the other three, the `default` returns 1 instead of falling through to `return 0`
(which is what turned a missing case into a segfault), and `Rotate.Quanty` now covers
complex functions, complex rotations and the `d` type. It was registered in neither
`All.Quanty` nor the CI matrix, so it had never run; it is in both now.

| symptom on a `d` | what happens | control | the function without a `d` case |
|---|---|---|---|
| `CalculateGbathfromGandSigma` | **SIGSEGV** | `t`, `n` → `And` | `ResponseFunctionEffectiveHamiltonian`, see below |
| `EffectiveHamiltonian` | returns an object whose **type is uninitialised**; any use of it fails and it leaks when collected | `t` → `Tri` | `ResponseFunctionEffectiveHamiltonian` |
| `FromEffectiveHamiltonian` | the same | `t` → `Tri` | `ResponseFunctionFromEffectiveHamiltonian` |
| `CalculateRho` | prints "unrecognized type d", no usable result | `a` gives rho | `ResponseFunctionToDensityMatrix` |
| `Eigensystem` | prints an error, returns nil | `t` fine | `LuaResponseFunctionEigensystem` |
| `CalculateGfromG0andSigma` | prints an error but computes the right answer, and returns a **`ListOfPoles`** instead of a `d` | — | `ResponseFunctionAddSelfEnergyDenseViaTri` |
| `#G` | clean lua error, "undifined type d" | `t` → 2 | `LuaResponseFunctionLen` |
| `ToTightBinding` | clean lua error | block `t` works | `LuaResponseFunctionToTightBinding` |

**Two of those are one bug.** `CalculateGBathFromGandSigmaDirectSubtract`
([BasicMath_GreenFunctionArithmetic.cpp:1845](../../Core/BasicMath/BasicMath_GreenFunctionArithmetic.cpp#L1845))
calls `ResponseFunctionEffectiveHamiltonian` and hands the result straight to
`ConvertResponsefunction` on the next line but one, so the uninitialised type becomes a
segfault there. Giving `EffectiveHamiltonian` a `d` case fixes the crash, the corrupt
object and `FromEffectiveHamiltonian` together.

**`Transpose` — DONE.** `BlockDoubleTriDiagonalMatrixTypeTranspose` follows the
`A_i → A_i^T`, `B_i → B_i^*` of every other transpose in the family:

    G   = A0   + B0^* ( w - Gval(w)   - Gcon(w)   )^-1 B0^T
    G^T = A0^T + B0   ( w - Gval(w)^T - Gcon(w)^T )^-1 B0^+

The two chains are added inside the denominator, so each transposes as the tri-diagonal
response function it is — `BlockTriDiagonalMatrixTypeTranspose` on both — `a0` transposes
and `b0`, standing in the head convention, conjugates. Four lines, the same shape as
`BlockNaturalImpurityOrbitalMatrixTypeTranspose`. No `QComplex` branch is needed: `a0` and
`b0` are `CompactMatrixType`s carrying their own flag, and `CompactMatrixConjugate` reads it
and leaves a real matrix alone. The single-valued half needs nothing — the function returns
early on `Block==false`, a 1×1 matrix being its own transpose. `ResponseFunctionTranspose`
has a `default` returning 1 now and `LuaResponseFunctionTranspose` raises on it.

That single rule covers both conjugation conventions, which is the thing that usually goes
wrong here: for the head, `(B^* X B^T)^T = B X^T B^+` and putting `B' = B^*` into the format
gives `B'^* X' B'^T = B X' B^+`; for a level below, `(B X B^+)^T = B^* X^T B^T` and the same
substitution gives `B'X'B'^+ = B^* X' B^T`. Verified against `EvaluateBlockTriDiagonalMatrix`,
which uses `TransformCompactMatrixInlinecT` (`b^* ( ) b^T`) for `iblock==1` and
`TransformCompactMatrixInlineNC` (`b ( ) b^+`) below.

`Transpose.Quanty` asserted only that the four formats agree with each other before and
after the transpose — which they would also do if all of them were wrong in the same way —
and otherwise drew plots. It now also checks `Transpose(G)(w,γ) == Transpose(G(w,γ))` for
all five formats, real and complex, that transposing twice returns `G`, and that a single
valued `G` of every format is its own transpose. **The head needs a hand built `G`**: every
response function that comes out of a conversion has a hermitian head, `A0` hermitian and
`B0` from `CompactMatrixSqrt` of a weight, and for a hermitian `B0` the conjugate and the
transpose are the same matrix — a build that transposes `b0` where it should conjugate it
passes every converted case. So the file ends with a `d`, a `t`, an `a` and an `l` built by
hand whose `A0` and `B0` are neither symmetric nor hermitian.

**`InvertEnergy` — DONE.** The map the function implements is `A0 - P(2 mu - w)^+` for
`G = A0 + P(w)`: the pole part is inverted and the constant term stays, so the `G(mu+w) ->
-G(mu-w)^*` of its header holds as it stands only for `A0 = 0`. For a `d` the two chains
stand in the *denominator*, so unlike the chains of a `n` — added, their heads a constant of
`G` that the map leaves alone — their heads follow: `val.a[0] + con.a[0]` is the onsite
energy of the observed site, half in each, and it maps to `2 mu` minus itself, which is
`mu - a[0]` for each of the two. Below the head `a[i] -> 2 mu - a[i]`, `a0` and `b0` do not
move, and the two chains are exchanged. Two static helpers,
`DoubleTriDiagonalMatrixChainInvertEnergy` and its block twin, do one chain each so that
`val`, `con`, real and complex are not four copies of the same loop. Both `default`
branches return 1 now and `LuaResponseFunctionInvertEnergy` raises on it.

Two repairs came with it, both in the single-valued `n` case. Its `val ↔ con` exchange
swapped the `a` and `b` pointers but not `N`, which corrupts the structure whenever the two
chains differ in length; it swaps the whole structs now. And it recomputed the seven cached
numbers of the observed head combination by hand — `acvB` and `acvA` are each invariant
under the exchange, `acvBA` changes sign, and the `+= 2 mu` silently assumed
`val.b[0]² + con.b[0]² == 1` — where it now calls
`SetBADoubleTriDiagonalMatrixConValFromChains`, exactly as the block case calls
`SetTMTFromValandConInBlockNaturalImpurityOrbitalMatrixType`.

Names survive it now as well. `ListOfPolesCopy` and `BADoubleTriDiagonalMatrixCopy` copy
the `Name`, and every block copy starts with `*B = A`, which copies the array with the rest
of the struct; `TriDiagonalMatrixCopy`, `AndersonMatrixCopy` and `DoubleTriDiagonalMatrixCopy`
did not, so a single-valued `t`, `a` or `d` came back from `InvertEnergy` — and from `Chop`,
and from a `ChangeType` to its own type, all of which copy — carrying whatever
`InitTriDiagonalMatrix` and friends had written into it. One `strcpy` each.

Still open, and a bigger surface: the block **conversions** do not carry the name across.
A block `l → t` gives "Matrix", `l → a` and `l → d` give "Block Anderson Matrix". The copies
are fine, so this sits in the conversion table itself rather than in one function.

`InvertEnergy.Quanty` gained a hand written `d` in all six of its sections (single-valued
and matrix, real and complex, `mu = 0` and `mu ≠ 0`), both directions; the identity
`Ginv(mu+w) = A0 - [G(mu-w) - A0]^+` across all five formats, which is the only place the
convention is written down; and a check on the stored numbers that the two chains really
are exchanged — no comparison of two response functions can see that, a `d` and the same
`d` with its chains swapped being equal as functions — and a check that the name of the
response function survives.

**`×` and `÷` — DONE.** `ResponsefunctionMul` has a `d` case on both halves: `a0 *= x`,
`b0 *= sqrt(x)` and the two chains untouched, which is the whole of it as `b0` appears
twice in `G = a0 + b0^* ( w - Gval - Gcon )^-1 b0^T`. Its `default` returns 1 now and both
lua wrappers raise on that, so a format nobody handled can no longer come back unscaled.
`__mul.Quanty` and `__div.Quanty` compare a hand written `d` with a hand written `2 d`,
single-valued and block — the table constructor for a `d` is
`{ {a0, b0}, val=Gval, con=Gcon, mu=0, type="DoubleTri" }`, documented at
[TableToArray.cpp:470](../../LuaInterface/lua_extensions/TableToArray.cpp#L470), with `val`
and `con` tri-diagonal where a `n` takes lists of poles. Both files also check the identity
`x G(w,γ) == (G*x)(w,γ)` for all five formats, single-valued and block, real and complex;
before that nothing in the tree compared the values of a scaled response function with the
scaled values of the original, which is why the no-op survived.

**Found by scanning, not reachable from lua in these probes.** `FreeBlockMatrix` has no
`d` case and leaks; it is the legacy one, its own comment says to use
`FreeBlockMatrixFromLua`, and it has two callers, `ssp.cpp:3991` and
`TableToArray.cpp:4537`. `ResponseFunctionAddSelfEnergyOnGridViaG0Tri` and
`…ViaAndLP`, the `Method` variants of the self-energy calls, also lack `d`.

**Not `d` gaps — do not chase these.** `-G` is refused for every type. Single-valued
`Rotate` wants a block representation whatever the type. `add`, `sub`, `ReducePoles`,
`RemovePolesInInterval`, `MergePolesInInterval` and `PlaceOnMesh` return a `ListOfPoles`
for every type, by design. Single-valued `ToTightBinding` is unimplemented for `a` too.
`CalculateHybridizationFunction` was **not** exercised: a one-argument call raises for `a`
as well, so it needs a signature that was not worked out.

**Suggested order.** `CalculateRho`, `Eigensystem`, the `CalculateGfromG0andSigma` return
type and `#G`; then `FreeBlockMatrix`, the `Method` variants and `ToTightBinding`, which are
worth least.

### Still open, non-block

**§6 is done.** `DoubleTriDiagonalMatrixToImpurityHamiltonian`,
`BlockDoubleTriDiagonalMatrixToImpurityHamiltonian`, the `d` case of
`LuaResponseFunctionToOperator` with `Nv`/`Nc`, and the `d` case of
`ResponseFunctionToLuaBitMask` all work and are tested — see the audit above and
`ToOperator.Quanty`. What is left of the `d` row is the list of broken entry points there.

**Zero-length chains — DONE.** `InitDoubleTriDiagonalMatrix` refused
`val.N == 0` or `con.N == 0` and returned 2 without allocating, so
`DoubleTriDiagonalMatrixCopy` and `DoubleTriDiagonalMatrixPart` — which ignore
the return value — then wrote through uninitialised pointers. `ChangeType` of a
`d` with an empty chain to `d` died on it, verified against a build of
`815f1e04`. The reject is gone; `InitTriDiagonalMatrix` already handled `N = 0`
(it allocates `a[0]`, no `b`, and guards `b[0] = 1`), as does
`InitBlockTridiagonalMatrix` on the block side, and `EvaluateTriDiagonalMatrix`
returns `a[0]` for an empty chain, which is the right answer with no special
case. `AndersonMatrixToDoubleTriDiagonalMatrix` no longer needs to allocate its
two chains by hand and calls `InitDoubleTriDiagonalMatrix` again.
`InitIdentityDoubleTriDiagonalMatrix` and `InitIdentityTriDiagonalMatrix` accept
`N = 0` as well; the latter now guards `b[0] = 1` exactly as `InitTriDiagonalMatrix`
does. Neither is reachable from lua — `InitIdentityTriDiagonalMatrix` is called
only by `InitIdentityBADoubleTriDiagonalMatrix`, which has no callers at all — so
those two are verified by compilation and inspection, not by a test.
**Both `BADoubleTriDiagonalMatrixType` constructors keep their `N = 0` reject, and
it is not an oversight.** Maurits: the `n` Green's function is defined as a
*linear combination of the first bath site of the val and the con chain*, and that
combination is only well defined if both sites exist. So an empty chain is
meaningless for `n` — not merely awkward to allocate — and both
`InitBADoubleTriDiagonalMatrix` and `InitIdentityBADoubleTriDiagonalMatrix` are
right to refuse it. (I removed the reject from the identity one and put it back;
the reason it looked safe — that it sets `acv`/`bv`/`bc` to constants inline
rather than through `SetBADoubleTriDiagonalMatrixConVal` — was an argument about
memory safety, not about what the type means.) The implementation agrees:
`SetBADoubleTriDiagonalMatrixConVal`
([BasicMath_BADoubleTriDiagonalMatrixVector.cpp:96](../../Core/BasicMath/BasicMath_BADoubleTriDiagonalMatrixVector.cpp#L96))
reads `val.b[0]`, `con.b[0]`, `val.a[1]` and `con.a[1]` unconditionally to build
`acvB`/`acvA`/`acvBA`, and is called from fourteen sites in
`BasicMath_MatrixTypeConversion.cpp`.

**This is a real difference between `d` and `n`, not an inconsistency** — see the
end of §1.

**Dimension zero, `b0 == 0` — DONE by Maurits in `8b112855`.** Distinct from an
empty chain: no bath *and* no observed site, so `G(w) = a0`, a constant with no
matrix behind it. Other types write it as `{{a0}, {}}`; in `d` it is `b0 = 0` with
two empty chains, `b0` being the hopping onto the observed site. So the guards
test `b0` and not `val.N + con.N` — an impurity with no bath also has two empty
chains but keeps `b0` nonzero and stays a 1×1 matrix.
`AndersonMatrixToDoubleTriDiagonalMatrix` accepts `M.N == 0` and emits it;
`…ToTriDiagonalMatrix`, `…ToAndersonMatrix` and `…ToListOfPoles` short-circuit on
it. `malloc(0)` returning `NULL` is no longer treated as failure in
`InitListOfPoles`, `InitTriDiagonalMatrix`, `InitIdentityTriDiagonalMatrix`.
Tested both ways; `a0` comes back bit for bit through every conversion.

**Three `N = 0` heap-buffer-overflows found with `-fsanitize=address`, none of
them introduced by that commit, two of them now on the `d` row. 1 and 2 are
fixed; 3 is open.**

1. `EvaluateAndersonMatrix`
   ([BasicMath_GreenFunctionArithmetic.cpp:4046](../../Core/BasicMath/BasicMath_GreenFunctionArithmetic.cpp#L4046))
   has no `N == 0` branch: it reads `A->a[1]` and `A->b[0]` unconditionally, and
   for `N == 0` those are `a` of length 1 and `b` of length 0. It is the **only**
   evaluator missing that branch — `EvaluateTriDiagonalMatrix` has it explicitly,
   `EvaluateListOfPoles` and `EvaluateBADoubleTriDiagonalMatrix` are safe by
   construction. Newly reachable because `d → a` now produces `H->N = 0`.
   **Fixed** with the same early return `EvaluateTriDiagonalMatrix` has.
2. `SetBADoubleTriDiagonalMatrixConVal`
   ([BasicMath_BADoubleTriDiagonalMatrixVector.cpp:96](../../Core/BasicMath/BasicMath_BADoubleTriDiagonalMatrixVector.cpp#L96))
   reads `val.b[0]`, `con.b[0]`, `val.a[1]`, `con.a[1]`.
   `ListOfPolesToBADoubleTridiagonalMatrix`
   ([BasicMath_MatrixTypeConversion.cpp:5218](../../Core/BasicMath/BasicMath_MatrixTypeConversion.cpp#L5218))
   pads an empty chain to `N = 1` with zeros when *one* chain is empty — exactly
   so both heads exist — but not when **both** are, and then calls it. Reachable
   as `ChangeType(G, "Nat")` on a dimension-zero `ListOfPoles`, with no `d`
   involved at all; `d → n` hits it because it routes through `l`. **Fixed**: with
   no poles the function returns early with `b0 = 0`, `a0` from the pole list, and
   two chains of **one** site carrying only zeros — the existing padding extended
   to the both-empty case. `SetBADoubleTriDiagonalMatrixConVal` is deliberately
   not called there, since it divides by `b0`. So `d → n` of a dimension-zero
   matrix has chains of length 1 where `d → d` has length 0: the §1 asymmetry,
   visible in the result.
3. `luaL_checkNatrualImpurityOrbitalMatrixType`
   ([TableToArray.cpp:1877](../../LuaInterface/lua_extensions/TableToArray.cpp#L1877))
   reads `buff->val.b[1]` and `buff->con.b[1]` when a chain has `N == 1`.
   **`LengthZero.Quanty` trips this today** on its "size 1+1 natural" case,
   unrelated to `d`. Those seven lines duplicate the second half of
   `SetBADoubleTriDiagonalMatrixConVal` **without** its `N == 1` guards. Proposed
   fix: split that function — leave the `b0 = sqrt(val.b[0]² + con.b[0]²)`
   normalization and head rescaling in `SetBADoubleTriDiagonalMatrixConVal`, move
   the `acv`/`bv`/`bc` derivation and its guards into a new
   `SetBADoubleTriDiagonalMatrixConValFromChains`, and have the reader call the
   new one. The reader cannot call the existing one as-is: it has already read
   `b0` from the table's `{a0, b0}` pair, and the normalization would overwrite it
   with `sqrt(val.b[0]² + con.b[0]²)`, which is 1 for a correctly normalized
   table. **Fixed** that way. `…FromChains` accepts a chain of length zero and
   gives it one site of zero energy and zero coupling — the padding factored out
   as `AddSiteToEmptyTriDiagonalChain`, which `SetBADoubleTriDiagonalMatrixConVal`
   also calls before its normalization, since that reads the two head hoppings.
   That normalization now skips the division when `b0 == 0`, which would otherwise
   make the two heads NaN. The reader warns when it is handed a chain of length
   zero; in the current suite that fires 3 times in `LengthZero.Quanty` and once
   in `Chop.Quanty`, all on `Nat` tables written with an empty chain.

**Fixing 3 uncovered two more of the same one-line kind, both fixed:**
`AndersonMatrixToCompactMatrix` and `AndersonMatrixToMatrix`
([BasicMath_MatrixTypeConversion.cpp:196](../../Core/BasicMath/BasicMath_MatrixTypeConversion.cpp#L196))
read `b[0]` of a dimension-zero Anderson matrix, and `ResponseFunctionGetB0`
([BasicMath_CompactMatrixVector.cpp:4150](../../Core/BasicMath/BasicMath_CompactMatrixVector.cpp#L4150))
does the same for both `t` and `a`. All now yield `b0 = 0`, which is the same
convention `d` uses: nothing couples to the matrix. ASan halts at the first
error, so each was hidden behind the previous one — the five are one bug repeated,
not five independent ones.

**A different class, also fixed: 14 self-overlapping `memcpy` calls.** They are
all one idiom — delete an element and shift the tail down over it — and `memcpy`
is undefined when the ranges overlap. The same shift is already written with
`memmove` three times in `BasicMath_BlockListOfPoles.cpp` (:400, :443, :1068), so
these are slips rather than a convention. ASan found two of them one run apart;
the rest came from scanning every `memcpy` in the tree whose destination and
source index the same array, which is what makes this a complete sweep rather
than another peel:

| file | sites |
|---|---|
| `BasicMath_AndersonMatrixVector.cpp` | 2, in `AndersonMatrixChopReduceDimension` |
| `BasicMath_ListOfPoles.cpp` | 2, in `MakeListOfPolesPhysical` |
| `BasicMath_BlockAndersonMatrix.cpp` | 2, in `BlockAndersonMatrixChopReduceDimension` |
| `BasicMath_BlockListOfPoles.cpp` | 8, in `ChopBlockListOfPolesReduceDimension` and the merge routines |

`memmove` is valid wherever `memcpy` is, so the conversion cannot change a
result; four of the eight in `BasicMath_BlockListOfPoles.cpp` move only the `a0`
block, where the gap provably exceeds the length, and were converted for
uniformity rather than necessity. One apparent hit is a false positive and was
left alone: `OperatorArithmetic.cpp:300` copies `Oset.Res.Determinant[i]` into
`O->Res.Determinant[i]`, two different objects, into memory just `malloc`ed.

All the `N = 0` reads landed on heap that happened to be zero, so the affected
tests passed and gave the right answer; only ASan saw them. `ChangeTypeDoubleTri`, `ChangeTypeFromDoubleTri`,
`LengthZero`, `ChangeType` and `NaturalOrbitals` all now run to completion under
ASan with zero errors, and the whole suite is unchanged number for number by the
five fixes. The **whole** suite now runs to completion
under ASan with zero errors.
Tested in both directions — an Anderson matrix with no bath at all, the `d → d`
copy with either or both chains empty, and hand-written `d` matrices of shape
`(0,7)`, `(7,0)` and `(0,0)` through the whole `d →` row and back.

**Eight helpers that do not exist** (§3): `DoubleTriDiagonalMatrixTypeToEigenSystemCompact`,
`DoubleTriDiagonalMatrixToEigenvalues`, `SymmetrizeDoubleTriDiagonalMatrix`,
`ShiftDoubleTriDiagonalMatrix`, `DoubleTriDiagonalMatrixTranspose`,
`DoubleTriDiagonalMatrixEffectiveHamiltonian`,
`DoubleTriDiagonalMatrixFromEffectiveHamiltonian`,
`SetDoubleTriDiagonalMatrixConVal`. Note `…ToEigenvalues` and
`…ToEigenSystemCompact` are now nearly free: the `d → l` reordering makes the
matrix plain tri-diagonal, so they are one `dstevd` each.

**§2(c)**, still there: `DoubleTriDiagonalMatrixAdd` loops `i < M->val.N - 1`, so
the last element of each chain is never added, and with `val.N == 0` the
`uint32_t` underflows into a runaway loop. `DoubleTriDiagonalMatrixScalarMultiply`
scales the matrix and then resets `a0 = 0`, `b0 = 1`; which of the two it should
do needs deciding before `ResponsefunctionMul` reaches it.

**21 dispatch sites with no `d` case**, regenerated by scanning for functions that
branch on `'n'` and not on `'d'`:

`ResponseFunctionInvertEnergy`, `ResponseFunctionTranspose`,
`ResponseFunctionEffectiveHamiltonian`, `ResponseFunctionFromEffectiveHamiltonian`,
`ResponsefunctionChop`, `ResponsefunctionMul`
(`BasicMath_GreenFunctionArithmetic.cpp`); `ResponseFunctionToDensityMatrix`,
`ResponseFunctionAddSelfEnergyDenseViaTri`,
`ResponseFunctionAddSelfEnergyOnGridViaG0Tri`,
`ResponseFunctionAddSelfEnergyViaAndLP` (`Core/ResponseFunction.cpp`);
`LuaResponseFunctionEigensystem`, `ResponseFunctionToLuaBitMask`,
`LuaResponseFunctionToTightBinding`, `LuaResponseFunctionToOperator`,
`LuaResponseFunctionLen`, `LuaResponseFunctionIndexMu`,
`LuaResponseFunctionIndexType`, `LuaResponseFunctionIndexName`,
`LuaResponseFunctionNewIndexMu`, `LuaResponseFunctionNewIndexName`
(`LuaInterface/ResponseFunction.cpp`); `LuaReadSingleSelfenergy`
(`TableToArray.cpp`).

The six index metamethods are the cheap ones and are user visible: `G.type`,
`G.mu` and `G.name` do not work on a `d` userdata. Several of these `switch`es
end in a `default:` that prints and continues, so a half-finished `d` gives wrong
numbers rather than stopping — add `d` to the guard before adding it to the body.

**§3.1 beyond `d`:** `ListOfPolesCopyVal` / `…Con` still use the linear ramp, so
the non-block `→ n` path still disagrees with `EffectivekBT`'s documented meaning,
and five structs still carry `Epsilon`.

### Still open, block

`ConvertBlockMatrix` has no `d` row or column at all
([BasicMathBlock.cpp:321](../../Core/BasicMathBlock/BasicMathBlock.cpp#L321) is
its "not implemented" default). Nine stubs, all in
`BasicMath_BlockTypeConversion.cpp`:

`BlockListOfPolesToBlockDoubleTriDiagonalMatrix`,
`BlockAndersonMatrixToBlockDoubleTriDiagonalMatrix`,
`BlockTriDiagonalMatrixToBlockDoubleTriDiagonalMatrix`,
`BlockNaturalImpurityOrbitalMatrixToBlockDoubleTridiagonalMatrix`,
`BlockDoubleTriDiagonalMatrixToBlockTriDiagonalMatrix`,
`BlockDoubleTriDiagonalMatrixToBlockAndersonMatrix`,
`BlockDoubleTriDiagonalMatrixToBlockNaturalImpurityMatrix`,
`BlockDoubleTriDiagonalMatrixToBlockListOfPoles`, and
`BlockDoubleTriDiagonalMatrixToCompactMatrix` **by value** — a second overload of
a name that already has a working by-pointer version in
`BasicMath_CompactMatrixVector.cpp`. Delete the by-value declaration and stub
rather than implementing it; a caller who passes by value binds to the stub
silently.

Exists for the block type: `Init`, `Free`, `Copy`, `Evaluate`,
`…ToCompactMatrix[WithDeflation]`, `…ToSpectra[WithBroadening]`, and the two
`Chop` routines — but both `Chop`s are `//TODO:` bodies that return 0 without
doing anything. Missing: `…Part` (the `Nv`/`Nc` primitive), `…Add`,
`…ToImpurityHamiltonian`.

**`a → d` is written and works.** See "Decided: how `a → d` is built in block
form" in §8 for the algorithm and the reasoning.
`BlockAndersonMatrixInverseToBlockListOfPoles` and
`BlockAndersonMatrixToBlockDoubleTriDiagonalMatrix(A, B, kbT, SingularValue, MergeTolerance)`
exist, the `a → d` cell of `ConvertBlockMatrix` is wired, and
`ResponseFunction.ChangeType(blockAnd, "DoubleTri", {{"EffectivekBT", …}})`
reproduces the original `G` to 1e-13 on a 120-point grid at `Gamma` 0.3 and 0.02,
including `mu` outside the bath in both directions:

| case | val.N | con.N | max relative error |
|---|---|---|---|
| kBT finite, 20 bath blocks of 2 | 23 | 23 | 6.7e-14 |
| kBT zero | 20 | 20 | 1.8e-13 |
| mu off centre | 31 | 40 | 8.7e-14 |
| valence chain empty | 0 | 40 | 4.3e-14 |
| conduction chain empty | 40 | 0 | 4.2e-14 |
| 3 bath blocks | 6 | 6 | 6.4e-14 |
| 1 bath block | 1 | 2 | 3.7e-14 |

Also `d → d` and the lua table round trip, both bit exact.

Two supporting fixes were needed. `BlockListOfPolesToBlockTridiagonalMatrix` now
returns early for a pole list with no poles — it used to hand the block Lanczos a
starting space of dimension zero, which produced a matrix claiming one more block
than it owned. And `ConvertResponsefunction` now **passes on**
`ConvertBlockMatrix`'s return value instead of discarding it, so a `d →` cell that
is not implemented reaches lua as an error rather than as a response function whose
`data` was never allocated. One new helper,
`BlockListOfPolesRemoveZeroWeightPoles`, exact and with no tolerance.

**Chain lengths exceed `NBath` blocks, and the cap can bind.** With 20 bath blocks
of size 2 the chains come out 23 and up to 40 blocks long, so `NBath` would have
been the wrong bound and `NBath·BlockSize` is right — the rank-1 weights make the
starting block deflate, blocks shrink, and the chain lengthens. But `con.N = 40` in
the `mu` off centre case is **exactly** the cap, so there the chain was truncated
rather than exhausted. `G` is still reproduced to 8.7e-14, which is precisely the
situation [TODO_ConversionAccuracy.md](TODO_ConversionAccuracy.md) warns about: the
tail of a chain barely shows on a uniform grid. Worth measuring before trusting long
block chains.

**Both `DoubleTri` readers now convert on all four branches**, block and non-block:
`Tri`, `And`, `Nat` and `ListOfPoles`. They used to `malloc` a buffer, convert
nothing and return it uninitialised without a word. Each is now one call to the
corresponding conversion, with the `malloc` checked and a `luaL_error` on failure.
`mu` comes from the source table and `kBT` is 0 — none of those table formats
carries a temperature, so the split is a sharp step at `mu`.

**All eight are unreachable from any current call site.** Every caller of either
reader — `luaL_checkResponsefunction` and the `"DoubleTri"` branches of the `Tri`,
`And` and `Nat` readers — passes a table that is already of type `DoubleTri`, so
only that branch ever runs. The value is removing the uninitialised-memory trap,
not fixing a live bug.

**The non-block reader `luaL_checkDoubleTriDiagonalMatrixType` had the identical four
dead branches, and all four are now written**, since every non-block conversion into
`d` exists: `TriDiagonalMatrixToDoubleTriDiagonalMatrix`,
`AndersonMatrixToDoubleTriDiagonalMatrix`,
`BADoubleTriDiagonalMatrixToDoubleTriDiagonalMatrix` and
`ListOfPolesToDoubleTridiagonalMatrix`. `mu` comes from the source table and `kBT` is
0 — none of those table formats carries a temperature, so the split is a sharp step at
`mu`, the same choice the block reader's `And` branch makes.

**How to reach them for testing.** A temporary hook at the top of the type dispatch
in `luaL_checkResponsefunction`:

```c
lua_getfield(L, narg, "AsDoubleTri");
if(lua_isnil(L,-1)==0)
{
  lua_pop(L, 2);
  buff->type='D';
  if(buff->Block==true) buff->data=(void*)luaL_checkBlockDoubleTriDiagonalMatrixType(L, narg);
  else                  buff->data=(void*)luaL_checkDoubleTriDiagonalMatrixType(L, narg);
  return buff;
}
lua_pop(L, 1);
```

Then a lua table carrying `AsDoubleTri=1` is read by the `DoubleTri` reader whatever
its own `type` says. Both sets were checked this way and the hook reverted. `G` after
the round trip, against the source:

| source | non-block | block |
|---|---|---|
| `And` | 4.2e-16 | 5.2e-16 |
| `Tri` | 2.5e-16 | 1.1e-15 |
| `ListOfPoles` | 2.2e-15 | 1.6e-15 |
| `Nat` | 4.7e-16 | 5.4e-15 |

`mu` came across in all eight.

**Tested by `TestScripts/Responsefunction/ChangeTypeDoubleTriBlock.Quanty`**,
registered in `All.Quanty` and in `.gitlab-ci.yml`. Beyond the grid comparison it
pins down three things:

- `val.a[0] + con.a[0] == A1` **exactly**, element for element, and the two chain
  heads carry the whole hybridization weight `Σ_k B_k B_k^T` to 1e-12 relative.
- **Rotating the bath by an exactly orthogonal reflection.** `(c,s) = (0.6,0.8)`
  satisfies `c*c + s*s - 1 == 0` in double, so the test is not limited by the
  precision of `U`. For a **degenerate** bath block the two `d` matrices come out
  **bit identical**, every element — that is the property the rank-1 weights and the
  merge exist for, since the eigenvectors inside a degenerate subspace are otherwise
  arbitrary. For a **non-degenerate** bath the statement is weaker and the test says
  so: `v v^+` is sign-independent so the pole list is invariant in exact arithmetic,
  but `dsyev` on `U^T A U` perturbs the eigenvalues at 1e-16, the deflation inside
  the block Lanczos then cuts elsewhere, and the chain comes out a **different
  length** — 4 blocks against 6 in the case tested. `G` still agrees to 1.8e-14. So
  the chain is not a basis-invariant object for a non-degenerate bath, only `G` is.
  That is the "chain tail is not determined" behaviour of
  [TODO_ConversionAccuracy.md](TODO_ConversionAccuracy.md) appearing in the length
  rather than in the coefficients.
- the four unimplemented `d →` cells each raise a lua error, which is the
  regression test for the `ConvertBlockMatrix` → `ConvertResponsefunction` →
  `ChangeType` failure path.

**`t → d` is written too, through the Anderson form**, the same route the scalar
type takes. Both steps need the bath diagonalized, so composing costs one dense
`O((N·BlockSize)³)` diagonalization — what a direct route would cost as well — plus
the block-local ones of the second step. It reproduces `G` to 5e-13 for 20 bath
blocks of size 2, and `mu` outside the bath empties a chain as it should.

`BlockTriDiagonalMatrixToBlockAndersonMatrix`'s `!!Does not work as intended yet!!`
header turns out to be about the *shape* of its output and not about `G`: it
reproduces `G` to 5e-13, and its two complaints — hoppings mixed between the
orbitals of a block, bath blocks reordered because the eigenvalues are sorted — are
both irrelevant on this path, since `a → d` diagonalizes every bath block again on
its way to the hybridization poles.

**But it dropped `mu`.** The `N == 0` and `N == 1` branches set `B->mu = ATri.mu`
and the general branch did not, so anything converted through it came out with the
chemical potential at zero and a Fermi split further downstream then happened in the
wrong place: `t → d` with `mu = ±6` gave the same chain lengths as `mu = 0` instead
of emptying a branch. Fixed, one line. Nothing in the suite moved by a digit, which
means no existing test converts `t → a` or `t → n` with a nonzero `mu` — the bug was
invisible to it.

**The six index metamethods are done.** `G.mu`, `G.type` and `G.name`, read and
written, now work on a `d` userdata, block and non-block — `LuaResponseFunctionIndexMu`,
`…IndexType`, `…IndexName`, `…NewIndexMu`, `…NewIndexName`. Exercised in
`__index.Quanty`, which already loops over one response function of every type; the
`Chop` loop there skips `d`, since `ResponsefunctionChop` still has no `d` case.

That work started because a probe died reading `G.mu`, and it turned up a second
defect worth knowing about. **Nine `luaL_error` calls in
`LuaInterface/ResponseFunction.cpp` formatted a `char` with `%s`** —
`luaL_error(L, "ResponseFunction of undifined type %s\n", a->type)` and two variants
on `M->type`. A `char` promoted into a `%s` is read as a pointer, so the error path
itself segfaulted instead of raising. That is why the probe *died* rather than
reporting. All nine are now `%c`. Worth remembering when adding the remaining `d`
cases: several of the 21 dispatch sites end in a `default:` like this, so the failure
mode is a crash and not a message.

**`l → d` and `n → d` go the same way**, `l → a → d` and `n → a → d`, with the
signature of the family: `(source, B, kbT, SingularValue, MergeTolerance)`. Both
stubs had different parameter lists — `l → d` had `NBlockMaxVal`/`NBlockMaxCon`,
which have no role once the chain lengths come from the bath size — and both were
changed to match.

For the pole list the detour is not avoidable: `l` holds the poles of `G`, while
the two chains of `d` are the two halves of the hybridization and sit in the
**denominator**, so getting from one to the other means inverting, and that is what
the Anderson form is. Worth keeping straight from
`BlockListOfPolesToBlockNaturalImpurityOrbitalMatrix`, which splits the poles of `G`
directly and never inverts — it looks like the same function and is a different
calculation. `n → d` is not a change of basis either, for the same §1 reason.

Measured, against the source, on the 120-point grid at `Gamma` 0.3 and 0.02:

| route | worst error |
|---|---|
| `a → d`, `t → d` | 5e-13 |
| `l → d` | 4.6e-13 |
| `n → d` | 1.3e-12 |

`n → d` composes four conversions (`n → l → t → a → d`) and shows it, so the test
allows 1e-10 on those two rows rather than the 1e-11 used elsewhere, with the
measured numbers in the printout so the headroom is visible. `mu` outside the bath
empties a chain on every route.

Both intermediates were measured before being built on: block `l → a` and `n → a`
reproduce `G` to 5e-13, like `t → a`.

**The whole inbound block row is now done: `l`, `t`, `a`, `n` and the copy all reach
`d`.** `Chop` for both `d` types followed, so `__index.Quanty` no longer skips them.
Still to do: the `d →` cells — `d → a` first, see below — and
`Part`/`Add`/`ToImpurityHamiltonian`.

**The non-block insights carry over directly**, and are the reason the block row
should be cheap: a block `d` read as `val_N … val_1, imp, con_1 … con_N` is
**block tri-diagonal**, and interleaved behind the impurity it is a band matrix of
half bandwidth `2·BlockSize`. `dsbtrd` is already used on block tri-diagonal
matrices in
[BasicMath_BlockTypeConversion.cpp:1425](../../Core/BasicMathBlock/BasicMath_BlockTypeConversion.cpp#L1425).
`BADoubleTriDiagonalMatrixToTriDiagonalMatrix` is still on the dense O(N³) route
and the same trick fixes it.

### Not testable yet

The deflated block path. Nothing can build a deflated block `d`: the lua reader
derives `BlockSize` from the matrices and takes no per-block list, and both
`BlockListOfPolesToBlockDoubleTriDiagonalMatrix` and
`BlockDoubleTriDiagonalMatrixChopReduceDimension` are stubs. `WithDeflation` is
verified only to agree with the plain version when there is nothing to deflate.

### Tooling notes, learned the hard way

- **Quanty lua options are a list of `{"name", value}` pairs.** A table with named
  keys (`{Emin=-4}`) has an empty array part, passes `LuaTestOptionsList` and is
  then silently ignored. Write `{{"Emin",-4},{"GammaL",0.3}}`.
- Piped stdout is block buffered and the buffer is lost when the process dies, so
  a crash looks silent. Re-run under `script -q /dev/null` to see the last
  `printf`.
- In a lua `ListOfPoles` table `a[1]` is the offset `a0` and the poles start at
  `a[2]`, while `b[1]` is the first weight. A `Tri` or `And` table has `#a = N+1`
  and `#b = N`, `a[1]` and `b[1]` being `a0` and `b0`.
- `Matrix.Sub(M,n)` returns the leading `n × n` **matrix**, not an element, and is
  1-based. `Complex.Abs` and `Complex.Conjugate` do not exist, but `Complex.Re`,
  `Complex.Im`, `Complex.New`, `Complex.I` and `Complex.Exp` do.
  `Spectra.ToTable(sp)` gives `sp[ispectrum][ipoint] = {energy, value}`.
- **Read a number with `Complex.Re(x)`, never `x.Re`.** A method modifies its
  object in Quanty and a function returns one, so `x.Re` is the wrong spelling;
  and it fails outright on a plain lua number, which is what an element of a block
  matrix becomes when its imaginary part is exactly zero. `Complex.Re` and
  `Complex.Im` take either kind. See the entry in [TODO.md](TODO.md).
- `ToSpectra` broadening is `GammaL`, not `Gamma`.
- **Build outside the source tree.** `g++ … -o` takes any path and nothing in the
  build or the tests depends on where the binary lands, so put it in a scratch
  directory (or `build/`, which `.gitignore` already covers) rather than leaving a
  `Quanty` in the repo root. CI builds `Quanty.exe` in the root and runs
  `./Quanty.exe "TestScripts/Responsefunction/NAME.Quanty"` from there; both that
  and running from inside the test directory work, since `All.Quanty` finds its own
  directory with `debug.getinfo`.
- Running the test suite does **not** normally rewrite the `.svg` files under
  `TestScripts/*/OutputPlots/`; they come out byte-identical, from the repo root or
  from inside the test directory, at `-O2` or at `-O1` with ASan. One did change
  once, during the window when `AndersonMatrixToCompactMatrix` and
  `ResponseFunctionGetB0` were still reading uninitialised memory for a
  dimension-zero matrix — an uninitialised read is the one thing that would make a
  plot change intermittently, and it has been stable since those were fixed. So
  check `git status` after a run rather than reflexively `git checkout --`ing the
  plots: a changed `.svg` is now a signal, not noise.

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

### An empty branch is legal for `d` and meaningless for `n`

A direct consequence of the above, and the reason the two types differ on
`Init`. In `d` the two branches are nested in **one** denominator, and each
enters through `Gval` / `Gcon`, which for an empty chain is just its offset
`a[0]`: `EvaluateTriDiagonalMatrix` returns `a[0]` when `N == 0` and never touches
`b`. So `val.N = 0` means "this branch is absent, and contributes only its share
of the impurity energy" — a well-defined matrix of dimension `1 + con.N`, which is
exactly what the Fermi split produces whenever `mu` lies outside the bath.

In `n` the observed orbital is a **linear combination of the two chain heads** —
that is what `T` in the struct comment does and what the cached `acvB`/`acvA`/
`acvBA`/`bvB0`/`bcB0`/`bvA0`/`bcA0` fields hold. That combination is only defined
if both heads exist, so a chain of length zero is not a degenerate case of the `n`
form; it is outside it. Both `BADouble…` constructors reject it, correctly, and
that is not a guard to be tidied away.

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

### 3.0 Convention sweep: the shared impurity site — DONE

`avc0 = val.a[0] + con.a[0]`, the plain sum, everywhere. Every consumer that read
`(val.a[0] + con.a[0]) / 2` has had the `/ 2` deleted, and every producer writes half of the
impurity energy into each head slot:

| file | function |
|---|---|
| `BasicMath_DoubleTriDiagonalMatrixVector.cpp` | the file header comment, `DoubleTriDiagonalMatrixInverse`, `Real…MatrixVector`, `Complex…MatrixVector`, `RealVector…MatrixVector`, `ComplexVector…MatrixVector` |
| `MatrixOperator.cpp` | `DoubleTriDiagonalMatrixToAndersonImpurityOperator` and `SecondDoubleTriDiagonalMatrixToAndersonImpurityOperator`, four sites each |
| `StructDef.cpp` | the `SecondDoubleTriDiagonalMatrixType` comment, which said `/2` about the `DoubleTriDiagonalMatrixType Dt` it wraps |
| `BasicMath_MatrixTypeConversion.cpp` | `AndersonMatrixToDoubleTriDiagonalMatrixPlusT` now writes `0.5 * M.a[1]` into each head, as `AndersonMatrixToDoubleTriDiagonalMatrix` already did |
| `Test/LanczosDyson.cpp` | the one site outside `Core/` |

`DoubleTriDiagonalMatrixToTriDiagonalMatrix` needed no edit of its own: it goes through
`RealDoubleTriDiagonalMatrixVector`, so its first Lanczos coefficient is now the same
impurity energy that `DoubleTriDiagonalMatrixToCompactMatrix` puts in `VRe[0]`.

**Deliberately not touched**, and still on the average: `MapSecondoubleTridiagonalMatrixToBA`
([BasicMath_BADoubleTriDiagonalMatrixVector.cpp:136](../../Core/BasicMath/BasicMath_BADoubleTriDiagonalMatrixVector.cpp#L136)),
which reads `A.Dt.val.a[0]` — a `DoubleTriDiagonalMatrixType` head — but lives in the
`BADouble…` file and writes a `BADouble…` field. It is one token, and it belongs with the
`BADouble…` pass rather than this one. `DoubleTriDiagonalMatrixToSpectrumOld` averages both
the head and `b0^2`; `…Old` functions are kept for reference and are out of scope.

**Not verified by any test**, because none of these functions has a live caller yet: the four
matrix–vector products are reached only through `DoubleTriDiagonalMatrixToTriDiagonalMatrix`,
which `ConvertResponsefunction` does not call, and the two operator builders only from
`Test/`, which is not in the build. The sweep is verified by inspection and by now agreeing
with `Evaluate` and `…ToCompactMatrix`, which are tested. The first thing that will exercise
it is the `d →` row.

Note that deleting the `/ 2` does **not** make `DoubleTriDiagonalMatrixInverse` correct: it
still builds the conduction block from `M.val` instead of `M.con`, which is §2(b).

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
- `BADoubleTriDiagonalMatrixToDoubleTriDiagonalMatrix` — routes via the Anderson form
- `BADoubleTriDiagonalMatrixToAndersonMatrix` — rewritten. The bath of an `n` object, seen
  from the observed orbital, is **one tri-diagonal eigenvalue problem**: drop the observed
  orbital B and the extra natural orbital A is left with exactly two neighbours, the head of
  each chain, so the whole minor is a single chain read `av_N … av2 A ac2 … ac_N`. It used to
  build the dense matrix and call `CompactMatrixToAndersonMatrix`, which is O(N^3) with
  `2*N^2` of *stack* workspace — it segfaults above about N = 600 — and which also lost `a0`
  and `b0`, so the result did not reproduce `G` at all (that half is fixed in
  `CompactMatrixCopy`). `n → d` is
  about twice as fast as the old route through the poles, which needed three tri-diagonal
  eigenvalue problems and a tri-diagonalization where this needs one.
- `DoubleTriDiagonalMatrixToAndersonMatrix` — the first of the `d →` row. Read the matrix as
  `val_N … val_1, impurity, con_1 … con_N` and it is one ordinary tri-diagonal chain with the
  impurity in the middle, so deleting the impurity leaves the two chains **disconnected**: the
  Anderson bath is the two chain spectra side by side, and the impurity reaches eigenvector `k`
  of a chain only through that chain's head, with hopping `b[0]·c_k`. Two tri-diagonal
  eigenvalue problems of size `val.N` and `con.N`, against roughly four times that for a route
  through `l`, which would form `G`'s poles only to undo them. Dimension preserving, bath
  ordered, zero-hopping levels kept for `Chop` to remove. It was the stub called
  `…ToAndersonMatrixMatrix`, with `Matrix` twice; renamed, and its one caller — the lua reader
  that reads a `DoubleTri` table as an `And` — follows
- `DoubleTriDiagonalMatrixToTriDiagonalMatrix` — rewritten as a **band reduction**. Behind the
  impurity, interleave the chains as `impurity, val_1, con_1, val_2, con_2, …` and every
  coupling is at distance 1 or 2, so the matrix is a symmetric band matrix of half bandwidth 2
  and the starting vector is `e_1`, with no rotation needed. `dsbtrd` annihilates only in
  planes ≥ 2, so it leaves `e_1` alone and returns the Lanczos chain from the impurity —
  computed by orthogonal similarity rather than by the three-term recurrence, hence backward
  stable, with no re-orthogonalization and `4N` doubles of memory. Measured against
  `d → a → t`: 4 to 8 times faster and slightly more accurate (5e-14 against 7e-14…2e-13 on
  `G`). `NTriDiag` is gone, the full chain always comes back and `Chop` trims it. `val.N == 0`
  or `con.N == 0` is a real physical case and is handled separately, the minor then being a
  single chain that is already tri-diagonal. The five `d` matrix–vector routines this used to
  go through are now uncalled and renamed `…Old`
- **`DoubleTriDiagonalMatrixToCompactMatrix[WithDeflation]`** and
  **`BlockDoubleTriDiagonalMatrixToCompactMatrix[WithDeflation]`** — done, see §0
- the `d → d` copy cell in `ConvertResponsefunction` and `ConvertBlockMatrix`, which is what
  `ChangeType(…,"Table"/"Userdata")` uses

To write:

| # | function | note |
|---|---|---|
| 1 | `ListOfPolesToDoubleTriDiagonalMatrixOccupation` | only if the occupation-constrained variant is wanted, cf. `ListOfPolesToBADoubleTridiagonalMatrixOccupation` — note that one stores `Abs(mu − muOccupation)` in the field, a third meaning again (§11) |
| ~~2~~ | ~~`DoubleTriDiagonalMatrixToListOfPoles`~~ | **done**. The poles of `G` itself, not of the hybridization, so this is **not** the two chain spectra — it is the spectrum of the whole matrix. Read as `val_N … val_1, impurity, con_1 … con_N` that matrix is plain tri-diagonal, half bandwidth 1, so it is one eigenvalue problem of size `1+val.N+con.N`, with the weights read off the impurity component, which sits in the middle rather than at the end. The interleaved ordering `d → t` needs is only for putting the observed site first, which a Lanczos chain requires and an eigendecomposition does not |
| ~~3~~ | ~~`DoubleTriDiagonalMatrixToTriDiagonalMatrix`~~ | **done**, see the entry above |
| ~~4~~ | ~~`DoubleTriDiagonalMatrixToBADoubleTriDiagonalMatrix`~~ | **done**, `d → l → n`. The two forms combine their chains differently (§1), so there is no rotation between them; the poles of `G` are the only thing they agree on. Going through them also means `d → n` splits at `mu` exactly as `l → n` does, by construction |
| 5 | `DoubleTriDiagonalMatrixCopy` | exists |

Each of these wants its own thought. `d → a` is not a template for the rest: it works because
the *hybridization* splits into the two chains, which is a statement about the minor, not
about `G`. `d → t` is not one either: it works because the *whole* matrix, impurity included,
is narrow-banded in an interleaved ordering.

### A property of long chains, found while checking `d → t`

The map from a measure to its recurrence coefficients is **exponentially ill conditioned in
the chain index**. The band reduction and the route through `a` are different, both backward
stable, algorithms for the same chain, and they agree to 3e-14 for a chain of 41 sites — but
part company at site 31 of 61, and at site 27 of 115, where they end up differing by 5. Both
still reproduce `G` to 1e-13, because the deep tail of a chain has almost no influence on it.

So a long chain reproduces `G`, but its tail coefficients are not determined at double
precision, and two correct routines will not agree on them. Compare `G`, never the deep
coefficients — and note that "reproduces `G`" here means on a uniform grid at two broadenings,
which is the whole of what any of these conversions has been measured against; see
[TODO_ConversionAccuracy.md](TODO_ConversionAccuracy.md). It is also why `Chop` matters: the tail hoppings of the 121-site chain in that
test are ~1e-4 and carry nothing.

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

## 8. The block type — the struct and the evaluator are done, the conversions are not

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

### Decided: how `a → d` is built in block form

Settled with Maurits before writing any of it. Three corrections to earlier
readings of this document are folded in.

**The formula is `G = a0 + b0^*·(w − Gval − Gcon)⁻¹·b0^T`** — `b0^*`, the
elementwise conjugate, on the left and `b0^T` on the right. That is what
[StructDef.cpp:709](../../Core/StructDef.cpp#L709) says for
`BlockTriDiagonalMatrixType`, what the `BlockAndersonMatrixType` comment says
with `B0^*`/`B0^T`, and what `EvaluateBlockDoubleTriDiagonalMatrix` implements as
`CompactMatrixMultiply3('*','N','T', b0, det, b0)`. It costs nothing on the real
path. The same `^*`/`^T` question runs through `DiagonalMatrixBlockTridiagonalize`
(contract written `psi . 1/(w−A) . psi^C`) and the `CompactMatrixTranspose` that
`BlockListOfPolesToBlockTridiagonalMatrix` applies to its starting block, so the
`QComplex` path is to be **verified against `Evaluate`**, not reasoned about.

**Failures return 1 and the caller listens.** Every new function returns 1 on
failure, `ConvertBlockMatrix`'s `default:` returns 1 instead of 0 with
`Mout->data` unassigned, and the four dead branches of
`luaL_checkBlockDoubleTriDiagonalMatrixType` — which today `malloc` and return
*uninitialised memory silently* — `luaL_error` until they are written. A failed
conversion has to reach the lua side as an error.

**The rank-1 outer product is the mechanism, not a cost.** For a bath block with
a degenerate group `D`, `Σ_{m∈D} (b·u_m)(b·u_m)^† = b·P_D·b^†` with `P_D` the
spectral projector — basis-independent, so the arbitrary rotation LAPACK picks
inside the degenerate subspace cancels exactly, while each individual rank-1 term
depends on it. Spectral weight is the only additive object, which is why a
`BlockListOfPoles` stores `W = B B^†` and not `B`. Maurits' case:
`{{1,0,-0.811957,0.583717},{0,1,-0.583717,-0.811957},{-0.811957,-0.583717,2,0},{0.583717,-0.811957,0,2}}`
with `BlockSize = 2`, one doubly degenerate bath block `2·I` and an orthogonal
`b`, so `P_D = I`, the merged weight is `b b^† = I`, `sqrt(I) = I`, and the chain
starts from the symmetric coupling of the rotated form
`{{1,0,1,0},{0,1,0,1},{1,0,2,0},{0,1,0,2}}` whatever the eigensolver returned.
This matters at the µeV level: magnetic splittings on a band of tens of eV.

**One new primitive, then existing code.** The non-block precedent is
`AndersonMatrixInverseToListOfPoles`
([BasicMath_MatrixTypeConversion.cpp:809](../../Core/BasicMath/BasicMath_MatrixTypeConversion.cpp#L809)),
which sets `a0 = M.a[1]` and lists the bath as `b_i²` at `a_{i+2}` — the pole
list of *the thing in the denominator*. The block twin is what is missing:

```
BlockAndersonMatrixInverseToBlockListOfPoles(A, Delta)
  for each bath block k = 1 … A->N-1:
    CompactMatrixEigensystem(a[k+1]) -> U_k, e_k       // block-local, O(BlockSize^3)
    V_k = b[k] . U_k
    emit BlockSize poles at e_{k,m} with weight v_{k,m} v_{k,m}^dag
  Delta.a0 = A->a[1];  Delta.mu = A->mu
```

Cost `O(NBath·BlockSize³)`, no dense assembly. Deliberately **not**
`BlockAndersonMatrixToBlockListOfPoles`, which builds the whole matrix dense and
diagonalizes it, `O((N·BlockSize)³)`, and returns the poles of `G` — the right
source for `n`, the wrong one for `d`, whose chains are the hybridisation seen
from the impurity and sit in the denominator.

**Written and verified.** `CompactMatrixEigensystem` returns `fun` with
`fun^*·A·fun^T` diagonal, so the eigenvectors are the **rows** of `fun` and
`U_k = fun^T`; the code forms `V_k = B_k·U_k` as
`CompactMatrixMultiply('N','T', &Bk, &fun, &Vk)` and reads `v_{k,m}` off as column
`m`. The poles are ordered on the way out, which is what the merge expects.
`Abath` and `Bk` are views on the Anderson matrix's own storage and are never
freed. `A->N == 0` gives an empty pole list with `a0 = 0`, per the `b0 == 0`
convention.

It has no caller yet, so it was verified by routing the block `a → l` cell of
`ConvertBlockMatrix` through it temporarily, checking, and reverting. Against
`A1 + Σ_k B_k (w − A_{k+1})⁻¹ B_k^+` computed independently in lua with
`Matrix.Inverse`, on the real axis away from every pole:

| case | max difference |
|---|---|
| `BlockSize` 1, 4 bath blocks — must match the scalar convention | 8.9e-16 |
| `BlockSize` 2, 3 bath blocks, none diagonal | 4.4e-16 |
| `BlockSize` 3, one bath block | 5.6e-17 |
| a doubly degenerate 2×2 bath block | 1.1e-16 |

And the property the rank-1 form exists for. Maurits' example to 64 digits, the
bath block `2·I` with hopping `{{c,s},{s,-c}}`:

```
c = 0.5830044176323191546188243313891656343050509643544691478828829605
s = 0.8124689834210290580594052941697557903194370438951662443780806694
```

The two rank-1 weights come out `{{c²,cs},{cs,s²}}` and `{{s²,-cs},{-cs,c²}}` —
basis-dependent, as expected — and their sum has an off-diagonal that is
**exactly** `0` and a diagonal that equals `c*c + s*s` **bit for bit**. So the
merge recovers `b·b^+` with no error at all; the `±cs` cancel exactly.

Those doubles satisfy `c² + s² = 1 - 1 ulp` (`-1.11e-16`), so the merged weight
differs from the unrotated form's exact identity by exactly that — the input's own
precision, not the conversion's. **Two consequences for the test:** assert against
`b·b^+`, never against the identity; and note that
`(c,s) = (0.6, 0.8)` and `(0.28, 0.96)` *are* exactly orthogonal in double
(`c*c + s*s - 1 == 0`), so with a Pythagorean reflection the whole merged pole
list is **bit identical** between the two bath bases, with no tolerance at all.
Verified for `(0.6,0.8)`, `(0.28,0.96)` and `(-0.8,0.6)`. That is the form the
bath-rotation test should take.

Only the real path is exercised. The `QComplex` branch is written but nothing
builds a complex block Anderson matrix in these tests, so it is unverified — and
per the `b0^*`/`b0^T` note above that is exactly where a convention slip would
hide.

```
BlockAndersonMatrixToBlockDoubleTriDiagonalMatrix(A, B, kbT, SingularValue, MergeTolerance)
  BlockAndersonMatrixInverseToBlockListOfPoles(A, &Delta)
  BlockListOfPolesMergeDegeneratePoles(&Delta, MergeTolerance)   // restores the symmetry
  CopyBlockListOfPolesVal(Delta, &Dval, mu, kbT)                 // weights times f
  CopyBlockListOfPolesCon(Delta, &Dcon, mu, kbT)                 // weights times 1-f
  Dval.a0 *= 0.5;  Dcon.a0 *= 0.5                                // the shared site, §1
  BlockListOfPolesToBlockTridiagonalMatrix(Dval, &B->val, NBath*BlockSize, SingularValue)
  BlockListOfPolesToBlockTridiagonalMatrix(Dcon, &B->con, NBath*BlockSize, SingularValue)
  B->a0 = A->a[0];  B->b0 = A->b[0];  + mu, kBT, BlockSize, QComplex, Name
```

Two things fall out. `CopyBlockListOfPoles*` copies `a0` and scales only the
weights, and `BlockListOfPolesToBlockTridiagonalMatrix` writes the list's `a0`
into `B->a[0]`, so halving `a0` on the two lists puts `0.5·A1` into `val.a[0]`
and `con.a[0]` with no special-casing. And the Fermi factors are applied to poles
that have already been merged, so degenerate partners get an identical `f` by
construction and the split cannot break the symmetry the merge just restored.

**Signature.** `SingularValue` and `MergeTolerance` are parameters, not
constants. `ConvertBlockMatrix` carries only `T` and `SingularValue`, so
`MergeTolerance <= 0` means "use the default" and `ChangeType` keeps working
until an explicit lua option exists.

**`MergeTolerance` — the existing default would erase real physics.**
`BlockAndersonMatrixToBlockListOfPoles` passes
`G->BlockSize * G->N * DBL_EPSILON * BandWidth`. At `BlockSize = 10`, a 50-block
bath and a 60 eV band that is `5100 · 2.2e-16 · 60 ≈ 7e-11 eV ≈ 70 µeV`; even
`BlockSize = 2` with 100 poles gives ≈ 3 µeV. That is the scale of the Zeeman
splittings this type has to resolve. The error we are hiding under is
`O(ε·‖a_k‖)` from *independent* `BlockSize × BlockSize` diagonalizations — it does
not grow with the number of poles — so the default should be a few·`ε`·bandwidth,
about `1e-13 eV`, nine orders below µeV. Belongs in
[TODO_ConversionAccuracy.md](TODO_ConversionAccuracy.md) as a measured quantity.
Keeping it tight also defuses the `//TODO: do not shift spectral weight through
the chemical potential` in `BlockListOfPolesMergeDegeneratePoles`: we merge
*before* splitting, and at `1e-13 eV` two merged poles straddling `mu` differ in
`f` by `1e-13/kBT`.

**`NBlockMax = NBath · BlockSize`, and why the obvious `NBath` is wrong.**
`CompactMatrixBlockBandDiagonalizeReal`
([BasicMath_CompactMatrixVector.cpp:5417](../../Core/BasicMath/BasicMath_CompactMatrixVector.cpp#L5417))
*does* implement deflation — `CompactMatrixSVDOrthogonalizeRow` shrinks
`veclist[i+1].Ni` and writes `B->BlockSizeList` (the
`"deflation needed but not implemented"` printf is in the routine ending at
:5409, a different function). So blocks shrink and a chain can run longer than
`NBath` blocks, up to `NBath·BlockSize` of them in the limit of every block
deflating to dimension 1. And the loop breaks on Krylov exhaustion —
`if(veclist[i+1].Ni==0) { B->N=i; break; }` — so an over-large `NTri` is an
allocation bound, not a truncation cap filled with rubbish. Over-allocating costs
`NBath·BlockSize³` doubles, about 4 MB at `NBath = 50`, `BlockSize = 10`; the
`veclist` entries are allocated inside the loop and unused ones cost nothing.

**An empty pole list must not reach the primitive.** With `A->Ni == 0`,
`CompactMatrixBlockBandDiagonalizeReal` sets `B->N = Min(0,NTri) = 0`, allocates
`a` for one block and `b` for none, then runs one iteration and hits the
exhaustion break, which sets `B->N = i = 1`. Every dgemm in that iteration has a
zero dimension so nothing is written, but the result claims one more block than it
owns and `b[0]` does not exist. So the conversion gets an explicit early return
for the empty chain — `val.N = 0`, `InitBlockTridiagonalMatrix`,
`val.a[0] = 0.5·A1` — the idiom `BlockListOfPolesToBlockAndersonMatrix` already
uses for `A.N == 0`.

**Drop exactly-zero-weight poles from each half.** `FermiFunction` returns exact
`1.0` below `mu − 40·kbT` and exact `0.0` above, so after the split a pole outside
the window has an exactly zero weight matrix in one of the two lists. Dropping
those is exact, needs no tolerance, shrinks the working space, and makes "`mu`
outside the bath" produce a genuinely empty list rather than a Lanczos that
deflates to nothing on step one. Test for exact zero only — the tiny-weight poles
near the window edge are kept, as the non-block code keeps them, so the two
implementations stay comparable.

**This path uses full re-orthogonalization.** `DiagonalMatrixBlockTridiagonalize`
passes `ReOrthogonalize = true`, so block `a → d` is Lanczos with full
re-orthogonalization while non-block `a → d` is RKPW. Per
[TODO_ConversionAccuracy.md](TODO_ConversionAccuracy.md) that is the more accurate
and more expensive of the two. **The two will not agree to the last digit on the
same input, and the block one should be the better.** Expected, not a bug.

**Memory, stated plainly.** After merging, the working space is
`Npoles·BlockSize` and the dense diagonal matrix is its square. With a fully
non-degenerate bath `Npoles = NBath·BlockSize`, so that is `BlockSize²` more than
the bath dimension warrants — 200 MB at `NBath = 50`, `BlockSize = 10`. With
degeneracy, which is the physically interesting case, merging collapses it toward
`NBath·BlockSize` and 2 MB. Write it the reuse way and measure; if it bites, the
fix is a rank-revealing square root emitting only `rank(W)` columns per pole,
which is local to one helper and does not change the answer.

**The test the merge asks for.** Build a block Anderson matrix, rotate its bath by
an arbitrary orthogonal `U` — which leaves `G` invariant — convert both to `d`,
and require the two results to agree **element by element** to machine precision,
not just their spectra. Nothing in `TestScripts/Responsefunction/` tests that
today, and it is the property the merge exists to protect.

### Decided: how `d → a` is built in block form

Settled with Maurits. The construction is fixed; only the merge tolerance is still open,
and he is working on that separately.

**One Anderson block per merged pole.** For a merged pole at energy `e_k` with weight
`W_k`:

```
A_{k+1} = e_k · I          (BlockSize × BlockSize)
B_k     = sqrt(W_k)        (full, dense, no columns dropped)
```

which is exact, since `Σ_k sqrt(W_k)·(w − e_k·I)⁻¹·sqrt(W_k)^† = Σ_k W_k/(w − e_k) = Δ(w)`.
Plus `A0 = a0`, `B0 = b0`, `A1 = val.a[0] + con.a[0]`, and the `b0 == 0` dimension-zero
early return the non-block version has.

**Why not the minimal-dimension form.** Maurits: the dimension of the matrix
representing a response function is an efficiency question, not a correctness one, while
breaking symmetry is not acceptable at all. Shifting every bath level of a `t2g` irep by
the same 100 meV is a poor but usable approximation; shifting one of the three by 1 meV
lifts the degeneracy and changes magnetic susceptibilities. So dimension loses.

**`A = e·I` makes the degeneracy structural rather than numerical.** All `BlockSize`
states of a block share *one stored number*, so nothing downstream can split them —
not roundoff, not a later shift, not a chop. The rejected alternative, packing states of
different energies onto one block's diagonal, gives every state its own independent
number and hence its own opportunity to drift. That is the real argument for this form,
beyond dimension bookkeeping.

**Why `sqrt(W)` is canonical although it is computed from an arbitrary eigenbasis.**
`CompactMatrixSqrt`
([BasicMath_CompactMatrixVector.cpp:3453](../../Core/BasicMath/BasicMath_CompactMatrixVector.cpp#L3453))
eigendecomposes `W`, takes the square root of the eigenvalues (clamping negatives to zero
and printing if one is beyond `1000·ε·Σ|λ|`) and reassembles `U·sqrt(Λ)·U^†`. Any
`f(W) = U f(Λ) U^†` depends only on the **spectral projectors**, so the arbitrary basis
inside a degenerate eigenvalue group of `W` cancels exactly.

That is the line, and it is worth stating because a plausible-looking construction
crosses it: **functions of `W` are basis-independent, individual eigenvectors of `W` are
not.** An earlier proposal here was to rank-reveal `W` and keep its `r = rank(W)`
eigenvectors as amplitudes. Maurits' counterexample kills it —
`v = {(1,1,0,0,0),(1,0,1,0,0),(0,1,1,0,0)}` gives
`W = {{2,1,1,0,0},{1,2,1,0,0},{1,1,2,0,0},0,0}` with eigenvalues `4, 1, 1, 0, 0`, and the
two `λ = 1` eigenvectors are an arbitrary orthonormal pair. That reintroduces exactly the
arbitrariness the merge removed. `sqrt(W) = (1/3){{4,1,1},{1,4,1},{1,1,4}}` on that
3-space is symmetric and unique.

**And there is no way to be cleverer.** If `A A^† = A' A'^† = W` with both
`BlockSize × r` of rank `r`, then `A' = A·U` for a unitary `U(r)`. So a factorization
into exactly `rank(W)` vectors is only ever defined up to `U(r)`: **minimal dimension
and canonicality are in genuine conflict**, and any future attempt to compact must fix a
gauge. The canonical objects are `W` and `sqrt(W)`, both `BlockSize × BlockSize`. Do not
try to drop columns either — dropping *exactly zero* columns of `sqrt(W)` is canonical,
but it only compacts when `null(W)` happens to be spanned by basis vectors, and with
crystal field plus spin-orbit there is no pre-determined eigenbasis and the `sqrt(W)` are
dense. Zero couplings are kept.

**Dimension, and where compression actually belongs.** The bath comes out
`N1 · BlockSize` states for `N1` merged poles. Compression is a separate and existing
concern: Maurits compresses by **placing poles on a grid**, which merges them onto shared
energies *before* this conversion, so each grid point becomes one block. So the dimension
is controlled by the grid, not by `d → a`, and the two compose cleanly — gridding is just
another way of producing merged poles with a shared energy, which is the form this
construction wants.

**The pipeline**, all existing primitives except the last step:

```
BlockTriDiagonalMatrixToBlockListOfPoles(&M.val, &Lval)   /* poles of the val hybridization, a0 = val.a[0] */
BlockTriDiagonalMatrixToBlockListOfPoles(&M.con, &Lcon)
AddBlockListOfPoles(&Lval, Lcon)                          /* Lval.a0 becomes A1; poles concatenate */
OrderBlockListOfPoles(&Lval)
BlockListOfPolesMergeDegeneratePoles(&Lval, tol)
/* then one Anderson block per pole, as above */
```

`BlockTriDiagonalMatrixToBlockListOfPoles` builds the chain's minor with
`BlockTriDiagonalMatrixToCompactMatrix` (dimension `N·BlockSize`, `a[0]` excluded as the
convention requires), couples through `b[0]`, and copies `a[0]` into the list's `a0`; its
weights are rank 1. `AddBlockListOfPoles` adds the two `a0` blocks, which is exactly the
§1 rule that the shared site is the sum.

**The caveat that used to stand here is out of date, checked against the source.**
`AddBlockListOfPoles` does not merge — it reallocs, adds the two `a0` blocks and
concatenates, and ends `A->N+=B.N; return 0;`. The hard-coded
`BlockSize · N · DBL_EPSILON · BandWidth` it was accused of lived in
`BlockAndersonMatrixToBlockListOfPoles` and `BlockTriDiagonalMatrixToBlockListOfPoles`,
and both took `MinimalPoleDistance` as an argument in `371eab0d` / `bd2eb8d9`. So the
caller's tolerance now reaches every merge on this path and nothing arrives
pre-merged. The bullet on it in
[TODO_ConversionAccuracy.md](TODO_ConversionAccuracy.md) is stale for the same reason.

`OrderBlockListOfPoles` in the pipeline above is redundant and the implementation
leaves it out: `BlockListOfPolesMergeDegeneratePoles` orders first, both in its
default-epsilon branch and in the real and complex workers it dispatches to.

### What writing `d → a` added to the above

Four things the decision did not cover, all found against the source.

**Where the conjugates go, and the only test that can see them.** A block matrix couples
its head as `b0^* ( ) b0^T` and every level below it as `b ( ) b^+` — the struct comments
in [StructDef.cpp](../../Core/StructDef.cpp) are explicit, and `EvaluateBlockTriDiagonalMatrix`
implements exactly that, `TransformCompactMatrixInlinecT` for `iblock == 1` and
`TransformCompactMatrixInlineNC` below. Both forms share the head, so `a0` and `b0` copy
across untouched. One level down they differ, and the conversion crosses that boundary:
a chain's pole list comes out in the **head** convention — `HermitianCompactMatrixToBlockListOfPoles`
computes `vec^* (w−A)^-1 vec^T` and stores `W = x^* x^T` — while an Anderson bath block
sits **below** and wants `B B^+`. The two read the same Hermitian matrix with the index
order swapped, and the Hermitian square root is what crosses between them:
`CompactMatrixSqrt` returns `B = B^+` with `B B^+ = W`, so `sqrt(W)` is already what the
bath wants. The same `B` in a head slot would give `W^T`.

Note the two pole-list producers disagree on this: `BlockAndersonMatrixInverseToBlockListOfPoles`
stores `v v^+` (below), `HermitianCompactMatrixToBlockListOfPoles` stores `x^* x^T` (head).
Each is right for the function it represents, but they are transposes of one another.

**A real test cannot see any of it.** With the conjugate deliberately put on the wrong
side, every real check in `ChangeTypeDoubleTriBlock.Quanty` still passes to the last
digit, and the complex ones go from `1e-13` to `5.5`. That was measured, not argued, and
it is why the complex section exists. It also gives `a → d` its first complex test, which
crosses the same boundary in the other direction.

**Four `QComplex` flags in, one out.** A `BlockDoubleTriDiagonalMatrixType` carries a flag
per part and `SetBlockDoubleTriDiagonalMatrixQComplex` sets its own to **2** for a mixed
one; a `BlockAndersonMatrixType` has a single flag. So the result is complex when any part
is, and the pole list is promoted with `MakeRealBlockListOfPolesComplex` before the two
halves are added. `H->QComplex = M->QComplex` would be wrong.

**The dimension, measured.** `BlockTriDiagonalMatrixToCompactMatrix` builds the full
`N·BlockSize` minor and ignores `BlockSizeList`, so a chain yields exactly `N·BlockSize`
poles and a deflated slot arrives as a zero-weight pole at energy 0; those merge together
and cost one bath block carrying nothing. For the `BlockSize = 2` test bath, 20 Anderson
blocks in come back as 42 — a factor `BlockSize`, which is the rank-1 weight of a
non-degenerate pole being stored as a full `BlockSize × BlockSize` `sqrt(W)`. That is the
accepted price of §8; `Chop` is what removes the decoupled orbitals, and it matters for §6,
where they become real bath sites in an ED calculation.

**Two defects in `AddBlockListOfPoles`, found here and now fixed.** Its
`else if(B.QComplex)` branch — `A` real, `B` complex, so `A` is made complex first —
realloced `(A->N+B.N+1)·BlockSize²` `Complex16` where the other two branches use
`(A->N+B.N+1)·(1+BlockSize²)−1`. That drops one energy slot per pole, so `B`'s poles were
written past the end of the block: `82 + 82` poles of `BlockSize 2` aborted in malloc,
while short lists survived on what malloc rounded up. `d → a` never reached it — the
promotion above happens first — but `G1 + G2` on a real and a complex block response
function did. Second, the move of `a0` to the end of the grown list is a self copy when
`B.N == 0`, which an empty `con` chain produces, so both are `memmove` now.
`__add.Quanty` has the regression test, real list plus complex list at three lengths;
it aborts with SIGABRT on the old code.

### Done: `d → t` through `a`, and what `SingularValue` does there

Maurits' choice, not a derivation: `d → a → t`. The two forms split `G` differently and
agree only on `G`, so something in between has to hold `G` plainly, and the Anderson form
does — the two chains of a `d` are the two halves of the hybridization and sit in the
denominator, while a chain of the `t` form is the whole of `G` expanded from the observed
block. `BlockAndersonMatrixToBlockTridiagonalMatrix` then takes the poles of `G` and runs
a block Lanczos from that block. This is **not** the interleaved `dsbtrd` band reduction
the scalar `d → t` uses, and the chain is not the scalar version's `1 + val.N + con.N`.

What survives exactly is the head: `a[0]` is carried through untouched, and the first
hopping carries the whole spectral weight, `b[0]^* b[0]^T = B0^* B0^T`. `a[1]` comes back
as `val.a[0] + con.a[0]` only because every `b0` in the tests is a multiple of the
identity, so the observed block is not rotated on the way; with a general `b0` it is that
block whitened by `b0`, and the test says so.

**`SingularValue` decides the dimension on this route, and the standard value leaves it
at its worst.** Two separate places expand a pole into `BlockSize` states coupled through
`sqrt(W)`: `d → a` by construction (§8), and `BlockListOfPolesToBlockTridiagonalMatrix`,
whose `alpha` array carries `BlockSize` copies of each pole energy with the
`BlockSize × BlockSize` `sqrt(W)` as the starting block — the same construction, reached
independently. For a rank-1 weight `BlockSize − 1` of those states are decoupled, and
deflation is the only thing that takes them back out. Measured, `BlockSize = 2`:

| source `a` | `d → a` | `d → a → t`, `1000·ε` | the same, `SingularValue 1e-6` |
|---|---|---|---|
| 7 blocks, 14 states | 13 blocks, 26 | 26 blocks, **54 states** | 7 blocks, **16 states**, same `G` |
| 21 blocks, 42 states | 43 blocks, 86 | 83 blocks, **168 states** | no clean threshold found |

So at the `ChangeType` default the answer is a factor `BlockSize²` wider than the Anderson
matrix it came from, and one loose enough threshold gives the original dimension back at
the same accuracy — `2.48e-13` against `2.45e-13`. But where that threshold sits is a
property of the problem, not a constant: for the six-block bath it is between `1e-8` and
`1e-7`; for the twenty-block one nothing tried deflates cleanly, and `1e-4` makes the
chain *longer* in blocks while losing four digits of `G`. The test asserts only that the
option is honoured and that deflating costs no accuracy, never a number.

That is the thing to attack if this route is to be cheap: deflation on this path should
key on the **rank of `W`**, which is known exactly at the point `sqrt(W)` is formed, and
not on a singular value threshold discovered later by the Lanczos.

### Done: `d → l`, by band-reducing the whole matrix

Maurits' choice again, and deliberately not through `a`: diagonalize the matrix whose
resolvent is `G`, merge the degenerate poles, keep `W`.

The poles of `G` are the eigenvalues of the **whole** matrix, with weights what the
observed block sees of each eigenvector — not the two chain spectra, which are the poles
of the hybridization and are what `d → a` wants. Read the matrix as

```
val_N … val_2, val_1, observed, con_1, con_2 … con_N
```

and it is an ordinary block tri-diagonal matrix of `1 + val.N + con.N` blocks, half
bandwidth `2·BlockSize−1`, which is exactly what `BlockTriDiagonalMatrixEigenSystem`
already band-reduces with `dsbtrd`/`zhbtrd`. The observed block sits in the middle rather
than at an end, which costs nothing: the weights are read from the middle block of the
eigenvectors instead of the first. This is the block twin of the scalar
`DoubleTriDiagonalMatrixToListOfPoles`, which does the same with `dstevd` at `kd = 1`.

**The list stores `W`, not `sqrt(W)`.** `B_i = B'_i B'^C_i` is what the struct comment
defines and what `EvaluateBlockListOfPoles` sums. `sqrt(W)` is what turns a pole back
into a matrix — `d → a`'s business, and `BlockListOfPolesToBlockTridiagonalMatrix`'s.

**Where the conjugates go, and which of them anything can see.** The two chains enter the
denominator through their heads, `b0^* ( ) b0^T`, so the block coupling the observed block
to `val_1` is `conj(val.b[0])` and not `val.b[0]`; eliminating the chain then gives
`conj(b0) ( ) conj(b0)^+ = b0^* ( ) b0^T`. Reversing the valence half exchanges upper and
lower on top of that:

| `T.b[j]` | | |
|---|---|---|
| `j = 1 … val.N−1` | `val.b[val.N−j]^+` | inner val bonds, reversed |
| `j = val.N` | `val.b[0]^T` | `val_1` to the observed block |
| `j = val.N+1` | `conj(con.b[0])` | observed block to `con_1` |
| `j = val.N+1+k` | `con.b[k]` | inner con bonds, as stored |

Measured by injecting each one wrong in turn: the two head conjugates are caught only by
the **complex** tests, where `G` goes from `1e-13` to O(10) while every real check still
passes to the last digit. The reversal `^+` of the inner bonds is caught by **nothing**,
and provably so: the block Lanczos returns Hermitian `b` blocks — verified, `|B − B^+|`
is `1e-14` for both a real and a complex bath, while `|B − B^T|` is `0.57` for the complex
one — so `b^+ = b` for any `d` that came out of a conversion. It is still written `^+`
because a hand-written `d` need not have Hermitian bonds.

**The eigenvector convention.** `BlockTriDiagonalMatrixEigenSystem` returns `fun` with
`fun^* H fun^T` diagonal — the rows are the eigenvectors, the same convention
`CompactMatrixEigensystem` uses. Nothing in the tree consumed it in a way that documented
this, so it was checked against a dense matrix built by hand: off-diagonal `8e-16`,
`|diag − val|` `7e-16`.

**Dimension and cost.** One band reduction of a matrix of `(1 + val.N + con.N)·BlockSize`,
giving that many poles before the merge and no deflation anywhere — 94 poles merging to 82
for the twenty-block test bath, where `G` has 42. The surplus carry no weight and are for
`Chop` and `ReducePoles`. Note `BlockTriDiagonalMatrixEigenSystem` finishes with a full
`n×n×n` back-transform to form every eigenvector component, where only the `BlockSize`
columns of the observed block are used; that is the obvious thing to improve if this ever
becomes hot, and it is in shared code.

`T.N == 1`, both chains empty, asks the band reduction for a half bandwidth wider than the
matrix itself. `dsbtrd` accepts it and gives the right answer, so there is no special case.

### Done: `d → n` through `l`

The natural orbital form **adds** its two chains where the double tri-diagonal form
**nests** them — §1 — so there is no rotation from one to the other; the two split the
same `G` differently. What they agree on is `G`, and the form that holds `G` plainly is
the list of poles. So `d → l → n`, which is what the scalar
`DoubleTriDiagonalMatrixToBADoubleTriDiagonalMatrix` does, and it means `d → n` splits at
`mu` exactly as `l → n` does — both through `CopyBlockListOfPolesVal`/`Con` — so the two
agree by construction rather than by accident.

The poles here are the poles of `G`, not of the hybridization: the observed level is a pole
of `G` and not one of the bath. That is the whole difference from `d → a`.

**Three parameters, none decided here.** `MinimalPoleDistance` merges in the first step,
`kBT` is the temperature of the Fermi split, `SingularValue` is what the two block Lanczos
runs deflate on. That last one needed threading:
`BlockListOfPolesToBlockNaturalImpurityOrbitalMatrix` hard-coded `1000·DBL_EPSILON` for
both of its `BlockListOfPolesToBlockTridiagonalMatrix` calls. It now takes the value; its
three other callers pass `1000·DBL_EPSILON`, so their behaviour is unchanged, and the
`l → n` cell of `ConvertBlockMatrix` forwards the option it already had in hand — which is
also `1000·DBL_EPSILON` when unset, so that changes nothing by default either and only
makes the option live. **Since threaded everywhere.**
`BlockAndersonMatrixToBlockNaturalImpurityOrbitalMatrix`,
`BlockTriDiagonalMatrixToBlockNaturalImpurityMatrix` and
`BlockListOfPolesToBlockNaturalImpurityOrbitalMatrixReduceDimension` had no
`SingularValue` in their signatures and hard-coded `1000·DBL_EPSILON`, so
`ChangeType(…, "Nat")` on an `a` or a `t` source silently ignored an option it had itself
defaulted. They take it now and the `ConvertBlockMatrix` cells forward it.

One more of the same kind turned up on the way and is fixed with them:
**`BlockListOfPolesToBlockAndersonMatrix` took `SingularValue` and ignored it**,
hard-coding `1000·DBL_EPSILON` over the value `ConvertBlockMatrix`'s `l → a` cell was
already handing it. That reached `l → a` and, through
`BlockListOfPolesToBlockDoubleTriDiagonalMatrix`, `l → d`.

All of it is behaviour-preserving at the default — `ChangeType` fills `SingularValue` with
`1000·DBL_EPSILON` when unset, which is what was hard-coded — and that was measured, not
argued: 692 coefficients over twenty conversions between the four block sources come out
bit identical to `f16811b4`. `ChangeType.Quanty` now asserts the option is honoured for all
six cells that take it, and that leaving it out is bit for bit the same as passing
`1000·DBL_EPSILON`.

A loose threshold does not simply shorten the chain: at `1e-4` all of these get **longer** —
24/26 blocks becoming 36/39 — and lose four to five digits of `G`, deflating a direction
that carries weight and leaving the rest of the chain to make up for it. So the tests check
that the option changes the answer and that a tight value keeps `G`, never a length.

**An empty branch is legal for `d` and outside the `n` form** (§1), and the route respects
that on its own. With `mu` below the bath the `d` has `val.N = 0`, but the poles of `G` are
split by a Fermi function and the observed level is a pole of `G` wherever `mu` sits, so
the valence side keeps something — measured `val.N = 1`, never zero. The test asserts both
chains are non-empty for exactly that reason.

The chains come out long, since the split is of the poles of `G` and the `±40 kBT` window
keeps most poles in both halves: 38 and 49 blocks for the twenty-block test bath, from a
list of 82 poles. `SingularValue` is the lever there, as in `d → t`.

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
| ~~1~~ | ~~the `§3.0` head sweep~~ | **done**, see §3.0 |
| ~~2a~~ | ~~`DoubleTriDiagonalMatrixToAndersonMatrix`~~ | **done**, see §4 |
| ~~2b~~ | ~~`DoubleTriDiagonalMatrixToListOfPoles`~~ | **done**, and with it `d → n`. The `d →` row is complete |
| **next** | **the block conversions, §8** | the nine stubs in §0. The non-block insights carry over: a block `d` is block tri-diagonal in the reversed-val ordering and a band matrix of half bandwidth `2·BlockSize` interleaved behind the impurity, and `dsbtrd` is already used on block tri-diagonal matrices next door |
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
