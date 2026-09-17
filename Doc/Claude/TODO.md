# Quanty — open TODOs

Standing items on the codebase: things known to be wrong, missing, or deferred,
that are not bugs with a ticket. Kept apart from the architecture notes so those
stay a description of how things *are*.

Each entry records what the problem is, why it has not simply been fixed, and
what to do in the meantime.

> **The order carries no meaning** — entries are added as they come up, not
> ranked. This list is also incomplete: the items below are what surfaced while
> writing the orientation documentation, which is not the same as what matters
> most. Work will likely start elsewhere.


## Coherent guards and asserts across the codebase

**Status:** deferred, deliberately.

The C routines are broadly unguarded and do no type checking — junk in, junk
out, as `README.md` puts it. Quanty does not warn or error when asked a
physically meaningless or self-inconsistent question; it proceeds and produces a
confusing failure somewhere further downstream.

The concrete case that motivated recording this: a starting state containing
determinants outside the requested restrictions makes `BlockGroundState` cycle
without converging — basis size static, energy variance stuck and repeating with
a short period — until it aborts at the 100-outer-loop limit, with no diagnostic
pointing at the real cause. The user-side workaround is to set `H.Restrictions`,
so that products like `H * psi` respect the restrictions too.

**Why it is deferred rather than fixed:** the intent is to add asserts and
guards *coherently across the whole codebase*, as one designed pass with a
consistent policy on what is checked and what it costs. Guards added ad hoc
while fixing unrelated bugs make that harder, not easier, so they should not be
introduced in passing.

**Meanwhile:** when a calculation fails to converge, suspect inconsistent input
before suspecting the numerics.


## Parallel scaling on many-core machines

**Status:** open, and likely to need a rethink rather than a fix.

Quanty is parallelised with **OpenMP on shared memory only** — no MPI, no
distributed memory. The user controls it with `OMP_NUM_THREADS`. There are 379
`#pragma omp` directives in the tree, and they are heavily concentrated in the
operator-times-wavefunction kernel:

| file | pragmas |
|---|---|
| `Core/OperatorWaveFunction.cpp` | 178 |
| `Core/BlockWaveFunctionArithmetic.cpp` | 39 |
| `Core/ComplexWaveFunctionArithmetic.cpp` | 36 |
| `Core/RealWaveFunctionArithmetic.cpp` | 24 |
| `Core/Spectral.cpp` | 12 |

That is loop-level parallelism added where the time was going, which worked
well on the core counts of the time.

**The problem:** using all cores of a 256-core machine does not scale well.
Maurits expects this may require rethinking the parallelisation setup rather
than tuning the existing pragmas.

**Constraint on any solution:** backward compatibility is guaranteed, so a
redesign cannot change Lua-visible behaviour, and it has to be affordable for a
maintainer who also runs a group.

**User documentation: already updated.**
`https://www.quanty.org/documentation/parallel` no longer promises near-linear
scaling unconditionally, and the stale "macOS waits on clang regaining OpenMP"
caveat has been removed. It now says scaling holds for few cores, that "for many
cores this scaling is no longer valid", and that users should test the optimum
on their own machine. That is an honest description of the current state, so
nothing further is needed there until the underlying behaviour changes.


## Migrating the remaining `MatrixType` uses to `CompactMatrixType`

**Status:** in progress indefinitely; new code is already expected to comply.

`MatrixType` is the original dense-matrix struct. It allocates a row-pointer
array plus one block per row, and keeps real and imaginary parts separately, so
a complex `Ni × Nj` matrix costs on the order of `2*(Ni+1)` allocations. This
became the dominant cost — at one point the code spent **99%+ of its time
allocating memory**.

`CompactMatrixType` is the answer: a single flat allocation, complex values
interleaved as `Complex16` in a union, element `(i,j)` at `[i*Nj+j]`. It is what
the Lua `Matrix` object uses.

**Write new code against `CompactMatrixType`.**

**Why it is not finished:** converting every remaining `MatrixType` call site is
expensive, so the old struct is kept and conversion is allowed via
`MatrixToCompactMatrix` (`Core/BasicMath/BasicMath_MatrixTypeConversion.h:60`).
There is deliberately no reverse conversion.

**Meanwhile:** do not add new `MatrixType` interfaces. When touching code that
already uses it, converting is welcome but is not required as part of an
unrelated change.


## Interpolating functions — overhaul pending, C or C++ undecided

**Status:** waiting on the next student; direction deliberately open.

`Core/BasicMath/Interpolation/` is the one genuinely C++ corner of the code: a
class hierarchy in `namespace Interpolation` with virtual dispatch. It was not
planned that way. Maurits wrote the first version in C — the commented-out
`InterpolatingFunctionType` at `Core/StructDef.cpp:1796` is that original — and
a student subsequently rewrote it in C++.

The next student to work on it will overhaul it again, and whether the result
stays C++ or moves back to C is an open question rather than a decided one.

**Meanwhile:** do not treat this directory as the model for new code, and do not
"fix" it in either direction as a side effect of unrelated work. If new code
must interact with it, keep the coupling thin so the eventual rewrite stays
cheap.


## Splitting struct definitions out of the monolith

**Status:** started, mostly not done.

`Core/StructDef.cpp` holds essentially every data structure in the program in
one 1800-line file, including several that have nothing to do with the numerical
core.

`Plot/StructDef.cpp` (163 lines) is the first piece of the intended split:
`GraphicsType`, `LineType`, `CurveType`, `ArcType`, `TextType`, `BitMapType` and
the plot option structs now live with the plotting code, independent of the
physics. More of this separation is intended.


## Tight binding models with `Units = "Relative"`

**Status:** open. One confirmed bug with a known one-token fix, held back
because it sits inside a larger unresolved convention question. Needs
attention; normal use is `Absolute`, which is the default and is unaffected.

`TightBindingDefHopping.ConnectingVector` (`Core/StructDef.cpp:1005`) is the
vector from atom 1 to atom 2 of a hopping, atom 2 taken in its translated cell:

    ConnectingVector = TransVecs·(a,b,c) + Pos2 - Pos1

That definition is fixed by the absolute-units branch of `LuaReadHopping`
(`LuaInterface/TightBinding.cpp:2289` and `:2296`), and `CreateClusterHamiltonian`
depends on it — it matches hoppings by testing `ConnectingVector == Position[j] -
Position[i]` (`LuaInterface/TightBinding.cpp:1704`) against a cluster it has
forced to absolute positions (`:1572`). So `ConnectingVector` is meant to be
**Cartesian, always**, whatever the model's `Units`.

### The confirmed bug

`LuaInterface/TightBinding.cpp:2327`, the relative-units branch, reads

    Cell->Atom[H->Atompos2].Position[j] - Cell->Atom[H->Atompos2].Position[j]

`Atompos2` twice. The difference is identically zero, so the intra-cell atom
offset drops out and `ConnectingVector` collapses to `TransVecs·(a,b,c)`. The
fix is to make the second one `Atompos1`.

Measured on one physical model written both ways — simple cubic `a = 4`, atom A
at the origin, atom B at `(2,0,0)` = `(0.5,0,0)` in units of `a,b,c`, hopping
A→B within the cell and to the cell at `+a`:

| | bond within cell | bond to next cell |
|---|---|---|
| `Absolute` | `2.0, 0, 0` | `6.0, 0, 0` |
| `Relative` | `0.0, 0, 0` | `4.0, 0, 0` |

Consequences for relative-units models: `ChopTightBindingDefDistance`
(`Tight_Binding/Tight_Binding_lib.cpp:203`) sees wrong bond lengths and chops the
wrong terms; `TightBindingDerivative` (`:391`) forms `dH/dk` from
`ConnectingVector` and is wrong; `CreateClusterHamiltonian` silently matches no
hopping between atoms at different positions in the cell and drops those terms
from the cluster Hamiltonian without a diagnostic. Onsite terms, where
`Pos1 == Pos2`, are unaffected.

### Why it has not simply been fixed

Two sibling routines write `ConnectingVector` in **lattice** coordinates in their
relative branches, with no multiplication by `rUnitCell`:
`RemoveCrystalInformationFromTightBinding`
(`Tight_Binding/Tight_Binding_lib.cpp:2915`) and
`AddCrystalInformationFromTightBinding` (`:2990`). That is the opposite
convention to the one above. Repairing line 2327 alone is a strict improvement —
it makes the two unit conventions agree for the same physical model — but it
does not resolve which convention is canonical, and the struct comment at
`Core/StructDef.cpp:1005` ("in units according to Units") reads as endorsing the
lattice-coordinate one. That decision is yours to make before anyone edits the
siblings.

### Unverified, related

`CreateHamiltonianFromTightBindingDef` builds its phase vector as
`Pos1 - Pos2 + TransVecs·(a,b,c)` (`Tight_Binding/Tight_Binding_lib.cpp:587`),
deliberately, per the comment at `:572` giving the phase as `e^{i k (x_1 - x_2)}`.
`ConnectingVector` carries `Pos2 - Pos1 + TransVecs·(a,b,c)` — the atomic part
has the opposite sign. `TightBindingDerivative` forms `dH/dk` from
`ConnectingVector`, so for cells with more than one atom position the two do not
match. This was found by reading, not by testing, and is independent of the
`Units` question above: it would apply to absolute-units models too. It needs
checking before it is either fixed or dismissed. `TestScripts/TightBindingDef/PhaseTests.Quanty`
is the natural place to start.

**Meanwhile:** define tight binding models with `Units = "Absolute"`. Treat any
existing relative-units model as suspect, in particular if it has more than one
atom in the unit cell and has been through `Chop` with a distance,
`TightBinding.Derivative`, or `CreateClusterHamiltonian`.


## Who sizes the output matrix — the preallocation convention for `Get…` routines

**Status:** deferred — needs one codebase-wide decision, not a local fix.

Routines that hand back a `CompactMatrixType` follow the original convention that
the **caller** sets `Ni`, `Nj` and `QComplex` before the call and the routine only
fills in the data. As a consequence many of them choose which union member to read
or write from the *destination's* `QComplex`, not from the source's.

The concrete case that surfaced this is `BlockResponseFunctionGetB0`
([Core/BasicMath/BasicMath_CompactMatrixVector.cpp](../../Core/BasicMath/BasicMath_CompactMatrixVector.cpp))
for response functions of type `d`. It sets `B0->QComplex` from the enclosing
`BlockDoubleTriDiagonalMatrixType` and then `memcpy`s from `b0.V` or `b0.VRe`
according to that. For `d` the two can disagree:
`luaL_checkBlockDoubleTriDiagonalMatrixType` assigns `buff->a0 = *tempM` and
`buff->b0 = *tempM` wholesale from a user-supplied `Matrix_Type`, so
`b0.QComplex` is whatever the user's matrix happened to be, independent of the
`QComplex` of the struct that holds it. A real `b0` inside a complex struct, or
the reverse, then reads the wrong union member — silently.

The `d` branch of `BlockResponseFunctionGetA0` avoids it by using
`CompactMatrixCopy`, which takes `Ni`, `Nj` and `QComplex` from the source and
allocates accordingly. That is very likely the behaviour wanted everywhere.

**Why it is deferred rather than fixed:** adopting it means *relaxing* the
preallocation requirement, and that requirement is relied on both by callers and
by a long list of sibling implementations — the `t`/`a`/`l`/`n` branches of these
same two functions, `BlockResponseFunctionGetB0`, `BlockTriDiagonalMatrixToCompactMatrix`
and its friends, and the `Evaluate…` family. Relaxing it one function at a time
produces a codebase where some routines size their own output and some require it
presized, with nothing at the call site to say which — strictly worse than either
convention applied consistently. It wants one pass with a stated policy, and a
decision on whether `InitCompactMatrix` on an already-initialised matrix is a
no-op, an error, or a reallocation.

There is a smaller alternative worth weighing in that pass: instead of changing
the accessors, make the readers normalise on the way in, so
`luaL_check…DoubleTriDiagonalMatrixType` converts `a0` and `b0` to the struct's
`QComplex` and the invariant "a member matrix has the same QComplex as its
struct" holds by construction. That fixes the observed case without touching any
accessor, but only for the types whose readers are updated.

**Decided since, and no longer part of this item.** `QComplex` on an object whose
sub-parts can each be real or complex is `0` all real, `1` all complex, `2`
mixed — the convention `SetQComplexOperator` already uses for `OperatorType`.
Sub-parts carry their own `0`/`1`. A representation with no mixed form (a
`CompactMatrixType`, a block tri-diagonal) narrows a `2` to `1` on conversion.
The block `DoubleTri` and `NaturalImpurityOrbital` readers and the `a0`/`b0`
accessors follow this; the rest of the codebase has not been audited against it.

That decision does **not** settle this item: it makes the mismatch *legitimate*
rather than accidental, so reading a member matrix by the enclosing struct's flag
is now definitely wrong rather than merely fragile.

**Meanwhile:** set `Ni`, `Nj` and `QComplex` on the destination before calling
these routines, and where a struct holds a `CompactMatrixType` member, read that
member by its own `QComplex` — for a mixed object it will not match the
enclosing struct's.


## The mixed real/complex case for response functions

**Status:** open, and only partly handled. Flagged after the convention was
settled, not before — most of the code predates it.

`QComplex` on an object whose sub-parts can each be real or complex is `0` all
real, `1` all complex, `2` mixed, as `SetQComplexOperator` does for
`OperatorType`. Sub-parts carry their own `0`/`1`. A representation with no mixed
form narrows a `2` to `1` on conversion.

Of the response function types only two are compound and can therefore ever be
`2`:

| type | sub-parts | can be mixed |
|---|---|---|
| `BlockDoubleTriDiagonalMatrixType` | `a0`, `b0`, `val`, `con` | **yes** |
| `BlockNaturalImpurityOrbitalMatrixType` | `a0`, `TMT`, `val`, `con` | **yes** |
| `BlockTriDiagonalMatrixType` | — | no, always 0/1 |
| `BlockAndersonMatrixType` | — | no |
| `BlockListOfPoles` | — | no |
| every non-block type | — | no `QComplex` at all |

So the audit is bounded: everything that touches those two types.

### Already done

The two lua readers set `0`/`1`/`2` by counting. `ResponseFunctionGetA0`,
`BlockResponseFunctionGetA0` and `BlockResponseFunctionGetB0` read each member
matrix by that member's own flag and narrow to `0`/`1` on the way out.
`BlockDoubleTriDiagonalMatrixToCompactMatrix` decides from `val` and `con` only,
promotes whichever chain is real, and leaves the matrix real when only `a0` or
`b0` is complex. `SetTMTFromValandConInBlockNaturalImpurityOrbitalMatrixType`
already narrowed correctly.

### Known to be wrong for a mixed object

- **`InitBlockDoubleTriDiagonalMatrix` and `InitBlockNaturalImpurityOrbitalMatrix`
  push the compound flag down into every sub-part** (`A->val.QComplex=A->QComplex;`
  and so on for `con`, `a0`, `b0`/`TMT`). Handed a `2` they write `2` into
  sub-objects that have no mixed state, and they cannot construct a mixed object
  at all. This is where the convention breaks first, and it is the reason a mixed
  object can currently only arrive from lua.
- **`BlockNaturalImpurityOrbitalMatrixChop`** does `if(A->QComplex)` and then reads
  `A->val.a[i]` as `Complex16`. That is the compound flag selecting a sub-part's
  union member: with `val` real it reads past the end of the allocation.
- **`BlockNaturalImpurityOrbitalMatrixChopReduceDimension`** does
  `if(A->QComplex==true)`. `true` converts to `1`, so `2` compares false: a mixed
  object silently takes the real branch, then reads `A->a0.VRe` when `a0` may be
  complex.

### Not yet checked

Everything else reading `QComplex` of those two types, in particular
`ComplexBlockNaturalImpurityOrbitalMatrixToReal` and its inverse,
`CopyBlockNaturalImpurityOrbitalMatrix` / `CopyBlockDoubleTriDiagonalMatrix`
(which copy the struct wholesale, so they are probably fine),
`ResponsefunctionMul`, `ResponseFunctionInvertEnergy`, `ResponsefunctionChop`,
`BlockResponseFunctionRotate`, `ResponseFunctionEq`, and every `... To ...`
conversion, each of which has to narrow to `1` because no other representation
has a mixed form.

**Why it is not simply fixed:** the two `Init` functions define how every
compound object is built, so changing them means deciding whether `Init` takes a
per-part specification or whether mixed objects are only ever assembled by the
readers and the conversions. That is a small design decision, but it wants making
once rather than per call site. `Core/Lanczos.cpp` and
`Core/BlockWaveFunctionArithmetic.cpp` show the shape to copy: `if(==0) … else
if(==1) … else` for the mixed case, on the operator flag.

**Meanwhile:** a mixed response function can only be created from a lua table.
Read every member matrix by its own `QComplex`, never by the enclosing struct's,
and never compare `QComplex` against `1` or `true` when the object is compound —
use `!=0` for "has any complex part". See the entry above on who sizes the output matrix.


## `ReducePoles{Nval, Ncon}` corrupts memory on a non-block response function

**Status:** open, confirmed by reading. A live crash with a mechanical fix; held
out of the `d` work only because it is unrelated to it.

The four-argument overload `ReduceListOfPolesNumberFermi`
([Core/BasicMath/BasicMath_ListOfPoles.cpp:888](../../Core/BasicMath/BasicMath_ListOfPoles.cpp#L888))
reads:

```c
void ReduceListOfPolesNumberFermi(ListOfPoles * M, uint32_t Nval, uint32_t Ncon, double Ef)
{
  ListOfPoles Mval, Mcon;
  ListOfPolesCopyVal(M, Mval, Ef, DBL_MIN);
  ListOfPolesCopyCon(M, Mcon, Ef, DBL_MIN);
  FreeListOfPoles(M);
  *M = Mval;
  ReduceListOfPolesNumber(&Mval,Nval);
  ReduceListOfPolesNumber(&Mcon,Ncon);
  AddListOfPoles(M, Mcon);
}
```

Three separate faults:

1. **Destination and source are swapped.** `ListOfPolesCopyVal(ListOfPoles *M,
   ListOfPoles B, …)` takes the destination as the pointer and the source by
   value, so this copies the *uninitialised stack variable* `Mval` into `M`.
   `ListOfPolesCopy` ([:1507](../../Core/BasicMath/BasicMath_ListOfPoles.cpp#L1507))
   sets `M->N = B.N` from garbage, calls `InitListOfPoles` on that length —
   leaking the caller's two arrays — and then reads `B.a[i]` / `B.b[i]` through
   garbage pointers. Undefined behaviour, in practice a crash or silent heap
   corruption. `Mval` and `Mcon` are never written by anything.
2. **`*M = Mval` happens before the reduction.** Even with the arguments the
   right way round, `M` would take a snapshot of `Mval` and then
   `ReduceListOfPolesNumber(&Mval, Nval)` would lower `Mval.N` without `M`
   seeing it, so `M->N` would be the *pre*-reduction count over the same
   arrays.
3. **`Mcon` is leaked.** `AddListOfPoles`
   ([:285](../../Core/BasicMath/BasicMath_ListOfPoles.cpp#L285)) appends and does
   not take ownership.

The three-argument overload
([:675](../../Core/BasicMath/BasicMath_ListOfPoles.cpp#L675)) is a separate,
self-contained implementation and is fine.

**It is reachable from Lua.** `G.ReducePoles{Nval, Ncon}` on a **non-block**
response function lands here: `LuaReducePolesResponseFunction`
([LuaInterface/ResponseFunction.cpp:318](../../LuaInterface/ResponseFunction.cpp#L318)).
`G.ReducePoles{N}` takes the sound three-argument overload, and the block case
takes `ReduceBlockListOfPoles`, so only the two-number non-block call is
affected.

The block sibling `ReduceBlockListOfPoles(BlockListOfPoles*, uint32_t, uint32_t)`
([Core/BasicMathBlock/BasicMath_BlockListOfPoles.cpp:1210](../../Core/BasicMathBlock/BasicMath_BlockListOfPoles.cpp#L1210))
is the correct model and is very likely the code this was transcribed from —
note that `CopyBlockListOfPolesVal(*A, &Aval, …)` takes source first and
destination second, the **opposite** order to the non-block
`ListOfPolesCopyVal(&Mval, *M, …)`. That inconsistency is presumably how the
swap happened, and is worth removing at the same time.

**Why it is not simply fixed:** it is not a design question — the fix is to swap
the two argument pairs, move the reduction ahead of `FreeListOfPoles(M)`, and
free `Mcon`. It is recorded rather than done because it was found while working
on the `d` conversions and does not belong in that commit. Two things do want
deciding when someone picks it up: the overload is **not declared in
[BasicMath_ListOfPoles.h](../../Core/BasicMath/BasicMath_ListOfPoles.h)** (only
the three-argument one is) and resolves only through the unity build; and it
splits with `ListOfPolesCopyVal/Con` at `Epsilon = DBL_MIN`, i.e. the linear
ramp degenerated to a step, which the Fermi-function decision below replaces.

**Meanwhile:** call `G.ReducePoles{N}` with a single number, or convert to a
block response function first.

## The dormant single-valued DoubleTriDiagonal helpers

**Status:** open, cheap, and nothing calls them. The `d` type is complete and
tested through `ChangeType`, `ToOperator`, `ToMatrix`, `Eigensystem`, `Rotate`,
`InvertEnergy`, `Transpose`, `CalculateRho`, `EffectiveHamiltonian`,
`FromEffectiveHamiltonian`, `×`, `÷` and `CalculateGfromG0andSigma`. What is
listed here is the remainder of the family that no code path reaches, kept
because the defects are real and will bite whoever wires them up.

- `DoubleTriDiagonalMatrixInverseOld`
  ([BasicMath_DoubleTriDiagonalMatrixVector.cpp:137](../../Core/BasicMath/BasicMath_DoubleTriDiagonalMatrixVector.cpp#L137))
  builds the conduction block from `M.val` — `MI->VRe[0][i] = M.val.b[0]` and
  `M.val.a[i - M.val.N]` where `con` is meant. A copy-paste slip. It also uses
  `MatrixType`, which should be `CompactMatrixType`.
- `DoubleTriDiagonalMatrixAdd` (same file) loops `i < M->val.N - 1`, so
  `val.a[val.N-1]` and `val.b[val.N-1]` are never added, and with `val.N == 0`
  the `uint32_t` underflows into a runaway loop. It also ignores `a0`/`b0`.
- `DoubleTriDiagonalMatrixScalarMultiply` scales every `a` and `b` and then
  *resets* `a0 = 0`, `b0 = 1` — scaling the matrix rather than the Green's
  function, and discarding the offset. Which of the two is intended needs
  deciding before it is used; `ResponsefunctionMul` does not call it and scales
  `a0` and `b0` itself, which is the Green's-function reading.

`AndersonMatrixToDoubleTriDiagonalMatrixNonOrthogonal` is a `/* To DO */` stub
returning 1 with no callers; leave it.

A deflated block `d` can be built only through a conversion that deflates, not
from a hand-written table: the lua reader derives `BlockSize` from the matrices
and takes no per-block size list. So `…ChopReduceDimension` and the `WithDeflation`
matrix routine are exercised, but only along that one route.

## Reserved-zero slots in the block natural impurity orbital type

**Status:** open, latent, in live lua-exposed code.

`BlockNaturalImpurityOrbitalMatrixType` has an outer `CompactMatrixType a0` and
no `b0`; the scaling lives in `val.b[0]` / `con.b[0]`. Its convention is additive,
`G(w) = a0 + G.val(w) + G.con(w)`, so the total constant offset is
`a0 + val.a[0] + con.a[0]`. The head slots are *intended* to be zero, nothing
enforces that except the constructors, and the code does not agree on what a
non-zero one means:

| path | `val.a[0]` non-zero → |
|---|---|
| `EvaluateBlockNaturalImpurityOrbital` | **honoured** — delegates to `EvaluateBlockTriDiagonalMatrix`, which adds its own `a[0]` |
| `…MatrixToSpectra` | **honoured** — `Gval[i] + Gcon[i] + A->a0` |
| `…MatrixToBlockListOfPoles` | **honoured** — and hence `…ToBlockAndersonMatrix` / `…ToBlockTridiagonalMatrix`, which route through it |
| `BlockResponseFunctionGetA0` | **honoured** — returns `a0 + val.a[0] + con.a[0]`, carrying a `//TODO: test if this is correct` |
| `…MatrixChop` | touched — loops all `N+1` a-blocks, block 0 included |
| `…MatrixToCompactMatrix`, `…WithDeflation`, `…ToImpurityHamiltonian` | **dropped** — built from `TMT` and the chain blocks from index 2 on |
| `GetTFromValandCon…`, `SetTMTFromValandCon…` | **dropped** — read only `b[0]` and `a[1]`, `a[2]`, … |
| `ResponsefunctionMul` (block `n`) | **wrong** — scales `a0` by `x` and `val.b[0]`/`con.b[0]` by `√x`, but not `val.a[0]`/`con.a[0]`, so `x·G ≠ Mul(G,x)` |
| `ResponseFunctionInvertEnergy` (block `n`) | skips block 0, i.e. treats it as an offset |

The "dropped" rows are not wrong in themselves — an additive constant on `G` has
no home in a matrix or a Hamiltonian. The problem is that `Evaluate` and
`ToMatrix` can disagree about the same object.

Two places do enforce zero, so an `n` built from lua or from a list of poles is
fine and one built through the C API is not:
`BlockListOfPolesToBlockNaturalImpurityOrbitalMatrix` `memset`s both head blocks
after folding their offset into `B->a0`, and
`luaL_checkBlockNatrualImpurityOrbitalMatrixType` does the full normalisation,
transforming both head blocks with the table's `B0`, adding them into `a0`,
zeroing them and folding `B0` into the two head couplings.

Worth doing, in this order:

1. Fix block-`n` `ResponsefunctionMul` to scale `val.a[0]` / `con.a[0]` by `x` as
   well. It is the one operation that is outright wrong rather than merely
   inconsistent. `__mul.Quanty` has the identity that catches it,
   `x G(w,γ) == (G*x)(w,γ)`, but only exercises it on an `n` whose head slots are
   zero.
2. Resolve the `BlockResponseFunctionGetA0` TODO — it *is* correct, in that it
   agrees with `Evaluate` and `ToSpectra`. Delete the TODO and say so.
3. State the reserved-zero invariant in the struct comment, since two
   constructors already depend on it.

The `d` type deliberately avoids this design: its shared impurity energy is the
sum `val.a[0] + con.a[0]` with every slot live. See the architecture notes.

## The Epsilon split parameter and the three `…PlusT` variants

**Status:** convention settled, sweep unfinished.

`Epsilon` became `kBT` in `DoubleTriDiagonalMatrixType` and
`BlockDoubleTriDiagonalMatrixType` and the convention is a true Fermi function
with the temperature stored on the object. The other five structs still carry
`Epsilon`, and `AndersonMatrixToDoubleTriDiagonalMatrixPlusT` still splits with
the old linear ramp over `[mu-ε, mu+ε]`.

`AndersonMatrixToDoubleTriDiagonalMatrixPlusT`,
`TriDiagonalMatrixToDoubleTriDiagonalMatrixPlusT` and
`TriDiagonalMatrixToSecondDoubleTriDiagonalMatrixPlusT` hold the original
implementation and are kept to be rewritten rather than replaced. Nothing outside
`Test/` calls them and `Test/` is not in the build. They no longer agree with the
routine `ChangeType` uses: they write the whole impurity energy into both chains
instead of half in each, split with the ramp, invent a zero-coupling bath site for
an empty branch, and tri-diagonalize with re-orthogonalized Lanczos, which costs
the O(N²) memory that `PolesToTriDiagonalCoefficients` exists to avoid. A comment
on each says so, so that a caller who binds to one by accident finds out.

Rewriting them on top of `PolesToTriDiagonalCoefficients` means deciding how `T`
is produced: the recurrence gives it in O(NPoles·N) but only to 3e-8
orthogonality, so machine precision still costs one Gram–Schmidt sweep and the
O(NPoles·N) basis. `T` is not exposed by any `ChangeType` on the lua level, so
that is a project of its own.

## Two leftovers from the DoubleTriDiagonal work

**The wiki page for `ResponseFunction`** does not mention the `DoubleTri` type or
that `Nv`/`Nc` in `ToOperator` now apply to it as well as to `Nat`. Both are
user-facing. See [WIKI_EDITING_WORKFLOW.md](WIKI_EDITING_WORKFLOW.md).

**The `Blok…` family is dead code awaiting retirement.**
`BlokTriDiagonalMatrixType` and `BlokDoubleTriDiagonalMatrixType` are the
predecessors of the `Block…` types, spelling and all, and they hold `MatrixType`
where the `Block…` family holds a flat array with a `BlockSize`. Seventeen
functions each live in `BasicMath_BlokTriDiagonalMatrixVector.cpp` and
`BasicMath_BlokDoubleTriDiagonalMatrixVector.cpp`, three more mentions sit in
`BasicMath_MatrixTypeConversion.cpp`, and nothing outside those files calls any
of them. Removing them is a separate decision from anything about `d`; the
codebase convention is to rename rather than delete, so they may want a `…Old`
suffix and a comment instead.

## The self-energy `Method` variants are deliberately partial

**Status:** by design, not a defect — recorded so it is not re-opened.

`ResponseFunctionAddSelfEnergyOnGridViaG0Tri` and `…ViaAndLP`, reached through the
`Method` option of `CalculateGfromG0andSigma`, support a varying subset of the
type and option combinations. These functions exist to compare algorithms rather
than for production use, and completing the matrix would be neither possible nor
useful.

## Accuracy of the response function conversions — its own file

**Status:** open, and a programme rather than a defect. See
[TODO_ConversionAccuracy.md](TODO_ConversionAccuracy.md).

Every conversion is checked against the same probe: the largest relative
difference in `G` over a uniform grid across the band, at `Gamma` of 0.3 and
0.02. So "agrees to 1e-13" means that and nothing more. Small `Gamma`, `w`
between poles, and above all `w` near the chemical potential on a logarithmic
mesh are untested, and that last one is where the physics is.

The separate file records what has been measured, the three things already known
that bear on it — the chain tail not being determined at double precision, full
re-orthogonalization being a measured factor of five to seven better on `G` than
RKPW at a hundredfold cost, and two tolerances that sit exactly where a log mesh
near `mu` would go — and the shape of the `Algorithm` option that should come out
of it, so that a user who needs the last digit can pay for it.

## A real number and a complex number are different lua types

**Status:** open, **very low priority, and no good solution is known.** Recorded
because it costs real time repeatedly, not because it is actionable.

Quanty's complex numbers are userdata with a `Complex_Type` metatable; real
numbers are lua's own `number`. To a user both are just numbers, so people write
a complex number where a real one is expected, or read `a.Re` where `a` is real —
and get an error, or worse, `nil`.

`LuaPushTable_CompactMatrixType`
([TableToArray.cpp](../../LuaInterface/lua_extensions/TableToArray.cpp)) makes
this visible in a way that surprises people. It pushes **element by element**,
and an element whose imaginary part is **exactly** zero is pushed as a plain
number rather than as a `Complex`:

```c
if(M->QComplex)
{
  if(M->V[i*M->Nj+j].Im==0)  lua_pushnumber(L, ...Re);
  else                       /* Complex userdata */
}
```

So one lua table can hold both kinds, and routinely does. An ordinary
block-diagonal `G` at finite `omega` and `Gamma`:

```
G[1][*] : userdata , number
G[2][*] : number   , userdata
G = { { (-0.7797 - 0.0949 I) , 0 } ,
      { 0 , (-1.1458 - 0.2067 I) } }
```

The diagonal elements are `Complex`, the exactly-zero off-diagonal ones are plain
`0`. Nothing has gone wrong: the matrix itself is complex — `Gw.QComplex` is 1 —
and a response function evaluated at finite `Gamma` never returns a real matrix.
It is only the per-element push. But it means `G(w,g)[i][j].Re` is relying on that
element having a nonzero imaginary part, which for a block-diagonal `G` is false
off the diagonal.

### `a.Re` was the wrong spelling to begin with

Maurits: the convention in Quanty is that **a method modifies its object and a
function returns a value**. So `a.Re` reads as "modify `a`", which is not what
anyone means by it, and it only works at all when `a` happens to be a complex
userdata. `Complex.Re(a)` is a function, returns the real part, and **works
whether `a` is a complex userdata or a plain lua number**:

```
  number    Complex.Re -> 2.5   Complex.Im -> 0
  userdata  Complex.Re -> 1     Complex.Im -> -2
```

So on the reading side there is nothing to fix in the implementation — the right
spelling already exists and already handles both types. What is needed is to stop
using `a.Re` and `a.Im`. `Complex` also provides `New`, `I`, `Exp` and the
lowercase aliases `re`, `im`; note `Complex.Abs` and `Complex.Conjugate` do
**not** exist.

The four `DoubleTri` test files were converted as the first instalment:
`ChangeTypeDoubleTri`, `ChangeTypeFromDoubleTri`, `ChangeTypeDoubleTriBlock` and
the `DoubleTri` sections of `Chop`. Converting the rest of `TestScripts/` and the
documentation is a mechanical sweep nobody has done.

**What this does not fix.** Only reading. Going the other way — *writing* a
number into a lua table that C will read — the two formats are still distinct and
the reader still has to be told which it is getting. `luaL_checkarray_complex16`
accepts a plain number or a `{re,im}` pair per element but rejects a table mixing
them, which is a real trap: `{{w, gamma/2}, {0, 0}}` works and
`{{w, gamma/2}, 0}` does not. Nor does it fix the per-element push above; a table
coming back from `Matrix.ToTable` still mixes numbers and userdata, it is just
that `Complex.Re` no longer cares.

**Why the rest is not simply fixed:** it would mean changing the lua number
implementation itself so one type carries both, and that is vendored lua nobody
has wanted to touch. Uniformly pushing `Complex` userdata for every element of a
complex matrix would remove the mixing, but it allocates a userdata per element
and changes what every existing script sees back from `Matrix.ToTable`, so it is
not obviously an improvement either.

**Meanwhile:** use `Complex.Re(x)` and `Complex.Im(x)`, never `x.Re`. In C the
question does not arise: every routine reads the `QComplex` of the operand it is
handed.

## `PyramidalHarmonicsP` has no rotation matrix

`NewOperator("...", NF, IndexUp, IndexDn, {{"Basis","P"}})` is accepted by the
option parser and then fails. Every other basis is now handled by building the
operator in the spherical harmonics basis and rotating with the matrix
`C[i][m] = <new_i|Y_m>` that `YToZMatrix`, `YToKMatrix` and `YToJJzMatrix`
return, but there is no `YToPMatrix`: the pyramidal harmonics exist only as an
option string and a comment saying "alternating dn, up order rotating around z
axis". Nothing in the source says which linear combinations they are.

Until someone writes that matrix, `CreateOperatorInShellBasis` prints
`The PyramidalHarmonicsP basis is not implemented, the operator is returned in
the SphericalHarmonicsY basis` and returns non zero, and
`AngularMomentumOperators` raises a lua error. That is an improvement on what it
did before, which was to leave the operator uninitialized and segfault, but it
is still a basis that is advertised and does not work. Either write
`YToPMatrix` or remove the option.

## `CreateOperatorCFjjz` for a single j shell

The seven argument `CreateOperatorCFjjz`, which takes two j shells, is now built
in the spherical harmonics basis and rotated. The five argument one, reached
from lua as `NewOperator("CF", NF, Jindex, Akm)` with a single relativistic
shell, is still not implemented. It cannot be done the same way: a single j
shell is a *projection* of an l shell, not a rotation of one, and the caller
hands over no indices for the other j shell to put the discarded part in. It
needs its own reduced matrix elements.

## `NewOperator` repeats its basis option parsing thirty one times

Every `LuaOperatorNewXxx` in `LuaInterface/Operator.cpp` carries its own
identical twenty five line block turning the `basis` string into `ibasis`, plus
its own copy of the index length validation. That is around a thousand lines of
duplication, and it is the reason a basis that does not work has to be rejected
in thirty one places instead of one. Factoring it into a single helper — the way
`LuaCoulombShellAngularMomentum` now does for the three coulomb branches — would
shrink the file considerably. Not done here because it touches thirty one
working functions at once.

## The jjz rotation is Hermitian only to round off

`Rotate` accumulates the `(i,j)` and `(j,i)` matrix elements of an operator in
different orders, so a rotated operator can differ from its own adjoint in the
last bit. With `YToZMatrix` and `YToKMatrix`, whose entries are all
`±1/sqrt(2)`, the roundings happen to pair up and the result is exactly
Hermitian; with `YToJJzMatrix`, whose entries are square roots of ratios, it is
not. `Qzz`, `Tz` and the coulomb operator on the jjz basis therefore differ from
their adjoint by about `1E-17`, and `Eigensystem` prints `called with an
operator that is not fully Hermitian, I now will try to repair this`. The
repair is correct and the warning is harmless, but it is noise, and it would
disappear if `RotateOperatorNormalOrder` summed the two orders symmetrically.
Explicitly symmetrizing after the rotation is *not* the fix: `Lplus`, `Lmin` and
the transition operators are not Hermitian to begin with.
