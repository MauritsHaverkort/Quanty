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

**Meanwhile:** set `Ni`, `Nj` and `QComplex` on the destination before calling
these routines, and where a struct holds a `CompactMatrixType` member, do not
assume that member's `QComplex` matches the enclosing struct's.
