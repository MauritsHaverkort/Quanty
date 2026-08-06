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
