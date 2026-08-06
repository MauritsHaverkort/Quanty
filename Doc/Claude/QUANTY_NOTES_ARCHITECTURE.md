# Quanty — architecture orientation notes

Working notes assembled while getting oriented in the codebase (August 2026).
Intended as a starting point for a session that continues exploring the general
structure of Quanty, not as authoritative documentation. Line numbers refer to
the state of the `Claude` branch at commit `5124437d`.

Everything below was checked against the source; where something is a
convention stated by Maurits rather than enforced by the code, it says so.


## Build

The whole program is one translation unit. `CMakeLists.txt` registers exactly
one source file:

    set(SOURCES LuaInterface/main.cpp)

`LuaInterface/main.cpp` (127 lines) `#include`s everything else — the Core
library, the Lua interpreter, and the Lua binding layer. There are no object
libraries and no per-directory builds. C++11, options `USE_INTEL`, `USE_OPENMP`
(default ON), `USE_MKL`, falling back to `find_package(BLAS/LAPACK)`.

`README.md` gives compile commands for Intel and GNU on Linux. A working
command on macOS, using the Accelerate framework for BLAS/LAPACK:

    g++ LuaInterface/main.cpp -D LUA_USE_MACOSX -DBLASunderscore -O2 \
        -o Quanty -framework Accelerate -lreadline

`-lreadline` is required because `LUA_USE_MACOSX` enables Lua's readline
support. `-DBLASunderscore` selects the underscore-suffixed BLAS/LAPACK symbol
convention, which Accelerate uses. Note that `g++ -fsyntax-only` on this
translation unit reports a spurious `luai_hashnum` error from the bundled Lua;
a real compile succeeds. Build takes roughly 80 s at `-O2`.


## The four core objects

Defined in `Core/StructDef.cpp` (1802 lines). Note that `Core/StructDef.h`
contains only a single function prototype — the struct definitions live in the
`.cpp`, which is consistent with the single-translation-unit build.

| Struct | Line | Role |
|---|---|---|
| `WaveFunctionType` | 172 | Slater-determinant expansion |
| `OperatorType` | 222 | Sum of ladder-operator products |
| `OperatorFixedLengthType` | 207 | One length-`l` block of an operator |
| `SpectraType` | 918 | Discretised spectrum |
| `ResponsefunctionType` | 869 | Green's function dispatcher |
| `DeterminantRestrictionType` | 195 | Occupation restrictions |
| `ConvergenceSettingsType` | 236 | Lanczos / convergence knobs |

### WaveFunctionType

A hash table of determinants. `Determinant` holds bitfields of `NWord` bytes,
blocked into chunks of `2^NBitsBlok`; element `i` lives at
`ValueRe[i >> NBitsBlok][i & BlokMask]`. Lookup goes through
`KeyList[hash(det, NWord, NBitsKey)]`, with `Next` chaining collisions.
`QComplex` decides whether `ValueIm` is allocated at all.

### OperatorType

Terms are grouped by *string length*: `l[0]` is the scalar constant, `l[2]`
holds all one-body terms (C A), `l[4]` all two-body terms (C C A A). Each
`l[length]` is an independent hash table (`o`, `Next`, `KeyList`, `NBitsKey`)
over ladder strings, with real (`VRe`) or complex (`V`) values.

`QComplex` is tri-state: 0 = all real, 1 = all complex, 2 = mixed.

### ResponsefunctionType

Deliberately thin — a `type` tag (`l,a,t,n` / `L,A,T,N`), a `Block` flag, and
`void *data`. The tag plus `Block` select which concrete representation `data`
points at: `TriDiagonalMatrixType`, `AndersonMatrixType`,
`BlockAndersonMatrixType`, `BlockNaturalImpurityOrbitalMatrixType`, and so on,
all defined further down `StructDef.cpp`.


## Hashing

`Core/Hash.cpp` is Bob Jenkins' `lookup3`, adapted: no `initval` argument, and
the result is masked with `hashmask(NBitsKey)` so the table has `2^NBitsKey`
slots. One function serves determinants, operator strings, and sparse-matrix
indices.

An important asymmetry: an **operator's** hash tables are used only while the
operator is being *built* — `RealOperatorAddLader`
(`Core/OperatorArithmetic.cpp:322`) uses them to merge duplicate terms, and
`OperatorRemoveSmallValues` (`:1309`) rebuilds them after a `Chop`. Every
routine that *applies* an operator to a wave function iterates `l[length].N`
linearly over `o[]`/`VRe[]` and never touches `KeyList`. Consequently two
operators holding identical terms may legitimately report different `NBitsKey`
values with no effect on results — `NBitsKey` is inherited from whichever
operand was copied first in a sum.


## Restrictions

`DeterminantRestrictionType` is a set of `N` filters: for each, a bitmask
`Determinant[i]` selecting orbitals, plus `Min[i]`/`Max[i]` allowed particle
counts. `QDetNotRestricted` (`Core/OperatorWaveFunction.cpp:2260`) returns 1
if a determinant passes all filters.

**Invariant worth knowing:** that function uses the *operator's* `NWord` as its
loop bound. `NWord` must be

    Min(NFermion, 1 + (NFermion - 1) / 8) + NBoson

If it is left at 0 the loops run zero times, the counted occupation stays 0,
and every determinant passes — restrictions are silently ignored. This caused
the order-dependent `Eigensystem` bug fixed in commit `5124437d`, where the
block matrix-to-operator converters in
`Core/BasicMathBlock/BasicMath_BlockMatrixToOperator.cpp` set `NFermion` and
`NBoson` but not `NWord`.

Restrictions reach a calculation by two routes, and both matter: the
`{"Restrictions", ...}` option to `Eigensystem`, and the operator's own
`O.Res` field (settable from Lua as `H.Restrictions = ...`). The latter also
governs plain operator-times-wavefunction products, so `H * psi` is
unrestricted unless `H.Restrictions` is set.


## Ground-state solver

`BlockGroundState` (`Core/BlockLanczos.cpp:5782`) alternates: apply `H` to the
current basis (`BlockOperatorPsiSerialRestricted`), grow the basis to the
result, then converge in-basis with `BlockLanczosGroundStateConserveBasis`.
Convergence is judged on the energy variance against `Settings->zero`. Small
bases (`N <= DenseBorder`) go to a dense LAPACK diagonalisation instead.

Two behaviours to be aware of when reading traces:

- The shift scheme at `:5913` calls `OperatorAddLader(H, NULL, 0, ...)`, which
  mutates `H`'s scalar term in place on every outer loop.
- The outer loop aborts with an error after 100 iterations (`:5843`).

`NBitsKey` for *wave functions* is re-tuned each loop
(`BlockWaveFunctionResetKeylist`) to track the growing basis — unlike operator
`NBitsKey`, this one is actively managed.


## Lua binding layer

`LuaInterface/` holds one `.cpp`/`.h` pair per exposed type — `Wavefunction`,
`Operator`, `Spectra`, `ResponseFunction`, `Matrix`, `Tensor`, `BraKet`,
`TightBinding`, `Plot`, `PostScript` — each implementing the metamethods
(`__add`, `__sub`, `__mul`, `__index`, `__newindex`, `__tostring`, `__gc`) that
make the C struct behave as a Lua object. `ssp.cpp` carries the large
free-function surface (`Eigensystem`, `CreateSpectra`, ...) and its option
parsing via `LuaOptionsType`/`LuaCheckOptions`. `eval.cpp` drives evaluation;
`lua/` is the vendored interpreter and `lua_extensions/` the additions
(including `get_put.cpp`, which serialises objects back to Lua source).

`Doc/Claude/LUA_API_OVERVIEW.md` is an existing survey of everything exported
to Lua.

Defaults that surprise: `Chop()` with no argument uses `10*DBL_EPSILON`
(`Operator.cpp:5024`), i.e. essentially a bit-exactness cutoff, applied as an
*absolute* threshold.


## Sub-folder structure

`Core/` is the numerical heart: wave-function and operator arithmetic split by
real/complex/block variants (`RealWaveFunctionArithmetic`,
`ComplexWaveFunctionArithmetic`, `BlockWaveFunctionArithmetic`, and the
mixed-type bridges), Lanczos variants, density matrices, spectra,
`BasicMath/` and `BasicMathBlock/` for matrix and Green's-function algebra.

Physics-model folders each still contain a legacy `main()` from the era when
input files *were* C main functions — `DMFT/`, `Hubbard_Model/`,
`Tight_Binding/`, `Model/`. They are not part of the CMake build.

Per-student folders (`Martin/`, `MartinDresden/`, `Robert/`, `Yi/`, `Simon/`,
`Charles/`-era code, ...) hold experimental work exposed in Lua as
`Name.something`. `README.md` marks these explicitly unsafe and unsupported
outside collaboration with their author.

`Plot/` is the first example of pulling structure definitions out of the
monolith: `Plot/StructDef.cpp` (163 lines) holds `GraphicsType`, `LineType`,
`CurveType`, `ArcType`, `TextType`, `BitMapType` and the plot option structs,
independent of the physics core. More of this splitting is intended.


## Coding conventions

From `README.md` and confirmed directly: Quanty began in C and moved to C++ on
request, but the core is deliberately written in the **common subset of C and
C++** despite the `.cpp` extension. Use arrays with `malloc`/`calloc` rather
than `vector`/`map`/`string`, and `struct` rather than `class`. The stated
reason is optimisability — inefficient C can be rewritten and optimised;
inefficient C++ abstraction cannot, easily. Exceptions are tolerated only in
genuinely separated code that does not touch the core.

Portability targets are deliberately wide: no external dependencies beyond
BLAS and LAPACK (no boost, HDF5, or FFT libraries), so it builds on very old
cluster toolchains.


## Open TODO noted during this work

The C routines are broadly unguarded, and the code does not warn or error when
asked a physically meaningless question. The concrete example that came up: a
starting state containing determinants outside the requested restrictions makes
`BlockGroundState` cycle without converging until it hits the 100-loop limit,
with no diagnostic pointing at the real cause. Maurits intends to add asserts
and guards **coherently across the whole codebase** rather than piecemeal, so
individual guards should not be added ad hoc while fixing unrelated bugs.
