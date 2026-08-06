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


## The central structures

Quanty is organised around a small set of data structures and the functions that
act on them. The style is functional rather than object-oriented — there are no
methods on the C side, only free functions taking struct pointers — but the
structs are still the spine of the program, and the natural way in.

Defined in `Core/StructDef.cpp` (1802 lines). Note that `Core/StructDef.h`
contains only a single function prototype — the struct definitions live in the
`.cpp`, which is consistent with the single-translation-unit build.

**Five structures carry the physics.** Each is exposed to Lua as an object, and
the user-facing behaviour of that object is documented on quanty.org:

| # | Struct | Line | Lua object | Role | User documentation |
|---|---|---|---|---|---|
| 1 | `WaveFunctionType` | 172 | `Wavefunction` | Slater-determinant expansion | [objects/wavefunction](https://www.quanty.org/documentation/language_reference/objects/wavefunction/start) |
| 2 | `OperatorType` | 222 | `Operator` | Sum of ladder-operator products | [objects/operator](https://www.quanty.org/documentation/language_reference/objects/operator/start) |
| 3 | `SpectraType` | 918 | `Spectra` | Discretised spectrum | [objects/spectra](https://www.quanty.org/documentation/language_reference/objects/spectra/start) |
| 4 | `ResponsefunctionType` | 869 | `ResponseFunction` | Green's function dispatcher | [objects/responsefunction](https://www.quanty.org/documentation/language_reference/objects/responsefunction/start) |
| 5 | `TightBindingDef` | 1005 | `TightBinding` | Periodic hopping model | [objects/tightbinding](https://www.quanty.org/documentation/language_reference/objects/tightbinding/start) |

Two naming mismatches to be aware of when moving between the Lua side and the C
side, both verified against the `luaL_checkudata` calls:

- The Lua `Spectra` object wraps **`SpectraTypeArray`** (`:939`), a bundle of
  `SpectraType` columns sharing one energy axis — not a bare `SpectraType`.
  That is why one `CreateSpectra` call can return several polarisations at once.
- `TightBindingDef` is the only one of the five not named `...Type`.

### Supporting structures

These exist to serve the calculations above rather than to represent physics in
their own right. The first three are also exposed to Lua as objects:

| Struct | Line | Lua object | Role | User documentation |
|---|---|---|---|---|
| `Complex16` | 53 | `Complex` | `{Re, Im}` pair; the scalar type throughout | [objects/complex](https://www.quanty.org/documentation/language_reference/objects/complex/start) |
| `CompactMatrixType` | 353 | `Matrix` | Dense matrix as exposed to Lua | [objects/matrix](https://www.quanty.org/documentation/language_reference/objects/matrix/start) |
| `CompactTensorType` | 376 | `Tensor` | Rank-`n` tensor, with per-rank co/contravariance | *(not documented on the site)* |
| `Interpolation::InterpolatingFunction` | — | `InterpolatingFunction` | Splines; see below | [objects/interpolatingfunction](https://www.quanty.org/documentation/language_reference/objects/interpolatingfunction/start) |
| `MatrixType` | 325 | — | Internal dense matrix, `VRe[i][j]` | |
| `CompactVectorType` | 397 | — | Vector, real or complex | |
| `SparseMatrixType` | 415 | — | Sparse matrix | |
| `OperatorFixedLengthType` | 207 | — | One length-`l` block of an `OperatorType` | |
| `DeterminantRestrictionType` | 195 | — | Occupation restrictions | |
| `ConvergenceSettingsType` | 236 | — | Lanczos / convergence knobs | |
| `LuaOptionsType` | 884 | — | Option-table parsing for the bindings | |
| `AtomType`, `CrystalStructureType` | 949, 965 | — | Geometry, used by `TightBindingDef` | |
| `BandstructurePointType`, `BandstructurePathType` | 1016, 1024 | — | k-paths for band structures | |

#### MatrixType versus CompactMatrixType — use CompactMatrixType in new code

Both are live, and the distinction is historical rather than functional.
`MatrixType` came first. It stores elements as `double **VRe, **VIm`: a row
pointer array plus one allocation per row, and real and imaginary parts in two
separate blocks — so a complex `Ni × Nj` matrix costs on the order of
`2*(Ni+1)` allocations. At some point the code was spending **99%+ of its time
allocating memory**. Removing allocations was the fix, and `CompactMatrixType`
is the result:

| | `MatrixType` (325) | `CompactMatrixType` (353) |
|---|---|---|
| storage | `double **VRe, **VIm` | `union { double *VRe; Complex16 *V; }` |
| allocations | one per row, per part | **one**, flat |
| element `(i,j)` | `VRe[i][j]` | `VRe[i*Nj+j]`, or `V[i*Nj+j]` if complex |
| complex layout | two separate blocks | interleaved `Complex16` |
| `Name` | inline `char[MAXStrLength]` | `char *` |
| Lua ownership | — | `OnLuaStack` suppresses the free |

**New code should use `CompactMatrixType`.** A full rewrite of the remaining
`MatrixType` call sites is expensive and has not been done, so the old struct
stays and conversion is allowed: `MatrixToCompactMatrix`
(`Core/BasicMath/BasicMath_MatrixTypeConversion.h:60`). Note there is no reverse
conversion — the flow is one-way, which matches the intended direction of
travel. `MatrixType` is still what several internals pass around; it is the
hopping-matrix member of `TightBindingDefHopping`, for instance.

Both carry `a0`, `b0` and `mu`, used when the matrix represents a continued
fraction. See `TODO.md` for the migration.

`InterpolatingFunction` is the one place among these that is genuine C++: a
class hierarchy in `namespace Interpolation` under
`Core/BasicMath/Interpolation/` (`InterpolatingFunction.h`, `CubicSpline`,
`CubicHermiteSpline`), used through virtual dispatch.

It got there by drift rather than by design. Maurits started it in C — the
commented-out `InterpolatingFunctionType` at `StructDef.cpp:1796` is that
original — and a student then took it over and rewrote it as real C++. It sits
within what `README.md` tolerates, being separated from the core, but it should
not be read as a worked example of how to write new Quanty code. **An overhaul
is expected when the next student picks it up, and whether it ends up in C or in
C++ is genuinely open.** Treat the current shape as provisional; see `TODO.md`.

A further large family of structs — `TriDiagonalMatrixType`,
`AndersonMatrixType`, their `Double`, `Second` and `Block` variants, and
`BlockNaturalImpurityOrbitalMatrixType` — are the concrete payloads behind
`ResponsefunctionType`; see that entry below.

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

### SpectraType

A discretised spectrum. The Lua object is the array form, `SpectraTypeArray`
(`:939`) — several spectra over a shared energy grid — which is what makes
`Spectra` support `+` and scalar `*` for Boltzmann averaging over a list of
eigenstates.

### ResponsefunctionType

Deliberately thin — a `type` tag (`l,a,t,n` / `L,A,T,N`), a `Block` flag, and
`void *data`. The tag plus `Block` select which concrete representation `data`
points at: `TriDiagonalMatrixType`, `AndersonMatrixType`,
`BlockAndersonMatrixType`, `BlockNaturalImpurityOrbitalMatrixType`, and so on,
all defined further down `StructDef.cpp`.

This is the one place where the codebase uses a tagged-union dispatch instead of
a concrete struct, which is why response functions have their own conversion
routines rather than plain arithmetic.

### TightBindingDef

A periodic hopping model: a `CrystalStructureType Cell` for the geometry, and
`nH` hopping entries `TightBindingDefHopping` (`:993`). Each hopping carries the
two atom/shell endpoints, the integer translation between unit cells
(`TransVecs[3]`), the resulting `ConnectingVector[3]` in the cell's units, and a
`MatrixType Hop` holding the hopping matrix itself. `NFermion` and
`FermionIndex` map the model onto spin-orbital indices, which is how a tight-
binding model becomes an `OperatorType`.

This is the entry point for external electronic-structure data: `IO/` reads
Wannier90, Dresden FPLO and Stuttgart NMTO output into this struct
(`TightBindingDefFromDresdenFPLO`, `TightBindingDefFromStuttgartCTRLHAMR`).


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
independent of the physics core. More of this splitting is intended — see
`TODO.md`.


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

Note that the C routines are broadly unguarded and do no type checking. This is
a known state of affairs rather than an oversight; see `TODO.md`.
