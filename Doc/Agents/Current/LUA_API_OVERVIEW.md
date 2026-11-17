# The Lua interface — where bindings live and how they work

How the C code is exposed to Lua: where each kind of binding is registered, the
binding idioms, the method and property names that are hidden in hand-written
`__index` functions, and the calling conventions that trap. For what a function
does, read its page in the language reference on <https://www.quanty.org>. For the
complete list of functions, read the registration tables named below: they are
the source of truth, and this note does not copy them.


## 1. Where a binding is registered

Everything is compiled as one translation unit through `LuaInterface/main.cpp`, and
`luaopen_ssp` in `LuaInterface/ssp.cpp` wires it up.

| What | Registration table | File |
|---|---|---|
| global functions (`Eigensystem`, `CreateSpectra`, `NewOperator`, ...) | `global_functions[]` | `LuaInterface/ssp.cpp` |
| global constants `I`, `Pi`, `pi`, `E`, `e`, `inf` | `luaconstants` | `LuaInterface/ssp.cpp` |
| `PhysicalConstants`, `EnergyUnits` (plain tables; `EnergyUnits.Set("eV"\|"cm-1"\|...)` re-derives the conversion table) | `luaopen_units` | `LuaInterface/units.cpp` |
| `math.*` | `mathlib[]` | `LuaInterface/lua/lmathlib.c` |
| `io.*` | `iolib[]` | `LuaInterface/lua/liolib.c` |
| `Complex.*` | `Complex_func[]` | `LuaInterface/Complex.cpp` |
| `Matrix.*` | `Matrix_func[]` | `LuaInterface/Matrix.cpp` |
| `Tensor.*` | `Tensor_func[]` | `LuaInterface/Tensor.cpp` |
| `Spectra.*` | `Spectra_func[]` | `LuaInterface/Spectra.cpp` |
| `ResponseFunction.*` | `ResponseFunction_func[]` | `LuaInterface/ResponseFunction.cpp` |
| `TightBinding.*` | `TightBinding_func[]` | `LuaInterface/TightBinding.cpp` |
| `InterpolatingFunction.*` | `InterpolatingFunction_func[]` | `LuaInterface/Interpolate.cpp` |
| `Graphics.*`, `BitMap.*` | `Graphics_func[]`, `BitMap_func[]` | `Plot/LuaBinding.cpp` |
| `Yi.*`, `Martin.*`, `Simon.*` | `*_func[]` | `Yi/functions_yi.cpp`, `Martin/functions_martin.cpp`, `Simon/functions_simon.cpp` |

`Operator` and `Wavefunction` have no static table: their constructors
`NewOperator` and `NewWavefunction` are global functions, and everything else is
reached through `__index` (§3). The same holds for the constructors of `Spectra`
(`CreateSpectra`, `CreateResonantSpectra`, `CreateFluorescenceYield`).

**`math` is Quanty's mathematical library, and it is global too.** Besides Lua's
own functions, `lmathlib.c` holds Quanty's: factorials and binomials, `ThreeJ`,
`Wignerd`/`WignerD`, spherical, tesseral and cubic harmonics, `SphericalBesselJ`,
`LaguerreL`, Gauss–Legendre quadrature, `erf`/`erfc`, `FaddeevaW`, a Mersenne
twister, and small vector and matrix helpers. `luaopen_ssp` registers the whole
`mathlib` table as globals as well, so `ThreeJ(...)` and `math.ThreeJ(...)` are the
same function, and most names exist in a lower-case and a capitalised spelling.
Search here before adding a mathematical function anywhere else.

**`io` has three additions:** `io.put` and `io.get` write and read Quanty objects as
Lua source (the serialisation is in `LuaInterface/lua_extensions/get_put.cpp`; what
it cannot save is in `TODO.md`), and `io.mkdir` creates a folder portably.

**The author namespaces** `Yi`, `Martin` and `Simon` are research code by
individual group members. `CONTRIBUTING.md` marks them unsafe and unsupported
outside collaboration with their author. `Simon/Simon.cpp` also includes
`Simon/Gaussian.cpp` (the `GaussianRadialFunction` object) and
`Simon/SparseMatrix.cpp` (the `SparseMatrix` object, built with
`Simon.MatrixToSparseMatrix`).

Commented-out entries in the tables are dead but kept as intent — in
`global_functions[]` for example `CalculateFermi`, `EvolveToGroundstate`, `Help`
and `ImpurityGreensFunction` (which therefore has no wiki page).


## 2. Binding idioms

- **`X_func` tables** (`static const luaL_Reg X_func[]`) are registered as a global
  table `X`, called as `X.Function(...)`.
- **`X_meta` tables** are the metatable `"X_Type"` of a userdata object. Most hold
  only the operator overloads (`__add`, `__mul`, `__tostring`, `__gc`, ...) plus
  `__index` and `__newindex`.
- **Hand-written `__index`.** For `Operator`, `Wavefunction`, `Spectra`,
  `TightBinding`, `Matrix`, `InterpolatingFunction`, `Graphics`, `Complex` and
  others, `__index` is a C function that compares the key with `strcmp`, so those
  method and property names are not visible in any `luaL_Reg` table. They are
  listed in §3. `ResponseFunction` uses small `luaL_Reg` tables per property
  instead (`ResponseFunction_Index`, `ResponseFunction_NewIndex`,
  `ResponseFunction_Method`).
- **A method is a C closure over its object.** `__index` returns the function with
  the object as upvalue (`lua_pushvalue`, `lua_pushcclosure`, read back with
  `lua_upvalueindex(1)`). That is why the dot call `TB.ShiftOnsiteEnergy(dE)` works
  without a `self` argument. `Operator`, `Wavefunction`, `Spectra` and
  `TightBinding` all use this idiom.


## 3. Names hidden behind `__index`

Read from the `__index` and `__newindex` functions; check there before relying on
a name.

| Object | `__index` function | Methods | Properties (read) | Settable |
|---|---|---|---|---|
| `Operator` | `LuaOperatorIndex` | `Print`, `Chop`, `Conjugate`/`C`/`Dag`, `ConjugateTranspose`/`CT`, `Transpose`/`T` | `Name`, `NF`, `NB`, `Complex`, `N`, `Restrictions`, `Offset` | `Name`, `Restrictions`; setting `NF`, `NB`, `Complex` or `N` raises an error, as they follow from the content. `IndexShift` raises "not yet implemented". |
| `Wavefunction` | `LuaWavefunctionIndex` | `Print`, `Chop`, `Conjugate`/`C`, `Normalize`, `Randomize` | `Name`, `N`, `NF`, `NB`, `Complex`, `Det`, `PrefactorDet`, `NBitsKey` | |
| `Spectra` | `LuaSpectraIndex` | `Print`, `Broaden`, `Element`, `MultiplyByEnergy`, `Shift`, `Sum`, `Conjugate`, `Transpose`, `ContractSpectraTensorWithPolarization`, `GetPowderTensor` | `Name`, `N`, `NE`, `dE`, `Emin`, `Emax`, `EminPole`, `EmaxPole`, `Gamma`, `Integral`, `ResonantEnergy` | |
| `TightBinding` | `LuaTightBindingIndex` | `ShiftOnsiteEnergy` | `Name`, `NF`, `NAtoms`, `Atoms`, `Cell`, `ReciprocalCell`, `Hopping`, `Units` | `Name`, `Atoms`, `Cell`, `Hopping`, `Units` (`Absolute`, `Relative`; length units `Bohr`/`a0`, `Angstrom`/`A`, `nm`/`nanometer`; `NoPi`, `2Pi`). `TB(k)` (`__call`, `LuaTightBindingCall`) gives `H(k)`, with options `AsOperator`, `ConstantPhaseInUnitCell`, `AddSpin`. |
| `Matrix` | `LuaMatrixIndex` | | `dim`/`Dim`, `elem`/`Elem`, `Nrow`/`NRow`, `Ncol`/`NCol`/`Ncolumn`/`NColumn`, `Ni`, `Nj`; a number gives a row object (`Matrix_Row_meta`) | |
| `ResponseFunction` | `ResponseFunction_Index` | `Chop` | `mu`/`Mu`, `type`/`Type`, `name`/`Name` | `mu`/`Mu`, `name`/`Name` |
| `InterpolatingFunction` | `LuaInterpolateIndex` | | `x`/`X`, `y`/`Y`, `type`/`Type`; for B-splines `degree`, `knots`, `coefficients` | |
| `Graphics` | `LuaGraphicsIndex` | `Add`, `Rotate` | `NCurve`, `NLine`, `NArc`, `NText`, `NBitMap`, `NLayer`, `OriginX`, `OriginY`, `SizeX`, `SizeY` | `OriginX`, `OriginY`, `SizeX`, `SizeY` |
| `Complex` | `LuaComplexIndex` | | any key starting with `r`/`R` gives the real part, `i`/`I` the imaginary part — but see §5 | |


## 4. Value versus reference semantics — the intended convention

For an operation that modifies an object, the static form takes the object as its
first argument, leaves it untouched, and returns a modified copy; the method form
reached through `__index` modifies the object itself and returns nothing:

    TB2 = TightBinding.ShiftOnsiteEnergy(TB, dE)   -- TB unchanged, TB2 is new
    TB.ShiftOnsiteEnergy(dE)                       -- TB itself is shifted

The C routine does the modifying and takes a pointer; the static wrapper calls
`CopyTightBindingDef` first, so the copying lives in the Lua layer and never in
the C library. `TightBinding.ShiftOnsiteEnergy` and `TightBinding.Chop` follow
this; not every existing binding does, so do not infer the convention from an
arbitrary neighbouring function. New bindings follow it.


## 5. Calling conventions that trap

Collected from use; each one has cost an hour at least once.

- **An options argument is a list of `{"name", value}` pairs**, not a table with
  named keys. `{{"Emin",-4},{"GammaL",0.3}}` is right; `{Emin=-4}` has an empty
  array part, passes `LuaTestOptionsList` and is then silently ignored.
- **Broadening in `ToSpectra` is `GammaL`**, not `Gamma`.
- **Read a number with `Complex.Re(x)`, never `x.Re`.** A method modifies its
  object in Quanty and a function returns one, so `x.Re` is the wrong spelling,
  and it fails outright on a plain Lua number — which is what an element of a
  block matrix becomes when its imaginary part is exactly zero. `Complex.Re` and
  `Complex.Im` take either kind. `Complex.Abs` and `Complex.Conjugate` do not
  exist; `Complex.New`, `Complex.I` and `Complex.Exp` do. The global `Conjugate`
  and `Transpose` work on matrices.
- **Response function tables are 1-based with the offset first.** In a
  `ListOfPoles` table `a[1]` is the offset `a0` and the poles start at `a[2]`,
  while `b[1]` is the first weight. A `Tri` or `And` table has `#a = N+1` and
  `#b = N`, with `a[1]` and `b[1]` being `a0` and `b0`. A `Nat` and a `DoubleTri`
  have one positional element, the `{a0,b0}` pair, their two chains being the
  named `val` and `con` — which is what `#G` counts.
- `Matrix.Sub(M,n)` returns the leading `n × n` **matrix**, not an element, and is
  1-based. `Spectra.ToTable(sp)` gives `sp[ispectrum][ipoint] = {energy, value}`.
- **A chopped complex matrix is still a complex matrix** and does not compare
  equal to a real `Matrix.Zero(n,{{"userdata",false}})`. Compare element by
  element with `Complex.Re`/`Complex.Im` when the operands may differ in kind.
- Running and writing tests has traps of its own (output buffering, paths, globals
  shared between scripts); see `TESTING.md`.


## 6. Files that look like bindings but are not compiled

- `LuaInterface/Plot.cpp` and `LuaInterface/PostScript.cpp` (with its own
  `PostScript_func[]`): their `#include` in `main.cpp` is commented out. Plotting
  goes through `Plot/Plot.cpp` and `Plot/LuaBinding.cpp`, which are compiled.
- `Plot/TestPlot.cpp` includes the same Yi/Martin/Simon chain as `main.cpp` but is
  not in `CMakeLists.txt`: an orphaned alternative entry point.
- `Core/`, `IO/`, `Tight_Binding/`, `Model/`, `Hubbard_Model/` and `DMFT/` contain no
  Lua code; they are C libraries used by the binding files. `DFT/Lua*.cpp` and
  `AtomicOperators/` do hold Lua wrapper functions, but these are registered
  through `global_functions[]` in `ssp.cpp`, like all globals.
