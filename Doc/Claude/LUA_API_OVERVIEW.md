# Quanty Lua API — exported functions (survey of `develop`, entry point `LuaInterface/main.cpp`)

This is a first-pass map of everything the C++ core exposes to Lua scripts. Quanty embeds
Lua 5.x and wires up bindings almost entirely in `LuaInterface/ssp.cpp` (function `luaopen_ssp`),
which is `#include`d — along with every other source file below — directly into
`LuaInterface/main.cpp`. There is a single translation unit; nothing here is a separate
compiled `.o`.

Two binding idioms are used throughout:

- **`X_func` tables** (`static const luaL_Reg X_func[] = {...}`) — registered as a global Lua
  table `X`, i.e. called as `X.Function(...)`. These are "static"/namespace functions.
- **`X_meta` tables** — registered as the metatable `"X_Type"` for a userdata object, i.e.
  called as `obj:Method(...)`. Most `_meta` tables only hold operator overloads
  (`__add`, `__mul`, `__tostring`, `__gc`, ...) plus `__index`/`__newindex`. For several classes
  (notably **Operator**, but also Wavefunction, Spectra, TightBinding, Matrix, ResponseFunction,
  GaussianRadialFunction, InterpolatingFunction), `__index` is a hand-written C function that does
  `strcmp` against method/property names rather than delegating to a lookup table — so those
  method names are *not* visible by grepping for `luaL_Reg` and are listed explicitly below.

---

## 1. Global functions (no namespace prefix)

Source: `global_functions[]` in `LuaInterface/ssp.cpp`, registered via
`luaMWH_RegisterGlobalFunctions(L, global_functions)` in `luaopen_ssp`.

| Lua name | C function |
|---|---|
| BlockBandDiagonalize | LuaBlockBandDiagonalize |
| Braket | LuaBraKet |
| BraketDiagonal | LuaBraKetDiagonal |
| CalculateDOS | LuaCalculateDOS |
| CalculateG | LuaCalculateG |
| CalculateRho | LuaCalculateRho |
| Chop | LuaChop |
| Clone | Luacopy |
| Conjugate | Luaconjugate |
| ConjugateTranspose | Luaconjugatetranspose |
| Copy | Luacopy |
| CreateAtomicIndicesDict | LuaCreateAtomicIndicesDict |
| CreateAtomicIndicesList | LuaCreateAtomicIndicesList |
| CreateSpectra | LuaCreateSpectra |
| CreateFluorescenceYield | LuaCreateFluorescenceYield |
| CreateResonantSpectra | LuaCreateResonantSpectra |
| CreateClusterHamiltonian | LuaCreateClusterHamiltonian |
| DensityMatrix | LuaDensityMatrix |
| DeterminantString | LuaDeterminantString |
| Eigensystem | LuaEigensystem |
| GeneralEigensystem | LuaGeneralEigensystem |
| Expand | LuaWavefunctionExpand |
| ExpandToBasis | LuaWaveFunctionToDeterminantBasis |
| FileReadStuttgartCTRL | LuaFileReadStuttgartCTRL |
| FileReadStuttgartHAMR | LuaFileReadStuttgartHAMR |
| FileReadDresdenFPLO | LuaFileReadDresdenFPLO |
| FindAllAtomsInsideSphere | LuaFindAllAtomsInsideSphere |
| GetSlaterIntegrals | LuaGetSlaterIntegrals |
| GetMultipoleBesselIntegral | LuaGetMultipoleBesselIntegral |
| GetMultipoleIntegral | LuaGetMultipoleIntegral |
| HybridizationFunctionFromFPLO | LuaHybridizationFunctionFromFPLO |
| Inverse | LuaInverse |
| MeanFieldGroundState | LuaMeanFieldGroundstate |
| MeanFieldOperator | LuaMeanFieldOperator |
| NewMatrix | LauMatrixNew |
| NewOperator | LuaNewOperator |
| NewResponseFunction | LuaResponseFunctionNew |
| NewTightBinding | LuaNewTightBinding |
| NewWavefunction | LuaNewWaveFunction |
| NonRelToRelOrbitals | LuaNonRelToRelOrbitals |
| OperatorToMatrix | LuaOperatorToMatrix |
| OperatorToHoles | LuaOperatorToHoles |
| OperatorSetOnsiteEnergy | LuaOperatorSetOnsiteEnergy |
| OperatorSetTrace | LuaOperatorSetOnsiteEnergy |
| OrbitalRotationMatrix | LuaOrbitalRotationMatrix |
| OrbToMultiplicity | LuaOrbToMultiplicity |
| Orthonormalize | LuaOrthonormalize |
| PartialMeanFieldOperator | LuaPartialMeanFieldOperator |
| PartialOperator | LuaPartialOperator |
| PotentialExpandedOnClm | LuaPotentialExpandedOnClm |
| PlotBands | LuaPlotBands |
| PlotBandsAndDOS | LuaPlotBandsAndDOS |
| PlotDOS | LuaPlotDOS |
| PlotFermi | LuaPlotFermi |
| ReadFPLO | LuaReadFPLO |
| PrintExpectationValues | LuaPrintExpectationValues |
| ReadFPLOBasisFunctions | LuaReadFPLOBasisFunctions |
| ReadLenaBasisFunctions | LuaReadLenaBasisFunctions |
| RelToNonRelOrbitals | LuaRelToNonRelOrbitals |
| Rotate | Luarotate |
| Statistics | LuaStatistics |
| TensorProduct | LuaTensorProduct |
| TimeEnd | LuaTimeEnd |
| TimePrint | LuaTimePrint |
| TimeSet | LuaTimeSet |
| TimeStart | LuaTimeStart |
| TightBindingDefFromStuttgartCTRLHAMR | LuaTightBindingDefFromStuttgartCTRLHAMR |
| TightBindingDefFromDresdenFPLO | LuaTightBindingDefFromDresdenFPLO |
| Transpose | Luatranspose |
| Verbosity | LuaVerbosity |
| YtojjzMatrix | LuaYtojjzMatrix |
| YtoKMatrix | LuaYtoKMatrix |
| YtoZMatrix | LuaYtoZMatrix |
| ZtoKMatrix | LuaZtoKMatrix |
| IndexToBitMask | LuaIndexToBitMask |
| BitMaskToIndex | LuaBitMaskToIndex |

A handful of entries are commented out in the source (dead but present as intent): `CalculateBands`,
`CalculateFermi`, `EvolveToGroundstate`, `Help`, `ImpurityGreensFunction`.

### Standard math functions promoted to globals

`luaopen_ssp` also calls `luaMWH_RegisterGlobalFunctions(L, mathlib)`, where `mathlib` is Lua's own
standard math library table (`LuaInterface/lua/lmathlib.c`). This makes `sin`, `cos`, `tan`, `exp`,
`log`, `sqrt`, `abs`, `floor`, `ceil`, `pow`, `random`, etc. callable directly without the `math.`
prefix, in addition to the usual `math.*` table Lua already provides.

### Global constants / data tables

From `luaconstants` (`ssp.cpp`): `I` (imaginary unit, a `Complex` userdata), `Pi`, `pi`, `E`, `e`, `inf`.

From `LuaInterface/units.cpp` (`luaopen_units`): `PhysicalConstants` and `EnergyUnits` — these are
plain Lua tables of values with a custom `__tostring`/`__index`/`__newindex`, not function
namespaces. `EnergyUnits.Set("eV"|"cm-1"|...)` re-derives the whole conversion table.

---

## 2. Classes (metatable-backed userdata types)

### `Complex` — `LuaInterface/Complex.cpp`

- Static (`Complex.X(...)`): `New→LuaComplexNew`, `I→LuaComplexI`, `Exp→LuaComplexExp`,
  `Re/re→LuaComplexRe`, `Im/im→LuaComplexIm`
- Instance operators: `__index→LuaComplexIndex`, `__newindex→LuaComplexNewIndex`, `__add`, `__sub`,
  `__mul`, `__div`, `__pow`, `__unm`, `__eq`, `__tostring`

### `Matrix` — `LuaInterface/Matrix.cpp`

- Static (`Matrix.X(...)`): `New→LauMatrixNew`, `Identity→LuaMatrixIdentity`, `Zero→LuaMatrixZero`,
  `Random→LuaMatrixRandom`, `Diagonal→LuaMatrixDiagonal`, `AddSpin→LuaMatrixAddSpin`,
  `RemoveSpin→LuaRemoveSpin`, `Enlarge→LuaMatrixEnlarge`, `Expand→LuaMatrixExpand`,
  `Flatten→LuaMatrixFlatten`, `Sub→LuaSubmatrix`, `Conjugate→LuaMatrixConjugate`,
  `ConjugateTranspose→LuaMatrixConjugateTranspose`, `Transpose→LuaMatrixTranspose`,
  `Rotate→LuaMatrixRotate`, `Inverse→LuaMatrixInverse`, `Sqrt→LuaMatrixSqrt`, `Exp→LuaMatrixExp`,
  `Eigensystem→LuaMatrixEigensystem`, `Tridiagonalize→LuaMatrixTridiagonalize`,
  `NullSpace→LuaMatrixNullSpace`, `ToOperator→LuaMatrixToOperator`, `ToTable→LuaMatrixToTable`,
  `ToUserdata→LuaMatrixToUserdata`
- Instance operators (`Matrix_meta`): `__unm __add __sub __mul __div __pow __eq __len __tostring
  __index→LuaMatrixIndex __newindex→LuaMatrixNewIndex __gc`
- Properties/rows read via `LuaMatrixIndex` string dispatch (not a table): `dim`/`Dim`, `elem`/`Elem`,
  `Ncol`/`NCol`/`Ncolumn`/`NColumn`, `Nrow`/`NRow`, `Ni`, `Nj`, plus row-object access (`__Row`,
  backed by `Matrix_Row_meta`: `__len __tostring __index __newindex`).
- `MatrixTable_meta` mirrors the arithmetic operators for tables-of-matrices.

### `Tensor` — `LuaInterface/Tensor.cpp`

- Static: `Zero→LuaTensorZero`, `Identity→LuaTensorIdentity`, `ToTable→LuaTensorToTable`,
  `ToUserdata→LuaTensorToUserdata`
- Instance operators (`Tensor_meta`/`TensorTable_meta`): `__unm __add __sub __mul __div __pow __eq
  __len __tostring __index→LuaTensorIndex __gc`

### `Operator` — `LuaInterface/Operator.cpp`

This is the central object of the whole program (many-body operators). Its only `luaL_Reg` table
is the operator-overload metatable; every named method/property is hand-dispatched inside
`LuaOperatorIndex`/`LuaOperatorNewIndex` (`Operator.cpp:5129` / `:5248`):

- Methods: `Print`, `Chop`, `Conjugate`/`C`/`Dag`, `ConjugateTranspose`/`CT`, `Transpose`/`T`
- Read-only properties: `Name`, `NF`, `NB`, `Complex`, `N`, `Restrictions`, `Offset`
- Settable properties: `Name`, `Restrictions` (setting `NF`/`NB`/`Complex`/`N` directly raises a
  Lua error — they must follow from operator content)
- `IndexShift` is a recognized key that explicitly raises "not yet implemented"
- Operator overloads (`Operator_meta`): `__add __sub __mul __div __unm __eq __index __newindex
  __tostring __gc`

`NewOperator` (constructor) lives in the global function table, not here.

### `Wavefunction` — `LuaInterface/Wavefunction.cpp`

- Only `luaL_Reg` table is `Wavefunction_meta` (operators: `__add __sub __mul __div __unm __eq
  __index→LuaWavefunctionIndex __newindex __tostring __gc`)
- Methods/properties via `LuaWavefunctionIndex` string dispatch: `Print`, `Chop`, `Conjugate`/`C`,
  `Normalize`, `Randomize`, `Name`, `N`, `NF`, `NB`, `Complex`, `Det`, `PrefactorDet`, `NBitsKey`
- `NewWavefunction` lives in the global function table.

### `Spectra` — `LuaInterface/Spectra.cpp`

- Static (`Spectra.X(...)`): `Broaden→LuaSpectraBroaden`, `Copy`/`Clone→LuaSpectraCopy`,
  `Element→LuaSpectraElement`, `MultiplyByEnergy→LuaSpectraMultiplyByEnergy`,
  `Shift→LuaSpectraShift`, `Sum→LuaSpectraSum`, `ToTable→LuaSpectraToTable`,
  `FromTable`/`ToUserdata→LuaTableToSpectra`
- Instance operators (`Spectra_meta`): `__add __sub __mul __div __pow __index→LuaSpectraIndex
  __newindex __tostring __gc`
- Methods/properties via `LuaSpectraIndex`: `Print`, `Broaden`, `Element`, `MultiplyByEnergy`,
  `Shift`, `Sum`, `Name`, `N`, `NE`, `dE`, `Emin`, `Emax`, `EminPole`, `EmaxPole`, `Gamma`,
  `Integral`, `ResonantEnergy`
- Constructors (`CreateSpectra`, `CreateFluorescenceYield`, `CreateResonantSpectra`) live in the
  global function table.

### `TightBinding` — `LuaInterface/TightBinding.cpp`

- Static: `New→LuaNewTightBinding`, `Chop→LuaTightBindingChop`, `Copy→LuaTightBindingCopy`,
  `AddSpin→LuaTightBindingAddSpin`, `Derivative→LuaTightBindingDerivative`,
  `ShiftOnsiteEnergy→LuaTightBindingShiftOnsiteEnergy`, `CalculateG→LuaCalculateG`,
  `CalculateMu→LuaCalculateMu`, `CalculateBands→LuaCalculateBands`,
  `SlaterKosterTable→LuaSlaterKosterTable`, `ExportCIF→LuaCrystalStructureToCIF`
- Instance operators (`TightBinding_meta`): `__index→LuaTightBindingIndex __newindex __tostring
  __call __gc`
- Properties via `LuaTightBindingIndex`: `Name`, `NF`, `NAtoms`, `Atoms`, `Cell`, `ReciprocalCell`,
  `Hopping`, `Hk`, `Units`, `Absolute`, `Relative`, `Bohr`, `Angstrom`, `nm`, `NoPi`
- Instance methods via `LuaTightBindingIndex`: `ShiftOnsiteEnergy`. Returned as a C closure over
  the object (`lua_pushvalue`/`lua_pushcclosure`, read back with `lua_upvalueindex(1)`), the same
  idiom `Operator`, `Wavefunction` and `Spectra` use. This is what makes the dot call
  `TB.ShiftOnsiteEnergy(dE)` work without a `self` argument — see the value/reference convention
  in section 5.
- Related file-reader helper types with (currently empty) metatables: `StuttgartCTRL_Type`,
  `StuttgartHAMR_Type`; `DresdenFPLO_Type` has `__gc`/`__tostring` only. Their populating functions
  (`FileReadStuttgartCTRL/HAMR`, `FileReadDresdenFPLO`, `TightBindingDefFrom...`) are global functions.

### `ResponseFunction` — `LuaInterface/ResponseFunction.cpp`

- Static: `New→LuaResponseFunctionNew`, `ChangeType→LuaChangeTypeResponseFunction`,
  `EffectiveHamiltonian→LuaResponseFunctionEffectiveHamiltonian`,
  `FromEffectiveHamiltonian→LuaResponseFunctionFromEffectiveHamiltonian`,
  `ReducePoles→LuaReducePolesResponseFunction`,
  `RemovePolesInInterval→LuaRemovePolesInIntervalResponseFunction`,
  `MergePolesInInterval→LuaMergePolesInIntervalResponseFunction`,
  `Rotate→LuaResponseFunctionRotate`, `PlaceOnMesh→LuaPlaceOnMeshResponseFunction`,
  `ToTable→LuaResponseFunctionToTable`, `ToMatrix→LuaResponseFunctionToCompactMatrix`,
  `ToSpectra→LuaResponseFunctionToSpectra`, `ToUserdata→LuaResponseFunctionToUserdata`,
  `ToTightBinding→LuaResponseFunctionToTightBinding`, `ToOperator→LuaResponseFunctionToOperator`,
  `CalculateSelfEnergy→LuaCalculateSelfEnergy`,
  `CalculateHybridizationFunction→LuaCalculateHybridizationFunction`,
  `CalculateGfromG0andSigma→LuaCalculateGfromG0andSigma`,
  `InvertEnergy→LuaResponseFunctionInvertEnergy`, `Transpose→LuaResponseFunctionTranspose`,
  `Elements→LuaResponseFunctionElement`, `CreateBitMasks→LuaResponseFunctionToBitMask`,
  `Chop→LuaChop`
- Instance operators (`ResponseFunction_meta`/`ResponseFunctionTable_meta`): `__unm __add __sub
  __mul __div __call __eq __tostring __gc __index __newindex __len`
- Properties via index/newindex helpers: get/set `mu`/`Mu`, `type`/`Type`(get only), `name`/`Name`;
  method `Chop` also reachable as `obj:Chop()` (`LuaResponseFunctionIndexChop` /
  `LuaResponseFunctionTableIndexChop`)

### `Graphics` and `BitMap` — `Plot/LuaBinding.cpp`

- `Graphics.X(...)`: `New→LuaNewGraphics`, `Add→LuaAddToGraphicsFunction`, `Plot→LuaPlotGraphics`,
  `ToPS→LuaGraphicsToPS`, `ToSVG→LuaGraphicsToSVG`, `PlotDOS`, `PlotBands`, `PlotBandsAndDOS`,
  `PlotFermi`, `SVGHTMLViewer→LuaSVGHTMLViewer`
- `Graphics_meta`: `__add __len __tostring __index __newindex __gc`
- `BitMap.X(...)`: `New→LuaNewBitMap`, `Overlay→LuaBitMapOverlay`, `FromMatrix→LuaBitMapFromMatrix`,
  `Save→LuaBitMapToFile`, `Pixel→LuaBitMapSetPixel`
- `BitMap_meta`: `__add __len __tostring __gc`

### `GaussianRadialFunction` — `Simon/Gaussian.cpp`

- Static: `NewGaussianRadialFunction→LuaNewGaussianRadial`,
  `ToInterpolationTimesR→LuaGaussianRadialToCubicSplineTimesR`
- Instance operators: `__unm __add __sub __mul __div __pow __concat __eq __lt __le __mode __call
  __tostring __gc __index __newindex`
- Properties via index dispatch: `l`, `N`, `widths`, `amplitudes`

### `InterpolatingFunction` — `LuaInterface/Interpolate.cpp`

- Static: `Spline→LuaSplineInterpolate`, `HermiteCSpline→LuaCubicHermiteSplineInterpolate`,
  `Integrate→LuaInterpolationIntegrate`, `BSpline→LuaBSplineInterpolate`
- Instance operators: `__unm __add __sub __mul __div __pow __concat __eq __lt __le __mode __call
  __tostring __index __gc`
- Properties via index dispatch: `x`/`X`, `y`/`Y`

### `SparseMatrix` — `Simon/SparseMatrix.cpp`

- No dedicated static table; constructed via `Simon.MatrixToSparseMatrix` (see below).
- Instance operators: `__unm __add __sub __mul __div __pow __eq __len __tostring __index __gc`

---

## 3. Author-specific extension namespaces

These are research add-ons by individual group members, each exposed as its own global table with
an (empty) metatable and a `_func` table of free functions — i.e. `Yi.X(...)`, `Martin.X(...)`,
`Simon.X(...)`. All three are wired into the real build (included transitively from
`LuaInterface/main.cpp` via `Simon/Simon.cpp` for the Simon files).

### `Yi` — `Yi/functions_yi.cpp`
`ReadWavefunction`/`readwavefunction → LuaFileReadWavefunction`

### `Martin` — `Martin/functions_martin.cpp`
`NewCoulomb→LuaCreate_Relativistic_Coulomb_Operator`, `NewDirac→LuaOperatorNewDirac`,
`ContinuumBase→LuaCreateBaseForContinuum`, `ReadFPLOBasis→LuaReadFPLOBasis`,
`RB→LuaReadFPLOBasisFunctions`, `Test→LuaMartinTest`, `Metric→LuaMetricWaveFuns`,
`ReadData→LuaReadDataFromFile`, `WriteData→LuaWriteDataToFile`,
`EvaluateSelfEnergy→LuaEvaluateSelfEnergy`,
`PlotSpectrumWithSelfEnergy→LuaPlotSpectrumWithSelfEnergy`,
`CalcSpectrumWithSelfEnergy→LuaCalcSpectrumWithSelfEnergy`,
`CreateEigenDifferentials→LuaCreateEigenDifferentials`, `EigenDifferential`/`ED→LuaEigenDifferential`,
`SparseOperator→LuaSparseOperator`, `SparseTri→LuaSparseTridiagonalize`,
`SparseBlockTri→LuaSparseBlockTridiagonalize`, `OperatorToCSRMatrix→LuaOperatorToCSRMatrix`,
`CSRTridiagonalize→LuaCSRTridiagonalize`, `Export→LuaExport`,
`ReadFPLOEnergies→LuaReadFPLOEnergies`, `RelativisticDipoleOperator→LuaRelativisticDipoleOperator`,
`ResonantSpectrum→LuaResSpectCNSEDInput`, `ListU→LuaListU`, `PlotSpectrum→LuaCDMPlot`,
`BraketList→LuaBraketList`, `AMSel→LuaAugerMeitnerSelfEnergy`, `ReduceOperator→LuaReduceOperator`,
`CreateContinuumSlaterIntegrals→LuaCreateContinuumSlaterIntegrals`,
`CreateAugerMeitnerU→LuaCreateAugerMeitnerU`, `CNSE→LuaCoolNewSelfEnergy`,
`CNSEDiag→LuaCoolNewSelfEnergyDiagonal`, `LinCombine→LuaLinCombine`, `BasisChange→LuaBasisChange`,
`AUcons→LuaOperatorNewAtomicUOccupationConserving`,
`AUconsOs→LuaOperatorNewAtomicUOccupationConservingOneShell`,
`OrthonormalizeContinuum→LuaOrthonormalizeContinuum`, `Thomas→LuaThomasAlgorithm`,
`BSpline→LuaBSpline`, `Breit→LuaBreit`, `Breit4S→LuaBreitFourShells`, `HF→LuaHartreeFock`,
`Slater→LuaSlater`, `HFMatrix→LuaHFMatrix2`, `HFSC→LuaHartreeFockSC`,
`AtomicHartreeFock→LuaAtomicHartreeFock`, `TestHF→LuaTestHF`

Commented out (dead but present): `NewCoulomb→LuaOperatorNewFullCoulomb` (a second, disabled
definition — the live `NewCoulomb` above wins), `Propagator`, `Propagator2`, a second
`ResonantSpectrum→LuaCreateResonantSpectraSE`/`LuaCreateResonantSpectrum`, `EffectiveHamiltonian`.

### `Simon` — `Simon/functions_simon.cpp`
`AtomicNumberOperators→LuaAtomicNumberOperators`, `AtomicToWannierMatrix→LuaAtomicToWannierMatrix`,
`CorthogonalShells→LuaCorthogonalShells`, `ExtendedLigandField→LuaExtendedLigandField`,
`Flatten→LuaFlatten`, `IsosurfaceFixedInt→LuaIsosurfaceFixedInt`,
`IsosurfaceFixedVal→LuaIsosurfaceFixedVal`, `MatrixToOperator→LuaMatrixToOperator`,
`MatrixToSparseMatrix→LuaMatrixToSparseMatrix`, `MatrixTrace→LuaMatrixTrace`,
`OrthogonalizeWithNonOrthoBasis→LuaOrthogonalizeWithNonOrthoBasis`,
`lsCouplingConstants→LuaLSCouplingConstants`, `FileReadOrca→LuaFileReadOrca`,
`PrepareOperatorForEmptyOrbitals→LuaPrepareOperatorForEmptyOrbitals`,
`RealspaceDensity→LuaRealspaceDensity`, `RemoveSpin→LuaRemoveSpin`, `SortByKappa→LuaSortByKappa`,
`SparseMatrixToMatrix→LuaSparseMatrixToMatrix`, `test→LuaSimontest`

(`Simon/Gaussian.cpp` and `Simon/SparseMatrix.cpp` are `#include`d via `Simon/Simon.cpp` too — their
bindings are listed under the `GaussianRadialFunction`/`SparseMatrix` classes above.)

---

## 4. Build-wiring notes

- `CMakeLists.txt` builds a single target from `LuaInterface/main.cpp` only; every file above is
  pulled in transitively through `#include`, so there is exactly one translation unit for the
  whole binary.
- `LuaInterface/main.cpp` has `//#include "Plot.cpp"` and `//#include "PostScript.cpp"` **commented
  out** — `Plot/Plot.cpp` and `LuaInterface/PostScript.cpp` (which also defines a `luaL_Reg` table)
  are dead code, not compiled into the real binary. Plotting instead goes through
  `Plot/LuaBinding.cpp` (`Graphics`/`BitMap`), which *is* included.
- `Plot/TestPlot.cpp` `#include`s the same Yi/Martin/Simon chain as `main.cpp` but is not listed in
  `CMakeLists.txt`'s sources — it looks like an orphaned/alternate entry point, not part of the
  current build.

---

## 5. Suggested next steps

- This survey only covers what's reachable from `#include`d `.cpp` files compiled into
  `LuaInterface/main.cpp`. `Core/`, `IO/`, `Tight_Binding/`, `AtomicOperators/`, `Model/`,
  `Hubbard_Model/`, `DMFT/` were checked and contain no direct Lua registration — they're pure
  C++ helpers consumed by the files above.
- The `Operator`/`Wavefunction`/`Spectra`/`TightBinding`/`ResponseFunction` `__index`-dispatch
  functions are hand-written `strcmp` chains rather than tables — worth keeping in mind for any
  refactor that wants a uniform binding mechanism.
- **Value versus reference semantics — the intended convention.** For an operation that modifies
  an object, the static form takes the object as its first argument, leaves it untouched, and
  returns a modified copy; the method form reached through `__index` modifies the object itself
  and returns nothing:

      TB2 = TightBinding.ShiftOnsiteEnergy(TB, dE)   -- TB unchanged, TB2 is new
      TB.ShiftOnsiteEnergy(dE)                       -- TB itself is shifted

  The C routine does the modifying and takes a pointer; it is the static wrapper that calls
  `CopyTightBindingDef` first, so the copying lives in the Lua layer and never in the C library.
  `TightBinding.ShiftOnsiteEnergy` and `TightBinding.Chop` follow this; not every existing
  binding does, so do not infer the convention from an arbitrary neighbouring function.
- Given the branch is called `develop`, it may be worth diffing this against `master`/tags if the
  goal is understanding recent additions vs. the stable API surface.
