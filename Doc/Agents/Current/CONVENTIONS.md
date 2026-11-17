# Conventions — notation, index order, storage

Conventions that hold across the code and that are **silent when broken**: a real,
symmetric or rank-one test passes either way, and only a complex, non-symmetric
case shows the error. Every test of code that touches one of them needs such a
case. Each convention below has cost at least one debugging cycle.


## 1. Notation

In comments, notes and discussion:

| symbol | meaning |
|---|---|
| `^*` | complex conjugate, element by element, **no** transpose |
| `^T` | transpose, **no** conjugation |
| `^dag` or `^+` | conjugate transpose, the two combined |

Use them exactly; never write `^*` loosely for "the adjoint". They are not
interchangeable: `^dag` leaves a Hermitian residue `w = v v^dag` alone, while `^T`
sends it to `w^T`. A comment that said `^*` in one line and `^T` in the next is
how a missing transpose in `ResponseFunction.InvertEnergy` went unnoticed.

For a matrix-valued function of a complex variable, say what happens to the
argument as well. State the law that always holds first, and a form that holds
only for Hermitian data as a separate, labelled corollary.


## 2. Matrices

- **Rows are eigenvectors.** Wherever a routine returns an eigenvector matrix
  `vec` next to eigenvalues, `vec^* H vec^T` is the diagonal matrix of the
  eigenvalues. For `CompactMatrixEigensystem(A, fun, val)`:
  `A[j][k] = sum_n val[n] fun[n][j] fun[n][k]^*`. `FunTrDotFunCC` and
  `Eigensystem` on a response function follow the same convention.
- `CompactMatrixType` is row-major, LAPACK column-major, so LAPACK sees the
  transpose. `CompactMatrixEigensystem` already conjugates `fun` at the end to undo
  this; do not conjugate again.
- **`CompactMatrixAddMultiply` with `'*'`** (conjugate, no transpose) conjugates
  the operand in place before and after `zgemm` — two extra passes over the
  matrix. On long Lanczos vectors this dominates; rewrite the product to use
  `'C'`, which `zgemm` does for free (`A = Q^dag Y` as `(Y Q^dag)^T`).
- **`CompactMatrixAddMultiply` with `beta = 0` and real `A`, `B`** makes a complex
  `C` real and reallocates it. Never hand it a view on memory it does not own,
  unless `A` or `B` is complex.


## 3. Operators

**Storage order of a ladder string.** In `OperatorFixedLengthType`, term `i` of
length `L` occupies `o[L*i … L*i+L-1]`, and **`o[0]` is the rightmost factor of
the written product**, `o[L-1]` the leftmost. A printed `C 0 A 0` is stored as
`{An(0), Cr(0)}`. Nothing in `Core/StructDef.cpp` says so; it follows from
`RealOperatorAddDiagonalOffset` and from `StringPrintOperatorReadable`, which
prints from `k = length-1` down to `0`.

- The `QDetNotOrdered*` comparators are called as
  `QDetNotOrdered(o[j], o[j-1], NFermion)`: "A is left of B" means A sits at the
  higher index. With `QDetNotOrderedFcFaBcBa` (normal order, creators left) a
  normal-ordered term reads, from high index to low: fermion creators, fermion
  annihilators, boson creators, boson annihilators, the larger index to the left
  within each group.
- Bit 15 (`CrAnMask`) set means creation; `NotCrAnMask` masks out the index. An
  index below `NFermion` is a fermion, from `NFermion` on a boson.
- `NormalOrderOperator(O, &Oo, &QDetNotOrderedFcFaBcBa)` needs `Oo` prepared with
  `InitOperatorCopySettings(&Oo, O)` first; it does not grow `Oo->l[il]`.

**Rotations.** `Rotate(Opp, U)` and `RotateOperatorNormalOrder` take
`U[i][m] = <new_i|Y_m>`, exactly what `YtoZMatrix`, `YtoKMatrix` and
`YtojjzMatrix` return. Pass it as it comes, **no transpose**. The substitution is
`a_m^dag(old) = sum_i U[i][m]^* a_i^dag(new)`, so a one-particle matrix becomes
`A' = U^* A U^T`. Passing `U^dag` rotates back instead, with no warning. The
distinguishing test: `Rotate(Jz, YtojjzMatrix(l))` must come out diagonal.


## 4. Response functions and Green's functions

**Character versus residue.** A density matrix is `rho_ij = <a_i^dag a_j>`, built
from `u_i^* u_j` (`ComplexEigenFunctionToCharacter` in
`Tight_Binding/Tight_Binding_lib.cpp`). A Green's function residue is `u_i u_j^*`
(`G_ij = <a_i 1/(z-H) a_j^dag>`, `ComplexEigenFunctionToResidue`). Using the
character as a residue gives `G^T` — every `CalculateG` route did that for complex
models until commit `e3518bdc`. The unreachable MultiSite/AddSpin routines in
`BasicMath_BlockTypeConversion.cpp` still use the character (no callers).
`TestScripts/TightBindingDef/CalculateGAgainstCreateSpectra.Quanty` pins
`CalculateG`, `CreateSpectra` and the tight-binding Hamiltonian together.

**The head is not like the levels below it.** A block representation couples its
head as `b0^* ( ) b0^T` and every level below as `b ( ) b^dag`.
`EvaluateBlockTriDiagonalMatrix` shows both: `TransformCompactMatrixInlinecT` for
the first block, `TransformCompactMatrixInlineNC` for the rest. `CompactMatrixSqrt`
returns the Hermitian `B` with `B B^dag = W`, the below form.

The transpose of a response function is the same substitution at every level —
`A_i → A_i^T`, `B_i → B_i^*` — which covers both forms:
`(B^* X B^T)^T = B X^T B^dag` and `(B X B^dag)^T = B^* X^T B^T`.

**`BlockTriDiagonalMatrixType`.** `b[k]` for `k >= 1` is the block at position
`(k, k+1)`: `b[k][i][j] = <v_k^i|H|v_{k+1}^j>`, the row being the earlier block.
`b[0]` is different in kind: the seed transform `b[0][i][j] = <v_new_j|v_old_i>`
from `LanczosBlockTriDiagonalize*Krylov*`, which is why `G = b0^* g b0^T`. If a
Gram–Schmidt gives `R_old[i] = sum_j Q_new[j] L[i][j]`, then `b[k+1] = L^*`, not
`L^T`.

**`BlockAndersonMatrixType`** has the same layout: `a[0]` the offset, `a[1]` the
impurity, `a[2…N]` the bath (diagonal), `b[0]` the seed, `b[i]` coupling `a[i+1]`
to the impurity.

**`BlockListOfPoles` stores residues.** It represents
`G(w) = a0 + sum_n W_n / (w - a_n)` with `W_n` the full residue matrix, already
`B B^dag`, not the coupling `B`. Code that needs the coupling takes the matrix
square root itself (`BlockAndersonMatrixFromListOfPolesAsSelfEnergy` calls
`CompactMatrixSqrt`); a contraction `Z^dag W Z` never needs it. Layout: `V` holds
`N` blocks of `1 + BS^2` (`BS` the block size); pole `n` has its energy at
`V[n*(1+BS^2)]` and its residue after it; `a0` points into `V[N*(1+BS^2)]`, set by
`InitBlockListOfPoles`.

**Degenerate poles are merged, with one default in one place.** Converting to a
block list of poles sums the rank-one weights of a degenerate group into
`b P_D b^dag`, `P_D` the spectral projector. That sum is independent of the
arbitrary rotation LAPACK picks inside a degenerate subspace, so merging is not
optional. The threshold is the argument `MinimalPoleDistance` (also a `ChangeType`
option), threaded through `ConvertResponsefunction` and `ConvertBlockMatrix` to
every conversion that merges. **`BlockListOfPolesMergeDegeneratePoles` is the
single place that turns an unset (negative) value into the default**
`10 * DBL_EPSILON * (BandWidth + DBL_EPSILON)`, the band width taken from the
outermost poles. Callers only forward what they were handed; do not add a second
default anywhere. `MinimalPoleDistance == 0` merges **nothing**, not even two
bit-identical poles: the test is a strict `<`.

Why so tight: the eigenvalue error the merge hides is `O(eps ||a_k||)` from
independent `BlockSize × BlockSize` diagonalisations and does not grow with the
number of poles. The old hard-coded `BlockSize * N * DBL_EPSILON * BandWidth` was
70 µeV for a 10 × 10 block bath of 50 blocks on a 60 eV band — large enough to
merge away the µeV splittings of a magnetic field, the physics the block types
exist to resolve.

**`ToMatrix` returns the `H` in the denominator and neglects `A0` and `B0`.** For a
`t` that is `a[1…N]` with its couplings, `b[0]` being left to `GetB0`; for a `d` it
is the impurity with both chains. So the eigenvalues of `ToMatrix(G)` are the
poles of `G`, and a format that stores redundant levels returns more of them — one
for an `n`, two per bath block for a `d` — the extra ones at zero and decoupled.
The block `n` has no `b0` member, so the weight of `G` lives in its head couplings,
and a basis transformation built from them has to be Löwdin-normalised:
`GetNormalizedTFromValandConInBlockNaturalImpurityOrbitalMatrixType`.


## 5. Lua bindings written in C

**Error and warning formats.** `lua_pushfstring`, behind `luaL_error`,
`luaL_warning` and `luaL_argerror`, knows only `%s %c %d %f %p %%`: no widths, no
precisions, no `%g %e %u %ld`. Anything else raises "invalid option" instead of
the message. Write a number as `%f` with a `(double)` cast (printed as `%.14g`)
and an integer as `%d` with an `(int)` cast. Do not add a `vsnprintf` wrapper to
get `printf` formats; the casts keep the code smaller.

**Colours** are `0xAARRGGBB`, `AA` the **transparency**: `00` opaque, `FF`
invisible, so `0xRRGGBB` is opaque. They are stored as `uint32_t` in exactly this
form (helpers `ColorRed`, `ColorGreen`, `ColorBlue`, `ColorTransparency`,
`ColorOpacity`, `ColorVisible`, `ColorCombine`); only the SVG, PostScript and BMP
writers convert to an opacity. Bitmap pixels use the same form, and a pixel table
`{r,g,b,t}` has `t` a transparency. The option names `Transparancy` and
`Fill Transparancy` mean a transparency; `alpha` and `Fill alpha` stay an opacity
(255 opaque), for `PlottingFunctions.Quanty` in the tutorials repository.

Converting a script written before release 0.9: `Graphics.Plot` `LineColor`,
`PointColor` and `Frame Color` were `0xAABBGGRR` opacity — swap `RR` and `BB` and
replace `AA` by `FF-AA`. `PlotBands` was `0xAARRGGBB` opacity — replace `AA` by
`FF-AA`. `Graphics.Add` ignored the top byte — drop it. `Transparancy v` becomes
`255-v`.
