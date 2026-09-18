# Wiki drafts

Drop-in DokuWiki source for <https://www.quanty.org>, one file per page, named
after the page id. See `../WIKI_EDITING_WORKFLOW.md` for the division of labour
and the house style.


## Batch 1 — the SingularValue option

Drafted 11 September 2026, for the six commits `d60799fa .. 0c5980ab` on the
`Claude` branch, which give the orthonormalization its own threshold instead of
borrowing `Settings.epsilon`.

Two kinds of change are needed. Two pages are short enough that the whole
source is given here as a drop-in replacement; four pages need one bullet
changed and one added, and for those the exact line-level edit is safer than a
full paste — their sources contain long verbatim output blocks that I can only
read back through a summarising fetch, and I will not risk retyping a table of
numbers I cannot verify character for character.


### Drop-in replacements

Paste the whole file into the page's edit box.

| file | page |
|---|---|
| `documentation.language_reference.functions.blockbanddiagonalize.txt` | `documentation:language_reference:functions:blockbanddiagonalize` |
| `documentation.language_reference.functions.orthonormalize.txt` | `documentation:language_reference:functions:orthonormalize` |

Both keep the pages' own existing markup and wording, including the heading
levels and list indentation they already use, so the diff in the preview should
be only the added material. I left the two typos on the BlockBandDiagonalize
page alone (`Gran-Schmidt`, and the broken sentence "only mixes with the second
/ The basis orbital") so that the diff stays about one thing; say the word and I
will fix them too.


### Line-level edits

#### `documentation:language_reference:functions:createspectra`

Replace this line

```
    * "epsilon" Positive real defining the smallest absolute value considered different than zero. (Default value 1.49E-8)
```

with these two

```
    * "epsilon" Positive real defining the smallest absolute prefactor of a determinant that is kept in a wave-function. (Default value 1.49E-10)
    * "SingularValue" Positive real defining the smallest singular value of a list of wave-functions that is considered different from zero. A direction whose singular value is smaller is removed from the basis, which reduces the size of the blocks of the block Lanczos. Set independently of "epsilon". (Default value 1.49E-8)
```

and then, directly after the last bullet of the `===== Input =====` list and
before `===== Output =====`, insert

```
###
The two thresholds answer different questions. //epsilon// controls the wave-functions: a determinant whose prefactor is smaller than //epsilon// is dropped, which is what keeps the Krylov vectors from filling memory. //SingularValue// controls the blocks: the starting block and every block of the block Lanczos are orthonormalized, and a direction whose singular value is smaller than //SingularValue// is removed, so the blocks can get smaller as the recursion proceeds. Neither follows the other.
###

###
The orthonormalization of a block of wave-functions is built from the overlap matrix $S_{ij} = \langle \psi_i | \psi_j \rangle$, in which a singular value $s$ enters as $s^2$. A singular value below $3 \sqrt{\epsilon_{\mathrm{machine}}} \, s_{\mathrm{max}} \approx 4.5 \cdot 10^{-8} s_{\mathrm{max}}$ can therefore not be distinguished from the rounding error of $S$, and is removed whatever //SingularValue// is set to. Setting //SingularValue// above that limit is meaningful and removes directions that could still have been resolved; setting it below does not buy accuracy.
###
```

#### `documentation:language_reference:functions:createresonantspectra`

Replace

```
    * "epsilon" Positive real defining the smallest absolute value considered different than zero. (Standard value 1.49E-8)
```

with

```
    * "epsilon" Positive real defining the smallest absolute prefactor of a determinant that is kept in a wave-function. (Standard value 1.49E-10)
    * "SingularValue" Positive real defining the smallest singular value of a list of wave-functions that is considered different from zero. A direction whose singular value is smaller is removed from the basis, which reduces the size of the blocks of the block Lanczos. Set independently of "epsilon"; see [[documentation:language_reference:functions:createspectra|]] for what the value means numerically. (Standard value 1.49E-8)
```

#### `documentation:language_reference:functions:createfluorescenceyield`

Replace

```
    * "epsilon" Positive real defining the smallest absolute value considered different than zero. (Standard value 1.49E-8)
```

with

```
    * "epsilon" Positive real defining the smallest absolute prefactor of a determinant that is kept in a wave-function. (Standard value 1.49E-10)
    * "SingularValue" Positive real defining the smallest singular value of a list of wave-functions that is considered different from zero. A direction whose singular value is smaller is removed from the basis, which reduces the size of the blocks of the block Lanczos. Set independently of "epsilon"; see [[documentation:language_reference:functions:createspectra|]] for what the value means numerically. (Standard value 1.49E-8)
```

#### `documentation:language_reference:functions:eigensystem`

This page indents its options with two spaces, not four. Replace

```
  * "Epsilon" smallest coefficient before a determinant that is kept in the basis. (Standard value $1.49012*10^{-8}$) 
```

with

```
  * "Epsilon" smallest coefficient before a determinant that is kept in the basis. (Standard value $1.49012*10^{-10}$) 
  * "SingularValue" smallest singular value of a list of wave-functions that is considered different from zero. A direction whose singular value is smaller is removed from the basis, which reduces the size of the blocks. Set independently of "Epsilon"; see [[documentation:language_reference:functions:createspectra|]] for what the value means numerically. (Standard value $1.49012*10^{-8}$) 
```

Note the trailing space after each bullet on this page — it is in the existing
source, so it is kept.


### A correction that is not about SingularValue

All four pages above document the standard `epsilon` as `1.49E-8`. The code has
used `sqrt(DBL_EPSILON)/100 = 1.49E-10` since commit `923c2948` of 28 April
2018, which changed it from `sqrt(DBL_EPSILON)`. The wiki has been a factor of
100 out ever since. The replacements above correct it. Worth a second pair of
eyes, since it is the one number in these edits that is not new work.


### Not touched

- `ImpurityGreensFunction` also takes the new option, but has no wiki page.
  Nothing to edit; worth knowing when a page is eventually written.
- No `===== Index =====` footers are affected: no page is added or renamed, and
  the function namespace uses the automatic `{{indexmenu>.#1}}` and
  `{{indexmenu>.#2}}` listings rather than manual footers.


## Batch 2 — the Basis option and AngularMomentumOperators

Drafted 17 September 2026, for the two commits `d626daec` and `92fff701` on the
`Claude` branch, which make the `Basis` option of `NewOperator` work for the
tesseral, cubic and jjz bases and add the `AngularMomentumOperators` function.

### A new page

| file | page |
|---|---|
| `documentation.language_reference.functions.angularmomentumoperators.txt` | `documentation:language_reference:functions:angularmomentumoperators` |

The page does not exist yet, so create it at that id and paste the file in. The
function namespace lists its children with `{{indexmenu>.#2}}` on
`documentation:language_reference:functions:start`, sorted automatically, so no
index footer and no `indexmenu_n` value needs setting — the page appears on its
own once it is saved.

### A drop-in replacement

| file | page |
|---|---|
| `documentation.language_reference.functions.newoperator.txt` | `documentation:language_reference:functions:newoperator` |

Built from the page's own raw source, so the diff in the preview should show
only two things: one added sentence in the first paragraph, and a new
`===== Basis =====` section between `===== Output =====` and
`===== Example =====`. Everything else, including the verbatim output block, is
the existing text unchanged.

### Checked

Every Quanty snippet on both pages was run against the branch build, and every
numerical claim was verified rather than asserted: that
`NewOperator("Qzz", ..., {{"Basis","Z"}})` equals
`Rotate(NewOperator("Qzz", ...), YtoZMatrix(2))`, that
`L^2 = l(l+1)(2l+1) GG[1]`, that the Slater operators sum to the Coulomb
operator, and that the `Casimir` and `CoulombJJ` keys advertised in the tables
are the ones actually returned. The output block on the
`AngularMomentumOperators` page is the real output of the script above it.

### Not touched

- `documentation:standard_operators:start` and the thirty-odd pages beneath it
  (`Lz`, `Qzz`, `ldots`, `coulomb_repulsion`, ...) each show a bare
  `NewOperator("Lz", NF, IndexUp, IndexDn)` with no mention of a basis. They are
  still correct, since the default is the spherical harmonics, and repeating the
  basis table on each of them would be worse than one link. If you want a
  pointer there, the natural place is one sentence on
  `documentation:standard_operators:start` rather than on every child.
- The "P" (pyramidal harmonics) value has no rotation matrix anywhere in the
  source. The `NewOperator` draft says in one line that it is not implemented,
  rather than leaving it undocumented.


## Batch 3 — ShellAverage

Drafted 18 September 2026, for the commit on the `Claude` branch that adds the
`ShellAverage` function, the average of an operator over all unitary
transformations within one shell.

### A new page

| file | page |
|---|---|
| `documentation.language_reference.functions.shellaverage.txt` | `documentation:language_reference:functions:shellaverage` |

The page does not exist yet, so create it at that id and paste the file in. As
with `AngularMomentumOperators`, the function namespace lists its children with
`{{indexmenu>.#2}}`, so no index footer and no `indexmenu_n` value is needed.

### Checked

Both snippets were run against the branch build and the output block is the
real output of the first one. Every number claimed on the page was verified
numerically rather than taken from a book: the averages
$F^0-\frac{2}{25}F^2$ (p), $F^0-\frac{2}{63}(F^2+F^4)$ (d),
$F^0-\frac{4}{195}F^2-\frac{2}{143}F^4-\frac{100}{5577}F^6$ (f) and
$F^0-\frac{1}{15}G^1-\frac{3}{70}G^3$ (core to valence) all come out of
`ShellAverage` exactly, and are also what a direct trace over the two particle
sector gives. `TestScripts/ShellAverage.Quanty` contains the same checks.

### Not touched

- `documentation:language_reference:functions:operatorsetonsiteenergy` is the
  closest neighbour and is only linked to, not edited. If you want a pointer
  back, one sentence there saying that `ShellAverage` does the same for the
  two particle part and returns a new operator would be the natural place.
- `documentation:language_reference:functions:operatorsettrace` says only that
  it is deprecated in favour of `OperatorSetOnsiteEnergy`; nothing to add.
