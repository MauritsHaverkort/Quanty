# Wiki drafts — the SingularValue option

Drafted 11 September 2026, for the six commits `d60799fa .. 0c5980ab` on the
`Claude` branch, which give the orthonormalization its own threshold instead of
borrowing `Settings.epsilon`.

Two kinds of change are needed. Two pages are short enough that the whole
source is given here as a drop-in replacement; four pages need one bullet
changed and one added, and for those the exact line-level edit is safer than a
full paste — their sources contain long verbatim output blocks that I can only
read back through a summarising fetch, and I will not risk retyping a table of
numbers I cannot verify character for character.


## Drop-in replacements

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


## Line-level edits

### `documentation:language_reference:functions:createspectra`

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

### `documentation:language_reference:functions:createresonantspectra`

Replace

```
    * "epsilon" Positive real defining the smallest absolute value considered different than zero. (Standard value 1.49E-8)
```

with

```
    * "epsilon" Positive real defining the smallest absolute prefactor of a determinant that is kept in a wave-function. (Standard value 1.49E-10)
    * "SingularValue" Positive real defining the smallest singular value of a list of wave-functions that is considered different from zero. A direction whose singular value is smaller is removed from the basis, which reduces the size of the blocks of the block Lanczos. Set independently of "epsilon"; see [[documentation:language_reference:functions:createspectra|]] for what the value means numerically. (Standard value 1.49E-8)
```

### `documentation:language_reference:functions:createfluorescenceyield`

Replace

```
    * "epsilon" Positive real defining the smallest absolute value considered different than zero. (Standard value 1.49E-8)
```

with

```
    * "epsilon" Positive real defining the smallest absolute prefactor of a determinant that is kept in a wave-function. (Standard value 1.49E-10)
    * "SingularValue" Positive real defining the smallest singular value of a list of wave-functions that is considered different from zero. A direction whose singular value is smaller is removed from the basis, which reduces the size of the blocks of the block Lanczos. Set independently of "epsilon"; see [[documentation:language_reference:functions:createspectra|]] for what the value means numerically. (Standard value 1.49E-8)
```

### `documentation:language_reference:functions:eigensystem`

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


## A correction that is not about SingularValue

All four pages above document the standard `epsilon` as `1.49E-8`. The code has
used `sqrt(DBL_EPSILON)/100 = 1.49E-10` since commit `923c2948` of 28 April
2018, which changed it from `sqrt(DBL_EPSILON)`. The wiki has been a factor of
100 out ever since. The replacements above correct it. Worth a second pair of
eyes, since it is the one number in these edits that is not new work.


## Not touched

- `ImpurityGreensFunction` also takes the new option, but has no wiki page.
  Nothing to edit; worth knowing when a page is eventually written.
- No `===== Index =====` footers are affected: no page is added or renamed, and
  the function namespace uses the automatic `{{indexmenu>.#1}}` and
  `{{indexmenu>.#2}}` listings rather than manual footers.
