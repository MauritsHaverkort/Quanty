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


## Batch 4 — NTriMax on ChangeType

Drafted 18 September 2026, for the commit on the `Claude` branch that adds the
`NTriMax`, `NTriMaxVal` and `NTriMaxCon` options to
`ResponseFunction.ChangeType`.

### A drop-in replacement

| file | page |
|---|---|
| `documentation.language_reference.objects.responsefunction.functions.changetype.txt` | `documentation:language_reference:objects:responsefunction:functions:changetype` |

Built from the page's own raw source. The existing text, the
`{{indexmenu_n>999}}` header and the table-of-contents footer are unchanged; the
diff in the preview should be the `===== Input =====`, `===== Output =====` and
`===== Example =====` sections, which the page did not have, plus one line in
the list at the top.

### Two things in the diff that are not the new options

- **The type list gains //DoubleTri//.** It said there are 4 types and listed
  the four that predate the double tri-diagonal one. `NTriMaxVal` and
  `NTriMaxCon` apply to //Nat// and //DoubleTri// both, so leaving the list at
  four would have made the page contradict itself two sections later.
- **The options that were already there are now written down.** The page showed
  no options at all, so a reader had no way to know that the third argument
  exists or what it looks like. `EffectivekBT`, `SingularValue` and
  `MinimalPoleDistance` are listed alongside the three new ones. `Deflate`,
  `relativeDistance` and `relativeWeigth` are left out: the first belongs to the
  `Dense` conversion and the other two are refused by `ChangeType` as not
  implemented.

### Checked

The example was run against the branch build and the output block is its real
output, character for character. The claim that a chain cut to $N$ keeps the
moments $\mu_0 \ldots \mu_{2N-1}$ was measured rather than quoted: the poles of
the cut and of the full chain were summed into $\mu_k = \sum_i W_i a_i^k$, and
for $N = 2, 3, 4, 5$ the two agree to $10^{-8}$ of $\mu_0 3^k$ up to
$k = 3, 5, 7, 9$ and differ from the next one on.

### Not touched

- `documentation:language_reference:objects:responsefunction:functions:reducepoles`
  is linked to twice but not edited. It is the other way to make a response
  function smaller, and a sentence there pointing back at `NTriMax` would be
  the natural counterpart if you want the link to go both ways.
- The `Dense`, `Userdata` and `Table` targets are named in the new `Input`
  section but not explained, since the page never explained them and that is a
  separate piece of work.


## Batch 5 — Downfold

Drafted 19 September 2026, for the three commits `65a737b9 .. 7ee48b91` on the
`Claude` branch, which add `ResponseFunction.Downfold`, the symmetric improved
estimator of the self-energy on a discrete spectral representation.

### A new page

| file | page |
|---|---|
| `documentation.language_reference.objects.responsefunction.functions.downfold.txt` | `documentation:language_reference:objects:responsefunction:functions:downfold` |

Built on the structure of the `reducepoles` and `calculateselfenergy` pages in
the same namespace, which is where the `{{indexmenu_n>999}}` header, the
`###` blocks and the table-of-contents footer come from.

### One thing to fill in

**The paper is not cited.** The page describes the method and names nothing. I
have the draft (`DiscreteSelfEnergy-15.pdf`, dated 18 September 2026) but not a
DOI or a journal reference, and you are an author, so the citation is yours to
add rather than mine to guess. The natural place is the first sentence of the
`===== The discrete self-energy =====` section, next to "the symmetric improved
estimator of the self-energy". Kugler's original estimator would belong there
too; it is ref. [51] of the draft and I did not want to reproduce a reference I
could not check.

### Checked

The example was run against the branch build and the output block is its real
output. Three claims in it are measured rather than quoted:

- **The table is reproducible.** With `NTri = 14` it is not: the Lanczos block
  deflates from 4 to 3, fourteen steps then leave the augmented propagator short
  of spanning the 56-determinant sector, and the tails of $\Sigma$ move in the
  fourth digit between runs. Twenty two steps is where they stop moving; the
  example and the test both use twenty four. The two columns agree with each
  other in either case, since both are built from the same propagator — which is
  why this was worth finding.
- **The timings for `NTriMax`** (200 blocks / 13 s, 52 / 1.2 s / $10^{-16}$,
  22 / 0.25 s / $2 \cdot 10^{-9}$) are from a chain of 100 blocks of size 20
  downfolded onto 10, run on the same build.
- **$\Sigma^{\mathrm{HF}}$ two ways.** The page says the expectation value
  $\langle \{q,a^\dagger\} \rangle$ and the 12 block of $B_0^* B_0^T$ are the
  same thing; both were computed and agree.

The `Negative value in CompactMatrixSqrt` lines are described on the page rather
than hidden, because every user of this function will see them and their number
varies from run to run. They are the conversions taking the square root of a
pole weight that is numerically zero — see the entry in `../TODO.md` on the
conjugation conventions, which records the same underlying gap.

### Not touched

- `documentation:language_reference:objects:responsefunction:functions:calculateselfenergy`
  is linked to but not edited. It is the other way to get a self-energy and a
  sentence there pointing at `Downfold`, saying when the direct difference is
  the wrong tool, would be the natural counterpart.
- `documentation:language_reference:objects:responsefunction:start` still lists
  four internal representations and there are five. That is the same omission
  the `ChangeType` page had in batch 4 and it is a separate edit.


## Batch 6 — the command line, args and ScriptDirectory

Drafted 20 September 2026, for the three commits on the `Claude` branch that add
`ScriptDirectory`, stop the interactive prompt calling an empty stack when its
input ends, and test all of it.

### A new page

| file | page |
|---|---|
| `documentation.command_line.txt` | `documentation:command_line` |

The page does not exist yet, so create it at that id and paste the file in. It
covers the whole command line rather than only the new variable: the two ways
to start Quanty, the interactive prompt, `args`, and `ScriptDirectory`.

`args` had no documentation at all. A full text search of the wiki for "args"
returns one hit, a Java snippet in the DokuWiki syntax help page, so a script
author has so far had no way of learning that the arguments after the file name
are available, that they arrive as strings, or that the table is `nil` when
there are none.

### Where it sits

`documentation:` is a flat list ordered by `indexmenu_n`: getting_started 1,
install 2, lua 3, memory 4, parallel 5, download 6, basics 7, language_reference
10. The new page belongs directly after Getting started, which is where the
invocation `Quanty [Filename.Quanty]` is already described, but 2 through 7 are
taken and renumbering six pages for one insertion is not worth it. The draft
therefore carries `{{indexmenu_n>2}}` as well and relies on the tie being broken
by name: `command_line` sorts before `install`, so it lands in the right place
with no edit to any other page. Say the word if you would rather renumber
properly, it is a one character edit on each of install, lua, memory, parallel,
download and basics.

The alternative home would be the language reference, which is explicitly
"Constants, Global functions, and Objects". Neither `args` nor `ScriptDirectory`
is any of the three, so that would mean a fourth category and an edit to
`documentation:language_reference:start`. A page of its own next to Getting
started is the smaller change and is where a new user looks.

### One line to add elsewhere

#### `documentation:getting_started`

In the first paragraph, after

```
Once the program is installed one can run Quanty by typing Quanty [Filename.Quanty] in a terminal. [Filename.Quanty] can be any name of a text file that contains the commands you want to execute. 
```

add

```
The page [[documentation:command_line|]] describes the command line in full: an interactive session when no file name is given, the arguments that a script reads back from //args//, and the directory of the script itself in //ScriptDirectory//.
```

This is optional. The new page appears in the navigation menu and in the
`{{indexmenu>.#2|msort nsort}}` index of `documentation:start` on its own, but
the sentence about typing `Quanty [Filename.Quanty]` is the place a reader is
standing when the question comes up.

### Checked

Every snippet on the page was run against the branch build and the output blocks
are real output, with one edit: the directory in the `ScriptDirectory` example
is written as `/home/username/NiO/` rather than the scratch directory the
example actually ran in.

Four statements on the page are measured rather than assumed:

- **Arguments are strings and are not converted.** `Quanty Args.Quanty 3.5 NiO
  -x` gives `args[1] = "3.5"` as a string, `#args == 3`, `args[0] == nil`, and
  the standard Lua `arg` table is not set at all.
- **`args` is `nil`, not an empty table**, when nothing follows the file name.
  Same for `ScriptDirectory` in an interactive session.
- **The prompt.** A bare expression prints nothing; an unfinished statement is
  continued at `>>`; `= 1+1` prints `2`. The session ends on `os.exit()` or on
  end of file.
- **`ScriptDirectory` under `dofile`.** An included file sees the same value as
  the file that included it, since it belongs to the command line and not to the
  chunk being run. `TestScripts/ScriptDirectory.Quanty` pins this, along with
  the absolute form, the trailing separator, and the value being independent of
  the working directory.

### Not touched

- `documentation:language_reference:start` still says the reference has three
  parts. If `args` and `ScriptDirectory` should be findable from there, one
  sentence at the end of that list pointing at `documentation:command_line`
  would do it, without inventing a fourth category.
- `documentation:install:*` describes where to put the binary and says nothing
  about how it is called; the new page is the counterpart and neither needs to
  know about the other.
- The output block on `documentation:getting_started` shows the header of a 2017
  build, with three authors and three citations. It is out of date against every
  current build but that is a separate edit and not mine to decide.


## Batch 7 — CalculateG, Downfold and CreateSpectra on an energy mesh

Drafted 26 September 2026, for the commits on the `Claude` branch that add bath
excitations and `Method "Mesh"` to CreateSpectra (30a1dfdb, e5bba26f), the mesh
route of Downfold (195f2582), `Method "Mesh"` in CalculateG (b36a85aa) and the
growing grid extension `Lambda` (f24a1c90).

### Drop-in replacements

| file | page |
|---|---|
| `documentation.language_reference.functions.calculateg.txt` | `documentation:language_reference:functions:calculateg` |
| `documentation.language_reference.functions.createspectra.txt` | `documentation:language_reference:functions:createspectra` |

Both start from the live source fetched today and keep its wording, markup and
examples; the diff in the preview should be additions only, apart from the two
lines below.

- **CalculateG**: documents the self-energy argument (a response function or a
  table per atom and shell), `EnergyGrid`, `HybridizationFunction` and its second
  result, `Method "Mesh"` and `Lambda`, with a section on how the poles land on
  the grid and a second example (square lattice with an 81 pole self-energy,
  default against mesh, and the hybridization function against zeta - 1/G). The
  first sentence now says "local" Green's function and gives the formula.
- **CreateSpectra**: documents the hybridisation argument and the options
  `Method`, `DenseBorder`, `ForceDenseMethod`, `KrylovBasis`, `NBathExcitations`,
  `NonCrossing`, `NTriFull` and `GridPoints`, the two extra results with a
  hybridisation, and adds a second example (one orbital with an empty band,
  against the same band written out). **It also contains the batch 1 edit**
  (`epsilon` reworded, `SingularValue` added, the two paragraphs after the input
  list), which was never applied to the live page. Paste this file instead of
  doing the batch 1 line edit for createspectra.

### Updated draft, page still to be created

| file | page |
|---|---|
| `documentation.language_reference.objects.responsefunction.functions.downfold.txt` | `documentation:language_reference:objects:responsefunction:functions:downfold` |

The batch 5 draft, which was never pasted (the page does not exist yet), now with
`Method "Mesh"`, `Emin`/`Emax`/`NE`, `EnergyGrid`, `Lambda`, a section "On an
energy grid" and the mesh route added to the example. The batch 5 notes about
where the page sits and the sibling edits still apply.

### Checked

Every new example was run against the branch build (commit 69cd61da) and the
output blocks are that output, with the progress lines left out as the pages
say. The Downfold example reproduces the batch 5 numbers unchanged. Facts on the
pages that are measured rather than read off the code: the dense route refuses a
hybridisation, so a small system needs `Method "Lanczos"` (the CreateSpectra
example needed it); with a hybridisation CreateSpectra returns S, G, S0, G0, the
last two being the system without the baths; the default CalculateG route is the
faster one for a two pole self-energy and twenty times slower for 81 poles.

### Not touched

- The first CreateSpectra example still has "description text" as its
  description and loads `../definitions.Quanty`; both were there before.
- The first CalculateG example (the Emery model on 1000 x 1000 k-points) was not
  run; it is left as it was.


## Batch 8 — CalculateBands and PlotBands with a self-energy

Drafted 26 September 2026, for the commit on the `Claude` branch that adds a
self-energy to `TightBinding.CalculateBands` and spectral functions to
`Graphics.PlotBands`.

### New pages

| file | page |
|---|---|
| `documentation.language_reference.objects.tightbinding.functions.calculatebands.txt` | `documentation:language_reference:objects:tightbinding:functions:calculatebands` |
| `documentation.language_reference.objects.graphics.functions.plotbands.txt` | `documentation:language_reference:objects:graphics:functions:plotbands` |

Neither function had a page. `documentation:language_reference:objects:tightbinding:functions`
holds only its template `start` page, so CalculateBands is its first entry; the
`indexmenu` table of contents of the tight-binding pages picks it up.

**PlotBands has no obvious home.** There is no Graphics section on the site at
all (no object page, no functions namespace). The draft assumes a new namespace
`documentation:language_reference:objects:graphics`, which then also needs a
`start` page and a `functions:start` page (the template pages the other objects
have). Alternatives: put it next to CalculateBands under the tight-binding
functions, or under `documentation:language_reference:functions`. If it moves,
the two links to it (in the CalculateBands page, twice) have to follow.

### Updated draft

| file | page |
|---|---|
| `documentation.language_reference.functions.calculateg.txt` | `documentation:language_reference:functions:calculateg` |

The batch 7 draft gains one paragraph after the input list: in which basis the
self-energy is added to H(k) (Bloch sums whose phase depends only on the unit
cell), and what a self-energy between two atoms therefore means (a coupling
within the same unit cell). This was undocumented; CalculateBands uses the same
basis and refers to it.

### Checked

Both examples were run on the branch build and the output block is that output.
The PlotBands example prints nothing and writes four SVG files; the sizes quoted
on the page were measured. The costs quoted on the CalculateBands page (under a
millisecond per k-point for 50 poles, about 0.2 s for 400) were measured on this
laptop with one orbital.

### Fixed along the way

- The options `"Transparancy"` and `"SecondTransparancy"` of PlotBands. The old
  line `LineColor & 0x00FFFFFF + (T<<24)` parses as
  `LineColor & (0x00FFFFFF + (T<<24))`, which gives the alpha of `Color` AND T
  bit by bit: right whenever `Color` is opaque (the default), wrong when both
  are set. `SecondTransparancy` was read but never used. Both now multiply the
  alpha byte of their color, which is the old result whenever either is 255.
  Both are documented on the PlotBands page.
