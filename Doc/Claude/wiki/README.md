# Wiki drafts

DokuWiki source for <https://www.quanty.org>, one file per page, named after the
page id (`documentation.getting_started.txt` is the page
`documentation:getting_started`). Every file is the whole page, built on the
live source of 12 November 2026 with the change applied. Paste the whole file
into the page's edit box, and delete the file here once the page is saved.

See `../WIKI_EDITING_WORKFLOW.md` for the house style. The notes of the earlier
batches (what was checked and how) are in the git history of this file, up to
commit `ec5a6a39`.


## Corrections

The live page is wrong or incomplete without these.

| file | change |
|---|---|
| `documentation.standard_operators.atomic_coulomb_repulsion.txt` | The old page was pasted below the new one, so the page now has everything twice. The file is the new page alone. |
| `documentation.language_reference.objects.graphics.functions.plotbands.txt` | The defaults of `Color`, `SecondColor`, `Transparancy` and `SecondTransparancy` follow the old color convention (FF opaque). Now 0xAARRGGBB with 00 opaque, as in the code. |
| `documentation.language_reference.functions.createresonantspectra.txt` | `epsilon` is 1.49E-10 (not 1.49E-8, since 2018) and means the smallest prefactor of a determinant; the option `SingularValue` added. |
| `documentation.language_reference.functions.createfluorescenceyield.txt` | The same. |
| `documentation.language_reference.functions.eigensystem.txt` | The same. |
| `documentation.language_reference.objects.responsefunction.start.txt` | The list of representations gains the fifth, double tri-diagonal (`DoubleTri`), with its formula. Checked numerically: a 2 by 2 complex list of poles converted with `ChangeType(G, "DoubleTri")` and rebuilt from its `ToTable` parts by this formula agrees with G to 1E-14 at two energies; with `B_0 ... B_0` instead of `B_0^* ... B_0^T` it does not. |
| `documentation.language_reference.objects.responsefunction.functions.new.txt` | Five types, not four, and a new section for `DoubleTri` with a single valued and a matrix valued example. The matrix valued example of the natural impurity orbital representation was a placeholder ("some example code") and is filled in. All three examples were run and checked against the formula by hand; the order of the fields in the printed tables changes from run to run. |


## Links to new pages

New pages that are otherwise only reachable from the menu.

| file | change |
|---|---|
| `documentation.getting_started.txt` | One sentence pointing to the command line page, after "Quanty [Filename.Quanty]". |
| `documentation.language_reference.objects.start.txt` | Graphics added to the list of objects. |
| `documentation.language_reference.start.txt` | One sentence pointing to the command line page for `args` and `ScriptDirectory`. |
| `documentation.standard_operators.start.txt` | One sentence pointing to the `Basis` option of NewOperator. |
| `documentation.language_reference.functions.operatorsetonsiteenergy.txt` | One sentence pointing to ShellAverage. |
| `documentation.language_reference.objects.responsefunction.functions.reducepoles.txt` | One sentence pointing to `NTriMax` of ChangeType. |
| `documentation.language_reference.objects.responsefunction.functions.calculateselfenergy.txt` | Why the direct difference is ill-conditioned, and a link to Downfold. |
| `documentation.language_reference.functions.readfplobasisfunctions.txt` | One sentence pointing to CalculateAtomicBasisFunctions. |


## Open points without a file

- `documentation:standard_operators:atomic_coulomb_repulsion`: the first sentence
  writes `NewOperator("AtomicU",indices,...)` without the `NF` that the example
  passes. Left for you to decide.
- `documentation:language_reference:objects:responsefunction:functions:downfold`:
  Kugler's original estimator could be cited next to the PRB reference.
- `documentation:getting_started` shows the start-up banner of a 2017 build.
- `documentation:language_reference:functions:createspectra`: the first example
  has "description text" as its description and loads `../definitions.Quanty`.
- No pages yet for `ImpurityGreensFunction`, `NewOperator("Dirac")` and
  `NewOperator("NuclearMultipole")`.
- The images in `media/` are on the site; the folder (never committed) can be
  deleted.
