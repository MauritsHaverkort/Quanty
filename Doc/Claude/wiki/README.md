# Wiki drafts

DokuWiki source for <https://www.quanty.org>, one file per page, named after the
page id (`documentation.getting_started.txt` is the page
`documentation:getting_started`). Every file is the whole page, built on the
live source with the change applied. Paste the whole file
into the page's edit box, and delete the file here once the page is saved.

See `../WIKI_EDITING_WORKFLOW.md` for the house style. The notes of the earlier
batches (what was checked and how) are in the git history of this file, up to
commit `ec5a6a39`.


## Open pages

- `documentation:standard_operators:atomic_coulomb_repulsion`: the core form
  (kappas1, kappas2) documented, "1P-Operator" no longer tied to "conserve"
  (in the code it acts only in the core form, "conserve" only without it), and a
  second example, Fe [Ar] core on 3d and 4s, run on the present build. Built on
  the live page source of 2026-11-13.
- `documentation:language_reference:objects:interpolatingfunction:start`: a section
  on the fields x, y, knots, degree, coefficients and type, with an example run on the present
  build. The page still says interpolation is cubic splines only; the BSpline
  constructor is not yet documented. Built on the live page source of 2026-11-13.

## Open points without a file

- `documentation:language_reference:functions:createspectra`: the first example
  has "description text" as its description and loads `../definitions.Quanty`.
- No pages yet for `ImpurityGreensFunction`, `NewOperator("Dirac")` and
  `NewOperator("NuclearMultipole")`.
