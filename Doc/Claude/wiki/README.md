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
- `documentation:language_reference:objects:interpolatingfunction:start`: rewritten.
  What interpolation is, a table of the four types, a section per type (Spline,
  HermiteCSpline, BSpline, DKBSpline with its large and small components),
  evaluation, properties (x, y, knots, degree, coefficients, type), arithmetic
  (which results are exact), integration. Built on the live page source of
  2026-11-13; the old example had `FtimesG` made from `y2` instead of `y3`.
  Image `interpolation_types.png` (in `media/`, 2400 x 1697) to upload to the
  namespace `documentation:language_reference:objects:interpolatingfunction`.
- New pages in `documentation:language_reference:objects:interpolatingfunction:functions`:
  `start`, `spline`, `hermitecspline`, `bspline`, `integrate`. All examples run
  on the present build, outputs pasted from those runs.

## Open points without a file

- `documentation:language_reference:functions:createspectra`: the first example
  has "description text" as its description and loads `../definitions.Quanty`.
- No pages yet for `ImpurityGreensFunction`, `NewOperator("Dirac")` and
  `NewOperator("NuclearMultipole")`.
