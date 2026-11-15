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

- New pages `documentation:language_reference:functions:atomnumbertoname` and
  `documentation:language_reference:functions:atomnametonumber`: the element symbol of
  Z = 1 .. 118 and back (exported to Lua 2026-11-16); examples run on the present build.
- `documentation:language_reference:functions:calculateatomicbasisfunctions`: the
  potentials are now `rVN`, `rVH`, `rVxc` as BSplines (commit 50df4172): the Output
  entry, and a new section "The potentials" before "Accuracy" (what r V is, the
  derivative V' = (f' - f/r)/r, Kr LDA screening with a figure, the Ni 3d spin-orbit
  constant from dV/dr). The G entry links the interpolating function page. Built on
  the live page source of 2026-11-15; the four existing examples were rerun and
  their outputs are unchanged. Image `atomic_potentials.png` (in `media/`,
  2400 x 1697) to upload to the namespace `documentation:language_reference:functions`.
  Also the new call `G, F, Atom = CalculateAtomicBasisFunctions(...)`: `x` is gone
  (every returned function gives the quadrature points as `f.x`), and the fields
  `EtotEigenvalueSum`, `PotentialResidual`, `MixDensity`, `EavJJ` in Output; a
  configuration of empty orbitals only is the bare nucleus (Input, Configuration).
- `documentation:standard_operators:atomic_coulomb_repulsion` and
  `documentation:language_reference:functions:newoperator`: only the call
  `G, F, x, Atom = ...` becomes `G, F, Atom = ...` in the examples (outputs rerun,
  unchanged).
- `documentation:language_reference:objects:interpolatingfunction:start`: one
  sentence under Properties, f.x of the radial functions of
  CalculateAtomicBasisFunctions is now the quadrature grid also for the
  nonrelativistic BSpline; DiracRadialSolve (removed) is no longer named in the
  DKBSpline section; DKBSpline arithmetic (as a BSpline through f.x, result a BSpline)
  in the DKBSpline and Arithmetic sections, and the DKBSpline example now computes
  overlaps instead of showing the error (`G, F, Atom` call; rerun).

## Open points without a file

- `documentation:language_reference:functions:createspectra`: the first example
  has "description text" as its description and loads `../definitions.Quanty`.
- No pages yet for `ImpurityGreensFunction`, `NewOperator("Dirac")` and
  `NewOperator("NuclearMultipole")`.
