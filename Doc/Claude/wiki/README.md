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

Version tags (2026-11-16): every page with something that is not in release 0.9.0
(tag `0.9.0`, 2026-10-03) carries `New in 0.9.1` or `Modified in 0.9.1`, red and right
aligned under the page or section heading (`;;#` block, layout from the live
CalculateAtomicBasisFunctions page), or inline in red at the end of a list item.
Older tags stay, newest first: `Modified in 0.9.5 - New in 0.9.1`. Bug fixes get no tag.
Tagged in the drafts below; the pages without other changes:

- `documentation:language_reference:constants:physical_constants` and
  `...:constants:energy_units`: CODATA 2022 (commit 9ab94392), the outputs rerun
  (listing order of the live page kept, PhysicalConstants.cAU added at the end; Joule
  and Ha, cut off at the top of the live EnergyUnits listing, restored), one paragraph
  each on CODATA 2022 versus 2014. Built on the live sources of 2026-11-16.
- `documentation:language_reference:objects:interpolatingfunction:functions:bspline`,
  `...:spline`, `...:monotoniccspline` (HermiteCSpline): f.knots, f.degree,
  f.coefficients, f.type and f.x of a BSpline (commits 2e584f4b, 4caf4f00) in Output
  and in the examples that use them; B-spline arithmetic (a8a582b4). Live sources of
  2026-11-16, tags only.

- `documentation:language_reference:objects:graphics:functions:plot`: interpolating
  functions as data (PlotPoints + 1 equally spaced points and f.x, drawn as a list of
  points), in Data, PlotPoints and the color-map restriction. PlotPoints and NE count
  intervals (PlotPoints + 1 points; checked: 400 gives 401), the live text said
  "number of energies". Built on the live page source of 2026-11-16.
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
  The configuration is an option, no longer the second argument (2026-11-16):
  `CalculateAtomicBasisFunctions(28, {{"Configuration", "[Ar] 3d8"}})`; Input and the
  examples changed accordingly.
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
  `NewOperator("NuclearMultipole")`, `XCFunctional`, and no page for
  `math.SphericalBesselJ`; when written, all but ImpurityGreensFunction need a 0.9.1
  tag (NuclearMultipole and XCFunctional new; Dirac: Z real, c = cAU instead of 137.0;
  SphericalBesselJ: option Method).
