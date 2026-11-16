# Wiki drafts

DokuWiki source for <https://www.quanty.org>, one file per page, named after the
page id (`documentation.getting_started.txt` is the page
`documentation:getting_started`). Every file is the whole page, built on the
live source with the change applied. Paste the whole file
into the page's edit box, and delete the file here once the page is saved.

See `../WIKI_EDITING_WORKFLOW.md` for the house style. The notes of the earlier
batches (what was checked and how) are in the git history of this file, up to
commit `ec5a6a39`, and for the batch of 2026-11-16 up to `a6fa487a`.


## Open pages

None. The 12 drafts of 2026-11-16 (version tags, CODATA 2022, the interpolating functions,
CalculateAtomicBasisFunctions, AtomNumberToName, AtomNameToNumber, Graphics.Plot) were pasted
the same day (commit 2f697975). The live pages were compared with the drafts afterwards: they
are the same, apart from blank lines at the end and the HermiteCSpline link below.

## Version tags

Every page with something that is not in release 0.9.0 (tag `0.9.0`, 2026-10-03) carries
`New in 0.9.1` or `Modified in 0.9.1`, red and right aligned under the page or section heading
(`;;#` block, layout from the live CalculateAtomicBasisFunctions page), or inline in red at the
end of a list item. Older tags stay, newest first: `Modified in 0.9.5 - New in 0.9.1`. Bug fixes
get no tag.

## Open points without a file

- The image `atomic_potentials.png` (in `media/`, 2400 x 1697) is not on the site: the page
  `documentation:language_reference:functions:calculateatomicbasisfunctions` shows it, the file
  gives 404 (checked 2026-11-16). Upload it to the namespace
  `documentation:language_reference:functions`.
- The page `documentation:language_reference:objects:interpolatingfunction:functions:hermitecspline`
  does not exist; HermiteCSpline is documented on `...:functions:monotoniccspline`. Two links
  still point to the missing page: on `...:interpolatingfunction:functions:spline` (Example, "for
  data that must not overshoot") and in the text of `...:interpolatingfunction:start` (the
  paragraph on HermiteCSpline). The table on the start page now links the right page, but as a
  full web address; the internal form is
  `[[documentation:language_reference:objects:interpolatingfunction:functions:monotoniccspline|...]]`.
- `documentation:language_reference:functions:createspectra`: the first example has
  "description text" as its description and loads `../definitions.Quanty`.
- No pages yet for `ImpurityGreensFunction`, `NewOperator("Dirac")` and
  `NewOperator("NuclearMultipole")` (only the overview on the NewOperator page),
  `XCFunctional`, and `math.SphericalBesselJ`; when written, all but ImpurityGreensFunction
  need a 0.9.1 tag (NuclearMultipole and XCFunctional new; Dirac: Z real, c = cAU instead of
  137.0; SphericalBesselJ: option Method).
- Graphics: the pages Add, New, ToTable and Rotate under
  `documentation:language_reference:objects:graphics:functions`, and BitMap as its own object,
  are not written yet (the start, functions start, Plot, ToSVG, ToPS, SVGHTMLViewer, PlotBands
  and properties pages exist).
