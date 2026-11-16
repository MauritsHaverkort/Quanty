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

Batch of 2026-11-17, built on the live sources of that day; every example was run with the
build of commit `a6fa487a` and its output pasted. Images in `media/` (untracked), to upload to
the namespace given.

New pages:
- `documentation:language_reference:objects:graphics:functions:add`, `...:new`, `...:totable`,
  `...:rotate` (all already in 0.9.0, no tags). Images `add_elements.png`, `new_triangle.png`,
  `rotate_cards.png` to `documentation:language_reference:objects:graphics:functions`.
- The BitMap object: `documentation:language_reference:objects:bitmap:start`,
  `...:bitmap:functions:start`, `...:functions:new`, `...:functions:pixel`, `...:functions:save`
  (in 0.9.0, no tags; no properties page, a BitMap has no readable fields). Image
  `bitmap_drum.png` to `documentation:language_reference:objects:bitmap`, `bitmap_new.png` and
  `bitmap_save.png` to `documentation:language_reference:objects:bitmap:functions`.
  BitMap.Overlay and BitMap.FromMatrix get no page: they are empty in the code.
- `documentation:language_reference:functions:xcfunctional` (New in 0.9.1) and
  `documentation:language_reference:functions:sphericalbesselj`, title math.SphericalBesselJ
  (option Method New in 0.9.1).
- `documentation:standard_operators:dirac_operator` (Modified in 0.9.1; indexmenu_n 11) and
  `documentation:standard_operators:nuclear_multipole` (New in 0.9.1; indexmenu_n 12).

Changes to live pages:
- `documentation:standard_operators:atomic_coulomb_repulsion`: indexmenu_n 13 (after Dirac and
  NuclearMultipole), and one sentence linking the Dirac operator.
- `documentation:language_reference:functions:newoperator`: "Dirac", "AU" and
  "NuclearMultipole" in the first paragraph link their pages.
- `documentation:language_reference:functions:calculateatomicbasisfunctions`: the options
  DensityTolerance, NucleusLmax, NucleusNTheta, NucleusNPhi and the output
  NuclearMultipoleMoments added; links to the Dirac, NuclearMultipole and AU pages.
- `documentation:language_reference:objects:start`: BitMaps in the list of objects.
- `documentation:language_reference:objects:interpolatingfunction:start` and
  `...:interpolatingfunction:functions:spline`: the links to the missing page `hermitecspline`
  point to `monotoniccspline`, as internal links.
- `documentation:language_reference:functions:createspectra`: the first example has a
  description, its definitions written out instead of `dofile("../definitions.Quanty")`, and
  the output of the present version (the second output is a list of poles for this small
  space; the numbers differ in the 9th digit by the CODATA 2022 tesla).

The drafts of 2026-11-16 (version tags, CODATA 2022, the interpolating functions,
CalculateAtomicBasisFunctions, AtomNumberToName, AtomNameToNumber, Graphics.Plot) are live and
match the drafts (checked 2026-11-16, commit 2f697975); `atomic_potentials.png` is uploaded.

## Version tags

Every page with something that is not in release 0.9.0 (tag `0.9.0`, 2026-10-03) carries
`New in 0.9.1` or `Modified in 0.9.1`, red and right aligned under the page or section heading
(`;;#` block, layout from the live CalculateAtomicBasisFunctions page), or inline in red at the
end of a list item. Older tags stay, newest first: `Modified in 0.9.5 - New in 0.9.1`. Bug fixes
get no tag.

## Open points without a file

- `math.SphericalBesselJ` is wrong for negative x with the default method "Recurrence" (j0(-5)
  gives 0.018 instead of -0.192); 0.9.0 was right. The draft says x >= 0; drop that once the code
  is fixed.
- ImpurityGreensFunction is not registered in Lua (commented out in `LuaInterface/ssp.cpp`), so
  it gets no page.
