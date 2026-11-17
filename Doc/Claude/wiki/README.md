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

None. The batch of 2026-11-17 (Graphics Add/New/ToTable/Rotate, the BitMap object,
XCFunctional, math.SphericalBesselJ, the Dirac and NuclearMultipole operators, and the changes
to CalculateAtomicBasisFunctions, NewOperator, AU, objects start, the interpolating function
pages and CreateSpectra) was pasted the same day (commit fe35c5f2). The live pages were
compared with the drafts: the same, apart from blank lines, except
- math.SphericalBesselJ is the page `documentation:language_reference:functions:math.sphericalbesselj`
  (drafted as `...:sphericalbesselj`);
- the first paragraph of CalculateAtomicBasisFunctions was rewritten by Maurits (the occupation
  comes from "Configuration"; without options Hartree-Fock for the ground-state configuration and
  a Gaussian nucleus).

All six images of the batch are uploaded (`new_triangle.png` last, checked 2026-11-17).

## Version tags

Every page with something that is not in release 0.9.0 (tag `0.9.0`, 2026-10-03) carries
`New in 0.9.1` or `Modified in 0.9.1`, red and right aligned under the page or section heading
(`;;#` block, layout from the live CalculateAtomicBasisFunctions page), or inline in red at the
end of a list item. Older tags stay, newest first: `Modified in 0.9.5 - New in 0.9.1`. Bug fixes
get no tag.

## Open points without a file

- ImpurityGreensFunction is not registered in Lua (commented out in `LuaInterface/ssp.cpp`), so
  it gets no page.
