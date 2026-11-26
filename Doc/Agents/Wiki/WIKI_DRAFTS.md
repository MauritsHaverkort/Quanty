# Wiki drafts

DokuWiki source for <https://www.quanty.org>, one file per page, named after the
page id (`documentation.getting_started.txt` is the page
`documentation:getting_started`). Every file is the whole page, built on the
live source with the change applied. Paste the whole file
into the page's edit box, and delete the file here once the page is saved.

See `WIKI_EDITING_WORKFLOW.md` for the workflow and the house style. The notes of the earlier
batches (what was checked and how) are in the git history of this file, up to
commit `ec5a6a39`, for the batch of 2026-11-16 up to `a6fa487a`, and for the batch of
2026-11-17 up to `d18b9abb`.


## Open pages

Batch of 2026-11-26, the Stevens operators (branch `StevensOperators`). Seven new pages,
no existing page changes: both namespaces list their pages with `indexmenu`, so no
index footer needs an entry. No images.

| file | page |
|---|---|
| `documentation.standard_operators.stevens_operators.txt` | the overview: definitions of the five normalizations, low order forms, remarks, the Ce³⁺ example with Operators. `indexmenu_n>14`, after "Atomic coulomb repulsion" (13) |
| `documentation.language_reference.functions.stevenso.txt` | `StevensO`, with the zero field splitting of a spin 1 as example |
| `documentation.language_reference.functions.stevensd.txt` | `StevensD` |
| `documentation.language_reference.functions.stevensz.txt` | `StevensZ` |
| `documentation.language_reference.functions.stevensc.txt` | `StevensC` |
| `documentation.language_reference.functions.stevensy.txt` | `StevensY` |
| `documentation.language_reference.functions.stevensprefactor.txt` | `StevensPrefactor` |

What was checked:

- Every example was run with the build of commit `bf448269` (exit 0) and its output pasted.
  The matrix examples build the spin matrices in the script, ordered from m = -S to S;
  checked that [Sx, Sy] = i Sz and that S+ raises m.
- The low order forms on the overview page (O20, O21, O2-2, O22, O40, O44), D(1,m) = Jx, Jy,
  Jz and C(1,m) were compared with the matrices of J = 5/2 (deviations below 5e-14).
- The Result block of the Ce³⁺ example shows only the lines the script prints. Eigensystem
  also prints its progress and, at present, a warning that the operator "is not fully
  Hermitian" (see "Stevens operators are Hermitian only to rounding" in `TODO.md`). The
  energies are right: they are the analytic -2ζ - 240 B4, -2ζ + 120 B4 and 3ζ/2 +
  (-1080, 120, 840) B4.

## Version tags

The current tags are `New in 0.9.1` and `Modified in 0.9.1`, for everything that is
not in release 0.9.0 (tag `0.9.0`, 2026-10-03). How tags are written is in
`WIKI_EDITING_WORKFLOW.md`.

## Open points without a file

- ImpurityGreensFunction is not registered in Lua (commented out in `LuaInterface/ssp.cpp`), so
  it gets no page.
