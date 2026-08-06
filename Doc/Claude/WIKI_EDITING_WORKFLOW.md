# Editing quanty.org — workflow and style reference

How to prepare changes to the user documentation at <https://www.quanty.org>.
Established August 2026.


## The division of labour

The wiki requires a login to edit, and I cannot register an account or enter a
password — that holds regardless of who authorises it, and it is the *mechanism*
that is blocked, not the editing. Relaying a password from a confirmation
e-mail does not change this. So:

1. **I draft** complete, drop-in DokuWiki source and save it under
   `Doc/Claude/wiki/`, one file per wiki page, named after the page
   (`development.code_layout.txt` for `development:code_layout`).
2. **Maurits pastes** it into the page's edit box, previews, and saves.

A side benefit worth keeping: nothing reaches a public page without a human
having read it first.

Useful URLs, given a page id like `documentation:basics:basis`:

| purpose | URL |
|---|---|
| read rendered | `https://www.quanty.org/documentation/basics/basis` |
| **read exact source** | `https://www.quanty.org/documentation/basics/basis?do=export_raw` |
| open the edit box | `https://www.quanty.org/documentation/basics/basis?do=edit` |

`?do=export_raw` needs no login and works from WebFetch. **Always fetch the raw
source of a neighbouring page before drafting** — it is the only reliable way to
match the house markup, and it shows which plugins are actually in use.

When drafting a *replacement* for an existing page, fetch its raw source first
and edit that text, rather than writing a fresh page from the rendered version.
The rendered view silently drops the prose fences and the index footer.


## House style

Reconstructed from the raw source of `documentation:start`,
`documentation:basics:start` and `documentation:basics:basis`.

### Page skeleton

```
{{indexmenu_n>1}}
====== Basis ======

###
Prose paragraph.
###

<code Quanty Example.Quanty>
-- a Quanty snippet
NF = 6
</code>

###
More prose.
###

===== Index =====
  - Basis
  - [[documentation:basics:operators|]]
  - [[documentation:basics:wave_functions|]]
```

### The elements

- **`{{indexmenu_n>N}}`** — first line. Sets the page's sort position within
  its namespace in the navigation menu. Pick N to place the page where it
  belongs among its siblings; check the neighbours' values first.
- **`====== Title ======`** — page title, once. `===== Section =====` for
  sections. `== Subtitle ==` appears on `documentation:start` as a small
  subtitle directly under the title.
- **`###` … `###`** — every prose paragraph is wrapped in these fences. This is
  consistent across all pages inspected and is clearly deliberate; drafts should
  follow it even though the purpose is presentational.
- **`<code Quanty Example.Quanty>` … `</code>`** — Quanty script blocks. First
  argument is the syntax name, second the filename shown on the block. Use
  `Example.Quanty` unless the page refers to a specific named script.
- **LaTeX** — display maths in `$$ … $$`, taking full amsmath environments
  (`equation`, `eqnarray`, with `\nonumber`); inline maths in `$…$`. The basis
  page uses this heavily, so a MathJax-class plugin is present.
- **Links** — `[[namespace:page|]]`. The empty pipe makes DokuWiki substitute
  the target page's own title, which is why index lists stay correct when a
  page is renamed. Prefer it over hard-coded link text.
- **Emphasis** — `//italic//`.
- **Auto index** — `{{indexmenu>.#1|msort}}` on a namespace `start` page lists
  its children one level deep, sorted by title. `documentation:start` uses
  `{{indexmenu>.#2|msort nsort}}` for two levels.
- **Manual index footer** — content pages inside `documentation:basics` end
  with an `===== Index =====` list of every sibling, with the current page as
  plain text and the rest as links. If a page is added to such a namespace,
  **every sibling's footer needs the new entry**, so list the affected pages
  when handing over a draft.


## Site structure

Top-level namespaces in use: `documentation`, `physics_chemistry`, `forum`,
`download`, `mathematica`, `matlab`, `python`, `ctm4`, `crispy`, `dft`.

Everything currently published is user-facing — what functions do, never how
they are implemented. There is no developer or internals section, which is the
obvious home for the material in `QUANTY_NOTES_ARCHITECTURE.md` and
`QUANTY_OVERVIEW.md` if that is ever wanted on the site. Whether such a section
belongs at top level (`development:`) or under `documentation:development:` is
an open question for Maurits; DokuWiki can move pages afterwards, but the
`indexmenu_n` values and any `===== Index =====` footers would need revisiting.


## Checklist for a draft

- [ ] Raw source of a neighbouring page fetched and used as the template
- [ ] `{{indexmenu_n>N}}` set, with N checked against siblings
- [ ] Prose in `###` fences
- [ ] Code in `<code Quanty …>` blocks, and the code actually runs
- [ ] Links written as `[[target|]]`
- [ ] Index footers of sibling pages updated if a page was added
- [ ] Saved as one file per page under `Doc/Claude/wiki/`
- [ ] Sibling pages needing a follow-up edit listed explicitly in the handover
