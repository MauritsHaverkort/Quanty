# Editing quanty.org — workflow and style reference

How to prepare changes to the user documentation at <https://www.quanty.org>, a
DokuWiki.


## The division of labour

Editing the wiki needs a login. An agent does not register accounts or enter
passwords, whoever authorises it, so an agent never edits the wiki directly:

1. **The agent drafts** complete, drop-in DokuWiki source and saves it in
   `Doc/Agents/Wiki/`, one file per wiki page, named after the page id
   (`documentation.basics.basis.txt` for `documentation:basics:basis`). Every file
   is the whole page, built on the live source with the change applied.
2. **A maintainer pastes** it into the page's edit box, previews and saves, and
   deletes the file once the page is live.

Nothing reaches a public page without a human having read it first.

Drafts are merged in **batches**. The maintainer pastes a whole batch, deletes the
files and commits. The agent then compares every live page with its draft and
reports on all of them at once in `WIKI_DRAFTS.md`, briefly: what differs, and
whether that matters. There is no page-by-page loop.

Useful URLs, given a page id like `documentation:basics:basis`:

| purpose | URL |
|---|---|
| read rendered | `https://www.quanty.org/documentation/basics/basis` |
| **read exact source** | `https://www.quanty.org/documentation/basics/basis?do=export_raw` |
| open the edit box | `https://www.quanty.org/documentation/basics/basis?do=edit` |

`?do=export_raw` needs no login. **Fetch the raw source of the page, or of a
neighbouring page for a new one, before drafting** — it is the only reliable way to
match the house markup, and it shows which plugins are in use. Edit that text
rather than writing a page afresh from the rendered view, which silently drops the
prose fences and the index footer.


## Images

Pages show PNG images that the maintainer uploads with the media manager, into the
namespace of the page, with lower-case file names. Make them from the example's
SVG output (written with `{{"RelativeSize",true}}`) with headless Firefox, at
2400 × 1697:

    /Applications/Firefox.app/Contents/MacOS/firefox --headless \
        --window-size=2400,1697 --screenshot <file.png> <file.svg>

`qlmanage` crops to a square, and `rsvg-convert` drops the bitmaps of colour
maps; do not use them. The PNGs go in `Doc/Agents/Wiki/media/`, which is **not
committed**: stage drafts by file name, never with `git add Doc/Agents/Wiki`.


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
  consistent across all pages and deliberate; drafts follow it.
- **`<code Quanty Example.Quanty>` … `</code>`** — Quanty script blocks. First
  argument is the syntax name, second the file name shown on the block. Use
  `Example.Quanty` unless the page refers to a specific named script.
- **LaTeX** — display maths in `$$ … $$`, taking full amsmath environments
  (`equation`, `eqnarray`, with `\nonumber`); inline maths in `$…$`.
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
  **every sibling's footer needs the new entry**, so draft those pages too.
- **Spelling** — the wiki writes "color", as the option names do.

### Version tags

Every page with something that is not in the latest release carries `New in
<version>` or `Modified in <version>`, red and right-aligned under the page or
section heading (a `;;#` block; copy the layout from the live
`CalculateAtomicBasisFunctions` page), or inline in red at the end of a list item.
Older tags stay, newest first: `Modified in 0.9.5 - New in 0.9.1`.
`WIKI_DRAFTS.md` says which release the current tags refer to.


## Site structure

Top-level namespaces in use: `documentation`, `physics_chemistry`, `forum`,
`download`, `mathematica`, `matlab`, `python`, `ctm4`, `crispy`, `dft`.

The language reference is `documentation:language_reference`, with `functions`,
`constants` and `objects`. Each object has a namespace with `start`,
`functions:start` and one page per function, and `properties:start`; the
`Graphics` object (`documentation:language_reference:objects:graphics`) is the
most recent example of the layout.

Everything published is user-facing — what functions do, never how they are
implemented. There is no developer section. Whether the notes in `Doc/Agents`
should ever become one, at top level (`development:`) or under
`documentation:development:`, is an open question for the maintainer.


## Checklist for a draft

- [ ] Raw source of the page (or of a neighbouring page) fetched with
      `?do=export_raw` and used as the base
- [ ] `{{indexmenu_n>N}}` set, with N checked against siblings
- [ ] Prose in `###` fences
- [ ] Code in `<code Quanty …>` blocks, and the code actually runs
- [ ] Links written as `[[target|]]`
- [ ] Nested Lua tables in running text written `{ {`, not `{{`: outside a
      `<code>` block DokuWiki reads `{{` as the start of an image
- [ ] Version tag set where the page documents something new or changed
- [ ] Index footers of sibling pages updated if a page was added
- [ ] Saved as one file per page in `Doc/Agents/Wiki/`, the whole page ready to
      paste; a change to an existing page is a file of its own too, not a
      line-level edit described in `WIKI_DRAFTS.md`
- [ ] Images made and listed for upload, with their namespace
