# Notes for coding agents

Quanty is a Lua interpreter with quantum-mechanical objects (operators, wave
functions, spectra, response functions) written in C. The user's Lua script is the
calculation. This file holds the rules that apply to every change; the
implementation details are in `Doc/Agents/`, to be read when a task touches them.

Read `CONTRIBUTING.md` (branches, backward compatibility, tests, C vs. C++, code
style, code layout) and `BUILDING.md` (compilers, macOS, CMake, flags to avoid)
first. They are short and are not repeated here.

The user documentation is the wiki at <https://www.quanty.org>. It says what a
function does; `Doc/Agents/` says how it is implemented. Read the wiki page of a
function before reverse-engineering it from the C source: append `?do=export_raw`
to a page URL to get its source.


## Rules

- **C, not C++**, even though the files end in `.cpp`: arrays with `malloc`/`calloc`,
  `struct` rather than `class`, no STL containers or strings. Two-space indentation,
  braces on their own line.
- **Use `CompactMatrixType`, not `MatrixType`**, in new code. Allocation was once
  more than 99% of the run time.
- **Backward compatibility is a hard constraint.** Lua function names, option names,
  defaults and output formats do not change. Add an option or a function instead.
- **Keep changes small.** One maintainer reviews everything. A bug fix fixes the bug.
  Do not add guards, asserts or warnings along the way: they will be added
  consistently across the whole codebase, as a project of its own.
- **Never delete a function**, even one with no callers. Rename it `…Old` and add a
  comment saying what replaced it. Functions already named `…Old` are left out of
  sweeps and convention changes.
- **One implementation per function.** Search the whole tree first, including
  `LuaInterface/lua/lmathlib.c`, which registers several `math.*` functions. A second
  method goes next to the first in the same file, named by method, and is selected
  with a `{{"Method", ...}}` option.
- **Numerical algorithms are discussed before they are written**: soundness,
  stability, cost and alternatives. Then wait for agreement.
- **Branches.** Work on a branch and merge into `develop` through a merge request on
  GitLab. Never commit to `develop` or `master` directly. Agents may commit on their
  own, but a human reviews before anything is pushed: push, and open a merge request,
  only when asked.


## Building and testing

- Build outside the working tree; the binary is not in `.gitignore`. On macOS:

      g++ LuaInterface/main.cpp -D LUA_USE_MACOSX -DBLASunderscore -O2 \
          -o <outside the repo>/Quanty -framework Accelerate -lreadline

  About 80 s. `g++ -fsyntax-only` reports a spurious `luai_hashnum` error from the
  bundled Lua; a real compile succeeds.
- The GitLab CI runs every test script on each push (Intel compiler, Linux). Locally,
  run only the tests that cover the change. Register a new test in the matching
  `TEST_NAME` list in `.gitlab-ci.yml`; the `All.Quanty` runners are no longer used.
- A test passes on **exit code 0**. Grepping the output for errors misses
  segmentation faults (exit 139, no message). Run a test by its path from the
  repository root, as the CI does, never after `cd`-ing into its folder.
- Tests run on Linux: no `math.random` (its sequence differs between systems), and
  thresholds derived from the data, not from a run on one machine. More in
  `Doc/Agents/Current/TESTING.md`.


## `Doc/Agents/` — implementation notes

The folder says what kind of note a file is:

- `Current/` describes the code as it is and is kept up to date with it.
- `Open/` lists known problems and planned work.
- `Projects/` holds the design of a project from start to finish. The first line
  of each says its status: "active plan" while the work runs, then "finished
  <date>", after which it is a dated record that is not updated. A description of
  the resulting code, kept current, goes into `Current/`.
- `Wiki/` is the workflow for quanty.org and the drafts of wiki pages.

A change that makes a note in `Current/` or `Open/` wrong corrects that note in the
same commit: a changed default, a renamed option, a fixed TODO entry. The reviewer
then sees the code and its description change together.

A note moves to another folder when its status changes. So that a move breaks
nothing, every file name is unique within `Doc/Agents`, and code comments and notes
cite a note by its file name alone ("`TODO.md` in `Doc/Agents`"), never by its
folder, and do not link to it with a relative path. This table is the only place
that gives the folder.

| File | Read it when |
|---|---|
| `Current/QUANTY_OVERVIEW.md` | new to Quanty: what a calculation looks like, the basis, the five-step script, the named operators |
| `Current/QUANTY_NOTES_ARCHITECTURE.md` | touching the C core: one-file build, the structs in `Core/StructDef.cpp`, hashing, restrictions, the ground-state solver, the response-function representations, the folder layout |
| `Current/CONVENTIONS.md` | before touching matrices, operators, rotations, response functions, Lua error messages or colours: notation (`^*`, `^T`, `^dag`) and index and storage conventions that fail silently |
| `Current/LUA_API_OVERVIEW.md` | adding or changing a Lua binding: where each kind is registered, the binding idioms, method names hidden in `__index`, value versus reference semantics, calling conventions that trap |
| `Current/TESTING.md` | before writing or debugging a test: the CI, what differs from a laptop, running a test, thresholds, reference values that change, CI-only crashes |
| `Current/NUMERICS_TRUNCATION.md` | accuracy, `epsilon` versus `zero`, restrictions, and how the Lanczos truncation relaxes with `Gamma` |
| `Current/DMFT.md` | DMFT and Green's-function building blocks: the loop function by function, the Weiss field, `Downfold`, bath excitations, large `U`, converging a loop |
| `Current/GRAPHICS.md` | plotting: the `Graphics` and `BitMap` objects, `Graphics.Plot`, colour maps, band structures, plotting block elements, output formats, rules for scripts |
| `Open/TODO.md` | before working in an area: known problems, deferred work and why it was deferred |
| `Open/TODO_ConversionAccuracy.md` | response-function conversions (`ChangeType`): what is measured and what is not |
| `Current/QuantyDFT_StageA_Implementation.md` | the relativistic atom (`DFT/`, `CalculateAtomicBasisFunctions`): Lua interface, methods, code map, verification |
| `Projects/BSPLINE_DESIGN.md` | why the B-spline C core in `Core/BasicMath/Interpolation/` looks as it does (finished 2026-11-09) |
| `Projects/QuantyDFT_StageA/` | the design of the relativistic atom, its checkpoints and decisions, and the Python prototypes (finished 2026-11-11) |
| `Wiki/WIKI_EDITING_WORKFLOW.md` | drafting pages for quanty.org |
| `Wiki/WIKI_DRAFTS.md` | the state of the current batch of wiki drafts |
