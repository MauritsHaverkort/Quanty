# Testing — what the CI runs, and how to write a test that passes there

The tests are Lua scripts, run one by one by the GitLab CI on Linux. Most problems
with tests come from the difference between that machine and the one a change is
developed on. This note says how the tests are run, how to write one, and how to
read a failure.


## 1. The layers

- **`TestScripts/`** — the tests. One `.Quanty` script per test, in a folder per
  area (`Math/`, `Matrix/`, `Responsefunction/`, `DiracAtom/`, ...), plus
  `TestScripts/NoSinglePrecision.sh`, which checks that no single-precision type
  enters the DFT code.
- **`.gitlab-ci.yml`** — the CI. `build-job` compiles Quanty once; each
  `test-jobs-<area>` job runs every script of its `TEST_NAME` list as
  `./Quanty.exe "TestScripts/<Area>/${TEST_NAME}.Quanty"`, each in its own process.
  **A test that is not in a `TEST_NAME` list is never run.** The `All.Quanty`
  runner scripts in some folders are no longer used; do not add tests to them.
- **`Test/`** — an older C test program with its own entry point (`Test/Test.cpp`).
  It is not part of the Quanty binary and not run by the CI.
- **`Debug examples/Test_All.sh`** with its `Reference/` set — the old test
  runner; see the entry on `Debug examples/` in `TODO.md`.

A test **passes when Quanty exits with code 0**. A Lua error prints "Error while
executing the script:" and exits 1; a check in a script fails by calling `error`.


## 2. Where the CI differs from a laptop

`build-job` compiles with the Intel compiler, `icpx -O2 -ipo -qopenmp
-fp-model=precise -fimf-precision=high -march=core-avx2`, against MKL. A Mac build
uses clang, Accelerate and usually no OpenMP. So the CI differs in compiler,
optimiser, LAPACK and threading at once, and three kinds of failure appear only
there:

- **Undefined behaviour.** `icpx -ipo` inlines across the whole translation unit
  and exploits every rule of the language. The case that cost most: an array of
  `void const *` passed on cast to `BlockListOfPoles const **` breaks strict
  aliasing; `icpx` treated the array as never written and read a stale register
  (fixed in `0026a400`). ASan, UBSan and `clang -fpointer-tbaa` locally were all
  clean. **Never pass an array of `void *` cast to an array of `T *`**; use one
  typed array per case.
- **A different LAPACK.** Results limited by an eigensolver cancel differently in
  MKL and Accelerate.
- **A different C library.** `math.random` gives a different sequence.


## 3. Running a test locally

- Build outside the working tree (see `AGENTS.md`), and run only the tests that
  cover the change. The CI runs everything on every push.
- **Run a test by its path from the repository root**, exactly as the CI does:

      <build>/Quanty TestScripts/Matrix/Sqrt.Quanty > out.txt 2>&1; echo "exit=$? lines=$(wc -l < out.txt)"

  Several scripts locate their data and plot files relative to their own path
  (`debug.getinfo(1, "S").source:sub(2):match("(.*/)")`), which is `nil` when the
  path contains no `/`. Run after `cd`-ing into the folder, such a script dies on
  "attempt to concatenate global 'directory'".
- **Judge by the exit code, and compare the number of output lines with a good
  run.** A segmentation fault prints nothing and exits 139, so grepping the
  output for "Error" reports a crash as a pass — and a test that dies early can
  look like a large speed-up.
- **Piped output is block buffered**, and the buffer is lost when the process
  dies: the last line in a log is not where it died, and a C `printf` can appear
  in the middle of a Lua `io.write` line. Re-run under `script -q /dev/null`
  (macOS) or `stdbuf -o0 -e0` (Linux) to see the true last output.
- **Running the tests does not normally rewrite the `.svg` files** under
  `TestScripts/*/OutputPlots/`; they come out byte-identical. A few tests do write
  them, so check `git status` after a run rather than reflexively reverting: a
  changed plot is a signal, and the one time it happened it was an uninitialised
  read.


## 4. Writing a test

- **Register it** in the matching `TEST_NAME` list in `.gitlab-ci.yml`, or in a new
  `test-jobs-<area>` job for a new folder.
- **No `math.random`.** Its sequence depends on the C library. Use a generator in
  the script, for example Park–Miller, `x = 16807 x mod (2^31 - 1)`, which is exact
  in doubles. (`rand()` inside Quanty is seeded from the time in
  `LuaInterface/main.cpp`, so the starting vectors of the Lanczos routines differ
  from run to run as well; a test must pass for any of them.)
- **Thresholds come from a bound or from both platforms**, never from one local
  run. Take the accuracy that limits the quantity — an eigensolver's, a basis
  function's own error — or measure on the CI machine before tightening. When only
  the Mac number is known, keep a margin of at least the spread seen between
  platforms before, about 10 to 50 times. A threshold of 1e-10 set from 5.2e-12 on
  a Mac failed at 1.2e-10 in the CI (fixed in `a651f281`).
- **Derive expected ranges from the data the script builds**, not from numbers
  read off a run.
- **Ask for no more accuracy than the settings deliver.** A test comparing
  eigenvalues to 1e-5 while passing `Zero = 1.5e-6` failed now and then: a small
  variance does not imply a small energy error in a near-degenerate spectrum
  (fixed in `ef84d295` with `Zero = 1e-8`).
- **Never ask for a number of states that ends inside a degenerate multiplet.**
  The top states converge last and their mixture depends on where the block
  Lanczos stops; choose `Npsi` to fall between multiplets (`7d490bfc`).
- **Include a complex, non-symmetric case** for anything that touches the
  conventions of `CONVENTIONS.md`; a real symmetric case passes whichever way the
  convention is broken.
- **Keep everything `local`.** When one script runs others with `dofile`, a global
  in one is visible in the next, and a global named like a library function
  shadows it.


## 5. When a reference value changes

A reference table recorded at one convergence level breaks when something
converges further, and that is not evidence the change is wrong. Before calling it
a regression, check whether the new answer is better converged. In the case that
taught this (the multiplet ligand field section of `TestScripts/Eigensystem.Quanty`,
September 2026), tightening `Zero` on unmodified code moved the top four of ten
states by up to 3.5e-4 in energy. Three signs said the tighter answer was right:

- the energies moved **uniformly down** — a variational estimate is an upper
  bound;
- the converged states formed clean near-degenerate doublets, the references did
  not;
- the real and complex calculations, which agreed only to 1e-5 in the references,
  agreed to 1e-13 once converged.


## 6. Debugging a failure that happens only in the CI

Reproduce on a Linux machine with `icpx` (in the group: `lin76n`), building with
the command of `build-job` plus `-g`. There `gdb` may not attach to a running
process (`ptrace: Operation not permitted`), but core dumps go to
systemd-coredump:

    coredumpctl gdb Quanty.exe

Useful inside `gdb`, one command per line: `p $_siginfo._sifields._sigfault.si_addr`,
`x/4i $pc`, `info registers`, `frame N`. Run the script under `stdbuf -o0 -e0` so
the log shows the true last line.
