# Quanty — what it is and how a calculation is shaped

Orientation notes for someone who has to work on the code and needs the
user-facing picture first. Companion to `QUANTY_NOTES_ARCHITECTURE.md` (which
describes the internals) and `LUA_API_OVERVIEW.md` (which lists the exports).

Two sources are mixed here and kept distinguishable:

- **[site]** — the user documentation at <https://www.quanty.org>, i.e. Maurits'
  own account of the concepts.
- **[code]** — checked directly against the `Claude` branch at commit `e5b5f0d2`.

Where the two could drift, the code reference is given so the claim can be
re-checked.


## 1. What Quanty is

> "a script language which allows the user to program quantum mechanical
> problems in second quantization" — [site], front page

The program was written to compute **core-level spectroscopy of correlated
materials** — XAS, XMCD, RIXS, NIXS, photoemission — and grew from there into a
general second-quantisation toolkit: exact diagonalisation, configuration
interaction, ligand-field and Anderson-impurity models, DMFT impurity solving.

Practically: Quanty is a **Lua 5 interpreter with physics objects bolted on**.
The user writes a Lua script; the objects it manipulates (`Operator`,
`Wavefunction`, `Spectra`, `ResponseFunction`, ...) are C structs behind
userdata with metatables, so ordinary Lua arithmetic — `H = a*OppA + b*OppB`,
`val = psi * H * psi` — drives the numerics. There is no input-file format and
no separate driver; the script *is* the calculation.

This explains several things about the codebase that otherwise look odd:

- Control flow, loops, temperature averaging, file output and plotting are all
  the user's problem, done in Lua. The C side only has to be fast at the
  primitives. Compare the Boltzmann-average loop in
  `Debug examples/Tutorials/03_XAS/03_XAS.Quanty`, which is plain Lua around
  calls to `CreateSpectra`.
- The Lua binding layer is large relative to the numerics, because every
  physics concept needs a constructor, a metamethod set, and option parsing.
- Legacy `main()` functions survive in `DMFT/`, `Hubbard_Model/`,
  `Tight_Binding/`, `Model/` from the era before the Lua front end, when an
  input file *was* a C `main`.

### Where it came from, and what follows from that

Quanty has a 20+ year history [site]. It began as a tool for x-ray
spectroscopy, which is why most of the examples are still core-level
spectroscopy, and grew outward from there.

Its original claim on users, against the established CTM4/Cowan-derived
multiplet codes in the field, was that it was **faster and more flexible** —
enough so that it was adopted even inside the groups maintaining those codes.
That comparison is Maurits' own and is now 10+ years old, so it is history
rather than a current benchmark. It does explain the culture of the codebase:
performance is not a late-stage concern here, it is the reason the program was
adopted, and it is why the C-in-C++ restriction and the allocation discipline
(see `CompactMatrixType` in the architecture notes) are enforced as firmly as
they are.

Two constraints follow, and both should be treated as hard:

- **Backward compatibility is guaranteed.** Scripts people wrote years ago are
  expected to keep running. Any change to Lua-visible behaviour — function
  names, option names, defaults, output formats — needs to be considered in
  that light, not just the C-level refactor that motivated it.
- **Maintenance capacity is finite.** Maurits runs a group *and* is personally
  involved in code management. Proposals should be sized accordingly: a change
  that is correct but demands a large coordinated rewrite is, in practice, a
  change that will not happen. This is the same reasoning behind the deferred
  items in `TODO.md`.

The open performance question is the other side of the same coin: the
shared-memory parallelism that made Quanty fast on the machines of the time
does not scale well to a 256-core node. See `TODO.md`.


## 2. The physical model

### Basis  [site: `/documentation/basics/basis`]

A calculation is fixed by two numbers:

| | meaning | index range |
|---|---|---|
| `NF` (`NFermion`) | number of fermionic modes — spin-orbitals | `0 … NF-1` |
| `NB` (`NBoson`) | number of bosonic modes — e.g. phonons | `NF … NF+NB-1` |

A "spin-orbital" is whatever single-particle state the physics calls for:
atomic (n, l, m, σ) for a multiplet calculation, Wannier or molecular orbitals
taken from DFT or Hartree-Fock for a solid or a molecule.

Spin-orbitals are grouped into **shells** as plain Lua index lists, which is how
the named operator constructors learn what is a p-shell and what is a d-shell:

```lua
NFermion = 16
NBoson   = 0
IndexDn_2p = {0,2,4}                -- 2p, spin down
IndexUp_2p = {1,3,5}                -- 2p, spin up
IndexDn_3d = {6,8,10,12,14}         -- 3d, spin down
IndexUp_3d = {7,9,11,13,15}         -- 3d, spin up
```

**Convention worth internalising:** within a shell the orbitals run over
`m = -l … +l`, and the two spins are *interleaved* — even index spin-down, odd
index spin-up. Shells are laid out consecutively, core shell first by custom
(not by requirement). Nothing enforces this; it is a convention the index lists
carry, and the named constructors trust the lists they are handed.

### Occupation representation  [site + code]

A many-particle basis state is a Slater determinant, stored as a bit pattern:
one bit per fermionic mode, and **8 bits per bosonic mode** (so a boson mode can
hold occupation 0–255). Hence the word count for a determinant, computed
identically everywhere it appears:

```c
NWord = NBoson + Min(NFermion, 1 + (NFermion - 1) / 8);
```

`LuaInterface/ssp.cpp:600`; `NWord > 4096` is rejected there. This is the same
invariant flagged in `QUANTY_NOTES_ARCHITECTURE.md` — every code path that
builds an operator by hand must set it, or restrictions silently pass
everything.

A general state is a linear combination of determinants; that is exactly what
`WaveFunctionType` stores (a hash table keyed on the bit pattern).


## 3. The canonical script arc

Nearly every real Quanty script is the same five movements. Using
`Debug examples/Tutorials/03_XAS/03_XAS.Quanty` as the reference specimen — a
2p → 3d L-edge XAS calculation for a divalent transition-metal ion:

**1. Declare the basis.** `NFermion`, `NBoson`, and the shell index lists.

**2. Build operators.** Named constructors for the standard physics, combined
with ordinary arithmetic into a Hamiltonian:

```lua
OppF2   = NewOperator("U",     NFermion, IndexUp_3d, IndexDn_3d, {0,1,0})
Oppldots= NewOperator("ldots", NFermion, IndexUp_3d, IndexDn_3d)
Akm      = PotentialExpandedOnClm("Oh", 2, {0.6,-0.4})
OpptenDq = NewOperator("CF",   NFermion, IndexUp_3d, IndexDn_3d, Akm)

Hamiltonian = F0dd*OppF0 + F2dd*OppF2 + F4dd*OppF4
            + tenDq*OpptenDq + zeta_3d*Oppldots
```

Note the idiom of defining each Slater-integral component as its own operator
with a unit coefficient (`{1,0,0}`, `{0,1,0}`, `{0,0,1}`) and putting the
numbers in only at the summation step. That keeps the parameter scan cheap: the
operators are built once, the sum is rebuilt per parameter set.

A core-level calculation carries **two** Hamiltonians — initial state, and a
final state with the core hole (extra spin-orbit, core-valence Coulomb, and
screened valence Slater integrals).

**3. Find eigenstates.**

```lua
StartRestrictions = {NFermion, NBoson, {"111111 0000000000",6,6},
                                       {"000000 1111111111",Nelec,Nelec}}
psiList = Eigensystem(Hamiltonian, StartRestrictions, Npsi)
```

Each restriction is `{mask, min, max}` — the mask selects orbitals (spaces are
cosmetic), and the occupation of the selected set must lie in `[min,max]`. Here:
exactly 6 electrons in 2p, exactly `Nelec` in 3d. `Eigensystem` returns a Lua
list of the `Npsi` lowest states. It can equally be started from an explicit
wave function instead of restrictions, in which case it finds the lowest states
reachable from that starting vector [site: `/documentation/basics/eigen_states`].

**4. Extract observables.** Expectation values are literally a triple product,
with the bra conjugated automatically [site]:

```lua
E = psiList[1] * Hamiltonian * psiList[1]
```

Handing a *list* of operators gives a list of values in one call, which is the
usual way to print a level table:

```lua
oppList  = {Hamiltonian, OppSsqr, OppLsqr, OppJsqr, OppSz, OppLz}
expvalue = value * oppList * value
```

Spectra come from a Green's function [site: `/documentation/basics/spectra`]:

```
G(ω) = ⟨ψi| T† 1/(ω - H + iΓ/2) T |ψi⟩
```

```lua
Spectra_z = CreateSpectra(XASHamiltonian, TXASz, psiList[j],
              {{"Emin",-10}, {"Emax",20}, {"NE",3500}, {"Gamma",1.0}})
```

The choice of transition operator `T` is what selects the experiment —
annihilation for photoemission, creation for inverse photoemission, a
creation-annihilation product for absorption.

**5. Finite temperature, output, plotting.** A plain Lua loop over `psiList`
weighting by `exp(-(E_j - E_0)/T)`; `Spectra` objects support `+` and scalar
`*`, so Boltzmann averaging is arithmetic on spectra. Then `.Print{{"file",...}}`
and, typically, a gnuplot script written and `os.execute`d from Lua.

Two supporting conveniences visible in that script: `EnergyUnits.Kelvin.value`
and `EnergyUnits.Tesla.value` convert to the internal unit (eV) — the unit table
is built in `LuaInterface/units.cpp:233` — and `Verbosity(0)` silences the
progress reporting that longer runs emit.


## 4. Operators in more detail

### Named constructors  [code: `LuaInterface/ssp.cpp:659` onward]

`NewOperator("Name", ...)` dispatches on the leading string. The full set of
recognised names, in dispatch order:

    Coulomb/U · AtomicU/AU · ContinuumU/CU · Dirac · H0
    Sx Sy Sz Splus Smin Ssqr/S2
    Lx Ly Lz Lplus Lmin Lsqr/L2
    Jx Jy Jz Jplus Jmin Jsqr/J2
    Tx Ty Tz Tplus Tmin Tsqr
    Qx Qxx Qy Qyy Qz Qzz Qxy Qyx Qxz Qzx Qyz Qzy
    ldots · CF
    Create/Cr · Annihilate/An · Number/CrAn/n · Identity/I/1
    HarmonicOscillatorOnStandingWaves/HOSW · HarmonicOscillatorOnGrid/HOGrid

`CF` is the workhorse: expanded on spherical harmonics via an `Akm` list, it
serves both as crystal/ligand field *within* a shell and as a multipole
transition operator *between* two shells — the XAS dipole operators in the
tutorial are `CF` with `Akm = {{1,-1,t},{1,1,-t}}` and 3d/2p index lists.
`PotentialExpandedOnClm("Oh", 2, {0.6,-0.4})` generates the `Akm` for a
point-group potential.

### Raw form  [code: `LuaInterface/ssp.cpp:563`]

When the first argument is a *number* rather than a string, `NewOperator` builds
an arbitrary ladder-operator sum directly:

```lua
Op = NewOperator(NFermion, NBoson, { {i, -j, value}, ... }, options)
```

Each inner table is one term: a list of mode indices followed by the
coefficient. **Sign encodes the ladder direction** — positive creates, negative
annihilates. The test uses `signbit`, so `-0` correctly annihilates mode 0 while
`0` creates it; that is a genuine trap when indices are computed rather than
written literally. Indices are bounds-checked against `NF+NB` and must be
integral. The Lua list order is reversed when written into the internal string,
and the term is inserted through `OperatorAddLaderNormalOrder`.

The optional fourth argument accepts `Restrictions`/`restrictions`, `Name`, and
`NBitsKey`.

### Mapping to the internals

The five structures that carry the physics, and what the user sees of each:

| user sees | struct | notes |
|---|---|---|
| `Wavefunction` | `WaveFunctionType` | hash table of determinants |
| `Operator` | `OperatorType` | terms bucketed by string length: `l[0]` scalar, `l[2]` one-body, `l[4]` two-body |
| `Spectra` | `SpectraTypeArray` | several `SpectraType` on one energy grid; supports `+` and scalar `*` |
| `ResponseFunction` | `ResponsefunctionType` | thin tag + `void *data`; concrete form chosen by tag |
| `TightBinding` | `TightBindingDef` | periodic hopping model; the route in for DFT data |

Supporting objects — `Complex`, `Matrix`, `Tensor`, `InterpolatingFunction` —
plus the internal-only structs are tabulated in `QUANTY_NOTES_ARCHITECTURE.md`.
`DeterminantRestrictionType` is not an object: it is the `N` filters (mask, min,
max) built from the restriction list, and it is also reachable as an operator's
own `.Restrictions` field.

That `.Restrictions` field governs plain `H * psi` products — see the
architecture notes, since this is a recurring source of confusion.


## 5. Where the real documentation lives

<https://www.quanty.org> is a DokuWiki and fetches cleanly page by page. The
parts worth knowing about:

**`/documentation/start`** — the index.

- `/documentation/basics/start` — the conceptual chapter this section
  paraphrases: `basis`, `operators`, `wave_functions`, `expectation_values`,
  `eigen_states`, `spectra`, `resonant_spectra`, `fluorescence_yield`.
- `/documentation/standard_operators/start` — one page per named operator
  family (S, L, J, T, Q, spin-orbit, Coulomb, crystal field, ligand field,
  transition operators, atomic Coulomb). The place to check what an `Akm` list
  or a Slater-integral argument actually means.
- `/documentation/language_reference/start` — `functions`, `constants`,
  `objects`. The authoritative per-function reference; more complete than
  `LUA_API_OVERVIEW.md`, which was assembled from the binding tables.
- `/documentation/tutorials/start` — quick start, introduction, NiO crystal
  field, NiO ligand field, connection to DFT, model examples. These correspond
  to the scripts under `Debug examples/Tutorials/`.
- `/documentation/memory`, `/documentation/parallel` — memory management and
  parallelisation notes.

**`/physics_chemistry/start`** — the physics background, alphabetical:
one-particle orbitals, point groups, Coulomb repulsion and Slater integrals,
*U* and Δ in ligand-field / Anderson-impurity / DMFT, photon-matter interaction.

Also on the site: `/authors`, `/how_to_cite`, `/publications`, the
Mathematica / MATLAB / Python front ends (`Quanty.nb`, `Quanty.m`, `Quanty.py`),
and the external programs that drive Quanty (CTM4, Crispy).

The site is user-facing: it documents *what the functions do*, never how they
are implemented. For implementation, the repo docs are the only source.


## 6. Still to write

Identified while assembling this, not yet covered anywhere in `Doc/Claude`:

- **Testing and verification** — `Test/`, `TestScripts/*.Quanty`, and
  `Debug examples/Test_All.sh` with its `Reference/` set. What each layer
  covers and how to check a change.
- **Spectra and response functions** — the concrete `ResponsefunctionType`
  representations, and the Lanczos family (`FiniteTemperatureLanczos`,
  `ImpurityLanczos`, `LanczosDyson`).
- **Physics-model layers** — `AtomicOperators/` (Slater integrals, Coulomb,
  Dirac), `Ligand_Field_Theory/OperatorDeffinitions`, and `IO/` (Wannier90,
  FPLO, Stuttgart NMTO) as the route for external electronic-structure data.
- **Call path** — how one Lua call travels from `ssp.cpp` down into `Core/`.

Written since: `NUMERICS_TRUNCATION.md`, on `epsilon`, restrictions, and how the
Lanczos truncation relaxes with Γ.
