# Accuracy and truncation — `epsilon`, restrictions, and the Lanczos relaxation

Why Quanty has an accuracy knob at all, what `epsilon` actually does, and the
one non-obvious trick: during a spectrum calculation the truncation is allowed
to get **progressively coarser** with each Lanczos step, at no cost to the
accuracy of the spectrum, and the rate at which it coarsens is set by the
Lorentzian broadening Γ.

Line references are to the `Claude` branch at commit `e5b5f0d2`.


## 1. Why truncation is unavoidable

Many-body quantum mechanics is exponentially hard twice over:

- **In particle number.** With `N` spin-orbitals and `n` electrons there are
  `binom(N, n)` Slater determinants. The `documentation/basics/basis` page puts
  the target scale at systems with on the order of `10^100` possible
  determinants — the only reason that is tractable is that the determinants are
  never all written down.
- **In requested accuracy.** Getting a more accurate answer generally means a
  larger one-particle basis, which feeds straight back into the first point.

So "solve exactly" is not an option, and "solve to machine precision" is
usually not one either. Quanty's position is the one stated on the front page
of the documentation: **the user specifies an accuracy and the code tries to
reach it** — cheaply if the accuracy asked for is modest.

Everything below is the machinery behind that sentence.


## 2. The two ways to make the problem smaller

They are independent and compose; most real calculations use both.

### A priori — restrictions

Fix the subspace before starting: allowed occupations per group of orbitals.
`DeterminantRestrictionType` (`Core/StructDef.cpp:195`) is `N` filters, each a
bitmask plus a `Min`/`Max` occupation. A determinant is kept only if it passes
all of them (`QDetNotRestricted`, `Core/OperatorWaveFunction.cpp:2260`).

```lua
StartRestrictions = {NF, NB, {"111111 0000000000",6,6},
                             {"000000 1111111111",Nelec,Nelec}}
```

This is exact within the subspace — it is a statement about which physics you
want, not an approximation of the physics you asked for. It reaches the
calculation two ways: the `Restrictions` option of `Eigensystem` /
`CreateSpectra`, and an operator's own `H.Restrictions` field, which also
governs plain `H * psi`.

### A posteriori — `epsilon`, determinant selection by coefficient

Let the state pick its own basis, and throw away determinants whose expansion
coefficient is too small. This is what `epsilon` controls:

> `epsilon` — smallest absolute prefactor of a determinant in wave-functions
> considered different from 0 — `Core/StructDef.cpp:247`

The workhorse is `WaveFunctionTruncateBasis(psi, epsilon)`
(`Core/WaveFunctionArithmetic.cpp:953`), which dispatches to real, complex and
block variants. For complex states the test is on the modulus, done squared to
avoid the square root:

```c
const double Sqrepsilon = Sqr(epsilon);
if (Sqr(a->ValueRe[iBlok][iElement]) + Sqr(a->ValueIm[iBlok][iElement]) > Sqrepsilon)
```
`Core/ComplexWaveFunctionArithmetic.cpp:1072`

The determinant is not merely zeroed — it is removed from the hash table, so
the state genuinely shrinks and every later operator application is cheaper.
This is determinant selection driven by the state itself, rather than by a
subspace chosen in advance.


## 3. `epsilon` and `zero`

Two distinct knobs, easy to confuse, both settable from Lua on `Eigensystem`
and `CreateSpectra` (both spellings accepted, `LuaInterface/ssp.cpp:1855`,
`:3614`):

| option | meaning | default | default value |
|---|---|---|---|
| `epsilon` | smallest determinant coefficient kept | `sqrt(DBL_EPSILON)/100` | ≈ 1.5e-10 |
| `zero` | smallest number distinguishable from one in convergence tests | `100*sqrt(DBL_EPSILON)` | ≈ 1.5e-6 |

`InitConvergenceSettings()` (`Core/StructDef.cpp:269`) sets these, along with
the other numerical knobs on `ConvergenceSettingsType`: `DenseBorder` (1513 —
below this dimension a dense LAPACK diagonalisation is used instead of
Lanczos), `NKrylovStart/Max/Step/Min`, and `NTri` (200 — the number of Lanczos
steps in the single-shot tridiagonalisation used for spectra).

Note the four orders of magnitude between them: `epsilon` is a coefficient
threshold on the wave function, `zero` is a relative convergence tolerance on
scalars. Setting `epsilon` where you meant `zero` will usually still run, and
will quietly do something quite different.


## 4. Truncation when finding eigenstates

`Eigensystem` converges on **the state itself**, so it cannot afford to relax
its truncation: everything downstream reads that wave function.

In the ground-state Lanczos (`LanczosGroundStateRestrictedCC`,
`Core/Lanczos.cpp:1675`, and its siblings), `epsilon` is passed as a
`double *` and truncation is applied flat:

```c
ComplexWaveFunctionTruncateBasis(&v[0], *epsilon);
```

The only thing that changes it is running out of memory (§6). `zero` drives the
convergence test — on both the energy gain and the variance:

```c
while ((((energy - energy1) > (zero * (1.0 + fabs(energy)))) || minloop)
       && ((energysqr - Sqr(energy1)) > zero * (fabs(energy1) + DBL_EPSILON)))
```
`Core/Lanczos.cpp:1717`

So for eigenstates: `epsilon` sets how much of the wave function is carried,
`zero` sets when to stop iterating. Neither adapts.


## 5. Truncation when calculating spectra — the interesting case

This is the part worth understanding, because it is where the accuracy that is
actually *requested* is much weaker than the accuracy of the wave functions
involved, and the code exploits that.

### What a spectrum needs

`CreateSpectra` does not need the Krylov vectors. It needs the continued
fraction

```
G(ω) = ⟨ψ| T† ────────────────────────────── T |ψ⟩
                ω - a₁ - ─────b₁²─────
                          ω - a₂ - b₂²/(ω - a₃ - …)
```

evaluated at `ω + iΓ/2`. Only the scalars `a[i]`, `b[i]` survive; the vectors
are scaffolding. And crucially, **an error introduced at level `i` of the
continued fraction is suppressed in `G` by the levels above it** — the deeper
in the recursion, the less it matters, and the more broadening there is, the
faster that suppression sets in.

### How the code exploits it

Every `LanczosTriDiagonalize*` variant takes `Gamma` as its last argument and
carries a running `errorestimate`, initialised to 1.0. Once per Lanczos step:

```c
errorestimate *= M->b[i] * M->b[i] / (M->b[i] * M->b[i] + 4 * Gamma * Gamma);
```
`Core/Lanczos.cpp:2003` — and identically in the RC, CC, MC, all four
`Restricted` and all four `PlusDensity` variants, twelve sites in total.

The factor `b²/(b² + 4Γ²)` is at most 1, so `errorestimate` decreases
monotonically. Truncation then uses it as a *divisor*:

```c
RealWaveFunctionTruncateBasis(&(v[vi]), epsilon / errorestimate);
```
`Core/Lanczos.cpp:1992`

`errorestimate` shrinking means `epsilon / errorestimate` **grows**: the
threshold gets looser at every step, the Krylov vectors are held ever more
coarsely, and the calculation accelerates as it proceeds — while the accuracy
of `G(ω)` is preserved, because the levels being computed sloppily are exactly
the ones the continued fraction suppresses.

### The role of Γ

Γ is the Lorentzian broadening, equivalently the distance from the real axis at
which the Green's function is evaluated. It sets the relaxation rate:

| Γ | factor `b²/(b²+4Γ²)` | effect |
|---|---|---|
| Γ = 0 | 1 | no relaxation; `epsilon` stays flat throughout |
| Γ ≪ b | ≈ 1 − 4Γ²/b² | slow relaxation |
| Γ ≫ b | ≈ b²/4Γ² | fast relaxation, few steps stay expensive |

This is the formal version of a physically obvious statement: a heavily
broadened spectrum cannot resolve fine structure, so there is no point
computing the fine structure. Asking for a broad spectrum should be cheap, and
here it is.

The same argument applies to any evaluation away from the real axis — Matsubara
frequencies, or a self-energy evaluated at finite imaginary part — not only to
a Lorentzian lineshape.

### Two consequences that will bite

**Set `Gamma` explicitly if you want the speedup.** In `LuaCreateSpectra` the
tridiagonalisation is called as

```c
LanczosTriDiagonalize(*H, Tpsi[...], Options.Settings.epsilon, &GTri[...],
                      Max(Options.Settings.Gamma, Options.Settings.epsilon));
```
`LuaInterface/ssp.cpp:2429`

but the automatic choice of `Gamma`,

```c
if(Options.Settings.Gamma==0) Options.Settings.Gamma = 10*(Emax-Emin)/NE;
```
`LuaInterface/ssp.cpp:2490`

happens **after** the Lanczos has already run. So if `Gamma` is left
unspecified it is still 0 at tridiagonalisation time, `Max(0, epsilon)` is
`epsilon` ≈ 1.5e-10, and the relaxation never engages — you get the default
broadening in the output but pay the un-relaxed price to compute it. Passing
`{"Gamma", 1.0}` is not only a statement about the lineshape; it is what makes
the calculation cheap.

**Requesting the Krylov basis disables the relaxation.** The
`LanczosTriDiagonalizeKrylov*` variants (`Core/Lanczos.cpp:3536` onward) take
no `Gamma` argument at all — necessarily, since with `KrylovBasis` the vectors
are the deliverable and may not be degraded. Expect the cost of a spectrum to
rise when you ask to keep the basis.


## 6. The out-of-memory fallback

Distinct from everything above, and not an accuracy strategy — a survival
strategy. Wherever an allocation can fail, the pattern is:

```c
while (OperatorPsiRR(H, v[0], &(v[1])))
{
  printf("Allocation of new memory failed ...\n");
  printf("reduce accuracy of the wave-function by a factor of two from %22.15E to %22.15E\n\n",
         epsilon, epsilon * 2.0);
  epsilon *= 2.0;
  if (v[0].NBlok == 1) AllocFailed();
  RealWaveFunctionTruncateBasis(&(v[0]), epsilon);
}
```
`Core/Lanczos.cpp:1954`

`epsilon` doubles, the state is re-truncated, the operation is retried; when the
state is down to a single block and still will not fit, `AllocFailed()`. The
same loop appears throughout `Lanczos.cpp` and `BlockLanczos.cpp`.

**If you see "reduce accuracy of the wave-function by a factor of two" in a
run, the result is no longer at the accuracy you asked for.** It is printed
with plain `printf`, not through the verbosity filter, so it is always visible —
but it scrolls past easily in a long run and is worth grepping for.


## 7. Watching it work

The per-step diagnostics are behind verbosity bit `0x0100`:

```c
VerbosityPrintf(0x0100,"BasisSetSize: %12i ", v[0].N);
VerbosityPrintf(0x0100,"%12i ", v[vi].N);
VerbosityPrintf(0x0100,"| M.a[%3i] M.b[%3i]: %12.5E %12.5E | Err: %12.5E\n",
                i, i-1, M->a[i], M->b[i-1], errorestimate);
```
`Core/Lanczos.cpp:1951`, `:1993`, `:1994`

`VerbosityPrintf(mask, …)` prints when `mask & GlobalVerbosity`
(`IO/Print.cpp:6`), and `GlobalVerbosity` defaults to `0xF`
(`Headers/header.h:39`) — so this bit is **off** by default. To watch the
relaxation happen:

```lua
Verbosity(0x10F)
```

The two basis-set sizes per line are the state before and after truncation, and
`Err:` is `errorestimate`. What you want to see is the sizes rising for a while
and then falling away as `Err` drops. If `Err` stays at 1.0, `Gamma` did not
reach the tridiagonalisation — see §5.


## 8. Summary of the knobs

| knob | where | what it does |
|---|---|---|
| `Restrictions` | `Eigensystem`, `CreateSpectra`, `H.Restrictions` | fixes the subspace a priori; exact within it |
| `epsilon` | both | smallest determinant coefficient kept; default ≈ 1.5e-10 |
| `zero` | both | relative convergence tolerance; default ≈ 1.5e-6 |
| `Gamma` | `CreateSpectra` | broadening — *and* the relaxation rate of the truncation |
| `NTri` | `CreateSpectra` | Lanczos steps, i.e. continued-fraction depth; default 200 |
| `DenseBorder` | both | below this dimension use dense LAPACK; default 1513 |
| `Verbosity(0x10F)` | global | show basis sizes and `errorestimate` per step |

Rules of thumb that follow:

- Restrict first. It is exact and it shrinks the problem the most.
- Always pass `Gamma` to `CreateSpectra` explicitly.
- Loosen `epsilon` before reducing `NTri` — `epsilon` trades accuracy you
  probably do not need, `NTri` truncates the continued fraction and can
  distort lineshape.
- Grep the output for "reduce accuracy" before trusting a long run.
