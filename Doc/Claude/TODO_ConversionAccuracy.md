# TODO — accuracy of the response function conversions, and letting the user choose

**Status:** open. Nothing here is known to be wrong; the point is that what has
been *measured* is narrower than what is *claimed*, and that the right long-term
answer is a choice offered to the user rather than one algorithm fixed by us.

Kept apart from [TODO.md](TODO.md) because it is a programme rather than a defect:
it wants a proper accuracy harness first, then an `Algorithm` option on
`ChangeType`.


## What is actually measured today

Every conversion test in `TestScripts/Responsefunction/` uses the same probe:

```lua
max over a uniform grid of  |G_new(w) - G_old(w)| / |G_old(w)|
```

with the grid spanning the band in 300 to 400 points and **two** broadenings,
`Gamma = 0.3` and `Gamma = 0.02`. That is a coarse, uniform, real-axis probe at
two values of one parameter. Everything reported as "1e-13" or "agrees to the
last digit" in the plan and in the commit messages means exactly that and no
more.

What it does **not** probe:

- **small `Gamma`.** As `Gamma` falls the resolvent approaches its poles and the
  relative error can behave completely differently. Nothing below `0.02` has been
  looked at.
- **`w` between poles.** A uniform grid of 300 points over a band carrying 1000
  poles samples roughly every third pole and never lands in the deep minima
  between them, where `G` is small and cancellation is worst.
- **`w` near the chemical potential.** This is the one that matters physically and
  is the least tested. The valence/conduction split happens there, the Fermi
  weights of the two copies of a level are `f` and `1-f`, and the interesting
  spectral structure sits there. A **logarithmic mesh in `|w - mu|`**, down to
  whatever the model's smallest scale is, is the honest test and has never been
  run.
- **the imaginary axis**, where a DMFT loop actually lives.
- **sum rules and moments** rather than pointwise values.

So the first piece of work is a harness, not a change of algorithm: sweep `Gamma`
over decades, use a log mesh in `|w - mu|` as well as a uniform one, and compare
against `G` evaluated from the definition — a direct sum over the poles of the
*source* object — rather than against another conversion, so that two wrong
answers cannot agree.


## Three things already known that this has to look at

**The chain tail is not determined.** The map from a measure to its recurrence
coefficients is exponentially ill conditioned in the chain index. Two different,
both backward stable, routes to the same chain — the `dsbtrd` band reduction and
the route through the Anderson form — agree to 3e-14 for a chain of 41 sites, part
company at site 31 of 61, and at site 27 of 115 differ by 5 on a scale of 4. Both
still reproduce `G` to 1e-13 *on the grid above*. Whether they still agree on a
log mesh near `mu` at small `Gamma` is exactly the open question.

**Full re-orthogonalization is measurably better on that metric.** On `l → t`,
against the same grid:

| poles | Lanczos, full re-orthogonalization | RKPW |
|---|---|---|
| 500 | 1.8e-14 | 8.0e-14 |
| 1000 | 1.3e-14 | 8.9e-14 |
| 2000 | 1.0e-14 | 8.5e-14 |

a factor of five to seven, for 0.53 s against 0.007 s at 1000 poles and O(N²)
memory against O(N). That is the trade the user should be allowed to make. Plain
Lanczos without re-orthogonalization is 1e-2 at 1000 poles and is not a candidate.

**Three tolerances sit exactly where the log mesh would go**, one of them now
settled.

- **the block degenerate-pole merge — RESOLVED, and the reasoning is worth
  keeping.** `BlockAndersonMatrixToBlockListOfPoles` used to hard-code
  `BlockSize * N * DBL_EPSILON * BandWidth`. At `BlockSize = 10`, a 50-block bath
  and a 60 eV band that is `5100 · 2.2e-16 · 60 ≈ 70 µeV`; at `BlockSize = 2`
  with 100 poles it is still ≈ 3 µeV. Magnetic fields give µeV splittings on
  bands of tens of eV, so that threshold merged away exactly the physics the
  block types exist to resolve. The eigenvalue error it is meant to hide under is
  `O(ε·‖a_k‖)` from *independent* `BlockSize × BlockSize` diagonalizations and
  does **not** grow with the number of poles, so a few·`ε`·bandwidth — about
  `1e-13 eV` — is both sufficient and nine orders below µeV.

  It is now an argument, `MinimalPoleDistance`, threaded from `ChangeType`
  through `ConvertResponsefunction` / `ConvertBlockMatrix` to every conversion
  that merges, and **`BlockListOfPolesMergeDegeneratePoles` is the single place
  that turns an unset (negative) value into the standard**
  `10 · DBL_EPSILON · (BandWidth + DBL_EPSILON)`, with the band width taken from
  the outermost poles. Every caller only forwards what it was handed; do not add
  a second default anywhere. Note `MinimalPoleDistance == 0` merges **nothing**,
  not even two bit-identical poles: the test inside the merge is a strict `<`.

  Merging is not optional: summing the rank-1 weights over a degenerate group
  gives `b·P_D·b^†` with `P_D` the spectral projector, which is basis-independent,
  and that is what makes the result independent of the arbitrary rotation LAPACK
  picks inside a degenerate subspace. So the threshold has to be tight *and*
  nonzero. Where exactly it should sit is still unmeasured — that part is open,
  and it is what the harness above is for.
- `PolesToTriDiagonalCoefficients` merges poles closer than
  `4·DBL_EPSILON·(|E_i|+|E_j|)`. That is a relative tolerance, and a logarithmic
  mesh accumulating near `mu` is precisely where level spacings get small. It is
  needed — the method is not robust against a measure whose support is smaller
  than the number of poles handed to it — but the threshold has never been
  stress-tested against a mesh designed to approach it.
- The Fermi split keeps a level in both chains throughout `±40 kBT`, with weights
  down to about 4e-18 at the window edge. Those tiny-weight poles are kept, they
  are the ones nearest the edge of the window, and what they do to the chain has
  not been looked at.


## The option

`ChangeType` already takes an options list — `SingularValue`, `EffectivekBT`,
`Deflate` — so an `Algorithm` option fits where the others are, in
`LuaChangeTypeResponseFunction` and through `ConvertResponsefunction`.

Candidates worth offering, all of which exist or nearly exist:

| | cost | memory | notes |
|---|---|---|---|
| RKPW | O(N·k) | O(N) | today's default, `PolesToTriDiagonalCoefficients` |
| band reduction | O(N²·kd) | O(N·kd) | `dsbtrd`; used by `d → t`, applies wherever the source is narrow-banded |
| Lanczos, full re-orthogonalization | O(N·k²) | O(N·k) | the `…Old` routines are this, kept for exactly this reason |
| re-run / iterative refinement | ×2 | — | convert, evaluate, convert the residual; untested here |
| higher precision accumulation | small | — | only the inner products, where the cancellation is |

The default should stay the cheapest one that passes the harness. The point of the
option is that someone who needs the last digit, on a small problem, can pay for
it.


## Why it is not simply done

The harness has to come first, or the option is a choice between algorithms whose
relative accuracy nobody has measured where it matters. And the harness needs a
reference that is independent of every conversion, which means evaluating `G` from
the source object's own definition at each probe point — cheap for a list of poles,
and the natural common ground for all the types.

There is also a policy question hiding in it: whether `Algorithm` is a hint the
conversion may ignore when it does not apply (the band reduction only makes sense
for a banded source), or an instruction that must fail loudly when it cannot be
honoured. The second is friendlier to a user who asked for accuracy and needs to
know they did not get it.


## Meanwhile

Read "agrees to 1e-13" in the plan and the commit messages as "agrees to 1e-13 on
a uniform grid at `Gamma` of 0.3 and 0.02".

Note also that the two families do not use the same algorithm.
`DiagonalMatrixBlockTridiagonalize` passes `ReOrthogonalize = true`, so every
`Block…` conversion that goes through a block Lanczos is the
full-re-orthogonalization row of the table above, while the non-block `a → d` is
RKPW. Block and non-block conversions of the same input therefore should *not*
agree to the last digit, and the block one should be the more accurate. Do not
read a 1e-14 disagreement between them as a defect in either. Do not read it as a statement about
small `Gamma`, about `w` near `mu`, or about the coefficients of a long chain,
none of which have been tested.
