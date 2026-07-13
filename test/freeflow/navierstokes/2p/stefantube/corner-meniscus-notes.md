# Corner-meniscus notes (Stefan-tube 3D benchmark)

Working notes on the square-tube corner physics, kept alongside
`concus-finn-1969-on-the-behavior-of-a-capillary-surface-in-a-wedge_.pdf`
(Concus, P. and Finn, R., "On the behavior of a capillary surface in a
wedge," PNAS 63(2):292-299, 1969, doi:10.1073/pnas.63.2.292).

## Two different corner solutions -- don't conflate them

**1. The bulging cap (what `interfacePositionAt` builds).** A meniscus
surface `x = xi(y,z)` that *terminates* the liquid axially, of constant
mean curvature, meeting all four tube walls at contact angle theta.
Concus-Finn's Theorem (their eq. 4a,b and the existence statement in
Section 2): a bounded solution of this type exists near a wedge vertex of
half-angle alpha **iff alpha + theta >= pi/2**. For our square corner,
alpha = pi/4, so existence needs theta >= pi/4 (45 deg). Below that, no
bounded terminating cap exists near the vertex -- volume grows like O(r),
height like O(1/r) approaching it (their Sections 4, 7-8, with `kappa=0`
gravity-free specialization in Section 2 for our case).

**2. The corner filament (this note).** A *translation-invariant* rivulet
of liquid sitting in the corner, constant circular-arc cross-section
extruded along the corner's axial direction (a portion of a circular
cylinder, not a cap). This is the local, gravity-free, axially-symmetric
constant-mean-curvature solution for an infinite straight corner. It is
the closed form that *does* exist unconditionally (for any r), and is a
genuinely different object from (1): it does not terminate the liquid, it
runs along the corner.

For a corner of half-angle alpha, contact angle theta, transverse
arc radius `r = gamma/p_c` (gamma = surface tension, p_c = capillary/
Laplace pressure of the filament):

```
A(r) = r^2 [ cos(theta) cos(theta+alpha)/sin(alpha) - (pi/2 - alpha - theta) ]   # wetting-phase cross-sectional area
b    = r cos(theta+alpha)/sin(alpha)                                             # apex-to-contact-line length, each wall
```

the free-surface arc subtends angle `2(pi/2 - alpha - theta)`. Existence
of a genuine (non-degenerate, A>0) filament again requires the Concus-Finn
condition `theta + alpha < pi/2`, i.e. **theta < pi/4 (45 deg) for the
square duct** -- same threshold as (1), which makes sense: both are
local corner constructions gated by the same wedge condition, just
different global topologies (cap vs. filament).

Specializing to alpha=pi/4, theta=0 (perfect wetting, square corner):
`A = r^2 (1 - pi/4)`.

## How this actually resolves for the Stefan-tube benchmark

Concus-Finn's "no bounded solution" result is a statement about an
**infinite** wedge with **no other physics**. Neither holds exactly here:

- **The corner wall is finite** (length L, only for x in [0,L]): a filament
  advancing along the corner is capped by the physical end of the wall at
  the opening, not by any capillary equilibrium. That's a geometric bound,
  not a solution to (4a,b).
- **Evaporation is present and destabilizing to an advancing filament.**
  A filament that has advanced far toward the opening sits in a region of
  short remaining vapor diffusion path (small local `ell`); since the
  evaporation flux scales like `~1/ell` (the same tube law used elsewhere
  in this benchmark), a far-advanced thin filament evaporates fast and is
  pulled back. So there is a genuine dynamic balance -- wicking advance
  (Cahn-Hilliard/wetting driven) vs. evaporative recession (vapor
  diffusion driven, accelerating as the filament thins and its remaining
  gas column shortens) -- rather than either runaway growth or a static
  equilibrium.

This is consistent with what the coupled (evaporation-ON) 3D runs actually
show: `vapor_l2_rel` plateaus at a persistent nonzero level (~2%) rather
than decaying to ~0 (as the plain extruded-arc model gives) or diverging.
That's the fingerprint of a bounded dynamic balance, not the unbounded
static Concus-Finn limit -- which only manifests when evaporation is
artificially disabled (`Problem.EnableCornerRelaxation=true` in `main.cc`
forces evaporation off specifically to isolate and confirm the static,
theorem-predicted unbounded-advance tendency as a diagnostic; it is not
meant to produce a physical initial condition, see the comment block in
`main.cc`).

## Open modeling question

None of this is yet reflected in `interfacePositionAt`'s initial condition
(still the plain extruded spherical cap, flat-clamped beyond its inscribed
circle). A genuinely corner-consistent IC would need to blend the bulk cap
with a filament construction (1)+(2) above, itself terminated at some
evaporation-vs-wicking balance point -- not a simple closed form, and not
yet attempted. The current flat clamp remains a numerical placeholder; the
`EnableCornerRelaxation` diagnostic and the coupled evaporation-ON runs are
the two available ways to probe the real behavior in the meantime.
