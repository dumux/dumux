<!--
SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
SPDX-License-Identifier: GPL-3.0-or-later
-->
# Diffuse-interface CH–NS rising bubble: model formulation

This test reproduces the **Hysing et al. (2009)** two-dimensional rising-bubble benchmark with
a **diffuse-interface Cahn–Hilliard / Navier–Stokes (CH-NS)** model, discretized with the
DuMux `PQ1Bubble` control-volume finite-element (CVFE) velocity space and a monolithic
multidomain coupling (momentum ⊕ Cahn–Hilliard mass).

The two reference papers live next to this doc:
- `Benchmark_computations_of_diffuse_interf.pdf` — **Aland & Voigt, IJNMF (2011/2012)**: the
  diffuse-interface benchmark of *exactly this problem*; defines the three models below.
- `AGG.pdf` — **Abels, Garcke & Grün, M3AS (2012)**: the original thermodynamically consistent
  variable-density model (Aland–Voigt's "Model 3").

> **TL;DR of the formulation choice.** We solve **Aland–Voigt Model 1** — the convective-form
> Navier–Stokes–Cahn–Hilliard system with a solenoidal (volume-averaged) velocity `∇·u = 0`.
> This matches our energy-stable convective/skew-symmetric momentum machinery, and Aland & Voigt
> show it gives results **≈ identical** to the thermodynamically consistent AGG model (Model 3).
> The AGG "interface flux" is therefore **OFF** — it belongs only to the conservative Model 3
> (see "The AGG interface flux" below).

---

## 1. The benchmark (sharp interface)

Domain `Ω = [0,1]×[0,2]` (we exploit the left–right symmetry and solve the **left half**
`[0,0.5]×[0,2]` with a symmetry plane at `x=0.5` — the benchmark's free-slip side wall *is* the
symmetry condition). A circular bubble of radius `R=0.25` is centred at `(0.5,0.5)` and rises
under gravity `g=0.98`. No-slip top/bottom, free-slip/symmetry sides.

| | ρ (bubble / outer) | μ (bubble / outer) | σ | Re | regime |
|---|---|---|---|---|---|
| **Case 1** | 100 / 1000 | 1 / 10 | 24.5 | 35 | mild, bubble stays ellipsoidal |
| **Case 2** | 1 / 1000 | 0.1 / 10 | 1.96 | 35 | violent, thin skirt / possible break-up |

**Benchmark QoIs** (Hysing §, Aland–Voigt §5.2), all from the bubble region `c<0`:
- center of mass `y_c = ∫_{c<0} y dx / ∫_{c<0} dx`
- rise velocity `V_c = ∫_{c<0} u_y dx / ∫_{c<0} dx`
- circularity `c = 2√(πA)/P` (perimeter of area-equivalent circle / actual perimeter)

Case-1 reference targets (finest curves): **peak V_c = 0.2417 @ t≈0.92**,
**min circ = 0.9013 @ t≈1.9**, **y_c(t=3) = 1.0813**.
Case 2 reference groups **disagree** (TP2D does not break up, MooNMD does), so they bracket
rather than pin the truth.

Reference time series are stored (downsampled, with provenance headers) in
[`../references/`](../references) and read directly by the plot script.

---

## 2. The diffuse-interface model (Aland & Voigt)

A phase field `c ∈ [-1,1]` (`c≈+1` outer fluid, `c≈-1` bubble; **in the code the sign
convention is `φ=+1` inside the bubble** — see note in §4) replaces the sharp interface by a
smooth layer of thickness `ε`. Density and viscosity are affine in `c`:
`ρ(c)=ρ₁(c+1)/2+ρ₂(c−1)/2`, likewise `ν(c)`. The chemical potential and its scaled
surface tension are
```
μ = σ̃ ε⁻¹ W'(c) − σ̃ ε Δc ,   W(c)=¼(c²−1)² ,   σ̃ = σ · 3/(2√2).
```
Aland & Voigt present **three** variable-density Navier–Stokes–Cahn–Hilliard models. All three
share the Cahn–Hilliard transport `∂_t c + u·∇c = ∇·(M(c)∇μ)` and a **solenoidal velocity**
`∇·u = 0`; they differ only in the momentum equation:

### Model 1 — convective form (eq. 8–11) ← **what we solve**
```
ρ(c) (∂_t u + (u·∇)u) = −∇p̃ + ∇·(ν(c) D(u)) + F + μ∇c
∇·u = 0
```

### Model 2 — Boyer (eq. 12–16)
Adds `(ρ₂−ρ₁)(1−c²)/4 · ∇(μ/ρ)` and uses a `ρ`-weighted CH flux. *Differs* from 1 & 3;
no known energy law. We do **not** use it.

### Model 3 — Abels–Garcke–Grün, thermodynamically consistent (eq. 17–21)
```
∂_t(ρ(c) u) + ∇·(ρ(c) u⊗u) = −∇p̃ + ∇·(ν(c) D(u)) + F + μ∇c − ∇·( u ⊗ (ρ₁−ρ₂)/2 · M(c)∇μ )
∇·u = 0
```
i.e. **conservative** storage `∂_t(ρu)` + **conservative** advection `∇·(ρu⊗u)` + an **extra
diffusive-momentum flux** `−∇·(u⊗J)`, `J = (ρ₁−ρ₂)/2 · M∇μ`.

**Key fact (why Model 1 is legitimate):** with `ρ` affine in `c`, the Cahn–Hilliard equation
and `∇·u=0` *imply* the mass balance `∂_t ρ + ∇·(ρu + J) = 0`. This is **derived, not imposed**.
Expanding Model 3's three conservative terms and substituting that identity telescopes them
exactly into Model 1's convective `ρ(∂_t u + u·∇u)`. So Models 1 and 3 are the *same PDE* to
O(ε), which is why Aland & Voigt find them numerically ≈ identical. **Model 2 (Boyer) is the
genuinely different one.**

---

## 3. Why Model 1 (and not the AGG Model 3) in *our* discretization

The two forms are continuum-equivalent, but their **discrete** behaviour differs:

- **Model 3 (conservative)** is discretely energy-stable *only if* the discrete scheme also
  satisfies `∂_t ρ + ∇·(ρu + J) = 0` node-wise. Our monolithic CVFE scheme does **not** satisfy
  that compatibility exactly (it enforces `∇·u=0`, and `ρ` is transported by the Cahn–Hilliard
  subsystem on a different stencil). The mismatch injects kinetic energy at the density jump →
  a growing **parasitic interfacial jet** (this is exactly the failure we observed with plain
  conservative advection: `max|u|` ran away 1.4 → 9).

- **Model 1 (convective)** written in **skew-symmetric (Temam) form** is **unconditionally**
  energy-neutral under `∇·u=0`, independent of any discrete mass-balance compatibility. With it
  on, `max|u|` plateaus instead of diverging.

So in our scheme Model 1 is both the *simpler* and the *more robust* choice, and it reproduces
the AGG result. (Aland & Voigt use P2/P1 Taylor–Hood with operator splitting, where Model 3's
compatibility is easier to satisfy — a different trade-off.)

---

## 4. Continuum term → code map

| Continuum term (Model 1) | Code | Notes |
|---|---|---|
| storage `ρ(c) ∂_t u` | `momentum/cvfe/localresidual.hh` + `2p/cvfe/localresidual.hh`, gated by `FreeFlow.ConvectiveFormMomentumStorage=true` | density frozen at current time level (convective storage), *not* `∂_t(ρu)` |
| advection `ρ(c)(u·∇)u` | skew-symmetric correction in `momentum/cvfe/localresidual.hh`, `FreeFlow.SkewSymmetricMomentumAdvection=true` | converts `∇·(ρu⊗u)` → `Σ_f m_f (u_up − u_cv)`, energy-neutral under `∇·u=0` |
| advection face value | `Flux.UpwindWeight=0.5` | **0.5 = central** (energy-neutral); note 0.0 is full *downwind* here, 1.0 full upwind |
| viscous `∇·(ν(c) D(u))` | `diffusiveMomentumFlux` (symmetric `D=½(∇u+∇uᵀ)`) | ν(c) affine, `mixtureViscosity()` |
| force `F + μ∇c` | `problem.hh` `source()`, `SurfaceTensionForm=wellbalanced` → `−φ∇μ` | `−φ∇μ ≡ μ∇φ` up to `∇(μφ)` absorbed into pressure ("well-balanced", ~0 spurious currents); gravity is Boussinesq `−ρ₂ g` |
| `∇·u = 0` | `FreeFlow.IncompressibleContinuity=true` | volume-averaged (solenoidal) velocity — shared by all three models |
| CH transport `∂_t c + u·∇c = ∇·(M∇μ)` | `mass/2p/localresidual.hh` | phase-field advection `PhaseFieldUpwindWeight=0.5` (central) + skew correction |
| `μ = σ̃ε⁻¹W'(c) − σ̃εΔc` | `problem.hh` `source()` (`energyScale_·φ(φ²−1)`) + gradient-energy (`σ̃ε`) | verified: `energyScale_=σ̃/ε`, `σ̃=σ·3/(2√2)` — matches Aland–Voigt exactly |
| mobility `M(c)` | `problem.hh` `mobility()`; degenerate `M₀(c²−1)²` via `FreeFlow.DegenerateMobility` | Case 1: constant **`M=5e-6 (=1e-3·ε)`** — calibrated by a full mobility sweep (§8). Case 2: **degenerate** `M₀=5e-6` |

**Sign convention note.** The code uses `φ=+1` *inside* the bubble (opposite to Aland–Voigt's
`c≈−1` in the bubble). All formulas above are written in the paper's `c` convention; the code's
`mixtureDensity`, `mixtureDensityDerivative=(ρ₁−ρ₂)/2`, and the QoI weights (`φ>0` = bubble) are
consistent with `φ=+1` in the bubble.

---

## 5. The AGG interface flux — the subtle bit

`problem.hh::interfaceFlux()` assembles, per face, `−v (M∇μ·n)(dρ/dφ)`, i.e. summed over a
control volume it is `∇·(v⊗J)` with `J=(dρ/dφ)M∇μ` — **exactly Model 3's extra term**.

It is therefore **only consistent inside the Model 3 package** (conservative storage +
conservative advection). Enabling it *on top of* our Model 1 convective storage + skew advection
**double-counts** the density-transport term (the piece the convective form already handles via
the `∂_tρ` / `∇·(ρu)` telescoping), leaving a spurious `v·∇·J` force localized at the interface.
That error is negligible at Case 1's density ratio 10 (the flux itself changes Case-1 peak V_c by
< 0.1%) but grows at Case 2's ratio 1000.

**Decision:** `Problem.EnableInterfaceFlux = false` (Model 1). To run a **clean Model 3**
instead, set *all three* together:
```
Problem.EnableInterfaceFlux            = true
FreeFlow.ConvectiveFormMomentumStorage = false   # -> conservative storage d_t(rho u)
FreeFlow.SkewSymmetricMomentumAdvection = false   # -> conservative advection div(rho u x u)
```
(expect the parasitic-jet stability issue described in §3 unless the discrete compatibility is
also addressed — e.g. by folding `J` into the skew advecting mass flux `m_f=(ρu+J)·n`, which
would be the energy-stable *and* thermodynamically-consistent hybrid; not yet implemented).

---

## 6. Measurement: sharp vs diffuse QoIs

The benchmark integrates over the **sharp** bubble region `∫_{c<0}`. A diffuse-interface code
can weight either **sharp** (indicator `φ>0`) or **diffuse** (`(1+φ)/2`). The diffuse weighting
biases the peak rise velocity low by ~1–2%. The test prints **both** (`rise_velocity` /
`rise_velocity_sharp`, `centroid_y` / `centroid_y_sharp`); we compare the **sharp** ones to the
reference. Circularity uses the `φ=0` marching-triangles contour (an independent total-variation
perimeter `∫|∇φ|` is printed as a cross-check).

---

## 7. Reproduce

```sh
# run (from the build test dir; writes QoI lines to stdout)
./test_ff_stokes_2p_pq1bubble_simplex params.input        -Problem.Name c1 | tee c1.log   # Case 1
./test_ff_stokes_2p_pq1bubble_simplex params_case2.input  -Problem.Name c2 | tee c2.log   # Case 2

# plot vs reference (self-contained; reads ../references/)
python3 postprocessing/plot_hysing_benchmark.py --case 1 --out case1.png c1.log:Model1
python3 postprocessing/plot_hysing_benchmark.py --case 2 --out case2.png c2.log:Model1
```
See [`postprocessing/plot_hysing_benchmark.py`](../postprocessing/plot_hysing_benchmark.py)
`--help` for details.

---

## 8. Mobility calibration (Case 1)

A full-`t=3` sweep of the constant mobility `M` (all other settings = clean Model 1, `ε=0.005`)
shows every benchmark QoI converging monotonically toward the reference as `M` decreases:

| `M` | peak `V_c` (ref 0.2417) | min circ (ref 0.9013) | `y_c(3)` (ref 1.0813) | Σ\|err\| |
|---|---|---|---|---|
| 5.7e-5 | +1.0 % | +2.0 % | +1.1 % | 4.2 % |
| 2.0e-5 | +0.6 % | +0.8 % | +0.3 % | 1.7 % |
| 1.0e-5 | +0.4 % | +0.5 % | +0.1 % | 1.0 % |
| **5.0e-6** | **+0.1 %** | **+0.5 %** | **+0.0 %** | **0.7 %** |

Mechanism: smaller `M` → less Cahn–Hilliard interfacial relaxation → the bubble is allowed to
deform more (min circularity drops toward the reference) and, being less round, has more form
drag → the peak rise velocity and centroid stop over-shooting. The circularity residual (+0.5 %)
is an **ε-floor** (finite interface thickness): it stops improving below `M~1e-5` and does **not**
over-deform, so it cannot be removed by lowering `M` — only by a smaller-`ε` convergence study.
`M=5e-6 = 1e-3·ε` is exactly the Aland & Voigt mobility scaling, so this is literature-grounded,
not a free fit. Plots: `hysing_case1_mobility_sweep.png`, `hysing_case1.png`.

**Harmonic vs arithmetic mobility averaging.** For a *constant* `M` (Case 1) the face average is
`M` regardless of averaging rule — a no-op. For a *variable/degenerate* `M(c)` (Case 2) a harmonic
face average is dominated by the smaller (bulk-side) value, reducing the effective interfacial
flux — same direction as lowering `M`, but localized to the interface flanks and nonlinear.
