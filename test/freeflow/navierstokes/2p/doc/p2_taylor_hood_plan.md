<!--
SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
SPDX-License-Identifier: GPL-3.0-or-later
-->
# Plan: P2–P1–P2–P2 Taylor–Hood CH–NS scheme for the rising-bubble test

**Goal.** Replace the current `PQ1Bubble`(MINI)–P1–P1–P1 monolithic Cahn–Hilliard / Navier–Stokes
scheme with a higher-order **Taylor–Hood** scheme: **P2 velocity, P1 pressure, P2 phase-field `c`,
P2 chemical potential `μ`**, with `c` and `μ` co-located on the momentum balance (they share the
P2 degree). Motivation: the current scheme sits at a **~0.5–1 % accuracy/noise floor** on the
benchmark QoIs (established this session — see [`formulation.md`](formulation.md) §8 and the
ε-convergence plots), and that floor is a P1+bubble + AMR-remesh artefact, not a modelling error.
Aland & Voigt reach clean monotone convergence with exactly this P2–P1–P2–P2 element.

> **Hard constraint.** The current implementation must stay **byte-for-byte runnable and untouched**
> so we can do clean A/B comparisons. The whole plan is built as *new files + new TypeTags* on a
> *new branch*, with the baseline pinned in a **separate git worktree**.

This doc is produced from a multi-agent mapping of the codebase; the two load-bearing facts were
verified by hand (see §2). The automated adversarial critique did not complete, so §7 is my own
review.

---

## 1. Verified starting point

- **Current wiring**: 2 subdomains — `MultiDomainTraits<MOMENTUM, MASS>`.
  - subdomain 0 = **momentum** = `PQ1BubbleModel` (MINI: P1+bubble velocity), `numEq = dim`.
  - subdomain 1 = **mass** = `BoxModel` (P1), `numEq = 3` carrying **`p`, `c`(φ), `μ`** on one basis.
  - Both grid geometries are on the **same leaf grid** → coupling is element-local.
  - `c,μ,∇c,∇μ` reach the momentum force/interface-flux through `couplingmanager.hh`
    (`gradients()/values()` evaluate the mass P1 elemSol at momentum integration points).
- **A CVFE subdomain uses ONE basis for all its equations** → **mixed polynomial orders require
  separate subdomains**. This single fact drives the whole design.
- **The P2 infrastructure already exists locally**: `dumux/discretization/pq2.hh`
  (`PQ2HybridModel`, `PQ2FVGridGeometry`, hybrid caches, `PQ2QuadratureTraits`), and it is **already
  used as a Taylor–Hood momentum space** in `test/freeflow/navierstokes/donea`
  (`DoneaTestMomentumPQ2Hybrid = <DoneaTest, NavierStokesMomentumCVFE, PQ2HybridModel>`, with the
  required `Experimental::GridVariables` + `HybridCVFEGridVariablesCache` specialisations). **This is
  the copy-paste template.** *(verified: `properties_momentum.hh:81`.)*
- **The momentum 2p model is discretization-agnostic** (`NavierStokesMomentumTwoPCVFE` only overrides
  `LocalResidual`+`CouplingManager`), so it already composes with `PQ2HybridModel` the same way it
  does with `PQ1BubbleModel`.
- **Branch state**: `feature/ff-2p` is **35 ahead / 180 behind `origin/master`** (merge-base
  `e271b58a`). Master has since refactored the whole `pq2`/`cvfe/hybrid`/`fem` stack (new
  `pq2/dofhelper.hh`; the `pq2hierarchical*.hh` basis files relocated; a `boundaryFaces()` API rename
  `intersectionIndex→boundaryFaceIndex`; the momentum `felocalresidual.hh`/`localresidual.hh` slimmed
  by ~300 lines). *(verified: `git rev-list` counts + merge-base.)* **So the P2 stack must be pulled
  from master before new P2 code is written.**

---

## 2. Target architecture (recommended: Route B, monolithic)

Two subdomains, monolithic Newton (as today):

- **subdomain 0 = coupled momentum + Cahn–Hilliard**, on `PQ2HybridModel`, **`numEq = dim + 2`**:
  velocity (`dim`), `c` (idx `dim`), `μ` (idx `dim+1`) — all on one P2 hierarchical basis.
  The entire capillary coupling (`−c∇μ` well-balanced force, `c`-advection by `v`, CH diffusion
  `−M(c)∇μ`, the `μ`-definition `μ − γε∇c…`) assembles **element-locally with no coupling manager**.
- **subdomain 1 = pressure/continuity**, on `Box` (P1), **`numEq = 1`**.
- **P2 velocity / P1 pressure = inf-sup stable Taylor–Hood** (no pressure stabilization).
- The only remaining inter-domain coupling is **pressure ↔ momentum** (velocity → continuity `∇·v`;
  `∇p` → momentum) plus continuity reading `c` from the momentum domain for `ρ(c)`.
- **Assembly is hybrid**: vertex (control-volume) dofs use box-style FV storage/flux/source; edge
  (non-CV, P2) dofs use a Galerkin FE weak form over element quadrature — momentum FE terms already
  exist in `felocalresidual.hh`; the **CH FE terms are new**.
- **Monolithic coupling is strictly stronger than Aland & Voigt's operator split** and converges to
  the same benchmark QoIs.

Why Route B over the alternatives (this is a **decision to confirm** — see §8):
- **Route A2** (separate P2 CH subdomain → 3 subdomains) keeps every model clean but needs a
  **3-domain free-flow coupling manager that does not exist** — larger new infrastructure.
- **Operator split** (solve CH then NS each step, A&V-faithful) — DuMux has **no split driver**, and
  it is weaker coupling.
- **Route B** matches the task hint, keeps 2 subdomains, and **deletes the momentum↔(c,μ) coupling**;
  its cost is momentum-model surgery (`numEq: dim→dim+2`, `velocity()` returns only the first `dim`
  components) — all in **new forked files**, so the current impl is untouched.

---

## 3. ⚠️ Recommended pragmatic first target: P2-velocity-only (P2b–P1–P1–P1)

**Before** committing to the XL full-P2 surgery, do the smallest change that tests the core
hypothesis "*P2 velocity beats the noise floor*":

- Swap **only** the momentum space `PQ1BubbleModel → PQ2HybridModel` (copy donea's wiring). Keep
  `c,μ,p` exactly as today (Box-P1 mass, existing coupling manager, existing AMR).
- This is essentially a **`properties.hh` change** + the donea `GridVariables`/`GridGeometry`
  specialisations. No new models, no `numEq` surgery, no CH FE weak form.
- It yields **real benchmark numbers** and answers: does higher-order *velocity* alone (the main
  accuracy lever for rise/deformation, per A&V) close the gap? If yes, we may never need P2 `c,μ`.
- **AMR advantage**: the delicate *conservation-preserving* transfer stays on **P1 `c`** (works
  today); only the **velocity** gains edge dofs, and velocity is non-conserved → the existing
  bubble-dof L2-fit transfer pattern extends to P2 velocity edge dofs much more easily than a
  conservative P2-`c` transfer would. This sidesteps the hardest part (Phase 6) for a first result.

**Recommendation:** run P2-velocity-only (§Phase 2a below) and only proceed to full P2–P1–P2–P2 if it
proves insufficient. The user's stated end-goal is the full co-located scheme, so the phases below
build to it — but this intermediate de-risks the whole effort at ~1/5 the cost.

---

## 4. Phased plan

Ordered by dependency, smallest-viable-first. Effort: **S**mall / **M**edium / **L**arge / **XL**.

### Phase 0 — Baseline isolation + master reconciliation of the hybrid stack — **L**
Freeze the current impl as an immovable comparison baseline; put master's refactored PQ2/hybrid/FEM
under the feet of all P2 work without editing current-impl source in place.
1. `git tag ff2p-baseline HEAD`; `git worktree add ../ff2p-baseline ff2p-baseline` — run all baseline
   PQ1Bubble comparisons from that worktree (current source literally never edited).
2. Create `feature/ff-2p-p2` from HEAD for all new work.
3. Enumerate the refactor: `git log e271b58a..origin/master -- dumux/discretization/pq2
   dumux/discretization/cvfe/hybrid dumux/discretization/fem dumux/freeflow/navierstokes/momentum/cvfe`.
4. **Rebase `feature/ff-2p-p2` onto `origin/master`** (preferred) *or* cherry-pick just the
   pq2/hybrid/fem commits; re-apply the 2p physics + working-tree edits, resolving the
   `boundaryFaces` rename and the `felocalresidual` slimming.
5. Build+run `test/discretization/pq2` and the donea `PQ2Hybrid` test to confirm the P2 stack is
   healthy post-rebase.
6. **Build+run the current 2p test from BOTH the baseline worktree and the rebased branch; confirm
   identical settled Case-1 output** (proves the master pull did not perturb PQ1Bubble physics — it
   should not, the FE path is gated by `hasNonCVLocalDofsInterface`).

_Risk_: 180-commit conflict surface (slimmed momentum residual + `boundaryFaces` API). Mitigation:
baseline lives in its own worktree, so comparison capability is never at risk; reconcile *before* any
new code exists.

### Phase 1 — Single-phase P2/P1 Taylor–Hood NS sanity on the bubble geometry — **S**
Prove PQ2 velocity + Box-P1 pressure works for *this* geometry/BCs/gravity/time-stepping **before**
any Cahn–Hilliard.
- New parallel dir (or reuse a channel) with `…MomentumPQ2 = <base, NavierStokesMomentumCVFE,
  PQ2HybridModel>` + `…MassBox` (P1); copy donea's PQ2 wiring verbatim (incl. `BoundarySegmentIndexFlag`
  so gmsh flags survive, and the `HybridCVFEGridVariablesCache` + `Experimental::GridVariables`
  specialisation).
- Set ICs **nodally at dof positions** (PQ2 `localInterpolation` throws — no `interpolate()`).
- Run Stokes then NS single-phase; verify convergence order and stable pressure (no inf-sup failure).

### Phase 2a — (recommended intermediate) P2-velocity-only two-phase — **S–M**
Full §3: momentum `→ PQ2Hybrid`, keep P1 `c,μ,p` mass + existing coupling manager. Run Case 1 to t=3,
compare to baseline + Hysing. Decide whether full P2 `c,μ` is warranted. *(New files/TypeTags only.)*

### Phase 2 — New co-located momentum+CH P2 model (`numEq = dim+2`) — **XL**
The core new development (only if Phase 2a shows P2 `c,μ` is needed).
- New model dir `dumux/freeflow/navierstokes/momentum/chns/cvfe/`: `numEq()=dim+2`, new `Indices`
  appending `phaseFieldIdx=dim`, `chemicalPotentialIdx=dim+1`.
- New volume variables (fork of `momentum/cvfe/volumevariables.hh`): `velocity()` returns only the
  first `dim` components; add `phaseField()/chemicalPotential()`.
- New local residual (fork of `momentum/2p/cvfe/localresidual.hh`): restrict momentum to `dim`
  equations; **internalize** the interface terms (Korteweg `−c∇μ`, AGG flux) reading `c,μ` from the
  **same** `elemVolVars` instead of `couplingManager`.
- New **CH FE weak-form helper** (mirror `felocalresidual.hh`): per non-CV dof, integrate over
  `CVFE::quadratureRule(element)` the `c`-transport (`∫N ∂ₜc` + advection `u·∇c` + CH diffusion
  `M(c)∇μ·∇N`) and the `μ`-definition (`μN − γε∇c·∇N − f'(c)N`); box-CV form on vertices as `mass/2p`
  does today. Wire via `addToElementStorage/FluxAndSourceResidual`.
- **Raise PQ2 quadrature order** (custom `PQ2QuadratureTraits`) for the nonlinear `c`-products and the
  degenerate mobility `M(c)=M₀(c²−1)²` so the P2 CH terms are not under-integrated.

_Risk_: breaks the `numEq==dim` invariant and mixes FE-Galerkin (momentum) with CV/box (CH) in one
residual. **Mass-lumped storage and full-upwind `c`-advection are P1/MINI constructs** and may degrade
formal P2 accuracy on edge dofs — re-examine (see §7).

### Phase 3 — Pressure-only mass subdomain + minimal pressure↔momentum coupling — **M**
- New mass model `mass/chns_pressure/` with `numEq=1`: continuity `∇·v` (v via `faceVelocity`) and
  `ρ(c)` with `c` read from the momentum elemSol at the shared element.
- New coupling manager (fork of the test + `multidomain/freeflow/couplingmanager_cvfe.hh`): momentum
  reads `p`; mass reads `v` **and** `c` (new direction: continuity-reads-`c`-from-momentum), via
  `elementSolution + evalSolutionAtLocalPos` (basis-agnostic; P1 pressure ↔ P2 velocity is cross-basis).
- Audit `couplingmanager_cvfe.hh` for `scvs()`/`shapeValues[scv.localDofIndex()]` spots (e.g. density
  interpolation ~line 295) and switch to `localDofs()` loops where the coupled field is hybrid P2.

### Phase 4 — Assemble the full monolithic P2–P1–P2–P2 test — **L**
- New test dir `test/.../2p/p2/`: `MultiDomainTraits<MomentumCHNS(PQ2, dim+2), MassPressure(Box, 1)>`,
  new coupling manager, PQ2 grid-geometry + `GridVariables` specialisations.
- New `main.cc`: both geometries on the same leaf grid; `x[momentumIdx]` sized `numDofs*(dim+2)`; set
  analytic `c,μ` ICs **nodally**; `MultiDomainFVAssembler` + Newton + UMFPack.
- Re-point diagnostics (rise velocity, com, mass-balance) to read `c,μ` from `x[momentumIdx]`.
- **Statically refine an interface band** (A&V `h_int=h/8`) so P2 runs **without adaptive transfer**.
- Run Case 1 to first physical solution; compare qualitatively to the baseline worktree.

_Risk_: first full run; the capillary force form (`−c∇μ` well-balanced) must stay consistent or the
**parasitic interfacial jet returns**. `numEq=dim+2` on P2 substantially increases Jacobian bandwidth.

### Phase 5 — Validation vs baseline + Hysing + ε-convergence — **M**
- Regenerate PQ1Bubble reference quantities from the baseline worktree.
- Run P2 to t=3; compare rise velocity / circularity / com to Hysing **and** to the baseline.
- Run the ε-convergence sequence on P2; confirm it reaches the same sharp-interface limit **at coarser
  meshes** (the P2 payoff) and — the real test — with **less remesh jaggedness**.
- Tabulate cost (dofs, Jacobian nnz, assembly time, Newton iters) P2 vs PQ1Bubble.
- Store P2 references in a **new** subdir; never overwrite the PQ1Bubble references.

### Phase 6 (deferred/optional) — P2-capable adaptive mesh refinement — **XL**
The single largest AMR unknown. See §5. Do P2 on a static band first; only build this once the P2
physics is validated.

---

## 5. AMR: the hard, deferred part (§Phase 6)

The current custom AMR (`adapt.hh`) rests on an assumption **P2 breaks**: the only persistent dofs are
P1 vertices and *dof value == nodal value*. PQ2 adds **edge dofs** that are:
1. **non-persistent** under refinement,
2. **non-control-volume** — `scvs()` covers only vertices, so scvs-based transfer loops **silently
   skip them** (no compile error), and
3. **hierarchical** — an edge coefficient is `u(midpoint) − ½(u_a+u_b)`, **not** the field value, so
   the existing `sol[dof]=evalSolution(dofPosition)` pattern **corrupts every edge dof silently**.

A correct P2 transfer must: iterate `localDofs()` not `scvs()`; for new refined edge dofs set
`coeff = evalSolution(midpoint) − ½(evalSolution(A)+evalSolution(B))` using reference-element
endpoints on the same source element; for coarsening, reconstruct non-persistent edge dofs of non-leaf
elements via a **local L2 / mass-matrix fit**. **`c` is more sensitive than velocity**: it enters the
conserved CH balance and is not re-solved from scratch, so it needs a **consistent mass-matrix
projection** (not a diagonal L2 fit, not the lumped-scvVolume trick) to avoid integral-`c` drift per
remesh. The **indicator itself is basis-agnostic** (φ-jump via `evalSolution` at element centres) —
just re-point it to read `c` from the momentum PQ2 solution. Keep **`Dune::conforming` (bisection)**
refinement so the P2 space stays conforming (no hanging-edge constraints).

> This is why **Phase 2a (P2-velocity-only, P1 `c`)** is attractive: the conservative transfer stays
> on P1 `c` (works today), and only velocity gains edge dofs (non-conserved → tolerant L2 fit).

---

## 6. What "keep the current impl untouched" concretely means

- **Baseline worktree** pinned to tag `ff2p-baseline` = the comparison oracle; its source is never
  edited even though the P2 branch upgrades the shared hybrid infra.
- All P2 work = **new files, new TypeTags, new test dir** on `feature/ff-2p-p2`.
- The one shared-code change is the **master rebase** of `dumux/discretization/{pq2,cvfe/hybrid,fem}`;
  validated harmless to PQ1Bubble by re-running the current test in both trees (Phase 0 step 6).

---

## 6b. Refined architecture (from a line-level read of the residual/flux stack, 2026-07-04)

Locked design decisions for the co-located momentum+CH model, more precise than §2–§5:

1. **`velocity()` must return a `dim`-vector, not the full `dim+2` priVars.** `flux.hh:189` hard-asserts
   `NumEqVector::dimension == dimWorld`, and the flux context's `velocity()` is a `GlobalPosition`
   (`FieldVector<Scalar,dim>`). So the new volvars store `priVars_` (size `dim+2`) but expose
   `velocity()` = first `dim` components, plus `phaseField()=priVars_[dim]`,
   `chemicalPotential()=priVars_[dim+1]`.
2. **Instantiate the momentum flux/FE helpers with a `dim`-vector**, not the `dim+2` `NumEqVector`:
   `using MomVec = Dune::FieldVector<Scalar,dim>; NavierStokesMomentumFluxCVFE<GG, MomVec>`. They
   then produce `dim`-vector momentum contributions; the new residual **embeds** those into
   `residual[·][0..dim-1]` and adds CH scalars to `[dim]` (φ-transport) and `[dim+1]` (μ-def).
   The momentum FE helper (`felocalresidual.hh`) already loops `eqIdx < NumEqVector::dimension`;
   instantiated with `MomVec` it touches only `0..dim-1` of the real `dim+2` residual — correct.
3. **Density/viscosity MUST come from the (deflected) volvars, not `problem.density(ipData)`.** The
   monolithic numeric Jacobian deflects φ per dof and rebuilds the elemVolVars; only
   `volVars.density()=mixtureDensity(volVars.phaseField())` responds to that deflection.
   `problem.density(element,fvGeometry,ipData)` reads a *stored* (non-deflected) coupled solution →
   would give a **wrong Jacobian** w.r.t. φ. Therefore the chns momentum flux/FE/residual must be
   **forked** to read ρ(φ), ν(φ) from `volVars`, and φ has **no coupling manager** (it is local).
   Only pressure↔momentum coupling remains.
4. **CH terms port cleanly from `mass/2p/localresidual.hh`** with velocity taken locally (no
   `problem.faceVelocity` coupling): storage `[φ]=c, [μ]=0`; flux `[φ]=upwind(c)·volFlux − M(c)∇μ·n`,
   `[μ]=−σ̃∇c·n`; source `[μ]=μ − (σ̃/ε)f'(c)`. The φ-advection volume flux uses the local
   interpolated velocity `v·n·dA`. The skew-symmetric φ-advection and momentum-advection corrections
   port as-is.
5. **Continuity (`∇·v`) is NOT in this model** — it lives in the separate P1 pressure subdomain, which
   reads `v` and `c` (for ρ(c)) from the momentum+CH domain.

**File-by-file build order for Phase 2** (all new, under `dumux/freeflow/navierstokes/momentum/chns/cvfe/`):
`indices.hh` (done) → `volumevariables.hh` (velocity split + φ/μ + ρ/ν from problem material laws) →
`flux_chns.hh` (fork of `flux.hh`, ρ/ν from volvars) → `felocalresidual_chns.hh` (momentum FE via
`MomVec` + CH FE weak form) → `localresidual.hh` (embed momentum `dim`-vector + CH scalars; skew
corrections) → `model.hh` (ModelTraits `numEq=dim+2`, TypeTag, property wiring).

---

## 6c. Implementation status (2026-07-04) — model DONE, harness needs the "new interface"

**Written and compiling into the assembler** (`dumux/freeflow/navierstokes/momentum/chns/cvfe/`):
`indices.hh`, `volumevariables.hh`, `flux.hh` (volvars-ρ/ν fork), `felocalresidual.hh` (CH Galerkin
weak form + Korteweg stress + Boussinesq, sign-matched), `localresidual.hh`, `model.hh`. Plus the test
`properties.hh` (PQ2 momentum-CHNS + 1p-mass pressure + FreeFlowCouplingManager + hybrid GridVariables).
The capillary force is the **Korteweg stress** face-flux form (needs only ∇c, available at faces/qps;
the well-balanced −φ∇μ would need ∇μ at CVs → later refinement). Material laws and the double-well /
buoyancy are computed **inline in the residual** via problem accessors (`mixtureDensity/Viscosity`,
`surfaceTension`=σ̃ε, `energyScale`=σ̃/ε, `referenceDensity`, `mobility`).

**The blocker (plan critique #7, confirmed):** `PQ2HybridModel` requires the **experimental
"new" interface** (`Experimental::GridVariables` + `Experimental::MultiDomainAssembler` +
`Experimental::BoundaryTypes`), which is a *different* problem/main API. The classic
`MultiDomainFVAssembler` + `boundaryTypes(setDirichlet/setNeumann)` + `applyInitialSolution` I first
wrote do **not** compile against it. Template to follow: `test/.../donea/problem_newinterface.hh`
+ the `NEW_PROBLEM_INTERFACE`/`NEW_VARIABLES_INTERFACE` path in `donea/{properties,main}.cc`.

**COMPILE STATUS (2026-07-04): driven 20 → 2 root errors.** The model + properties + new-interface
problem.hh + main.cc (Experimental::MultiDomainAssembler, experimental GridVariables for BOTH
subdomains, constraint-based BCs) all compile; the assembler instantiates the CHNS residual. The two
remaining obstacles are **architectural, not typos**, and both were predicted by §7:
- **(A) PQ2 initial condition.** `assembleInitialSolution` (classic) sets only vertex dofs — it misses
  the P2 **edge** dofs — and the experimental CVFE problem base exposes no `applyInitialSolution`. Fix:
  a small IC routine that loops ALL local dofs and sets `sol[dof] = initialAtPos(dofPosition)` (works
  for vertex+edge), or find/borrow the experimental applyInitialSolution used by donea's GridVariables.
- **(B) THE Route-B obstacle — `FreeFlowCouplingManager` assumes momentum `numEq == dim`.** It reads the
  momentum solution AS the velocity (`VelocityVector = GlobalPosition`, dim=2) at **4 sites** (`velocity`,
  `faceVelocity`, the axpy interpolation, `evalSolutionAtLocalPos`), but our momentum solution is
  `FieldVector<dim+2>`. → **must fork/generalize the coupling manager** to slice the first `dim`
  components as the velocity (a no-op for existing numEq==dim tests, so a generic slice is safe). Make a
  `chns` coupling manager (fork of `dumux/multidomain/freeflow/couplingmanager_cvfe.hh`) or add a
  velocity-slice that reads `volVars.velocity()` (already dim-sized) instead of the raw solution.
  This is the single remaining substantial fork; everything else is in place.

**Remaining harness work (well-defined):**
1. `main.cc`: switch to `Experimental::MultiDomainAssembler<Traits, CM, numeric>` (has the transient
   `(problems, ggs, gvs, cm, timeLoop, prevSol)` ctor); initial solution is applied by the
   experimental grid-variables `init(x)` via `problem.applyInitialSolution` — provide `initialAtPos`
   and let the base handle it (drop the direct `applyInitialSolution` calls).
2. `problem.hh`: rewrite to the new interface —
   - `boundaryTypesAtPos` returns `setAllFluxBoundary()` on Neumann faces (free-slip y-momentum, all
     of the pressure domain) and default elsewhere;
   - build Dirichlet **constraints** in the ctor (`appendDirichletConstraints_`-style) for the
     pinned dofs: no-slip velocity on top/bottom (all velocity comps), no-penetration `velocityX` on
     the side walls (**per-component** `ConstraintInfo`, not `setAll` — the free-slip walls are mixed),
     and the single pressure-anchor dof (internal constraint on the Box mass domain);
   - `c`, `μ` get **no** constraints and **no** flux BC (natural/zero-flux) — they are interior-like
     everywhere;
   - `boundaryFlux(fvGeometry, elemVars, faceIpData)` returns 0 for the closed box (no through-flow);
   - keep `sourceAtPos` = 0 (all physics is inline in the residual) and the material/CH accessors +
     `initialAtPos` (tanh c + Gibbs-Thomson μ, velocity 0 / pressure 0).
3. Then compile-iterate (still expect: the hybrid assembler may want `storageIntegral`/`fluxIntegral`/
   `sourceIntegral` on the residual — add them mirroring the base momentum residual; and check the
   `FreeFlowCouplingManager` extracts velocity via `volVars.velocity()` given numEq=dim+2).

---

## 7. Critical review (the adversarial pass the workflow could not finish)

1. **Biggest lever vs biggest cost are different targets.** A&V themselves note P2 for `c,μ` is a
   *convenience* (equal basis with `u`), not a necessity. The accuracy win almost certainly comes from
   **P2 velocity**. Doing the XL Phase 2 (P2 `c,μ`) before proving Phase 2a (P2 velocity only) risks
   spending 80 % of the effort on the last 10 % of the benefit. **Sequence 2a before 2.**
2. **The rebase (Phase 0) is the critical path and the real schedule risk**, not the physics. 180
   commits, a public API rename, and a 300-line residual slimming, against 2p physics written for the
   frozen API. Budget real time here; consider cherry-pick over full rebase if conflicts explode.
3. **Monolithic `numEq=dim+2` P2 conditioning/cost.** UMFPack direct solve on a P2 velocity + 2 scalar
   fields per node **plus edge dofs** will be much heavier than PQ1Bubble; the coupled `c–μ–v` block
   may also condition worse. Expect to need an iterative solver + block preconditioner before this
   scales to ε-refinement — not in the current stack. Flag as a likely Phase 5 blocker.
4. **First-order constructs leaking into P2.** Mass-lumped CH storage and full-upwind `c`-advection are
   P1/MINI tricks; on P2 edge dofs they can *destroy* the formal order you paid for. The P2 CH weak
   form should use consistent (unlumped) mass and a higher-order-compatible advection (central + the
   existing skew-symmetric stabilization) — otherwise P2 may not beat P1+bubble and the whole exercise
   is moot.
5. **The parasitic-jet risk resurfaces.** The well-balanced `−c∇μ` force and the AGG-flux consistency
   (this session's Model-1 finding) must be re-derived for the co-located P2 assembly; a sign/edge-dof
   inconsistency reintroduces the interfacial jet with no obvious error.
6. **"Clean comparison" is only clean if Phase 0 step 6 passes.** If the master rebase changes the
   settled Case-1 numbers at all, every downstream P2 vs baseline comparison is confounded. This gate
   is non-negotiable.
7. **Coupling-manager, momentum-residual, and mass-CH reader agents did not complete** (session
   limits), so Phases 2–3 rest partly on the rich subdomain finding rather than a line-by-line read of
   `couplingmanager_cvfe.hh` and the momentum FE residual. Re-read those two files before starting
   Phase 2/3.

---

## 8. Decisions I need from you before Phase 0

1. **Route B vs Route A2 vs P2-velocity-first.** Full co-located P2–P1–P2–P2 (your stated goal, XL), or
   start with **P2-velocity-only** (§3, S–M) and escalate only if needed? *(My recommendation: 2a
   first.)*
2. **Must `c,μ` genuinely be P2?** If P1 `c,μ` is acceptable, a **P2(v)–P1(p)–P1(c)–P1(μ)** scheme
   stays entirely within the working stack and touches almost only `properties.hh`.
3. **Rebase vs cherry-pick** for the 180-commit gap.
4. **Do your uncommitted working-tree edits** (`felocalresidual.hh`, `cvfe/localresidual.hh`,
   `momentum/2p/cvfe/localresidual.hh`) carry physics to preserve, or are they superseded by master's
   refactor of the same files?
5. **Goal = A&V element fidelity, or beating the benchmark noise floor?** The latter is likely reachable
   with less than the full P2–P1–P2–P2 and relaxes several constraints.

---

## Appendix — key files

| Concern | File |
|---|---|
| Subdomain/TypeTag wiring (edit target) | `test/.../2p/properties.hh`, `main.cc` |
| **Taylor-Hood template to copy** | `test/freeflow/navierstokes/donea/properties_momentum.hh` |
| P2 infra | `dumux/discretization/pq2.hh`, `dumux/discretization/pq2/*`, `dumux/discretization/cvfe/hybrid/*` |
| Momentum model (discretization-agnostic) | `dumux/freeflow/navierstokes/momentum/2p/cvfe/model.hh` |
| Momentum FE (edge-dof) residual (mirror for CH) | `dumux/freeflow/navierstokes/momentum/cvfe/felocalresidual.hh` |
| CH physics to move/fork | `dumux/freeflow/navierstokes/mass/2p/localresidual.hh` |
| Coupling manager to fork | `test/.../2p/couplingmanager.hh`, `dumux/multidomain/freeflow/couplingmanager_cvfe.hh` |
| AMR to re-point (Phase 6) | `test/.../2p/adapt.hh` |
