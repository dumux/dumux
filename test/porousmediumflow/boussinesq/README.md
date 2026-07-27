# Boussinesq Convection: Parallel, Adaptive, Dynamically Rebalancing {#benchmark-boussinesq-convection}

## Onset of convection in a gravitationally unstable diffusive boundary layer

**Problem description**

This test simulates density-driven (Rayleigh-Bénard-type) convective fingering in a
porous medium: a solute diffuses down from the top boundary, forming an unstable
diffusive boundary layer that eventually breaks up into descending fingers once the
local Rayleigh number exceeds the onset threshold. The setup follows
Riaz et al. (2006) @cite Riaz2006 (see the PDF in this directory), which studies exactly
this problem — onset time, wavelength selection, and finger growth — as a model for
CO<sub>2</sub> dissolution trapping in saline aquifers.

The problem is implemented in two independent formulations that solve for the same
physics and are cross-checked against each other:

- **`pressure/`** — the standard @ref OnePNCModel with `BoussinesqCVFEDarcyLaw`, using a
  pressure/mass-fraction primary variable set. Equation 0 is replaced by the total
  (divergence-free) mass balance; the pressure level, otherwise undetermined for a
  pure-Neumann/pure-density-driven problem, is fixed weakly by a Robin penalty
  (Nitsche-type) boundary condition at the top.
- **`vorticity/`** — a custom quasi-static streamfunction/vorticity model
  (`vorticity/model/`) for the same equations, avoiding the pressure unknown
  entirely.

Both variants read a dimensionless Rayleigh number `DimensionlessNumbers.Ra` from
the input file (see `common/1p_boussinesq_fluidsystem.hh`) and expose the same
diagnostics: a Sherwood-number-style dissolutive mass flux (`writeSherwood`) and a
finger-tracking metric (`writeFinger`, dominant horizontal wavenumber `nHat`,
vorticity energy, tip depth). There are two perturbation styles for the initial
condition — `cosine-perturbations/` (a single deterministic mode, useful for
comparing measured onset/growth rates against linear stability theory) and
`white-noise-perturbations/` (random per-column perturbations, closer to
Riaz et al.'s numerical experiments and to the adaptive/parallel machinery below).

@note This is primarily an **engineering stress test for DuMux's parallel h-adaptive
and dynamic load-balancing infrastructure**, not (yet) a validated quantitative
benchmark with a reference solution and comparison plots in the style of e.g.
@ref benchmark-buckley-leverett. See `NOTES.md` in this directory for the full,
chronological development/debugging log this README is distilled from.

## Parallel execution, adaptivity & dynamic load balancing

`white-noise-perturbations/main_pressure.cc` and `main_vorticity.cc` combine three
things that are individually supported by DuMux but, at the time this test was
written, had never been exercised together end-to-end: **MPI parallelism**,
**h-adaptive refinement/coarsening** (`dumux/adaptive/adapt.hh`), and **dynamic load
rebalancing** (elements migrating between ranks *during* a run, not just once at
startup). This section documents the moving parts, since none of it is obvious from
reading a single file — `common/adaptive/` and `dumux/linear/istlsolvers.hh` are
touched by every step.

### Why rebalancing is needed here

Refinement is driven by a solution-gradient indicator
(`common/adaptive/gridadaptindicator.hh`), so as fingers form and grow, refined cells
concentrate wherever the fingers currently are — which is not where the *initial*
(pre-refinement) domain decomposition put the rank boundaries. Without rebalancing,
whichever rank happens to own the finger-active region ends up carrying most of the
mesh, while other ranks sit mostly idle. `Adaptive.EnableDynamicLoadBalancing true`
periodically (`Adaptive.RebalanceEvery`) or reactively
(`Adaptive.ImbalanceTolerance`) repartitions the mesh so refined work stays spread
across ranks as the simulation evolves.

### The moving parts, in the order they execute

1. **Criterion** — a purely local+collective check, no grid communication beyond two
   reductions:
   ```cpp
   // main_pressure.cc — currentImbalance()
   const Scalar local = localInteriorElements();      // count this rank's own cells
   const Scalar maxLocal = comm.max(local);            // MPI_Allreduce(MAX)
   const Scalar sumLocal = comm.sum(local);            // MPI_Allreduce(SUM)
   return sumLocal > 0 ? (maxLocal / (sumLocal/comm.size()) - 1.0) : 0.0;
   ```
   Triggered by `stepsSinceRebalance >= RebalanceEvery || currentImbalance() > ImbalanceTolerance`.

2. **Migration** (`common/adaptive/boxrebalance.hh`, `rebalanceBox()`) — snapshots the
   current/old solution into a `Dune::PersistentContainer` keyed by element, then calls
   the grid-specific balancer (`common/adaptive/loadbalancer.hh`,
   `GridLoadBalancer<Grid>::apply`):
   - **ALUGrid**: `grid.loadBalance(weights, dataHandle)`, weighted by leaf-descendant
     count per macro element (`LeafCountWeights`), or `grid.repartition(destinations,
     dataHandle)` for the optional deterministic vertical-strip partition
     (`VerticalStripDestinations` — keeps each rank's slice a full top-to-bottom column,
     which matters here because fingers descend vertically and the default
     space-filling-curve partition would otherwise cut through every finger repeatedly).
   - **UGGrid**: UGGrid's plain `loadBalance()` refuses to redistribute an
     already-parallel grid (`Dune::NotImplemented`); `computeSpaceFillingCurveTargetProcessors()`
     computes an explicit target-rank-per-element assignment (Morton order over gathered
     element centroids) and feeds it to the `loadBalance(targetProcessors, fromLevel,
     dataHandle)` overload instead.

   Either way, `dataHandle`'s `gather()`/`scatter()` (a `Dune::CommDataHandleIF`) is
   DUNE's generic per-entity MPI pack/unpack hook, invoked once per element that
   actually crosses a rank boundary. **Gotcha**: a `PersistentContainer` sized once
   before `loadBalance()` has no slot yet for an element that migrates *in* during the
   call — `scatter()` must call `container_.resize()` before writing, or the container's
   internal heap gets corrupted (confirmed via valgrind: `malloc(): invalid next size`).

3. **Rebuild everything indexed by DOF**, in this order:
   ```cpp
   gridGeometry->update(grid.leafGridView());   // new DOF count/order
   gridVariables->updateAfterGridAdaption(x);
   assembler->updateAfterGridAdaption();
   linearSolver->updateAfterGridAdaption(gridGeometry->gridView(), gridGeometry->dofMapper());
   ```
   plus, inside `rebalanceBox()`, a volume-weighted reconstruction of `(x, xOld)` from
   the migrated container and a `syncOwnerToGhost` pass
   (`common/adaptive/griddatacommunication.hh`) to make ghost/overlap copies consistent
   with their owner.

### Why the linear solver needs its own `updateAfterGridAdaption()`

The last call above (`linearSolver->updateAfterGridAdaption(...)`, added to
`Dumux::Detail::IstlIterativeLinearSolver` in `dumux/linear/istlsolvers.hh`) exists
because ownership, ghost status, and global DOF indices are stale for essentially
every degree of freedom after a rebalance — cells (and the DOFs on them) may now live
on a different rank entirely. It rebuilds, in place:

- `Dumux::ParallelISTLHelper` — per-DOF ownership/ghost bitmasks
  (`isOwned_`/`isGhost_`), recomputed by walking the *new* grid view with
  `Dune::CommDataHandleIF` gather/scatter handles.
- `Dune::OwnerOverlapCopyCommunication` (`communication_`) — including its
  `Dune::ParallelIndexSet` (per-DOF `(globalIndex, (localIndex, owner|overlap|copy))`
  triplets) and `Dune::RemoteIndices` (which neighbor rank sees this DOF as what).
- `scalarProduct_` — must match the (possibly changed) `Dune::SolverCategory`
  (sequential/overlapping/nonoverlapping).

It also **drops the cached operator and solver**:
```cpp
linearOperator_ = MatrixOperatorHolder{};
solver_ = nullptr;
```
This matters because DUNE-ISTL's `OverlappingSchwarzOperator`/
`NonoverlappingSchwarzOperator` (`dune/istl/schwarz.hh`) store the communication
object as a **bare reference member** (`const communication_type& communication;`),
not a `shared_ptr`. Reseating `communication_` to a freshly constructed `Comm` frees
the old one; any previously cached operator/solver still referencing it is left
holding a dangling reference. Not resetting these was the actual root cause of the
segfaults documented in `NOTES.md` under "ALUGrid: use-after-free across repeated
`ALUCommunication` constructions" and "out-of-bounds write in
`NonoverlappingSchwarzOperator::novlp_op_apply`" — both symptoms of the same
underlying stale-communication-reference problem, not independent library bugs.

**Consequence for callers**: this test's `NewtonSolver` always calls the 3-argument
`solve(A, x, b)`, which builds a fresh operator from whatever matrix is passed every
call — it never touches the cached `linearOperator_`/`solver_`, so it needs no extra
handling beyond the `updateAfterGridAdaption()` call itself. The cached-operator path
(`setMatrix()` + repeated `solve(x, b)`, used via `Dumux::LinearPDESolver::reuseMatrix()`
for problems whose matrix is known not to change between solves — see
`test/porousmediumflow/tracer/constvel/main.cc` or `examples/diffusion/main.cc` for
examples of that pattern) is different: after `updateAfterGridAdaption()`, `solver_`
is `nullptr`, so the next `solve(x, b)` throws `Dune::InvalidStateException` until
`setMatrix()` is called again with the freshly (re-)assembled Jacobian. This test
does not use that pattern, so it isn't exercised here, but it's the reason the
`updateAfterGridAdaption()` doc comment says *"A matrix set via setMatrix has to be
set again after this call."*

### Status

Per the development log in `NOTES.md`: parallel + h-adaptive without rebalancing has
been confirmed clean (`-np` 1/2/4, box + ALUGrid, ILU+GMRES). With the
`updateAfterGridAdaption` fix in place, dynamic rebalancing was confirmed running to
completion at `-np 4` with `-Adaptive.EnableDynamicLoadBalancing true
-Adaptive.RebalanceEvery 3` (box + ALUGrid, ILU+GMRES) — see the "Update 2026-07-21"
entry in `NOTES.md`. AMG (`AMGBiCGSTABIstlSolver`/`AMGRestartedGMResIstlSolver`) was
evaluated and rejected for this operator: with the stale-index-set bugs fixed, AMG's
coarse-hierarchy setup no longer crashes, but the aggregation coarsening itself
produces `NaN`/singular coarse blocks on this strongly nonsymmetric, mixed
elliptic/parabolic 2×2 block operator regardless of tuning — a structural mismatch,
not a parameter issue. ILU+GMRES (`ILURestartedGMResIstlSolver`) is therefore the
only supported solver for these tests.

@note `NOTES.md`'s own "Bottom line" section (written before the "Update 2026-07-21"
entry above it) still says rebalancing "does not work" — that line predates the fix
and is stale; the dated update entries are the current status.

**Parameters relevant to the parallel/adaptive machinery** (see `params.input`,
`params3d.input`):

| Parameter | Meaning |
|---|---|
| `Adaptive.EnableDynamicLoadBalancing` | opt-in for rebalancing during the run (default `false`) |
| `Adaptive.RebalanceEvery` | rebalance at least every N timesteps regardless of imbalance |
| `Adaptive.ImbalanceTolerance` | rebalance reactively once `max/mean - 1` exceeds this |
| `Adaptive.VerticalStripPartitioning` | ALUGrid-only: use weighted top-to-bottom column strips instead of the default space-filling-curve partition (keeps a descending finger on one rank) |

## References

See the PDF in this directory: Riaz, A., Hesse, M., Tchelepi, H.A., Orr, F.M. (2006).
"Onset of convection in a gravitationally unstable diffusive boundary layer in porous
media." *Journal of Fluid Mechanics*, 548, 87-111.