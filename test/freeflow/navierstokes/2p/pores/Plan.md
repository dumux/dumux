Step 1 — Initial conditions
Tanh profile centered at y = 0 (channel/box interface):


φ(y) = tanh(−y / (√2 · ε))    [−1 in channel, +1 in pores]
μ = 0 everywhere, p = 0 everywhere, u = 0 everywhere.

Step 2 — Inlet/outlet BCs for the gas channel
The .geo already has inlet (left wall, Line 5) and outlet_chan (right wall, Line 3) tags. Use BoundarySegmentIndexFlag (already in properties.hh) to identify them.

Momentum:

Inlet: Dirichlet, prescribed (u_in, 0) — parabolic profile optional
Outlet: Neumann (free outflow / zero normal stress)
Top wall, channel bottom-left/right: no-slip Dirichlet
Pillar surfaces: no-slip Dirichlet (see Step 3)
Mass: All Neumann at inlet/outlet is acceptable if the inlet flow is pure gas. The advective flux carries φ = −1 naturally once the velocity field is driven. Alternatively, Dirichlet φ = −1 at the inlet is cleaner but requires enabling Dirichlet for the mass problem.

The pillar box sides and cone exit remain Neumann = 0 (closed, no flow in/out of the pore region laterally).

Step 3 — No-slip and contact angle on pillar surfaces
(Unchanged from previous plan — geometric distance check for pillar identification)

Momentum: Dirichlet v = 0

Mass (chemical potential equation): Neumann flux for contact angle as before

Step 4 — Evaporation as a phase-change source term
Evaporation is not a boundary condition. It is a non-conservative bulk source in the phase field equation that drives φ → −1 (liquid → gas) at the diffuse interface. Add to computeSource in localresidual.hh:


∂φ/∂t + ∇·(uφ) − M∇²μ = −Γ(φ)
The simplest physically motivated choice for Γ, concentrated at the interface:


Γ(φ) = γ_evap · (1 − φ²) / ε
where (1 − φ²) is a smoothed indicator that peaks at the diffuse interface (φ = 0) and vanishes in the bulk phases (φ = ±1). Dividing by ε keeps the integral over the interface profile O(1) as ε → 0. γ_evap is the evaporation rate coefficient with units [1/s] (controls the speed of interface recession).

This goes into problem.source() in the mass branch:


source[Indices::phaseFieldEqIdx] -= evaporationRate_ * (1.0 - phi*phi) / interfaceThickness_;
A more physical model would tie γ_evap to the local vapor flux from the NS solution (vapor concentration gradient at the interface), but this simple form is a reasonable starting point.

Note: this source breaks conservation of φ by design — liquid mass is not conserved, it is removed. This is appropriate for evaporation.

Step 5 — Remove Hele-Shaw drag, clean up dead parameters
HeleShawDragCoeff, EvaporationThreshold, and the old top-wall evaporation Neumann logic all go away. Add ContactAngle and revise EvaporationRate to mean the rate coefficient γ_evap above.

Verification sequence
Pure flow, no evaporation: inlet velocity, no phase change, check velocity profile in channel around pillars
Static meniscus: zero velocity, no evaporation, check contact angle equilibrium at pillars
Evaporation without flow: no inlet velocity, evaporation source active, watch interface recede into pores
Full scenario: all together — vapor flow drives evaporation, interface recedes, meniscus shape around pillars controlled by contact angle