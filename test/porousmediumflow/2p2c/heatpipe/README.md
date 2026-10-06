# Heatpipe Effect {#benchmark-heatpipe}

## One-dimensional non-isothermal two-phase two-component flow

### Problem description

This benchmark models heat transport in a horizontal porous column containing liquid water and a gas mixture of water vapor and air. Heating the right boundary evaporates water. Vapor flows toward the cooler left end, condenses, and releases latent heat. Capillary forces return liquid water toward the hot end, creating countercurrent flow. At steady state, a two-phase heat-pipe region can coexist with a dry region near the heater.

The transient DuMux simulation uses the BOX discretization of the non-isothermal two-phase two-component model (`TwoPTwoCNI`). Its final profiles are compared with a modified semi-analytical reference based on the heat-pipe formulation of Udell and Fitch @cite Udell1985, as presented by Huang, Kolditz and Shao @cite Huang2015 (see also @cite Meng2022). The reference and the numerical model share selected material properties, but solve different equations and treat dry-out differently. Their profiles are therefore expected to be close, rather than identical.

![Schematic description](heatpipe_schematic_description.png){html: width=80%}

### Geometry and boundary conditions

The column is 2.4 m long. Although the computational grid is two-dimensional, the geometry, material properties and boundary conditions produce flow along the horizontal coordinate $x$. Gravity is disabled.

- **Left boundary:** gas pressure $p_g = 101300\ \mathrm{Pa}$, liquid saturation $S_w = 0.99$ and temperature $T = 341.75\ \mathrm{K}$ ($68.6^\circ\mathrm{C}$). The gas-phase air mole fraction follows from phase equilibrium. Neglecting dissolved air, $x_g^a \approx 1-p_\mathrm{vap}(T)/p_g \approx 0.71$.
- **Right boundary:** an inward heat flux of $100\ \mathrm{W/m^2}$ and zero flux of both mass components. The input parameter is `Problem.HeatFlux = -100`, since negative Neumann flux denotes injection.
- **Other boundaries:** zero component and heat fluxes.
- **Initial state:** $p_g = 101300\ \mathrm{Pa}$, $S_w = 0.5$ and $T = 343.15\ \mathrm{K}$ ($70^\circ\mathrm{C}$), with both phases present.

The left-boundary saturation is slightly below full saturation so that both phases are present. The initial state determines the transient evolution. The reference describes only steady state.

### Material properties

| Parameter | Symbol | Value | Unit |
|-----------|--------|-------|------|
| Intrinsic permeability | $K$ | $10^{-12}$ | m² |
| Porosity | $\phi$ | 0.4 | – |
| Residual liquid saturation | $S_{wr}$ | 0.15 | – |
| Solid thermal conductivity | $\lambda_s$ | 2.8 | W/(m·K) |
| Solid density | $\rho_s$ | 2600 | kg/m³ |
| Solid specific heat capacity | $c_s$ | 700 | J/(kg·K) |

The test uses `HeatPipeReferenceFluidSystem` (`referencefluidsystem.hh`), an adapter of `FluidSystems::H2OAir`. It aligns selected properties with the current reference:

- Gas density follows the ideal-gas mixture law. Liquid density is that of pure water at the local temperature and liquid pressure.
- A finite Henry constant of $10^{20}\ \mathrm{Pa}$ makes dissolved air negligible.
- All components have the common sensible enthalpy $4187(T-273.15)\ \mathrm{J/kg}$. Water vapor has an additional $2.258\cdot10^6\ \mathrm{J/kg}$, giving a fixed latent heat. The common heat capacity preserves heat storage during the transient simulation.
- Local water and air viscosities, Wilke's gas-mixture viscosity rule, thermal conductivities and binary diffusion coefficients are retained from `H2OAir`.

These are benchmark-specific approximations. They do not reproduce the full thermophysical behavior of water and air or the original constant-property Udell–Fitch setup.

#### Capillary pressure and relative permeability

Both models use `HeatPipeLaw`, with effective liquid saturation

```math
S_e = \frac{S_w-S_{wr}}{1-S_{wr}}.
```

The Fatt–Klikoff relative permeabilities @cite Fatt1959 are

```math
k_{rw}=S_e^3, \qquad k_{rg}=(1-S_e)^3.
```

`HeatPipeLaw` bounds these functions between zero and one and replaces the cubic relation with a spline when its argument exceeds 0.95. The Leverett capillary-pressure relation @cite Leverett1941 is

```math
p_c = \gamma\sqrt{\frac{\phi}{K}}
\left[1.417(1-S_e)-2.120(1-S_e)^2+1.263(1-S_e)^3\right],
```

with constant surface tension $\gamma=0.0588\ \mathrm{N/m}$. The implementation extends capillary pressure linearly outside the effective-saturation interval $[0,1]$.

#### Heat conduction and gas diffusion

Both models use the Somerton effective thermal conductivity:

```math
\lambda(S_w)=\lambda_\mathrm{dry}
+\sqrt{S_w}\left(\lambda_\mathrm{wet}-\lambda_\mathrm{dry}\right),
```

```math
\lambda_\mathrm{wet}=\lambda_s^{1-\phi}(\lambda_w^\mathrm{fluid})^\phi,
\qquad
\lambda_\mathrm{dry}=\lambda_s^{1-\phi}(\lambda_g^\mathrm{fluid})^\phi.
```

Liquid-water conductivity is evaluated at the local state. Gas conductivity uses the constant air value $0.0255535\ \mathrm{W/(m\,K)}$. The resulting endpoints are approximately $1.57$–$1.59$ and $0.428\ \mathrm{W/(m\,K)}$ for the wet and dry medium, respectively.

Gas diffusion uses `DiffusivityMillingtonQuirk` @cite MILLINGTON1961, with $S_g=1-S_w$:

```math
D_{pm}=\phi S_g^3\sqrt[3]{\phi S_g}\,D_g^{aw}(T,p_g),
```

```math
D_g^{aw}(T,p_g)=2.13\cdot10^{-5}\ \mathrm{m^2/s}
\frac{10^5\ \mathrm{Pa}}{p_g}
\left(\frac{T}{273.15\ \mathrm{K}}\right)^{1.8}.
```

### Semi-analytical reference

`test_heatpipe_odesolver.cc` integrates four coupled spatial ODEs for $(S_e,p_g,x_g^a,T)$ from the left boundary using explicit fourth-order Runge–Kutta integration. The ODE structure follows the literature formulation, but the material properties are adapted to this DuMux test. Liquid density and component viscosities are evaluated with local `H2OAir` properties, gas viscosity uses Wilke's mixing rule, and heat conductivity and gas diffusivity use the Somerton and Millington–Quirk laws described above.

These choices are intended to use the same material-property laws in the reference and numerical test. They are not fitted to the numerical profiles. However, they change the reference problem. The plotted curve is computed from these adapted ODEs, rather than taken from the published Udell–Fitch results. The comparison therefore assesses agreement with a literature-based, DuMux-specific reference. It is not an exact reproduction of the original Udell–Fitch solution, which assigned constant values to these properties.

The reference uses pure liquid water, ideal-gas mixture density and fixed latent heat $h_v^w=2.258\cdot10^6\ \mathrm{J/kg}$. Define the kinematic viscosities $\nu_g=\mu_g/\rho_g$, $\nu_w=\mu_w/\rho_w$ and their ratio $\beta=\nu_w/\nu_g$. The auxiliary quantities implemented in the reference are

```math
\alpha=1+\frac{p_c}{\rho_w h_v^w}, \qquad
\xi=\frac{1}{k_{rg}}\left(1+\frac{\rho_wRT}{p_gM^w(1-x_g^a)}\right)
+\frac{\beta}{k_{rw}},
```

```math
\delta=\frac{\rho_w(h_v^w)^2K\alpha}{\lambda\nu_gT}, \qquad
\zeta=\frac{K\rho_wRT}{M^w\rho_g\nu_gD_{pm}}
\frac{x_g^a}{1-x_g^a}
\left(\frac{p_gM^w}{\rho_wRT}+\frac{1}{1-x_g^a}\right),
\qquad
\eta=\frac{\delta}{\delta+\xi+\zeta}.
```

Here $R$ is the universal gas constant, $M^w$ is the molar mass of water and $\eta$ partitions the heat flux between heat-pipe transport and conduction. With the signed flux $q=-100\ \mathrm{W/m^2}$, the ODEs are

```math
\frac{\mathrm{d}S_e}{\mathrm{d}x}
=-\left(\frac{1}{1-x_g^a}+\beta\frac{k_{rg}}{k_{rw}}\right)
\frac{\eta q\nu_g}{K h_v^w k_{rg}}
\Big/\frac{\mathrm{d}p_c}{\mathrm{d}S_e},
\qquad
\frac{\mathrm{d}p_g}{\mathrm{d}x}
=-\frac{\eta q\nu_g}{K h_v^w k_{rg}(1-x_g^a)},
```

```math
\frac{\mathrm{d}x_g^a}{\mathrm{d}x}
=\frac{\eta qx_g^a}{h_v^wD_{pm}\rho_g(1-x_g^a)},
\qquad
\frac{\mathrm{d}T}{\mathrm{d}x}=-\frac{q(1-\eta)}{\lambda}.
```

Integration approaches $S_e=0$. Steps crossing that limit or producing an invalid state are rejected and retried with a smaller step. The final accepted wet position is used as the reference front. Beyond it, the output sets $S_w=0$ and $x_g^a=0$, keeps gas pressure constant, and continues temperature by dry-medium conduction:

```math
T(x)=T_\mathrm{front}-\frac{q}{\lambda_\mathrm{dry}}(x-x_\mathrm{front}).
```

This continuation is appended to the two-phase ODE solution without resolving liquid disappearance.

### Differences between the numerical model and the reference

Matching material properties does not make the reference an exact solution of the transient `TwoPTwoCNI` equations.

| Aspect | Semi-analytical reference | Numerical test |
|--------|-------------------------|----------------|
| Governing equations | Four reduced steady-state spatial ODEs | Transient component and energy balances, evaluated near steady state |
| Phase equilibrium | Temperature and gas composition evolve through the reduced ODE relations. IAPWS vapor pressure initializes the boundary composition | Local compositional equilibrium uses fugacity coefficients and IAPWS vapor pressure throughout the domain |
| Dissolved air | Exactly zero in the liquid state used for property evaluation | Negligible, but finite, with the large Henry constant |
| Heat storage | Absent from the steady-state equations | Fluid and solid energy storage affect the transient approach to steady state |
| Dry-out criterion | $S_e\to0$, corresponding to $S_w\to S_{wr}=0.15$ | Liquid can evaporate below residual saturation before its phase disappears |
| Dry region | Saturation jumps from approximately 0.15 to zero. The dry region has constant pressure and a linear temperature profile | Gas-only states and their profiles follow from the discretized balances and phase switching |

The reference uses a reduced closure for heat transport and composition gradients, while the numerical model assembles advective and diffusive component fluxes and the energy balance. Agreement of their property functions alone does not establish equivalence of those balances. No complete equivalence of the reference's thermodynamic reduction with the numerical phase-equilibrium model is assumed here.

Consequently, the discrepancy can contain both discretization error and differences between the two models. Refinement may move a numerical profile across the reference, and a persistent difference does not by itself indicate a solver error. Exact convergence of every field to this reference is not guaranteed.

### Running the benchmark

From the heatpipe source directory, run

```bash
python3 compile_run_plot.py
```

The script expects a configured `build-cmake` directory and requires NumPy, Matplotlib and PyVista. It builds `test_heatpipe_box` and `test_heatpipe_odesolver`, writes the reference to `heatpipe_reference.csv`, and runs the numerical test with 120, 240 and 480 cells along $x$ ($\Delta x=0.02$, $0.01$ and $0.005\ \mathrm{m}$).

The following figures are written to `build-cmake/test/porousmediumflow/2p2c/heatpipe`:

- `heatpipe_lineplot_comparison.svg`: saturation, temperature, gas pressure and gas-phase air mole fraction at 480 cells, compared with the reference.
- `heatpipe_saturation_comparison.svg`: the saturation comparison alone.
- `heatpipe_sw.png`: the numerical saturation field at 480 cells.
- `heatpipe_grid_convergence.svg`: saturation and temperature at all three resolutions, zoomed to the reference front.

### Results and grid refinement

![Line plot](heatpipe_lineplot_comparison.svg)

![Saturation field](heatpipe_sw.png){html: width=80%}

The default CTest simulation uses 120 cells. Its regression check identifies the numerical front as the first vertex containing only gas and requires its position to be between the configured reference position and 0.12 m downstream. This one-sided tolerance reflects the observed behavior of the default benchmark. Other grids or parameter choices may produce errors in either direction.

![Grid convergence](heatpipe_grid_convergence.svg)
