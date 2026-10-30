# Heatpipe Effect {#benchmark-heatpipe}

## One-dimensional non-isothermal two-phase two-component flow

### Problem description

This benchmark models heat transport in a horizontal porous column containing liquid water and a gas mixture of water vapor and air. Heating the right boundary evaporates water. Vapor flows toward the cooler left end, condenses, and releases latent heat. Capillary forces return liquid water toward the hot end, creating countercurrent flow. At steady state, a two-phase heat-pipe region can coexist with a dry region near the heater.

The transient DuMux simulation uses the BOX discretization of the non-isothermal two-phase two-component model (`TwoPTwoCNI`). Its final profiles are compared with the steady-state semi-analytical solution of Udell and Fitch @cite Udell1985, as extended by Huang, Kolditz and Shao @cite Huang2015 (see also @cite Meng2022) to an ideal gas mixture with a composition-dependent viscosity. Setup, material properties and closure relations follow the supplementary code of Huang et al., and the numerical model adopts them, so that the two models describe the same physical problem.

![Schematic description](heatpipe_schematic_description.png){html: width=80%}

### Geometry and boundary conditions

The column is 2 m long. Although the computational grid is two-dimensional, the geometry, material properties and boundary conditions produce flow along the horizontal coordinate $x$. Gravity is disabled.

- **Left boundary:** gas pressure $p_g = 101325\ \mathrm{Pa}$, liquid saturation $S_w = 0.99$ and gas-phase air mole fraction $x_g^a = 0.7$. The temperature in local equilibrium with this state, eq. (11) of Huang et al., is $T = 341.7588\ \mathrm{K}$ ($68.6^\circ\mathrm{C}$). It is prescribed as `Problem.LeftTemperature`, and `test_heatpipe_odesolver` checks that it agrees with the reference.
- **Right boundary:** an inward heat flux of $100\ \mathrm{W/m^2}$ and zero flux of both mass components. The input parameter is `Problem.HeatFlux = -100`, since negative Neumann flux denotes injection.
- **Other boundaries:** zero component and heat fluxes.
- **Initial state:** $p_g = 101325\ \mathrm{Pa}$, $S_w = 0.7$ and $T = 343.15\ \mathrm{K}$ ($70^\circ\mathrm{C}$), with both phases present, as in the numerical setup of Huang et al.

The left-boundary saturation is slightly below full saturation so that both phases are present. The initial state determines the transient evolution. The reference describes only steady state.

### Material properties

| Parameter | Symbol | Value | Unit |
|-----------|--------|-------|------|
| Intrinsic permeability | $K$ | $10^{-12}$ | m² |
| Porosity | $\phi$ | 0.35 | – |
| Residual liquid saturation | $S_{wr}$ | 0 | – |
| Solid density | $\rho_s$ | 2600 | kg/m³ |
| Solid specific heat capacity | $c_s$ | 700 | J/(kg·K) |
| Surface tension | $\gamma$ | 0.05878 | N/m |
| Liquid density | $\rho_w$ | 1000 | kg/m³ |
| Liquid viscosity | $\mu_w$ | $2.938\cdot10^{-4}$ | Pa·s |
| Air and water-vapor viscosity | $\mu_g^a$, $\mu_g^w$ | $2.094\cdot10^{-5}$, $1.2\cdot10^{-5}$ | Pa·s |
| Molar masses of water and air | $M^w$, $M^a$ | 0.018016, 0.028 | kg/mol |
| Universal gas constant | $R$ | 8.314 | J/(mol·K) |
| Macroscopic gas diffusion coefficient | $D_{pm}$ | $2.6\cdot10^{-6}$ | m²/s |
| Wet and dry thermal conductivity | $\lambda_\mathrm{wet}$, $\lambda_\mathrm{dry}$ | 1.13, 0.582 | W/(m·K) |
| Latent heat of vaporization | $h_v^w$ | $2.258\cdot10^6$ | J/kg |

All values are taken from the supplementary code of Huang et al. @cite Huang2015, which partly differs from the parameter table printed in their paper. The solid density and heat capacity only affect the transient approach to steady state.

#### Fluid system

`HeatPipeReferenceFluidSystem` (`referencefluidsystem.hh`) adapts `FluidSystems::H2OAir` to these assumptions:

- Constant liquid density and viscosity. A Henry constant of $10^{20}\ \mathrm{Pa}$ makes dissolved air negligible.
- Ideal gas mixture with the mole-fraction weighted viscosity $\mu_g = x_g^a\mu_g^a + (1-x_g^a)\mu_g^w$ and a constant diffusion coefficient. The molar masses and the gas constant are those listed above.
- Vapor pressure $p_v$ from the Clausius–Clapeyron equation with constant latent heat, corrected for capillary effects with the Kelvin equation in the form of Huang et al.:

```math
p_v = p_0\exp\left[\frac{h_v^wM^w}{R}\left(\frac{1}{T_0}-\frac{1}{T}\right)-\frac{(p_c-x_g^ap_g)M^w}{\rho_wRT}\right],
\qquad p_0 = 101325\ \mathrm{Pa},\ T_0 = 373.15\ \mathrm{K}.
```

Since $p_c-x_g^ap_g = p_v-p_w$, this is an implicit equation for $p_v$, which is solved by fixed-point iteration.

- All components have the common sensible enthalpy $4187(T-273.15)\ \mathrm{J/kg}$. Water vapor has an additional $h_v^w$. At steady state, the total component mass flux vanishes, so only the latent heat contributes to the advective heat flux, as in the reference energy balance.

#### Capillary pressure and relative permeability

The test uses `HeatPipeLaw`, with effective liquid saturation

```math
S_e = \frac{S_w-S_{wr}}{1-S_{wr}}.
```

The Fatt–Klikoff relative permeabilities @cite Fatt1959 are

```math
k_{rw}=S_e^3, \qquad k_{rg}=(1-S_e)^3,
```

and the Leverett capillary-pressure relation @cite Leverett1941 is

```math
p_c = \gamma\sqrt{\frac{\phi}{K}}
\left[1.417(1-S_e)-2.120(1-S_e)^2+1.263(1-S_e)^3\right].
```

`HeatPipeLaw` replaces the cubic relative permeabilities with a spline when their argument exceeds 0.95, which affects only a short region next to the left boundary.

#### Heat conduction and gas diffusion

The effective thermal conductivity (`thermalconductivity.hh`) interpolates between constant wet and dry values:

```math
\lambda(S_w)=\lambda_\mathrm{dry}
+\sqrt{S_w}\left(\lambda_\mathrm{wet}-\lambda_\mathrm{dry}\right).
```

Since $S_{wr}=0$, $S_w$ equals the effective saturation used by Huang et al. This has the form of `ThermalConductivitySomertonTwoP`, but with fixed end points instead of geometric means of the solid and fluid conductivities.

The macroscopic gas diffusion coefficient $D_{pm}$ is used without scaling by porosity, saturation or tortuosity (`diffusivity.hh`). As in the reference, the diffusive fluxes follow Fick's law for mole fractions (`ReferenceSystemFormulation::molarAveraged`).

### Semi-analytical reference

`test_heatpipe_odesolver.cc` integrates four coupled spatial ODEs for $(S_e,p_g,x_g^a,T)$ from the left boundary using explicit fourth-order Runge–Kutta integration. It implements the material properties listed above independently of the numerical model. Define the kinematic viscosities $\nu_g=\mu_g/\rho_g$, $\nu_w=\mu_w/\rho_w$ and their ratio $\beta=\nu_w/\nu_g$. The auxiliary quantities are

```math
\alpha=1+\frac{p_c-x_g^ap_g}{\rho_w h_v^w}, \qquad
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

The equation for $\zeta$ uses $p_g$ as in the derivation of Huang et al.; their code contains $\rho_g$ at this place, which changes the solution by less than $10^{-3}$. Otherwise, the reference reproduces the profiles of their code to within $10^{-3}$. Integration starts at the boundary state above and stops where $S_e\to0$, which defines the dry-out front at $x\approx1.869\ \mathrm{m}$. Steps that cross this limit or produce an invalid state are rejected and retried with a smaller step. Beyond the front, the output sets $S_w=0$ and $x_g^a=0$, keeps gas pressure constant, and continues temperature by dry-medium conduction, $T(x)=T_\mathrm{front}-q(x-x_\mathrm{front})/\lambda_\mathrm{dry}$.

### Remaining differences between the numerical model and the reference

| Aspect | Semi-analytical reference | Numerical test |
|--------|-------------------------|----------------|
| Governing equations | Steady-state spatial ODEs | Transient component and energy balances, run until steady state |
| Dissolved air | Zero | Negligible, but finite, with the large Henry constant |
| Diffusive flux | $\rho_g D_{pm}\,\partial x_g^a/\partial x$ (mass density, mole fraction) | $c_g D_{pm}\,\partial x_g^a/\partial x$ (molar density, mole fraction) |
| Relative permeabilities | Cubic for all $S_e$ | Spline regularization for arguments above 0.95 |
| Dry-out | Integration stops as $S_e\to0$; beyond the front, $S_w=0$ | The liquid phase disappears where it has evaporated completely |

### Running the benchmark

From the heatpipe source directory, run

```bash
python3 compile_run_plot.py
```

The script expects a configured `build-cmake` directory and requires NumPy, Matplotlib and PyVista. It builds `test_heatpipe_box` and `test_heatpipe_odesolver`, writes the reference to `heatpipe_reference.csv`, and runs the numerical test with 100, 200 and 400 cells along $x$ ($\Delta x=0.02$, $0.01$ and $0.005\ \mathrm{m}$).

The following figures are written to `build-cmake/test/porousmediumflow/2p2c/heatpipe`:

- `heatpipe_lineplot_comparison.svg`: saturation, temperature, gas pressure and gas-phase air mole fraction at 400 cells, compared with the reference.
- `heatpipe_saturation_comparison.svg`: the saturation comparison alone.
- `heatpipe_sw.png`: the numerical saturation field at 400 cells.
- `heatpipe_grid_convergence.svg`: saturation and temperature at all three resolutions, zoomed to the reference front.

### Results and grid refinement

![Line plot](heatpipe_lineplot_comparison.svg)

![Saturation field](heatpipe_sw.png){html: width=80%}

The numerical profiles agree closely with the reference. The numerical front lies downstream of the reference front and approaches it under grid refinement: the distance is about 0.09 m, 0.04 m and 0.02 m for 100, 200 and 400 cells (0.006 m for 800 cells). In the two-phase zone, the saturation and gas-pressure differences roughly halve with each refinement. The temperature and air mole fraction differences in the air-rich zone near the left boundary remain below 0.3 K and 0.01, respectively, but do not decrease further, which is attributed to the different diffusive-flux formulations listed above.

The default CTest simulation uses 100 cells. Its regression check identifies the numerical front as the first vertex containing only gas and requires its position to be between the reference position and 0.12 m downstream.

![Grid convergence](heatpipe_grid_convergence.svg)
