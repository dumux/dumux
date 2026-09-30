# Heatpipe Effect {#benchmark-heatpipe}

## One-dimensional non-isothermal two-phase two-component flow

**Problem Description**

The Heatpipe Effect can be observed in a nonisothermal water-gas system in a porous medium, in which the heat transfer processes convection, conduction, and diffusion, as well as capillary forces, play an essential role. Udell and Fitch @cite Udell:1985 provide a semi-analytical solution for this system, which is practical for comparing with numerical results (e.g., see @cite emmertpromo).

A one-dimensional horizontal porous column is considered. A constant heat flux is applied at the right boundary. Due to the heat flux, the system is heated until boiling temperature is reached and steam is produced at the right-hand boundary. This causes a pressure gradient in the gas phase and the steam flows away from the heat source. After reaching cooler regions of the column, the steam condenses and sets free its latent heat of vaporization. After a while, a non-uniform saturation profile is obtained with a gradient from the cooler to the hot end of the heatpipe.

According to the capillary pressure–saturation relationship a gradient of the capillary pressure into the same direction is produced. Hence, the pressure gradients of the phases have opposite directions and a circulation flow is created. After a stationary system state has been reached, three regions can be distinguished, each of them associated with a dominant heat-transport process.

![Schematic description](heatpipe_schematic_description.png){html: width=80%}


**Semi-analytical Reference Solution**

Udell and Fitch @cite Udell:1985 derive four coupled first-order differential equations for pressure, saturation, temperature and gas-phase mole fraction. These equations are solved by numerical integration by means of a fourth-order Runge–Kutta method. The numerical simulation of the heatpipe system was carried out with the BOX discretization method. Note that the choice of BOX or CVFE makes no difference in the present one-dimensional case.

Here, the formulation given by Huang, Kolditz and Shao @cite Huang:2015 (see also @cite ogs:heatpipe) is used. The system is integrated from the left (Dirichlet) boundary towards the heat source, using the effective wetting-phase saturation $S_e = (S_w - S_{wr})/(1-S_{wr})$ (see the $p_c$, $k_{rw}$, $k_{rg}$ relations given below in **Setup**) as the integrated state variable, together with the gas-phase pressure $p_g$, the gas-phase air mole fraction $x_g^a$ and the temperature $T$. With the gas-phase density $\rho_g = \rho_g^a + \rho_g^w$, the gas-phase viscosity $\mu_g$ (mixture of $\mu_g^a$ and $\mu_g^w$ according to Wilke's rule), the kinematic viscosities $\nu_g = \mu_g/\rho_g$ and $\nu_w = \mu_w/\rho_w$, the mobility ratio $\beta = \nu_w/\nu_g$, the saturation-dependent heat conductivity $\lambda(S_w) = \lambda_{pm}^{S_w=0} + \sqrt{S_w}(\lambda_{pm}^{S_w=1}-\lambda_{pm}^{S_w=0})$ and the diffusive pore conductance $D_{pm}$ (both matched to the numerical model's actual effective-property laws, see **Setup**), the following auxiliary quantities are introduced:
```math
\alpha = 1 + \frac{p_c}{\rho_w h_v^w}, \qquad
\xi = \frac{1}{k_{rg}}\left(1 + \frac{\rho_w R T}{p_g M^w}\frac{1}{1-x_g^a}\right) + \frac{\beta}{k_{rw}},
```
```math
\delta = \frac{\rho_w (h_v^w)^2 K \alpha}{\lambda \nu_g T}, \qquad
\zeta = \frac{K \rho_w R T}{M^w \rho_g \nu_g D_{pm}}\frac{x_g^a}{1-x_g^a}\left(\frac{p_g M^w}{\rho_w R T} + \frac{1}{1-x_g^a}\right),
```
```math
\eta = \frac{\delta}{\delta + \xi + \zeta} \, ,
```
where $\eta \in [0,1]$ partitions the imposed heat flux $q$ between the phase-change-driven heat-pipe processes (saturation, pressure and composition gradients) and direct conduction. The state vector $(S_e, p_g, x_g^a, T)$ then evolves according to
```math
\frac{\mathrm{d}S_e}{\mathrm{d}x} = -\left(\frac{1}{1-x_g^a} + \beta\frac{k_{rg}}{k_{rw}}\right)\frac{\eta\, q\, \nu_g}{K\, h_v^w\, k_{rg}} \Big/ \frac{\mathrm{d}p_c}{\mathrm{d}S_e},
\qquad
\frac{\mathrm{d}p_g}{\mathrm{d}x} = -\frac{\eta\, q\, \nu_g}{K\, h_v^w\, k_{rg}}\frac{1}{1-x_g^a},
```
```math
\frac{\mathrm{d}x_g^a}{\mathrm{d}x} = \frac{\eta\, q\, x_g^a}{h_v^w\, D_{pm}\, \rho_g\, (1-x_g^a)},
\qquad
\frac{\mathrm{d}T}{\mathrm{d}x} = -\frac{q\,(1-\eta)}{\lambda} \, .
```
Integration stops once the wetting phase dries out ($S_e \to 0$), after which the temperature is continued analytically assuming pure conduction through the dry medium ($\mathrm{d}T/\mathrm{d}x = q/\lambda_{pm}^{S_w=0}$) to cover the remainder of the domain. The fluid properties are evaluated at the local state along the column with the same fluid system as the numerical model, `FluidSystems::H2OAir`: $\rho_w$ and $\mu_w$ at the liquid pressure $p_g - p_c$ from the IAPWS formulations of `Components::H2O`, and $\mu_g$ from `Components::H2O` and `Components::Air` combined with Wilke's mixing rule. As the ODE system is derived for an ideal gas, $\rho_g$ follows from the ideal gas law. Likewise, $p_c$, $k_{rw}$ and $k_{rg}$ are evaluated with `HeatPipeLaw`, and $\lambda(S_w)$ and $D_{pm}$ with the effective-property laws of the numerical model, `ThermalConductivitySomertonTwoP` and `DiffusivityMillingtonQuirk` (see **Setup**). Only the latent heat of vaporization $h_v^w = 2.258\cdot 10^{6}\ \text{J/kg}$ (at the normal boiling point) is a fixed value of a reference state.


**Setup**

A one-dimensional horizontal porous column is considered:

A constant heat flux of $q = 100\ \mathrm{W/m^2}$ is imposed at the right boundary of the horizontal column. Zero-flux (Neumann) boundary conditions are prescribed for all mass components. The initial conditions in the entire domain are $p_g = 101300$ Pa, $S_w = 0.5$ and $T = 70^\circ$C.

At the left boundary, Dirichlet boundary conditions are applied for the gas-phase pressure $p_g = 101300\ \mathrm{Pa}$, the water saturation $S_w = 0.99$ (i.e. $S_e \approx 0.988$; slightly below full saturation, so that both phases are present, which improves the convergence behavior) and the temperature $T = 68.6^\circ\mathrm{C}$. Assuming local thermodynamic equilibrium, the air mole fraction in the gas phase follows from the vapor pressure of water, $x_g^a \approx 1 - p_\text{vap}(T)/p_g \approx 0.71$ (neglecting the small amount of air dissolved in the liquid phase).

The following model parameters were used for the simulation run:

| Parameter                                    | Symbol    | Value                   | Unit  |
|----------------------------------------------|-----------|-------------------------|-------|
| Permeability                                 | $K$       | $1.0\text{e-}12$        | m²    |
| Porosity                                     | $\phi$    | $0.4$                   | -     |
| Residual wetting-phase saturation            | $S_{wr}$  | $0.15$                   | -     |
| Solid (grain) thermal conductivity           | $\lambda_s$ | $2.8$                 | W/(m*K) |
| Soil grain density                           | $\varrho_s$   | $2600$              | kg/m³  |
| Specific heat capacity of the soil grains    | $c_s$     | $700$                   | J/(kg*K) |

The fluid properties of water and air are not input parameters; they follow from the default `FluidSystems::H2OAir` fluid system in both the numerical model and the semi-analytical solution (see above).

The effective heat conductivity $\lambda(S_w)$ is not an independent input but computed by DuMux's `ThermalConductivitySomertonTwoP` (the default for `TwoPTwoCNI`) as a porosity-weighted geometric mean of $\lambda_s$ and the phase heat conductivities, interpolated between the dry and fully saturated endpoints with $\sqrt{S_w}$:
```math
\lambda_{pm}^{S_w=1} = \lambda_s^{1-\phi}\left(\lambda_w^\text{fluid}\right)^\phi \approx 1.57\text{--}1.59\ \text{W/(m*K)}, \qquad
\lambda_{pm}^{S_w=0} = \lambda_s^{1-\phi}\left(\lambda_g^\text{fluid}\right)^\phi \approx 0.428\ \text{W/(m*K)}.
```
Here, the IAPWS liquid water heat conductivity $\lambda_w^\text{fluid}(T, p_w)\approx 0.66\text{--}0.68\ \text{W/(m*K)}$ of `Components::H2O`, evaluated at the local state in both models, and the constant air heat conductivity $\lambda_g^\text{fluid} = 0.0255535\ \text{W/(m*K)}$ of `Components::Air` are used. These endpoints differ noticeably from the values historically quoted for this benchmark (1.13 and 0.582 W/(m*K), for a different solid conductivity); with the historical values, the semi-analytical dry-out front would lie about $0.16\ \text{m}$ closer to the left boundary.

A function according to Fatt and Klikoff @cite Fatt:1959 is chosen for the relative permeability-saturation relationship:
```math
k_{rg} = (1 - S_e)^3 \quad \text{for steam (gas phase)} \nonumber 
```
```math
k_{rw} = S_e^3 \quad \text{for water} \, 
```
where `HeatPipeLaw` regularizes both functions with a spline for arguments above $0.95$.
with the effective water-phase saturation
```math
S_e = \frac{S_w - S_{wr}}{1-S_{wr}}\; .
```
For the capillary pressure-saturation relationship, the following function of Leverett @cite lev1 is used:
```math
p_c = p_0 \, \gamma \left[ 1.417(1-S_e) - 2.120(1-S_e)^2 + 1.263(1-S_e)^3 \right] .
```

The surface tension is $\gamma = 0.0588\ \text{N/m}$, based on the literature value of $0.05878\ \text{N/m}$ at $T = 100.5^\circ\text{C}$, which is close to the temperature in the heat-pipe zone. It is constant in both models, and $p_0 = \sqrt{\phi/K}$ applies for the scaling pressure.

The diffusive pore conductance $D_{pm}$ is likewise not a constant but computed with DuMux's default effective diffusivity model for `TwoPTwoC`, `DiffusivityMillingtonQuirk` @cite MILLINGTON1961, given by
```math
D_{pm} = \phi\, S_g^3 \sqrt[3]{\phi\, S_g}\; D_g^{aw}(T, p_g),
```
using the binary diffusion coefficient of the (unoverridden) `H2OAir` fluid system, `BinaryCoeff::H2O_Air::gasDiffCoeff`:
```math
D_g^{aw}(T, p_g) = 2.13\cdot 10^{-5}\ \text{m}^2\text{/s} \cdot \frac{10^5\ \text{Pa}}{p_g} \left(\frac{T}{273.15\ \text{K}}\right)^{1.8} .
```

The dimension of the model domain in $x$-direction is chosen at 2.4 m. However, this is not important for the length of the heatpipe after the stationary state has been reached as long as the domain is sufficiently large for the heatpipe to be built. The domain is discretized with 120, 240 and 480 cells, i.e. $\Delta x = 0.02$ m, $0.01$ m and $0.005$ m (see **Grid Resolution**).


**Result**

To run the test and produce the plots below, execute:
```bash
python3 compile_run_plot.py
```
The script builds once, then runs the simulation at three grid resolutions (120, 240 and 480 cells in $x$-direction; see **Grid Resolution** below), computes the semi-analytical solution described above, and produces four figures in the `build-cmake` directory:
- `heatpipe_lineplot_comparison.svg`: comparison of the finest-resolution (480 cells) numerical solution against the semi-analytical solution for wetting-phase saturation $S_w$, temperature $T$, gas-phase pressure $p_g$, and gas-phase air mole fraction $x_g^a$, all plotted along $x$
- `heatpipe_saturation_comparison.svg`: the wetting-phase saturation panel alone (finest resolution), used as the thumbnail on the benchmarks overview page
- `heatpipe_sw.png`: numerical wetting-phase saturation field (finest resolution)
- `heatpipe_grid_convergence.svg`: wetting-phase saturation and temperature at all three resolutions, zoomed to the dry-out front (see **Grid Resolution**)

The script expects PyVista, Matplotlib and NumPy to be available for post-processing.

![Line plot](heatpipe_lineplot_comparison.svg)

![Saturation field](heatpipe_sw.png){html: width=80%}


**Grid Resolution**

The comparison above uses the finest of three grid resolutions (120, 240 and 480 cells in $x$-direction) that `compile_run_plot.py` runs for a grid-convergence study. The default test runs at the coarsest, 120-cell resolution (`grids/heatpipe.dgf`), in order to reduce runtime. It fails if the simulated dry-out front does not lie between the semi-analytical position and $0.12\ \text{m}$ downstream of it, as numerical diffusion shifts the front towards the heat source (see below).

![Grid convergence](heatpipe_grid_convergence.svg)

The dry-out front position, where the wetting phase becomes immobile ($k_{rw}\to 0$) and fully evaporates, is $2.260\ \text{m}$, $2.210\ \text{m}$ and $2.185\ \text{m}$ at 120, 240 and 480 cells, respectively, versus the semi-analytical value of $2.157\ \text{m}$. The gap of $0.103\ \text{m}$, $0.053\ \text{m}$ and $0.028\ \text{m}$ roughly halves with each grid refinement, i.e. the dry-out front converges with first order towards the semi-analytical solution.
