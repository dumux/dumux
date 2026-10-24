# Udell's Heat Pipe {#benchmark-heatpipe-udell}

## Steady counter-current flow of liquid water and vapor in a horizontal porous column

**Problem Description**

A horizontal porous column filled with liquid water and its vapor is heated at one end and kept at a fixed state at the other end. The vapor generated at the heated end flows towards the cold end, and capillary forces drive the liquid back towards the heated end. This counter-current flow transports heat as latent heat, like a conventional heat pipe. The benchmark follows the steady one-dimensional analysis of Udell @cite Udell1985HeatPipe and its assumptions:

- no gravity (horizontal heating)
- a single component in a liquid and a gas phase, in local thermodynamic equilibrium
- constant densities and viscosities of both phases
- heat transported only as latent heat of the counter-current flow, i.e. no heat conduction and a vapor enthalpy that exceeds the liquid enthalpy by the constant latent heat $h_{fg}$
- homogeneous medium

For the phases $\alpha \in \{l, g\}$ with saturations $S_\alpha$, Darcy's law gives

$$\mathbf{v}_\alpha = -\frac{K k_{r\alpha}}{\mu_\alpha} \nabla p_\alpha,$$

and the balances of mass and energy read

$$\frac{\partial}{\partial t}\Big(\phi \sum_\alpha \varrho_\alpha S_\alpha\Big) + \nabla\cdot\sum_\alpha \varrho_\alpha \mathbf{v}_\alpha = 0, \qquad \frac{\partial}{\partial t}\Big(\phi \sum_\alpha \varrho_\alpha u_\alpha S_\alpha + (1-\phi)\varrho_s c_s T\Big) + \nabla\cdot\sum_\alpha \varrho_\alpha h_\alpha \mathbf{v}_\alpha = 0.$$

The capillary pressure follows Leverett's scaling with the correlation used by Udell (his eq. (20)),

$$p_c = p_g - p_l = \sigma \sqrt{\frac{\phi}{K}}\, f(s), \qquad f(s) = 1.417(1-s) - 2.120(1-s)^2 + 1.263(1-s)^3,$$

with the liquid saturation $s = S_l$, and the relative permeabilities are $k_{rl} = s^3$ and $k_{rg} = (1-s)^3$. The temperature of the two-phase state follows from the Clapeyron equation for an ideal-gas vapor with Kelvin's vapor-pressure lowering (Udell's eq. (19)),

$$T = \frac{T_0 \left(1 + \dfrac{p_c}{h_{fg} \varrho_l}\right)}{1 - \dfrac{T_0 R_w}{h_{fg}} \ln\dfrac{p_g}{p_0}},$$

where $R_w$ is the specific gas constant of water and $(T_0, p_0)$ is a reference point on the saturation curve.

**Semi-analytical Reference Solution**

In the steady state, the net mass flux vanishes, $\varrho_l v_l + \varrho_g v_g = 0$, and the energy flux $(h_g - h_l)\varrho_g v_g$ equals the heat flux $q$. The vapor mass flux $q/h_{fg}$ therefore flows towards the cold end at $x = 0$ and the same liquid mass flux flows back. Darcy's law for both phases gives the capillary pressure gradient (Udell's eq. (9))

$$\frac{\mathrm{d} p_c}{\mathrm{d} x} = \frac{q}{K h_{fg}} \left( \frac{\nu_g}{k_{rg}} + \frac{\nu_l}{k_{rl}} \right),$$

with the kinematic viscosities $\nu_\alpha = \mu_\alpha/\varrho_\alpha$. With the saturation $s_c$ at the cold end, the position is a quadrature in the saturation,

$$x(s) = \int_s^{s_c} \frac{-\sigma \sqrt{\phi/K}\, f'(u)}{\dfrac{q}{K h_{fg}} \left( \dfrac{\nu_g}{k_{rg}(u)} + \dfrac{\nu_l}{k_{rl}(u)} \right)} \,\mathrm{d}u,$$

and the gas pressure follows from

$$p_g(s) = p_{g,c} - \int_s^{s_c} \frac{\sigma \sqrt{\phi/K}\, f'(u)}{1 + \beta\, k_{rg}(u)/k_{rl}(u)} \,\mathrm{d}u, \qquad \beta = \frac{\nu_l}{\nu_g},$$

with $p_l = p_g - p_c$ and the temperature from eq. (19). The test evaluates both integrals with `Dumux::integrateScalarFunction` and inverts $x(s)$ by bisection. Over the whole two-phase zone from $s = 0$ to $s = 1$, the same quadrature gives Udell's dimensionless wick limit $\omega \delta_l = \int_0^1 -f'/(1/k_{rg} + \beta/k_{rl})\,\mathrm{d}s$, which is $0.0369$ for $\beta = 0.0146$; Udell reports $0.0368$.

At both ends of the two-phase zone, a relative permeability vanishes and the saturation varies with the fourth root of the distance to the end. The computational domain therefore lies inside the two-phase zone, with $s_c < 1$ at the cold end and a heated end before the dry-out point, so that the solution is smooth on the whole domain.

**Parameters**

| Parameter                          | Symbol     | Value              | Unit     |
|------------------------------------|------------|--------------------|----------|
| Length of the domain               | $L$        | $1$                | m        |
| Permeability                       | $K$        | $10^{-12}$         | m²       |
| Porosity                           | $\phi$     | $0.4$              | -        |
| Surface tension                    | $\sigma$   | $0.0589$           | N/m      |
| Liquid density                     | $\varrho_l$ | $958$             | kg/m³    |
| Liquid dynamic viscosity           | $\mu_l$    | $2.8\cdot10^{-4}$  | Pa·s     |
| Vapor density                      | $\varrho_g$ | $0.6$             | kg/m³    |
| Vapor dynamic viscosity            | $\mu_g$    | $1.2\cdot10^{-5}$  | Pa·s     |
| Kinematic viscosity ratio          | $\beta$    | $0.0146$           | -        |
| Latent heat                        | $h_{fg}$   | $2.257\cdot10^{6}$ | J/kg     |
| Molar mass of water                | $M_w$      | $0.018015$         | kg/mol   |
| Reference saturation temperature   | $T_0$      | $373.15$           | K        |
| Reference saturation pressure      | $p_0$      | $101325$           | Pa       |
| Heat flux at the heated end        | $q$        | $100$              | W/m²     |
| Gas pressure at the cold end       | $p_{g,c}$  | $101325$           | Pa       |
| Liquid saturation at the cold end  | $s_c$      | $0.7$              | -        |

The liquid heat capacity ($4200$ J/(kg K), also used for the vapor), the solid density ($2650$ kg/m³) and the solid heat capacity ($800$ J/(kg K)) only affect the approach to the steady state.

**Setup**

The test uses the two-phase one-component model (@ref TwoPOneCModel) on a one-dimensional `Dune::YaspGrid` with the box scheme (@ref BoxDiscretization) and with TPFA (@ref CCTpfaDiscretization). The coarsest grid has 20 cells. At $x = 0$, Dirichlet conditions fix the gas pressure $p_{g,c}$ and the liquid saturation $s_c$ in the two-phase state. The heated end at $x = L$ is closed for mass and receives the heat flux $q$.

A test-local fluid system on top of `FluidSystems::TwoPOneC` with a constant component implements Udell's assumptions: constant densities and viscosities, the enthalpies $h_l = c_l (T - T_\mathrm{ref})$ and $h_g = h_l + h_{fg}$, and eq. (19) as the temperature of the two-phase state. The `FluidMatrix::HeatPipeLaw` provides Leverett's capillary pressure and the cubic relative permeabilities; its regularization above relative permeability arguments of $0.95$ lies outside the saturation range of the solution. A test-local thermal conductivity model returns zero. Starting from the cold-end state everywhere, implicit Euler steps of growing size reach the steady state at $t = 10^{11}$ s.

**Results**

To run the convergence study and produce the figures, run in the build directory of the test:

```bash
python3 plot_results.py test_2p1cni_heatpipe_box test_2p1cni_heatpipe_tpfa
```

The script runs both schemes with full upwinding of the mobilities (`Flux.UpwindWeight = 1`, the default) and with central weighting (`Flux.UpwindWeight = 0.5`) and produces two figures:

- **`profile.png`**: liquid saturation, temperature, gas pressure and liquid pressure on the coarsest grid compared with the reference solution.
- **`convergence.png`**: discrete $L^2$ errors of the four fields versus the mesh size $h$ on a log-log scale for 20 to 1280 cells.

![Steady profiles](heatpipe2p1c_profile.png)

![Grid convergence](heatpipe2p1c_convergence.png)

With full upwinding, both schemes converge with first order in all four fields. The vapor flows towards the cold end and takes its mobility from the drier neighbor, the liquid flows towards the heated end and takes its mobility from the wetter neighbor. Both mobilities are too large, so the gas pressure rises too little and the liquid pressure drops too little. The capillary pressure is then too small, the saturation too large and the temperature too low. On the coarsest grid, the gas relative permeability doubles across the first cell at the cold end. The errors of the liquid saturation and the temperature are

| cells | saturation, box | rate | saturation, TPFA | rate | temperature (K), box | rate | temperature (K), TPFA | rate |
|------:|----------------:|-----:|-----------------:|-----:|---------------------:|-----:|----------------------:|-----:|
| 20    | 1.86e-2 | -    | 1.71e-2 | -    | 4.76e-2 | -    | 4.33e-2 | -    |
| 80    | 5.13e-3 | 0.93 | 4.98e-3 | 0.89 | 1.32e-2 | 0.93 | 1.27e-2 | 0.88 |
| 320   | 1.32e-3 | 0.98 | 1.31e-3 | 0.96 | 3.40e-3 | 0.98 | 3.37e-3 | 0.96 |
| 1280  | 3.33e-4 | 0.99 | 3.32e-4 | 0.99 | 8.57e-4 | 0.99 | 8.55e-4 | 0.99 |

With central weighting, both schemes converge with second order, so the upwinding is the only first-order error. The errors are

| cells | saturation, box | rate | saturation, TPFA | rate | temperature (K), box | rate | temperature (K), TPFA | rate |
|------:|----------------:|-----:|-----------------:|-----:|---------------------:|-----:|----------------------:|-----:|
| 20    | 1.64e-3 | -    | 1.19e-3 | -    | 4.38e-3 | -    | 3.05e-3 | -    |
| 80    | 1.10e-4 | 1.95 | 9.99e-5 | 1.79 | 2.99e-4 | 1.94 | 2.66e-4 | 1.76 |
| 320   | 6.94e-6 | 2.00 | 6.76e-6 | 1.94 | 1.88e-5 | 2.00 | 1.82e-5 | 1.93 |
| 1280  | 4.34e-7 | 2.00 | 4.31e-7 | 1.99 | 1.18e-6 | 2.00 | 1.17e-6 | 1.98 |

with rates computed between successive rows in both tables. The CTests `test_2p1cni_heatpipe_box` and `test_2p1cni_heatpipe_tpfa` run `convergencetest.py` on 20 to 160 cells with full upwinding and require a mean convergence rate within $1 \pm 0.2$ for the saturation, both pressures and the temperature. The CTests `test_2p1cni_heatpipe_box_central` and `test_2p1cni_heatpipe_tpfa_central` run 20 to 640 cells with central weighting and require a mean convergence rate within $2 \pm 0.2$.
