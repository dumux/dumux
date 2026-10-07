# Benchmark: Wellbore Heat Transport (1D–3D Embedded) {#benchmark-wellbore-heat-transport}

## Heat Transport in a Fluid-Filled Wellbore Embedded in Rock

**Problem Description**

The benchmark simulates non-isothermal single-phase flow in a vertical wellbore (a borehole
heat exchanger with a grouted pipe) that is embedded as a 1D line element in a 3D rock domain.
Cold water is injected at the top of the pipe, is heated on its way down by the surrounding
warm rock, and leaves the pipe at the bottom. The quantity of interest is the fluid temperature
along the pipe and, in particular, the outlet temperature as a function of time.

The setup follows a [wellbore/borehole-heat-exchanger benchmark](https://ogs.ogs.xyz/ogs/docs/benchmarks/heat-transport-bhe/pipe_flow_ebhe/) of OpenGeoSys (OGS) and is verified against the classical
analytical wellbore heat-transmission solution of Ramey (1962). The OGS results are shipped
as reference data inside the plotting script, so that DuMux, OGS and the analytical solution
can be compared in a single plot.

Two balance equations are solved in each subdomain (@ref OnePModel in its non-isothermal
variant `OnePNI`, i.e. one-phase flow with an energy balance). In the 1D pipe domain $\Lambda$
the mass and energy balances read

$$\frac{\partial (\phi \varrho)}{\partial t} - \nabla\cdot\left(\varrho \frac{k}{\mu}\nabla p\right) = q_m$$

$$\frac{\partial (\phi \varrho u)}{\partial t}
- \nabla\cdot\left(\varrho h \frac{k}{\mu}\nabla p + \lambda_\mathrm{eff}\nabla T\right) = q_e$$

where all terms are extruded with the pipe cross-section $\pi r_i^2$. The pipe is modelled as a
porous medium with porosity $\phi = 1$ and a Hagen–Poiseuille-equivalent permeability
$k = r_i^2 \varrho$ (`spatialparams_voids.hh`), so that the Darcy flux reproduces laminar pipe flow.

In the 3D rock domain $\Omega$ the same model is used with porosity $\phi = 0$ and
$k = 10^{-20}\,\mathrm{m^2}$, i.e. the rock is effectively impermeable and transports heat by
conduction only. Gravity is switched off.

**Coupling**

The two domains exchange energy through a source/sink term that is evaluated at the
integration points of the 1D–3D coupling (`exchangefluxcalculator.hh`):

$$q_h = 2\pi r_o\, U\, \left(T_\Lambda - T_\Omega\right) \qquad [\mathrm{W/m}]$$

where $T_\Lambda$ is the fluid temperature in the pipe and $T_\Omega$ the rock temperature at
the position of the point source. In the surface coupling mode the exchange term of each 1D
integration point is carried by `MixedDimension.NumCircleSegments` individual point sources
placed on the borehole circumference; each of them evaluates the local rock temperature of the
bulk element it falls into and carries the corresponding fraction of the exchange surface.
The effective heat transfer coefficient $U$ combines the
fluid-side convective resistance and the conductive resistance of the pipe wall and grout in
series, both referred to the outer radius $r_o$:

$$\frac{1}{U} = \frac{r_o}{r_i\, h} + \frac{r_o \ln(r_o/r_i)}{\lambda_{pg}},
\qquad h = \frac{\lambda_f\, \mathrm{Nu}}{2 r_i}$$

With the flow rate of this benchmark the pipe Reynolds number is $\mathrm{Re} \approx 865$, so
the flow is laminar and the Nusselt number of a pipe with constant wall heat flux,
$\mathrm{Nu} = 4.364$, is used. $\lambda_{pg} = 0.962\ \mathrm{W/(m\,K)}$
(`MixedDimension.PipeWallThermalConductivity`) is the equivalent conductivity of the pipe wall
and the grout annulus, chosen such that $\ln(r_o/r_i)/\lambda_{pg}$ equals the series
resistance $\ln(r_{po}/r_i)/\lambda_{pi} + \ln(r_o/r_{po})/\lambda_g$ of the OGS setup.

There is no mass exchange between pipe and rock; the coupling is purely thermal.

**Parameters**

| Parameter | Symbol | Value | Unit |
|-----------|--------|-------|------|
| Borehole length | $L$ | 30 | m |
| Pipe inner radius | $r_i$ | 0.12913 | m |
| Borehole (outer) radius | $r_o$ | 0.14 | m |
| Pipe wall thickness | $t_{pi}$ | 0.00587 | m |
| Volumetric flow rate | $q$ | $2\times10^{-4}$ | m³/s |
| Inlet temperature | $T_i$ | 20 | °C |
| Undisturbed rock temperature | $T_d$ | 55 | °C |
| Fluid density | $\varrho_f$ | 1000 | kg/m³ |
| Fluid heat capacity | $c_{p,f}$ | 4190 | J/(kg·K) |
| Fluid dynamic viscosity | $\mu_f$ | $1.14\times10^{-3}$ | Pa·s |
| Fluid thermal conductivity | $\lambda_f$ | 0.59 | W/(m·K) |
| Rock density | $\varrho_{re}$ | 1800 | kg/m³ |
| Rock heat capacity | $c_{p,re}$ | 1778 | J/(kg·K) |
| Rock thermal conductivity | $\lambda_{re}$ | 2.78018 | W/(m·K) |
| Pipe wall thermal conductivity | $\lambda_{pi}$ | 1.3 | W/(m·K) |
| Grout thermal conductivity | $\lambda_g$ | 0.73 | W/(m·K) |
| Equivalent conductivity of pipe wall and grout | $\lambda_{pg}$ | 0.962 | W/(m·K) |
| Simulation time | $t_\mathrm{end}$ | 5 | d |

The fluid properties are constant (`simpleh2o.hh`), matching the OGS reference setup.

**Analytical Solution (Ramey 1962)**

Ramey's wellbore heat-transmission solution assumes steady pipe flow, negligible axial
conduction in the fluid and radial heat conduction into an infinite formation described by a
line-source time function. The fluid temperature at distance $z$ from the inlet is

$$T_f(z, t) = T_d + (T_i - T_d)\,\exp\left(-\frac{z}{X}\right)$$

with the characteristic length

$$X(t) = \frac{q\,\varrho_f\, c_{p,f}\left(\lambda_{re} + r_i\, U\, f(t_D)\right)}
{2\pi\, r_i\, U\, \lambda_{re}}$$

where the dimensionless time and the line-source time function are

$$t_D = \frac{\lambda_{re}\, t}{\varrho_{re}\, c_{p,re}\, r_o^2}, \qquad
f(t_D) = \begin{cases}
1.1281\sqrt{t_D}\left(1 - 0.3\sqrt{t_D}\right) & t_D \leq 1.5 \\[4pt]
\left(0.4063 + 0.5\ln t_D\right)\left(1 + \dfrac{0.6}{t_D}\right) & t_D > 1.5
\end{cases}$$

Both $U$ and the characteristic length $X$ must be referred to the same radius for the
solution to be consistent. Here the inner radius $r_i$ is used, i.e.

$$\frac{1}{U} = \frac{1}{h} + r_i\left(\frac{\ln(r_{po}/r_i)}{\lambda_{pi}}
+ \frac{\ln(r_o/r_{po})}{\lambda_g}\right), \qquad r_{po} = r_i + t_{pi}$$

**Setup**

The wellbore is vertical and aligned with the $z$-axis: the inlet is at the surface
($z = 0$) and the bottom of the borehole, where the outlet temperature is recorded, is at
$z = -30\ \mathrm{m}$. The depth along the borehole used in the plots and in the analytical
solution is therefore $-z$.

The 3D rock domain is discretized on a YaspGrid with `TensorProductCoordinates` spanning
$10\ \mathrm{m} \times 10\ \mathrm{m} \times 30\ \mathrm{m}$
($[-5, 5] \times [-5, 5] \times [-30, 0]$) with $20 \times 20 \times 10$ cells that are graded
towards the borehole axis in the two lateral directions with a grading factor of $\pm 1.3$.
`Soil.Grid.Refinement = 2` refines this grid twice, so the simulation runs on
$80 \times 80 \times 40$ cells.

Two discretizations of the rock domain are available. They are built as separate executables from
the same `main.cc`, selected by the compile definition `BULKTYPETAG` (`properties.hh`):
`test_wellbore_heat_transport_tpfa` uses cell-centered two point flux approximation and
`test_wellbore_heat_transport_box` the box method. `convergence_test.py` runs both and draws them
in one figure.

The 1D pipe uses the box method on a FoamGrid specified in `grids/pipe.dgf` with 10 elements of
3 m length, refined twice by `Voids.Grid.Refinement = 2` to 40 elements of 0.75 m. The inner and
outer radii are read per element from the DGF parameters. Both domains are assembled and solved
monolithically with the AMG-preconditioned BiCGSTAB block-diagonal solver.

Boundary conditions:

| Domain | Boundary | Condition |
|--------|----------|-----------|
| Pipe | inlet ($z = 0$) | Neumann: mass influx $q\varrho/(\pi r_i^2)$ and the corresponding enthalpy influx at $T_i$ |
| Pipe | outlet ($z = -L$) | Dirichlet pressure $10^5$ Pa; energy: advective outflow plus Fourier conduction evaluated from the local gradients |
| Rock | top and bottom ($z = 0$, $z = -L$) | no-flow / no-heat-flux Neumann |
| Rock | lateral faces | Dirichlet at the initial state ($T_d$, $10^5$ Pa) |

Initial conditions are $T = T_i = 20\ ^\circ\mathrm{C}$, $p = 10^5\ \mathrm{Pa}$ in the pipe and
$T = T_d = 55\ ^\circ\mathrm{C}$, $p = 10^5\ \mathrm{Pa}$ in the rock. The simulation runs for
5 days with a time step of 2 h and output at every time step.

The coupling manager (`wellborecouplingmanager.hh`) extends @ref Dumux::Embedded1d3dCouplingManager by an
`innerRadius()` accessor. 

**Results**

![Outlet temperature over time and temperature profile along the borehole](wellbore_temperatures.png)

The figure shows both discretizations of the rock domain at refinement 3, i.e. on
$160 \times 160 \times 80$ rock cells and 80 pipe elements, together with Ramey's solution and the
OGS reference. The left panel shows the outlet temperature at a depth of 30 m over the 5 days of
simulated time, the right panel the fluid temperature along the borehole at the final time
$t = 5$ days; the dashed curves on the secondary axes are the absolute errors with respect to
Ramey. It is produced by `convergence_test.py --no-run --refinements 3` once the runs of that
refinement are there.

Starting from the initial fluid temperature of $T_i = 20\ ^\circ\mathrm{C}$, the outlet temperature
rises steeply while the first cold water travels down the borehole, peaks at about 26.6 °C within
the first half day and then decays again as the rock immediately around the borehole cools down and
the temperature difference driving the heat exchange shrinks. During the first day all results lie
above Ramey, by up to 0.38 °C (tpfa), TODO °C (box) and 0.59 °C (OGS), because the line-source time
function $f(t_D)$ is an early-time approximation and the analytical solution neglects the thermal
capacity of the fluid in the borehole. From 2.5 days on the deviation stays below 0.05 °C.

The temperature profile is the expected exponential approach of the fluid temperature towards the
undisturbed rock temperature; over the 30 m of this benchmark, with a characteristic length $X$ of
about 200 m, it is still practically linear and rises from 20 °C at the inlet to 24.69 °C (tpfa) and
TODO °C (box) at the outlet, against 24.70 °C for Ramey and 24.67 °C for OGS. Over the whole
borehole the deviation from Ramey is at most 0.011 °C for tpfa and TODO °C for box, against
0.032 °C for OGS.

<!-- TODO: fill in once box refinement 3 has finished, from
     python3 convergence_test.py --no-run --refinements 3
     (1) max deviation above Ramey during the first day, box
     (2) outlet temperature at t = 5 d, box refinement 3
     (3) max |deviation| from Ramey along the borehole, box
     and add the figure as wellbore_temperatures.png -->


To reproduce the results, build both executables and run the convergence test from the build
directory:

```bash
cd <build-dir>/test/multidomain/embedded/1d3d/nonisothermal
make test_wellbore_heat_transport_tpfa test_wellbore_heat_transport_box
python3 convergence_test.py --refinements 0 1 2
```

The script runs every discretization at every refinement, each in its own directory
`convergence/<variant>_refinement_<refinement>`, and writes the figure to
`convergence/convergence_temperatures.png`. It has to be started from the build directory, where
`params.input` and the executables live. `--reuse` skips runs whose results are already there,
`--no-run` only redraws the figure from existing results, and `python3 convergence_test.py --help`
lists all options.

A single run on the grid given in `params.input` is
`./test_wellbore_heat_transport_tpfa params.input`.
The 1D problem writes the outlet temperature at every checkpoint to
`test_wellbore_heat_transport_1d.csv` (columns: time in days, temperature in °C) and the full
solution to VTK files (`*_1d-*.vtp`, `*_3d-*.vtu`).

**References**

- Ramey, H. J. (1962): *Wellbore Heat Transmission*. Journal of Petroleum Technology 14(4),
  427–435. [doi:10.2118/96-PA](https://doi.org/10.2118/96-PA)
- Reference data: OpenGeoSys benchmark documentation, <https://www.opengeosys.org/6.5.9/docs/benchmarks/heat-transport-bhe/pipe_flow_ebhe/>
