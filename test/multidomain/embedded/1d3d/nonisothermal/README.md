# Benchmark: Wellbore Heat Transport (1D–3D Embedded) {#benchmark-wellbore-heat-transport}

## Heat Transport in a Fluid-Filled Wellbore Embedded in Rock

**Problem Description**

The benchmark simulates non-isothermal single-phase flow in a vertical wellbore (a borehole
heat exchanger with a grouted pipe) that is embedded as a 1D line in a 3D rock domain.
Cold water is injected at the top of the pipe, is heated on its way down by the surrounding
warm rock, and leaves the pipe at the bottom. The only driving force for the flow is the pressure gradient, gravity is neglected in this setup. The quantities of interest are the outlet
temperature as a function of time and the fluid temperature along the pipe at the final time.

The setup follows the
[wellbore heat transport benchmark](https://www.opengeosys.org/6.5.9/docs/benchmarks/heat-transport-bhe/pipe_flow_ebhe/)
of OpenGeoSys (OGS) and is verified against the analytical wellbore heat-transmission solution
of Ramey @cite Ramey1962, which serves as the reference solution. The numerical results of OGS
are shipped with the test (`ogs_solution.py`), so that DuMux, OGS and the analytical solution
can be compared in a single plot.

Two coupling methods for the 1D–3D exchange are compared, the surface coupling and the
kernel coupling (see @ref EmbeddedCoupling). Both are built from the same `main.cc`. At
refinement 3 they are converged (Figure 3) and agree with each other (Figures 1 and 2).

**Model**

The 1D pipe domain $\Lambda$ uses the non-isothermal one-phase model (@ref OnePModel with
@ref NIModel), i.e. a mass and an energy balance. With gravity switched off they read

$$\frac{\partial (\phi \varrho_f)}{\partial t}
- \nabla\cdot\left(\varrho_f \frac{\mathbf{K}}{\mu_f}\nabla p\right) = q_m$$

$$\frac{\partial (\phi \varrho_f u_f)}{\partial t}
- \nabla\cdot\left(\varrho_f h_f \frac{\mathbf{K}}{\mu_f}\nabla p\right)
- \nabla\cdot\left(\lambda_f\nabla T\right) = q_e$$

where $u_f$ and $h_f$ are the specific internal energy and enthalpy of the fluid. $q_m$ and $q_e$
are mass and energy sources; there is no mass exchange ($q_m = 0$), and $q_e$ holds the exchange
term $q_h$ of the coupling. All terms are extruded with the pipe cross-section $\pi r_i^2$.
The pipe is modelled as a porous medium with porosity $\phi = 1$, so only the fluid stores and
conducts heat, and a Hagen–Poiseuille-equivalent permeability $\mathbf{K} = r_i^2/8$
(`spatialparams_well.hh`), so that the Darcy flux reproduces laminar pipe flow.

The 3D rock domain $\Omega$ is impermeable, so heat is transported by conduction only. It uses
the energy balance of a solid (@ref SolidEnergyModel)

$$\frac{\partial (\varrho_s c_s T)}{\partial t} - \nabla\cdot\left(\lambda_s\nabla T\right) = -q_e$$

where the source $-q_e$ is the exchange term of the coupling with the opposite sign. The model
scales storage and conduction with the solid fraction $1 - \phi$, which is one as the rock has
the porosity $\phi = 0$ (`spatialparams_rock.hh`).

**Coupling**

The two domains exchange energy only, there is no mass exchange. The exchange term per unit
length of the wellbore is evaluated at the integration points of the 1D–3D coupling
(`exchangefluxcalculator.hh`):

$$q_h = 2\pi r_o\, U_o\, \left(T_\Lambda - T_\Omega\right)\, \Xi\ \qquad [\mathrm{W/m}]$$

where $T_\Lambda$ is the fluid temperature in the pipe, $T_\Omega$ the rock temperature
seen by the coupling and $\Xi$ a flux scaling factor that is only present in the kernel coupling method. The effective heat transfer coefficient $U_o$ combines the fluid-side
convective resistance and the conductive resistance of the pipe wall and grout in series, both
referred to the outer radius $r_o$:

$$\frac{1}{U_o} = \frac{r_o}{r_i\, h_c} + \frac{r_o \ln(r_o/r_i)}{\lambda_{pg}},
\qquad h_c = \frac{\lambda_f\, \mathrm{Nu}}{2 r_i}$$

where $h_c$ is the convective heat transfer coefficient of the pipe flow.

With the flow rate of this benchmark the pipe Reynolds number is $\mathrm{Re} \approx 865$, so
the flow is laminar and the Nusselt number of a pipe with constant wall heat flux,
$\mathrm{Nu} = 4.364$, is used. $\lambda_{pg} = 0.962\ \mathrm{W/(m\,K)}$
(`MixedDimension.PipeWallThermalConductivity`) is the equivalent conductivity of the pipe wall
and the grout annulus, chosen such that $\ln(r_o/r_i)/\lambda_{pg}$ equals the series
resistance $\ln(r_{po}/r_i)/\lambda_{pi} + \ln(r_o/r_{po})/\lambda_g$ of the OGS setup, with the
outer pipe radius $r_{po} = r_i + t_{pi}$.

The coupling method is selected at compile time by the definition `COUPLINGMODE`
(`properties.hh`); the coupling manager (`wellborecouplingmanager.hh`) extends
@ref Dumux::Embedded1d3dCouplingManager by an `innerRadius()` accessor.

- *Surface coupling* (`Embedded1d3dCouplingMode::Surface`): the exchange term of each 1D
  integration point is carried by `MixedDimension.NumCircleSegments` point sources placed on
  the borehole circumference of radius $r_o$. The rock temperature $T_\Omega$ is the average
  over these points, and each point source puts its share of the exchange term into the rock
  cell it falls into. The benchmark uses 400 segments.
- *Kernel coupling* (`Embedded1d3dCouplingMode::Kernel`): the exchange term is distributed
  in the rock over a kernel support of width $w = 2 r_o = 0.28\ \mathrm{m}$
  (`MixedDimension.KernelWidthFactor = 2`) around the wellbore (`problem_rock.hh`, `source()`);
  the kernel weights are integrated by the coupling manager with
  `MixedDimension.KernelIntegrationCRL = 0.001`.
  Since the rock temperature seen by the coupling is then not the temperature at the borehole
  wall, the exchange term is scaled such that it matches the radially symmetric steady
  solution around the wellbore @cite Koch2020 :

  $$\Xi = \begin{cases}
  \left[1 + \theta\left(\dfrac{d^2}{2w^2} + \ln\dfrac{w}{r_o} - \dfrac{1}{2}\right)\right]^{-1}
  & d \leq w \\[4pt]
  \left[1 + \theta \ln\dfrac{d}{r_o}\right]^{-1} & d > w
  \end{cases}, \qquad
  \theta = \frac{r_o U_o}{\lambda_s}$$

  where $d$ is the average distance between the wellbore segment and the rock cell
  (`ExchangeFluxCalculator::energyFluxScalingFactor`). This scaling assumes a steady radial
  temperature profile within the kernel support. The simulated period is short, so the rock
  around the wellbore does not reach a steady state. The kernel width is therefore chosen small
  to keep the error of this assumption minimal.

**Parameters**

| Parameter | Symbol | Value | Unit |
|-----------|--------|-------|------|
| Borehole length | $L$ | 30 | m |
| Pipe inner radius | $r_i$ | 0.12913 | m |
| Borehole (outer) radius | $r_o$ | 0.14 | m |
| Pipe wall thickness | $t_{pi}$ | 0.00587 | m |
| Volumetric flow rate | $q$ | $2\times10^{-4}$ | m³/s |
| Inlet temperature | $T_\mathrm{in}$ | 20 | °C |
| Undisturbed rock temperature | $T_{s,0}$ | 55 | °C |
| Fluid density | $\varrho_f$ | 1000 | kg/m³ |
| Fluid heat capacity | $c_{p,f}$ | 4190 | J/(kg·K) |
| Fluid dynamic viscosity | $\mu_f$ | $1.14\times10^{-3}$ | Pa·s |
| Fluid thermal conductivity | $\lambda_f$ | 0.59 | W/(m·K) |
| Rock density | $\varrho_s$ | 1800 | kg/m³ |
| Rock heat capacity | $c_s$ | 1778 | J/(kg·K) |
| Rock thermal conductivity | $\lambda_s$ | 2.78018 | W/(m·K) |
| Pipe wall thermal conductivity | $\lambda_{pi}$ | 1.3 | W/(m·K) |
| Grout thermal conductivity | $\lambda_g$ | 0.73 | W/(m·K) |
| Equivalent conductivity of pipe wall and grout | $\lambda_{pg}$ | 0.962 | W/(m·K) |
| Simulation time | $t_\mathrm{end}$ | 5 | d |

The fluid properties are constant (`[Component]` in `params.input`), matching the OGS benchmark
setup. The radii
are read per element from the parameters of `grids/pipe.dgf`.

**Analytical Solution (Ramey 1962)**

Ramey's wellbore heat-transmission solution (`ramey.py`) assumes steady pipe flow, negligible
axial conduction and heat storage in the fluid, and radial heat conduction into an infinite
formation described by a line-source time function. The fluid temperature at distance $z$ from
the inlet is

$$T_f(z, t) = T_{s,0} + (T_\mathrm{in} - T_{s,0})\,\exp\left(-\frac{z}{X}\right)$$

with the characteristic length

$$X(t) = \frac{q\,\varrho_f\, c_{p,f}\left(\lambda_s + r_i\, U_i\, f(t_D)\right)}
{2\pi\, r_i\, U_i\, \lambda_s}$$

where the dimensionless time and the line-source time function are

$$t_D = \frac{\lambda_s\, t}{\varrho_s\, c_s\, r_o^2}, \qquad
f(t_D) = \begin{cases}
1.1281\sqrt{t_D}\left(1 - 0.3\sqrt{t_D}\right) & t_D \leq 1.5 \\[4pt]
\left(0.4063 + 0.5\ln t_D\right)\left(1 + \dfrac{0.6}{t_D}\right) & t_D > 1.5
\end{cases}$$

$U_i$ and $X$ must be referred to the same radius, here the inner radius $r_i$:

$$\frac{1}{U_i} = \frac{1}{h_c} + \frac{r_i \ln(r_o/r_i)}{\lambda_{pg}}$$

Note: The analytical solution on the OGS benchmark page refers the heat transfer coefficient to the
outer pipe radius $r_{po}$ but uses $r_i$ in $X$. The consistent formulation is used here; it
brings the analytical solution closer to both numerical models.

**Setup**

The wellbore is vertical and lies on the $z$-axis: the inlet is at the surface ($z = 0$) and
the bottom of the borehole, where the outlet temperature is recorded, is at
$z = -30\ \mathrm{m}$. The distance along the borehole used in the plots and in the analytical
solution is therefore $-z$. The rock domain is the box
$[-5, 5] \times [-5, 5] \times [-30, 0]\ \mathrm{m}$ around it.

The rock domain uses cell-centered two-point flux approximation on a YaspGrid, the pipe the
box method on a FoamGrid. Both domains are assembled and solved monolithically. The rock grid depends on the coupling
method; `params.input` contains both grid specifications in `[Rock.Grid]`:

- *Surface coupling*: `TensorProductCoordinates` with $34 + 34$ cells in each lateral
  direction, graded towards the borehole axis with a factor of 1.2 (cell widths from 2.0 mm at
  the axis and 26 mm at the borehole wall to 0.84 m at the lateral boundary), and a single cell
  over the depth.
- *Kernel coupling*: `EquidistantOffsetCoordinates` (compile definition `COORDINATETYPE`) with
  $68 \times 68 \times 1$ cells, i.e. a lateral cell width of 0.147 m.

The pipe grid `grids/pipe.dgf` consists of a single 30 m element.

A plain run of the executables, as in the CTest, uses these grids (refinement 0). The
convergence script refines them by passing the refinement as command line arguments (see Grid
Convergence below).

Boundary conditions:

| Domain | Boundary | Condition |
|--------|----------|-----------|
| Pipe | inlet ($z = 0$) | Neumann: mass influx $q\varrho_f/(\pi r_i^2)$ and the corresponding enthalpy influx at $T_\mathrm{in}$ |
| Pipe | outlet ($z = -L$) | Dirichlet pressure $10^5$ Pa; Neumann energy: advective outflow plus Fourier conduction evaluated from the local gradients |
| Rock | top and bottom ($z = 0$, $z = -L$) | no-heat-flux Neumann |
| Rock | lateral faces | Dirichlet at the initial state ($T_{s,0}$) |

The lateral extent of the rock domain is chosen such that the cooling around the wellbore does
not reach the lateral boundaries within the simulated time. The rock temperature next to the
lateral faces stays within 0.001 °C of $T_{s,0}$. There is thus no temperature gradient at these
faces, and the Dirichlet condition does not bias the solution.

Initial conditions are $T = T_\mathrm{in} = 20\ ^\circ\mathrm{C}$ (`Well.Problem.InjectionTemperature`),
$p = 10^5\ \mathrm{Pa}$ in the pipe and $T = T_{s,0} = 55\ ^\circ\mathrm{C}$
(`Rock.Problem.InitialTemperature`) in the rock. The simulation runs for
5 days. It starts with a time step of 1 s (`TimeLoop.DtInitial`), so that the flushing of the
cold initial pipe content is resolved in time; the time step then grows up to the episode length
of 4 h (`TimeLoop.EpisodeLength`). The outlet temperature and the solutions of both domains are
written at the end of every episode.
The output of the 3D rock solution can be switched off with `Vtk.EnableBulkOutput = false`.

**Grid Convergence**

`convergence.py` runs both coupling methods on a sequence of grids and plots them together
with Ramey's solution and the OGS results. Refinement $n$ means $68 \cdot 2^n \times 68 \cdot 2^n
\times 2^n$ rock cells and $2^n$ pipe elements for both methods. The script starts from the
grids in `params.input` and refines them by overriding parameters on the command line; these
parameters are not set in `params.input`. For both methods it passes `-Well.Grid.Refinement n`.
For the rock grid of the surface coupling it passes `-Rock.Grid.Refinement n`. The kernel
coupling does not support refined rock grids, so for it the script instead passes
`-Rock.Grid.Cells` with the cells of `params.input` multiplied by $2^n$.

| Refinement | Rock cells | Pipe elements | Lateral rock cell width at the borehole wall, surface / kernel |
|------------|------------|---------------|-----------------------------------------------------------------|
| 0 | $68 \times 68 \times 1$ | 1 | 26.1 mm / 147.1 mm |
| 1 | $136 \times 136 \times 2$ | 2 | 13.1 mm / 73.5 mm |
| 2 | $272 \times 272 \times 4$ | 4 | 6.5 mm / 36.8 mm |
| 3 | $544 \times 544 \times 8$ | 8 | 3.3 mm / 18.4 mm |

**Results**

![Figure 1: Outlet temperature over time](wellbore_outlet_temperature.png)

![Figure 2: Temperature profile along the borehole at the final time](wellbore_temperature_profile.png)

Figures 1 and 2 show both coupling methods on the finest grid (refinement 3) together with
Ramey's solution (dotted) and the OGS results. Figure 1 shows the outlet temperature at a depth
of 30 m over the 5 days of simulated time, Figure 2 the fluid temperature along the borehole at
the final time $t = 5$ days. In both figures the upper panel shows the temperatures and the
lower panel their deviation from Ramey's solution. The outlet temperature and its deviation are
shown at all outputs except $t = 0$, where the pipe is still filled with the cold initial water
that Ramey's solution does not contain.

The kink in Ramey's outlet temperature at $t \approx 9.4$ h comes from the time function: there
$t_D = 1.5$, where $f(t_D)$ switches from the early-time to the late-time approximation. The two
branches do not match ($f$ drops from 0.874 to 0.853), so the outlet temperature jumps by about
0.04 °C.

*Grid convergence.*

![Figure 3: Change with regard to the previous refinement](wellbore_refinement_convergence.png){html: width=70%}

Figure 3 shows the change of the solution with regard to the previous refinement: for the
outlet temperature the maximum over all outputs in time, for the temperature profile at 5 d the
maximum along the borehole. Ramey's solution is only approximate, so the
convergence is measured between successive grids. The changes decrease with every refinement,
for the profile at about second order, except for the outlet temperature of the kernel coupling
(see below). Going from refinement 2 to 3 changes the outlet temperature by at most 0.009 °C
(surface, at 8 h) and 0.054 °C (kernel, at 4 h); from 8 h on the change is below 0.01 °C for both
methods, and after the first day below 0.001 °C. The profile changes by up to 0.006 °C for both
methods. Refinement 3
is therefore considered converged; on this grid the solutions of both methods lie on top of each
other (Figures 1 and 2).

At the first output ($t = 4$ h) the kernel coupling still gives an outlet temperature 0.07 °C
below the surface coupling. At this time the thermal front in the rock has only penetrated about
0.1 m, less than the kernel width of 0.28 m, so the steady radial profile assumed by the kernel
flux scaling does not hold yet. This early-time error sets the maximum change of the kernel outlet
temperature in Figure 3, which therefore does not decrease from refinement 2 to 3. It is also the
reason for the small kernel width, see the coupling section above.

*Comparison with Ramey and OGS.* The residence time of the water in the pipe is 2.2 h, so at the
first output ($t = 4$ h) the cold initial pipe content has left the pipe. The outlet temperature
is highest there (26.87 °C surface, 26.81 °C kernel) and then decays as the rock around the
borehole cools down and the temperature difference driving the heat exchange shrinks
(Figure 1). During the first day DuMux is warmer than Ramey's solution, by up to 0.44 °C
(surface) and 0.37 °C (kernel) at 4 h. OGS lies 0.76 °C below Ramey at 4 h and up to 0.6 °C
above it later in the first day. From day 2 on the deviation stays below 0.04 °C for DuMux and
below 0.07 °C for OGS.

The temperature profile (Figure 2) is the expected exponential approach of the fluid temperature towards
the undisturbed rock temperature. With a characteristic length $X \approx 208$ m it is still
practically linear over the 30 m of this benchmark and rises from 20 °C at the inlet to
24.674 °C (surface) and 24.675 °C (kernel) at the outlet, against 24.702 °C for Ramey and
24.670 °C for OGS. Over the whole borehole the deviation from Ramey is at most 0.028 °C for both
methods, against
0.032 °C for OGS.

**Running the Benchmark**

Navigate to the source directory of the benchmark and run the convergence test:

```bash
cd dumux/test/multidomain/embedded/1d3d/nonisothermal
python3 convergence.py
```

The script builds both executables and runs every coupling method at refinements 0 to 3. The
output is saved in the matching build directory,
`dumux/build-cmake/test/multidomain/embedded/1d3d/nonisothermal/convergence`: each run in its own
subdirectory `<method>_refinement_<n>`, and the figures as `outlet_temperature.png`,
`temperature_profile.png` and `refinement_convergence.png` (Figures 1 to 3). The 3D rock output
is not used and would take several GB per run on the fine grids, so the script switches it off
(`-Vtk.EnableBulkOutput false`).

The script needs the Python packages `numpy`, `matplotlib` and `vtk`, which can be installed
e.g. with `python3 -m pip install numpy matplotlib vtk`.

**References**

- Ramey, H. J. (1962): *Wellbore Heat Transmission*. Journal of Petroleum Technology 14(4),
  427–435. [doi:10.2118/96-PA](https://doi.org/10.2118/96-PA)
- Koch, T., Schneider, M., Helmig, R., Jenny, P. (2020): *Modeling tissue perfusion in terms of
  1d-3d embedded mixed-dimension coupled problems with distributed sources*. Journal of
  Computational Physics 410, 109370.
  [doi:10.1016/j.jcp.2020.109370](https://doi.org/10.1016/j.jcp.2020.109370)
- OGS results: OpenGeoSys benchmark documentation,
  <https://www.opengeosys.org/6.5.9/docs/benchmarks/heat-transport-bhe/pipe_flow_ebhe/>
