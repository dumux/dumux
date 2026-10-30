# Semi-infinite 1D Heat Conduction in a Porous Medium {#benchmark-1pni-conduction}

## Problem Description

This benchmark verifies the non-isothermal model for transient heat conduction in a fully saturated porous medium.
Physically the setup represents a semi-infinite, one dimensional porous medium that is suddenly exposed to an elevated surface temperature.

**Governing Equations**

Assumptions:

- local thermal equilibrium
- one-dimensional conduction in the $x-$direction,
- incompressible fluid $f$ and constant soil properties,
- single phase flow,
- no gravity
- negligible advective heat transport (constant pressure of 1 bar),
- no source terms in the domain,
- homogeneous domain (constant porosity)
- fixed temperatures at left and right boundary

The benchmark solves the transient energy conservation equation for a porous medium:

With those assumptions, we solve the following energy balance for the porous meidum:

$$ \phi \rho_f\frac{\partial u_f}{\partial t} + (1-\phi) \rho_s \frac{\partial u_s}{\partial t}
= \nabla \cdot \left(\lambda_{\text{eff}} \nabla T\right) \, ,
$$
where
- $\phi$ is porosity,
- $u$ is the internal energy,
- $T$ the overall porous medium temperature
- $\rho_f$, $\rho_s$ is the density of the fluid/solid,
- $c_f$, $c_s$ is the mass-specific heat capacities of the fluid and solid phase,
- $\lambda_{\text{eff}}$ the effective thermal conductivtiy of the porous medium.


The effective thermal conductivity of the porous medium is hereby obtained as a volume-fraction average of the liquid and solid thermal conductivities with
$$ \lambda_{\text{eff}} = \phi\,\lambda_f + (1-\phi)\,\lambda_s  \; .$$

Reformulating the energy equation in terms of temperature and simplifying it to one dimension (here $x$), this can be written as:
$$ \phi \rho_f\frac{\partial(c_f T)}{\partial t} + (1-\phi) \rho_s c_s \frac{\partial T}{\partial t}
= \frac{\partial}{\partial x} \left(\lambda_{\text{eff}} \frac{\partial T}{\partial x} \right) \; .
$$
For the analytical solution, we further assume constant material properties for the fluid, such as the thermal heat capacity and the fluid thermal conductivity. The heat conduction equation can then be reformulated to:
$$ \left(\phi \rho_f c_f + (1-\phi) \rho_s c_s \right)\frac{\partial T}{\partial t}
= \lambda_{\text{eff}} \,\frac{\partial^2 T}{\partial x^2} \; .
$$

The first term can hereby be combined into the total volumetric heat capacity

$$
C_{\text{tot}} = \phi\,\rho_f c_f + (1-\phi)\,\rho_s c_s \; .
$$

## Analytical Solution

For the semi-infinite diffusion problem with a constant surface temperature, the analytic temperature profile is (Eq. 9.62 in @cite Poirier2016)

$$
T(x,t) = T_{\text{high}} + (T_{\text{init}} - T_{\text{high}}) \operatorname{erf}\left(\frac{x}{2\sqrt{\alpha t}}\right),
$$

where the thermal diffusivity is

$$
\alpha = \frac{\lambda_{\text{eff}}}{C_{\text{tot}}}.
$$

The implementation in the benchmark uses the equivalent form

$$
T(x,t) = T_{\text{high}} + (T_{\text{init}} - T_{\text{high}})
\operatorname{erf}\left(0.5 \sqrt{\frac{x^2 C_{\text{tot}}}{t\,\lambda_{\text{eff}}}}\right).
$$

<!-- The exact solution is computed in `OnePNIConductionProblem::updateExactTemperature()` and exported as the VTK field `temperatureExact`. -->


## Parameters

| Parameter                                    | Symbol    | Value                   | Unit  |
|----------------------------------------------|-----------:|:------------------------|:------|
| Porosity                                     | $\phi$    | 0.4                     | -     |
| Initial temperature                          | $T_{\mathrm{init}}$ | 290                     | K     |
| Surface (left) temperature                   | $T_{\mathrm{high}}$ | 300                     | K     |
| Pressure (initial & boundary)                   | $p_{\mathrm{init}}$       | 1e5                     | Pa    |
| Solid density                                | $\rho_s$ | 2700                    | kg m$^{-3}$ |
| Solid thermal conductivity                    | $\lambda_s$ | 2.8                   | W m$^{-1}$ K$^{-1}$ |
| Solid heat capacity                           | $c_s$     | 790                     | J kg$^{-1}$ K$^{-1}$ |

Further, the fluid is water. The fluid properties (density, thermal conductivity and heat capacity) are temperature- and pressure-dependent and are evaluated using the IAPWS Industrial Formulation 1997 for the thermodynamic properties of water and steam @cite IAPWS1997.
In DuMux, it is implemented by the component `Dumux::Components::H2O` (@ref Components), which builds on the region formulations in `dumux/material/components/iapws/` (@ref IAPWS).

Although the analytical derivation assumes constant fluid properties, the analytical profile used for comparison is evaluated with fluid properties depending on the pressure and temperature of the numerical simulation for each evaluated time $t$.


##  Simulation Setup

We solve the heat conduction equation using the single phase, non-isothermal model in DuMux @ref OnePModel @ref NIModel. For the spatial discretization we choose a finite volume discretization with a cell-centered TPFA (@ref CCTpfaDiscretization) or MPFA (@ref CCMpfaDiscretization) scheme, or use the vertex-centered box scheme (@ref BoxDiscretization).

**Simulation Domain**
The computational domain is represented by a two-dimensional `Dune::YaspGrid` with $5 \times 1$ m with $200 \times 1$ discretization cell. Hence, this is a finite quasi-1D domain. For the chosen simulation times the thermal front remains sufficiently far from the right boundary, so the finite-domain is suitable to be compared to the semi-infinite analytical solution.

**Initial and Boundary Conditions**
Initially, the temperature and pressure are set to:
$$
T(x,0)=T_{\text{init}} = 290\ \text{K} \, , \\
p(x,0)=p_{\text{init}} = 10^5\ \text{Pa}.
$$

At the left boundary the tmeperature is set to
$$
T(0,t) = T_{\text{high}} = 300\ \text{K},
$$
while the pressure is fixed at initial pressure
$$
p(0,t) = p_{\text{left}} = 10^5\ \text{Pa},
$$

At the right boundary the temperature and pressure are both set to their initial values:

$$
T(L,t)=T_{\text{init}} = 290\ \text{K}
$$
$$
p(L,t)=p_{\text{init}} = 10^5\ \text{Pa}.
$$

The boundary conditions and the simulation domain are shown below.

![Simulation domain and boundary conditions](1pni_1d_conduction_benchmark_domain.svg){html: width=70%}

## Results

The numerical temperature profile is compared to the analytical error-function solution for each time $t$ and each discretized position $x_i$.

To reproduce the benchmark results, run one of the tests from the build directory:

```bash
cd dumux/build-cmake/test/porousmediumflow/1p/nonisothermal/conduction
make test_1pni_conduction_tpfa
./test_1pni_conduction_tpfa params.input
```

or

```bash
cd dumux/build-cmake/test/porousmediumflow/1p/nonisothermal/conduction
make test_1pni_conduction_box
./test_1pni_conduction_box params.input
```

The output includes the VTK field `temperatureExact` for direct comparison with the numerical temperature field `T`.
The temperature field `T` over the whole domain at the end of the simulation ($t = 10^5\ \mathrm{s} \approx 1.16\ \mathrm{d}$), rendered with ParaView from this output:

![Temperature field](1pni_1d_conduction_benchmark.png){html: width=80%}

To produce the figures below, build the tests and run the plot script in the build directory:

```bash
cd dumux/build-cmake/test/porousmediumflow/1p/nonisothermal/conduction
make test_1pni_conduction_box test_1pni_conduction_tpfa test_1pni_conduction_mpfa
python3 plot_results.py test_1pni_conduction_box test_1pni_conduction_tpfa test_1pni_conduction_mpfa
```

The script expects PyVista and Matplotlib to be available for post-processing.
It writes the figures to the current directory; the images of this page are regenerated with
`--output-dir <dumux source directory>/doc/doxygen/images`. The script produces:

- `1pni_1d_conduction_benchmark_lineplot.png`: 1D comparison of the exact solution with the numerical solution for five representative output times, along the first $1.5\ \mathrm{m}$ of the domain, which the thermal front does not leave
- `1pni_1d_conduction_convergence_box.png` and `1pni_1d_conduction_convergence_cc.png`: convergence of the temperature error at the end time $t = 10^5\ \mathrm{s}$ in space and in time, for the box scheme and the cell-centered schemes

![Line plot](1pni_1d_conduction_benchmark_lineplot.png)

**Convergence**

The temperature error is evaluated at the end time $t = 10^5\ \mathrm{s} \approx 1.16\ \mathrm{d}$ in two separate studies:

- *Space convergence*: the cell size $h$ is halved from 25 to 400 cells. The time step size is fixed in each run and decreases with $h^2$ ($\Delta t = 2000\ \mathrm{s} \cdot (25/N)^2$ for $N$ cells), so that the time discretization error decreases at the same rate as the expected spatial error and stays a small part of the total error.
- *Time convergence*: the time step size $\Delta t$ is halved from $2 \cdot 10^4\ \mathrm{s}$ to $625\ \mathrm{s}$ on a fine grid with 3200 cells ($h = 1.5625\ \mathrm{mm}$), on which the spatial error is negligible.

The box scheme and the cell-centered schemes evaluate the error in different norms and are therefore shown separately.
For the box scheme, the solution is piecewise linear and the error is measured in the (continuous) $L^2$ norm, normalized with the domain volume,
$\left( \int_\Omega (T_h - T_\mathrm{exact})^2 \, \mathrm{d}x \, / \, |\Omega| \right)^{1/2}$, integrated with a quadrature rule.
For the cell-centered schemes, the error is measured in the discrete $L^2$ norm
$\left( \sum_i |V_i| \, (T_i - T_\mathrm{exact}(x_i))^2 \, / \, \sum_i |V_i| \right)^{1/2}$,
with the cell-center values $T_i$ and the cell volumes $|V_i|$.
In all tables, the rates are computed between successive rows.

*Box scheme*

![Convergence of the box scheme](1pni_1d_conduction_convergence_box.png)

Space convergence of the box scheme ($L^2$ norm of the error):

| cells | $h$ (m) | $\Delta t$ (s) | error (K) | rate |
|------:|--------:|---------------:|----------:|-----:|
| 25    | 0.2     | 2000           | 5.07e-2   | -    |
| 50    | 0.1     | 500            | 1.26e-2   | 2.01 |
| 100   | 0.05    | 125            | 3.15e-3   | 2.00 |
| 200   | 0.025   | 31.25          | 7.88e-4   | 2.00 |
| 400   | 0.0125  | 7.8125         | 1.97e-4   | 2.00 |

Time convergence of the box scheme ($L^2$ norm of the error, 3200 cells):

| $\Delta t$ (s) | error (K) | rate |
|---------------:|----------:|-----:|
| 20000          | 6.68e-2   | -    |
| 10000          | 3.40e-2   | 0.97 |
| 5000           | 1.72e-2   | 0.99 |
| 2500           | 8.63e-3   | 0.99 |
| 1250           | 4.32e-3   | 1.00 |
| 625            | 2.16e-3   | 1.00 |

*Cell-centered schemes*

![Convergence of the cell-centered schemes](1pni_1d_conduction_convergence_cc.png)

On this quasi one-dimensional grid, the MPFA fluxes reduce to the TPFA fluxes, so both schemes give the same errors.

Space convergence of the cell-centered schemes (discrete $L^2$ norm of the error):

| cells | $h$ (m) | $\Delta t$ (s) | error (K), TPFA | rate | error (K), MPFA | rate |
|------:|--------:|---------------:|----------------:|-----:|----------------:|-----:|
| 25    | 0.2     | 2000           | 6.18e-2         | -    | 6.18e-2         | -    |
| 50    | 0.1     | 500            | 1.45e-2         | 2.09 | 1.45e-2         | 2.09 |
| 100   | 0.05    | 125            | 3.58e-3         | 2.02 | 3.58e-3         | 2.02 |
| 200   | 0.025   | 31.25          | 8.92e-4         | 2.00 | 8.92e-4         | 2.00 |
| 400   | 0.0125  | 7.8125         | 2.23e-4         | 2.00 | 2.23e-4         | 2.00 |

Time convergence of the cell-centered schemes (discrete $L^2$ norm of the error, 3200 cells):

| $\Delta t$ (s) | error (K), TPFA | rate | error (K), MPFA | rate |
|---------------:|----------------:|-----:|----------------:|-----:|
| 20000          | 6.68e-2         | -    | 6.68e-2         | -    |
| 10000          | 3.40e-2         | 0.97 | 3.40e-2         | 0.97 |
| 5000           | 1.72e-2         | 0.99 | 1.72e-2         | 0.99 |
| 2500           | 8.63e-3         | 0.99 | 8.63e-3         | 0.99 |
| 1250           | 4.33e-3         | 1.00 | 4.33e-3         | 1.00 |
| 625            | 2.17e-3         | 1.00 | 2.17e-3         | 1.00 |
