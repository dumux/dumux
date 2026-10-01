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
cd dumux/build-cmake/test/porousmediumflow/1p/nonisothermal
make test_1pni_conduction_tpfa
./test_1pni_conduction_tpfa params_conduction.input
```

or

```bash
cd dumux/build-cmake/test/porousmediumflow/1p/nonisothermal
make test_1pni_conduction_box
./test_1pni_conduction_box params_conduction.input
```

The output includes the VTK field `temperatureExact` for direct comparison with the numerical temperature field `T`.

To build and run the test and produce the plot below in one step, execute:

```bash
python3 plot_conduction.py
```

The script expects PyVista and Matplotlib to be available for post-processing.

The script produces one figure in the `build-cmake` directory:

- `1pni_1d_conduction_benchmark_lineplot.png`: 1D comparison of the exact solution with the numerical solution for five representative output times, along the first $1.5\ \mathrm{m}$ of the domain, which the thermal front does not leave

![Line plot](1pni_1d_conduction_benchmark_lineplot.png)

The temperature field over the whole domain at the end of the simulation ($t = 10^5\ \mathrm{s}$), rendered with ParaView:

![Temperature field](1pni_1d_conduction_benchmark.png){html: width=80%}
