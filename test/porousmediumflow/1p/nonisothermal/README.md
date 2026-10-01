# 1D Heat Convection with a Retarded Thermal Front {#benchmark-1pni-convection}

## Problem Description

This benchmark verifies the non-isothermal model for combined advective and conductive
heat transport in a fully saturated porous medium. Physically the setup represents a
semi-infinite, one dimensional tube that is initially at a uniform temperature and into
which water of an elevated temperature is injected at a constant rate from the left.

The injected heat travels as a front through the domain. Because part of the heat is
stored in the solid grains that the fluid passes, the front lags behind the water itself.
The benchmark therefore tests two things at once: whether the model reproduces the
correct *retarded front velocity*, and whether it reproduces the correct *front width*,
which is set by heat conduction.

### Governing Equations

Assumptions:

- local thermal equilibrium,
- one-dimensional flow and heat transport in the $x-$direction,
- single phase flow of an incompressible fluid $f$, and constant soil properties,
- no gravity,
- no source terms in the domain,
- homogeneous domain (constant porosity and permeability),
- a constant fluid influx of an elevated temperature at the left boundary.

With those assumptions, we solve the following energy balance for the porous medium, in
which heat is transported both by conduction and by advection with the flowing fluid:

$$ \phi \rho_f\frac{\partial u_f}{\partial t} + (1-\phi) \rho_s \frac{\partial u_s}{\partial t}
+ \nabla \cdot \left( \rho_f h_f \boldsymbol{v}_D \right)
= \nabla \cdot \left(\lambda_{\text{eff}} \nabla T\right) \, ,
$$

where

- $\phi$ is the porosity,
- $u$ is the internal energy,
- $T$ the overall porous medium temperature,
- $\rho_f$, $\rho_s$ is the density of the fluid/solid,
- $c_f$, $c_s$ is the mass-specific heat capacity of the fluid/solid phase,
- $\lambda_f$, $\lambda_s$ is the thermal conductivity of the fluid/solid phase,
- $\lambda_{\text{eff}}$ the effective thermal conductivity of the porous medium,
- $h_f$ is the specific enthalpy of the fluid,
- $\boldsymbol{v}_D$ is the Darcy velocity of the fluid, whose magnitude is written
  $v_D$ in the one-dimensional setting of this benchmark.

The effective thermal conductivity of the porous medium is obtained as a volume-fraction
average of the liquid and solid thermal conductivities,

$$ \lambda_{\text{eff}} = \phi\,\lambda_f + (1-\phi)\,\lambda_s  \; ,$$

and the total volumetric heat capacity of the porous medium is

$$ C_{\text{tot}} = \phi\,\rho_f c_f + (1-\phi)\,\rho_s c_s \; . $$

For the analytical solution we further assume constant material properties, a constant
Darcy velocity $v_D$, and $h_f \approx c_f T$, that is, we neglect the pressure work.
The energy balance then reduces to the one-dimensional advection-diffusion equation

$$ \frac{\partial T}{\partial t} + v_T \frac{\partial T}{\partial x}
= D_T \, \frac{\partial^2 T}{\partial x^2} \; ,$$

with the **retarded thermal front velocity**

$$ v_T = \frac{v_D\,\rho_f c_f}{C_{\text{tot}}}
       = \underbrace{\frac{v_D}{\phi}}_{\text{fluid velocity}} \cdot \frac{\phi\,\rho_f c_f}{C_{\text{tot}}}
       = \frac{1}{R} \, \frac{v_D}{\phi} \; ,
\qquad R = \frac{C_{\text{tot}}}{\phi\,\rho_f c_f} \; ,$$

and the **effective thermal diffusivity**

$$ D_T = \frac{\lambda_{\text{eff}}}{C_{\text{tot}}} \; . $$

The retardation factor $R$ is the ratio of the heat stored in the whole porous medium to
the heat stored in the fluid alone. For the parameters of this benchmark $R = 1.77$, so
the thermal front moves at $v_T = 1.42 \cdot 10^{-4}\ \mathrm{m\,s^{-1}}$ while the water
moves at $v_D/\phi = 2.5 \cdot 10^{-4}\ \mathrm{m\,s^{-1}}$.

## Analytical Solution

Two analytical solutions of the equation above are used, one for each of the two things
the benchmark checks.

**1. Moving step front (pure advection)**
In the limit of vanishing conduction, $D_T \rightarrow 0$, the equation becomes a pure
advection equation and its solution is a sharp front travelling at the retarded velocity,

$$ T(x,t) = \begin{cases}
T_{\text{high}} & x < v_T t \\
T_{\text{init}} & x > v_T t
\end{cases} \; . $$

It reproduces the position of the front exactly, but says nothing about its shape. This
is the solution the test problem itself provides.

**2. Advection-diffusion solution (Ogata-Banks)**
With conduction retained, on a semi-infinite domain with the initial and boundary
conditions

$$ T(x,0) = T_{\text{init}} \, , \qquad T(0,t) = T_{\text{high}} \, , \qquad T(\infty,t) = T_{\text{init}} \, , $$

the equation is solved by the classical solution of Ogata and Banks @cite OgataBanks1961:

$$
\frac{T(x,t) - T_{\text{init}}}{T_{\text{high}} - T_{\text{init}}}
= \frac{1}{2}\left[
\operatorname{erfc}\left(\frac{x - v_T t}{2\sqrt{D_T t}}\right)
+ \exp\left(\frac{v_T x}{D_T}\right)
\operatorname{erfc}\left(\frac{x + v_T t}{2\sqrt{D_T t}}\right)
\right] \; .
$$

It describes a front centred at $x = v_T t$, which has reached $4.25\ \mathrm{m}$ at the
end of the simulation, and which is broadened by conduction as $\sqrt{4 D_T t}$. This is
the solution of the equation the model actually solves, so it is the reference for the
shape of the front.

## Parameters

| Parameter                        | Symbol               | Value | Unit  |
|----------------------------------|---------------------:|:------|:------|
| Porosity                         | $\phi$               | 0.4   | -     |
| Permeability                     | $K$                  | 1e-10 | m$^2$ |
| Darcy velocity at the inlet      | $v_D$                 | 1e-4  | m s$^{-1}$ |
| Initial temperature              | $T_{\mathrm{init}}$  | 290   | K     |
| Injection temperature            | $T_{\mathrm{high}}$  | 291   | K     |
| Pressure (initial & outlet)      | $p_{\mathrm{init}}$  | 1e5   | Pa    |
| Solid density                    | $\rho_s$             | 2700  | kg m$^{-3}$ |
| Solid thermal conductivity       | $\lambda_s$          | 2.8   | W m$^{-1}$ K$^{-1}$ |
| Solid heat capacity              | $c_s$                | 790   | J kg$^{-1}$ K$^{-1}$ |

Further, the fluid is water. The fluid properties (density, thermal conductivity and heat
capacity) are temperature- and pressure-dependent and are evaluated using the IAPWS
Industrial Formulation 1997 for the thermodynamic properties of water and steam
@cite IAPWS1997. In DuMux, it is implemented by the component `Dumux::Components::H2O`
(@ref Components), which builds on the region formulations in
`dumux/material/components/iapws/` (@ref IAPWS). The
pressure drop that the prescribed influx induces is only about $0.22\ \mathrm{bar}$ over
the $20\ \mathrm{m}$ of the shipped test and $0.09\ \mathrm{bar}$ over the $8\ \mathrm{m}$
benchmark domain, so the fluid properties stay close to those at the initial state and
the assumption of constant material properties is well satisfied. The resulting derived
quantities are

| Quantity                        | Symbol                 | Value    | Unit  |
|---------------------------------|-----------------------:|:---------|:------|
| Effective thermal conductivity  | $\lambda_{\text{eff}}$ | 1.92     | W m$^{-1}$ K$^{-1}$ |
| Total volumetric heat capacity  | $C_{\text{tot}}$       | 2.95e6   | J m$^{-3}$ K$^{-1}$ |
| Retardation factor              | $R$                    | 1.77     | -     |
| Retarded front velocity         | $v_T$                  | 1.42e-4  | m s$^{-1}$ |
| Effective thermal diffusivity   | $D_T$                  | 6.50e-7  | m$^2$ s$^{-1}$ |

##  Simulation Setup

We solve the energy and mass balance using the single phase, non-isothermal model in
DuMux @ref OnePModel @ref NIModel. For the spatial discretization we choose a finite
volume discretization with a cell-centered TPFA (@ref CCTpfaDiscretization) or MPFA
(@ref CCMpfaDiscretization) scheme, or use the vertex-centered box scheme
(@ref BoxDiscretization).

### Simulation Domain
The computational domain is a one-dimensional `Dune::YaspGrid`. The default test uses a
$20\ \mathrm{m}$ domain with $80$ cells; the benchmark script below uses a shortened
$8\ \mathrm{m}$ domain, which still keeps the outlet far behind the front, so that the
refined grids stay affordable. For the simulated time of $3 \cdot 10^4\ \mathrm{s}$ the
thermal front reaches $x = 4.25\ \mathrm{m}$ and its conductive tail does not reach the
outlet, so the finite domain is suitable to be compared to the semi-infinite analytical
solution.

### Initial and Boundary Conditions
Initially, the temperature and pressure are set to

$$
T(x,0)=T_{\text{init}} = 290\ \text{K} \, , \qquad
p(x,0)=p_{\text{init}} = 10^5\ \text{Pa}.
$$

At the left boundary a Neumann condition prescribes the mass influx and the energy influx
that the injected water carries,

$$
(\rho_f \boldsymbol{v}_D) \cdot \boldsymbol{n} = -v_D\,\rho_f \, , \qquad
(\rho_f h_f \boldsymbol{v}_D) \cdot \boldsymbol{n} = -v_D\,\rho_f\,h_f(T_{\text{high}}, p) \, ,
\qquad T_{\text{high}} = 291\ \text{K},
$$

with $v_D = 10^{-4}\ \mathrm{m\,s^{-1}}$. At the right boundary temperature and pressure are
both fixed at their initial values,

$$
T(L,t)=T_{\text{init}} = 290\ \text{K} \, , \qquad
p(L,t)=p_{\text{init}} = 10^5\ \text{Pa}.
$$

@note The inlet condition is not the constant-temperature inlet that the Ogata-Banks
solution assumes. The solution of prescribing the *total* energy influx
is presented by van Genuchten and Alves
@cite vanGenuchtenAlves1982 in A3. At the Péclet number of this benchmark
($v_T L / D_T \approx 4.4\cdot10^{3}$) the two solutions differ by less than
$10^{-4}\ \mathrm{K}$, which is two orders of magnitude below the discretization error,
so the simpler Ogata-Banks form is used as the reference.

### Resolution
Four uniformly refined grids are tested, so that the front is resolved to a different
degree on each. Cell size and time step are both refined by a factor of four per level:

| Level | Cells | $\Delta x$ [m] | $\Delta t_{\max}$ [s] |
|------:|------:|---------------:|----------------------:|
| 1     | 100   | 0.08           | 160                   |
| 2     | 400   | 0.02           | 40                    |
| 3     | 1600  | 0.005          | 10                    |
| 4     | 6400  | 0.00125        | 2.5                   |

## Results

To reproduce the benchmark results, run one of the tests from the build directory:

```bash
cd dumux/build-cmake/test/porousmediumflow/1p/nonisothermal
make test_1pni_convection_tpfa
./test_1pni_convection_tpfa params_convection.input
```

or

```bash
cd dumux/build-cmake/test/porousmediumflow/1p/nonisothermal
make test_1pni_convection_box
./test_1pni_convection_box params_convection.input
```

Next to the numerical temperature field `T`, the output contains the VTK field
`temperatureExact`, which `OnePNIConvectionProblem::updateExactTemperature()` fills with
the moving step front described above.

To build and run the refinement sequence and produce the plot and the table below in one
step, execute:

```bash
python3 plot_convection.py
```

The script expects PyVista and Matplotlib to be available for post-processing. It
evaluates the Ogata-Banks solution itself, using the front velocity and the total heat
capacity that the test problem reports on standard output, so that the references use the
same water properties as the simulation and no change to the test problem is needed. It
produces one figure in the `build-cmake` directory:

- `1pni_1d_convection_benchmark_lineplot.png`: the temperature profile around the thermal
  front at the end of the simulation, for all four grid resolutions, together with the
  Ogata-Banks solution and the retarded step front.

![Line plot](1pni_1d_convection_benchmark_lineplot.png)

The plot shows how the numerical profiles approach the Ogata-Banks solution under
refinement, and how far both are from the pure-advection step.

### What is compared
For comparing  the numerical results against the analytical solutions, we do this for two separate quantities.

*1. Front position, against the step front.* Because the step carries only a position,
it is compared against the position of the numerical front, its **barycenter**

$$ x_{\text{bary}} = \frac{1}{T_{\text{high}} - T_{\text{init}}}
\int_0^L \left(T - T_{\text{init}}\right)\,\mathrm{d}x \; , $$

that is, the position of the step front that stores the same heat as the smeared profile. Inserting the two analytical
solutions into this definition gives the two reference positions

$$ x_{\text{bary}}^{\text{step}} = v_T t = 4.24896\ \mathrm{m}\; , \qquad
   x_{\text{bary}}^{\text{OB}} = v_T t + \frac{D_T}{v_T} = 4.25355\ \mathrm{m} \; . $$

*2. Front shape, against Ogata-Banks.* The whole profile is compared in the discrete
$L^2$ norm, normalized by the domain length,

$$ \lVert T - T^{\text{OB}} \rVert_{L^2}
= \sqrt{\frac{1}{L_x} \sum_i \left(T_i - T^{\text{OB}}(x_i, t)\right)^2 \Delta x } \; . $$

The two comparisons are independent of each other and are reported separately below. Note
that the $L^2$ norm is taken over the temperature profile and has nothing to do with the
barycenters, and that the barycenters are positions, not errors: there is no exact value
for them to converge to other than the two analytical positions quoted above.

### Front position 
Barycenters of the numerical profiles, next to the analytical
positions:

| Cells | $x_{\text{bary}}$ [m] | $x_{\text{bary}} - x_{\text{bary}}^{\text{step}}$ [m] | $x_{\text{bary}} - x_{\text{bary}}^{\text{OB}}$ [m] |
|------:|----------------------:|-----------------------------------------------------:|---------------------------------------------------:|
| 100   | 4.25417               | $+5.21\cdot10^{-3}$                                  | $+0.62\cdot10^{-3}$                                 |
| 400   | 4.25452               | $+5.56\cdot10^{-3}$                                  | $+0.97\cdot10^{-3}$                                 |
| 1600  | 4.25465               | $+5.69\cdot10^{-3}$                                  | $+1.10\cdot10^{-3}$                                 |
| 6400  | 4.25464               | $+5.68\cdot10^{-3}$                                  | $+1.09\cdot10^{-3}$                                 |

The position of the barycenters from the numerical and the analytical solutions agree upto $5.7$mm against the sharp front and upto $1.1$mm against the Ogata-Banks solution on a $8$m long domain.

The difference of the numerical solution against the step function splits into the two contributions. Most of it,
$D_T / v_T = 4.6\ \mathrm{mm}$, is the diffusive shift that the Ogata-Banks barycenter
already has over the step front. What is left is the $1.1\ \mathrm{mm}$ of the third
column, the distance to the Ogata-Banks barycenter, which is related to the constant-property
assumption of the analytical solutions.

### Front shape
$L^2$ error of the temperature profile against the Ogata-Banks solution, the advective-diffusion solution,
with the convergence rate between successive levels:

| Cells | $\lVert T - T^{\text{OB}} \rVert_{L^2}$ [K] | rate |
|------:|--------------------------------------------:|-----:|
| 100   | 0.0947                                      | –    |
| 400   | 0.0456                                      | 1.06 |
| 1600  | 0.0168                                      | 1.44 |
| 6400  | 0.0053                                      | 1.68 |

The shape of the front, unlike its position, is recovered by resolving it.
