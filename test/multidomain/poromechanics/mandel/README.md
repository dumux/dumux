## 1. Introduction of Mandel's Problem
Mandel's problem is a fundamental benchmark in the study of poromechanics, describing the behavior of a fluid-saturated porous medium subjected to loading. It was first introduced by J. Mandel in 1953 and is commonly used to validate theoretical models and numerical methods in poromechanics.

## 2. Problem Description
The typical setup for Mandel’s problem involves a rectangular (plane-strain) domain that is initially at equilibrium with uniform pore pressure. At time $ t = 0 $, an instantaneous load is applied (uniformly on the top and bottom surfaces). The resulting response of the medium is governed by:


- **Mechanical deformation** due to the applied load.
- **Fluid flow** induced by the pressure gradients set up by the deformation.

The interplay between these processes leads to a diffusion-like behavior in pore pressure evolution and the characteristic Mandel–Cryer effect. The analytical solution of this process is well documented in [Ref](#references)[1-3].

todo: insert the domain figure.

## 3. Governing Equations

We solve Mandel's 2D consolidation problem using a fully coupled poromechanics model implemented in Dumux. In this formulation, the primary unknowns are the solid matrix displacement ($\mathbf{u}$) and the fluid pressure ($p$). The solver employs a P1-CVFE (Box) scheme for the displacement field and a CC-TPFA method for the pressure field, ensuring robust coupling between the mechanical and hydraulic responses.

The governing equations consist of two main parts:

### 3.1. Momentum Balance for the Solid Matrix

The momentum balance (equilibrium) equation for the porous matrix is given by (see the [poro-elastic model description](https://dumux.org/docs/doxygen/master/group___poro_elastic.html)):

$$
\nabla \cdot \boldsymbol{\sigma}_{\mathrm{eff}} + \rho\, \mathbf{g} + \mathbf{f} = \rho\, \ddot{\mathbf{u}},
$$

where:
- $\boldsymbol{\sigma}_{\mathrm{eff}}$ is the effective stress tensor,
- $\rho$ is the density of the solid,
- $\mathbf{g}$ is the gravitational acceleration vector,
- $\mathbf{f} = 0$ represents additional body forces,
- $\ddot{\mathbf{u}} \approx 0$ denotes the acceleration of the solid matrix.

### 3.2. Mass Balance for the Fluid Phase

The mass conservation (fluid balance) equation is expressed as (see the [single phase model description](https://dumux.org/docs/doxygen/master/group___one_p_model.html)):

$$
\frac{\partial (\phi\, \varrho)}{\partial t} + \nabla \cdot \left\{ - \varrho\, \frac{\mathbf{K}}{\mu} \left( \nabla p - \varrho\, \mathbf{g} \right) \right\} = q,
$$

where:
- $\phi$ denotes the porosity,
- $\varrho$ is the fluid density,
- $\mathbf{K}$ is the permeability tensor,
- $\mu$ is the dynamic viscosity of the fluid,
- $q$ represents a source (or sink) term for the fluid.

### 3.3. Coupling Equations
The two equations are coupled through 
- porosity change (see the [ducumentation](https://dumux.org/docs/doxygen/master/class_dumux_1_1_porosity_deformation.html))
$$ \phi = \frac{\phi_0 + \nabla \cdot \boldsymbol u}{1 + \nabla \cdot \boldsymbol u}, $$ 
where $\phi_0$ is the initial porosity.

- effective stress.
$$\boldsymbol{\sigma}_{\mathrm{eff}} = \boldsymbol{\sigma} - \alpha p  \boldsymbol {\mathrm I},$$
where $\alpha$ is the Biot's coefficient.

## 4. Boundary and Initial Conditions

Due to the symmetry of the problem, we consider only a quarter of the domain. The boundary conditions are defined as follows:

#### Mechanical:
Initial condition: $$\boldsymbol u(t=0) = 0$$
Boundary conditions:
$$ u_y = \begin{cases}
u_y(t) & \text{on } \Gamma_{\text{top}} \\
0 & \text{on } \Gamma_{\text{bottom}}
\end{cases},$$
$$ u_x = 0 \quad \text{on } \Gamma_{\text{left}},$$
$$ \sigma_{xx} = 0 \quad \text{on } \Gamma_{\text{right}},$$
$$  \tau_{xy} = \tau_{yx} = 0 \quad \text{on } \Gamma ,$$
where $u_y(t)$ is the analytical solution (see section 6.2 in [Ref](#references)[2]). The displacement is applied so that the rigid plate condition is satisfied.


#### Hydraulic:
Initial condition:
$$p(t=0) = p^+,$$
where $p^+$ is the initial pore pressure after Skempton effect in [Ref](#references)[1].

Boundary conditions:
$$ \begin{align*}
p &= 0 \quad \text{on } \Gamma_{\text{right}}, \\
q &= 0 \quad \text{on } \Gamma_{\text{else}}.
\end{align*}$$

## 5.  Material parameters
In the input parameters, besides the force intensity F $[N/m]$, 
additional material parameters are listed as they are needed in the analytical solution.

An overview of the relation among poroelastic constants is offered in appendix B in [Ref](#references)[3].
- $K_f$, bulk modulus [Pa] of water, which is needed for $\rho = \rho_{\text {ref}} + \frac{p - p_{\text {ref}}}{K_f}$
- $\nu$, drained Poisson's ratio [-]
- $\nu_u$, undrained Poisson's ratio [-]
- $M$, Biot's modulus [Pa]
- $\alpha$, Biot's coefficient

## References:
1. "Mandel’s problem revisited." [DOI: 10.1680/geot.1996.46.2.187](https://www.icevirtuallibrary.com/doi/abs/10.1680/geot.1996.46.2.187)
2. "A coupling of mixed and continuous Galerkin finite element methods for poroelasticity I: the continuous in time case" [DOI:10.1007/s10596-007-9045-y](https://link.springer.com/article/10.1007/s10596-007-9045-y)
3. Cheng, Alexander H-D. Poroelasticity. Vol. 27. Berlin: Springer, 2016.[DOI:10.1007/978-3-319-25202-5](
https://doi.org/10.1007/978-3-319-25202-5)



## Problem Description

The problem considers a porous elastic slab confined laterally and subjected to a uniform compressive load on its top surface. The system is initially undrained, leading to a buildup of pore pressure due to fluid incompressibility. Over time, fluid flow occurs due to pressure gradients, causing the system to transition towards a drained equilibrium state.

TODO:
1) Drawing of problem setup
2) Analytical solution
3) Material parameters

Results
--------

To compile execute the code and display the results, you can use the
provided Python script

TODO: add script
```bash
python3 run_and_plot_mandel_benchmark.py
```

TODO: add convergence test
