Mandel's problem
--------------------

We solve Mandel's 2D consolidation problem with the poromechanics model.
The solver solves the fully coupled system with primary variables displacement and fluid pressure
using the P1-CVFE (Box) scheme for the displacement and CC-TPFA for the pressure. The governing equations
are given by

For displacement, we use the momentum balance equation, see [model description](https://dumux.org/docs/doxygen/master/group___poro_elastic.html):
$$
 \nabla\cdot\boldsymbol{\sigma_{\mathrm{eff}}} + \rho \mathbf{g} + \mathbf{f} = \rho\ddot{\mathbf{u}},
$$


 For pressure, we use the mass balance equation, see [model description](https://dumux.org/docs/doxygen/master/group___one_p_model.html):
 $$
 \frac{\partial (\phi \varrho) }{\partial t} + \nabla \cdot \left\lbrace - \varrho \frac{\textbf K}{\mu} \left( \nabla p -\varrho {\textbf g} \right) \right\rbrace = q
 $$

Description of the benchmark case
-----------------------------------
For details, please using the following papers as reference.
- paper1: "Mandel’s problem revisited." [DOI: 10.1680/geot.1996.46.2.187](https://www.icevirtuallibrary.com/doi/abs/10.1680/geot.1996.46.2.187)
- paper2: "A coupling of mixed and continuous Galerkin finite element methods for poroelasticity I: the continuous in time case" [DOI:10.1007/s10596-007-9045-y](https://link.springer.com/article/10.1007/s10596-007-9045-y)

## Problem setup
# Mandel's Problem in Poromechanics

Mandel's problem is a fundamental benchmark in the study of poromechanics, describing the behavior of a fluid-saturated porous medium subjected to loading. It was first introduced by J. Mandel in 1953 and is commonly used to validate theoretical models and numerical methods in poromechanics.

## Problem Description

The problem considers a porous elastic slab confined laterally and subjected to a uniform compressive load on its top surface. The system is initially undrained, leading to a buildup of pore pressure due to fluid incompressibility. Over time, fluid flow occurs due to pressure gradients, causing the system to transition towards a drained equilibrium state.

## Key Features

- **Material**: A porous, elastic medium saturated with an incompressible fluid.
- **Boundary Conditions**:
  - Lateral boundaries are impermeable and fixed.
  - Vertical surfaces allow for fluid flow and stress application.
- **Loading**: Uniform vertical stress is applied on the top surface.

## Solution and Significance

Mandel's problem highlights:
- **Coupled Behavior**: The interaction between mechanical deformation and fluid diffusion.
- **Pore Pressure Evolution**: An initial rise in pore pressure followed by a decay as fluid diffuses out.
- **Consolidation Effects**: The time-dependent settlement of the porous medium due to fluid flow.

The analytical solution to Mandel's problem provides insights into:
- Stress and pore pressure distributions.
- Consolidation rates and time scales.
- Validation benchmarks for numerical models in poromechanics.

Mandel's problem remains a cornerstone in the study of coupled solid-fluid interactions and is extensively used in fields such as geomechanics, petroleum engineering, and hydrology.

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
