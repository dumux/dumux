# Embedded network 1D-3D model (tissue perfusion)

We solve a tracer transport problem in a domain that consist of a porous medium block
that has an embedded transport network. In this example the porous medium is brain tissue
and the embedded network is the vasculature (blood vessels).

__Table of contents__. This description is structured as follows:

[[_TOC_]]

## Problem set-up
In this example we simulate clearance of a substance present in the tissue through the blood. The tissue
cube is assigned no-flux/symmetry boundary conditions assuming that identical cubes are mirrored on all sides.
Therefore, the tracer has to cross the vessel wall into the network (vessel lumen). It then gets transported
in the blood stream by advection and diffusion. In the network, the tracer mole fraction is zero at the inlet
and at the outlet the mole fraction gradient is zero. Thus, the tracer is transported out of the domain by advection only.
VTK output is written in every time step, and the total tracer concentration in the tissue is written into a text file
along the simulation.

## Network data and blood flow
The domain consists of a small blood vessel network embedded
in a porous tissue (pore space is the pore space of the extra-cellular matrix).
The network is extracted based on the mouse cortical brain data from [Blinder (2013)](https://doi.org/10.1038/nn.3426).
With the boundary condition estimated by [Schmid (2017a)](https://doi.org/10.1371/journal.pcbi.1005392)
and the pressure and network data available from [Schmid (2017b)](https://doi.org/10.5281/zenodo.269650), blood flow
simulations in DuMu<sup>x</sup> give a pressure field in the entire network. A small part (200µm)³ has been extracted
for this example. The example also contains a blood flow solver which uses the pressure boundary data to compute
blood volume flux for every facet (vertex) in the network.
Note: The blood flow solver is currently not documented in detail. This might be added in the future.

<figure>
    <center>
        <img src="img/network.png" alt="Free-flow setup" width="50%"/>
        <figcaption> <b> Fig.1 </b> - Capillary blood vessel network used in this example.</figcaption>
    </center>
</figure>

## Mathematical model
We solve the following coupled, mixed-dimensional PDE system:

```math
\begin{align}
	 \frac{ \partial (\phi_\mathrm{T} \varrho_\mathrm{T} x_\mathrm{T})}{\partial t} - \nabla\cdot \left( \phi_\mathrm{T} D_{\text{app},\mathrm{T}} \varrho_\mathrm{T} \nabla x_\mathrm{T} \right) &= \hat{q} \Phi_\Lambda & \text{in} \quad \Omega, \\
	 \frac{\partial (A_\mathrm{B} \varrho_\mathrm{B} x_\mathrm{B})}{\partial t} + \frac{\partial}{\partial s}  \left(A_\mathrm{B}v_\mathrm{B}\varrho_\mathrm{B}x_\mathrm{B} - A_\mathrm{B} D_{\text{app},A_\mathrm{B}} \varrho_\mathrm{B} \frac{ \partial x_\mathrm{B}}{\partial s} \right) &= -\hat{q} & \text{on} \quad \Lambda, \\
     \hat{q} &= - \int_P C_M D \bar{\varrho} ( x_\mathrm{T} - \Pi x_\mathrm{B}) \mathrm{d}\zeta,
\end{align}
```

where the subscript T and B denote the tissue and the network (blood flow) compartment,
$`x`$ is the tracer mole fraction, $`\varrho`$ the molar density of the mixture, $`\phi`$ is the porosity,
$`A_\mathrm{B}`$ denotes the network (vessel lumen) cross-sectional area, $`P`$ denotes the cross-sectional perimeter,
$`D`$ is the free diffusion coefficient, $`D_{\text{app}}`$ apparent diffusion coefficients and $C_M$ a membrane diffusivity factor.
$`Q_\mathrm{B} := A_\mathrm{B}v_\mathrm{B}`$ is a given blood flow field transporting the tracer by advection.
Furthermore, isothermal conditions with a homogeneous temperature distribution of constant $`T=37^\circ C`$ are assumed.
The 1D network PDE is formulated in terms of the local axial coordinate $`s`$.

For the coupling, we use a formulation with line source and a perimeter integral operator, that is $`\Phi_\Lambda`$ is chosen to be a delta distribution centered on vessel centerline, and $`x_\mathrm{B}`$, which is assumed to be constant on each cross-section (well-mixed),
and every point on the surface is assumed to be uniquely identified with a position $`s`$ in the network such that we can extend
$`x_\mathrm{B}`$ with $`\Pi`$ to the surface. The integral is evaluated numerically. Each integration point is represented in
the code as a point source. The point source values (the integrand) are implemented in the problem class function `pointSource(...)`.

All equations are discretized with a cell-centered finite-volume scheme (with two-point flux approximation)
as spatial discretization method with tracer mole fraction as primary variables. The equations are based on DuMu<sup>x</sup>'s
`Tracer` model, an advection-diffusion-reaction model for porous media. In time, we use an implicit Euler scheme.
The arising linear system is solved with a stabilized bi-conjugate gradient solver (`BiCGSTAB`) with a block-diagonal
zero-fill incomplete LU factorization (`ILU0`) preconditioner. For details on the spatial discretization scheme,
we recommend the [DuMu<sup>x</sup> documentation](https://dumux.org/docs/doxygen/master/group___discretization.html).
or the [DuMu<sup>x</sup> paper](https://doi.org/10.1016/j.camwa.2020.02.012).

# Implementation

Below is an overview over the files in the example folder and links to the more detailed
descriptions of some individual files. The main program flow can be found in the `main`
function implemented in `main.cc`. This is the only source file that includes all the
supplementary header files which contain important parts of the implementation. The
`CMakeLists.txt` file instructs `CMake` to create an executable from `main.cc` and
`CMake` will configure and figure out the necessary compiler command for compilation.
(The executable is built in the build folder by typing `make example_embedded_network_1d3d`.
See [dumux-course basic exercise](https://git.iws.uni-stuttgart.de/dumux-repositories/dumux-course/-/tree/master/exercises/exercise-basic)
for more detailed instructions.)

Key classes that are implemented by the user are the problems (`problem.hh`), one for the tissue
and one for the network. The problems define initial and boundary conditions. They also implement source terms.
As the coupling between tissue and network is realized through point sources at integration points,
the coupling condition is implemented in the problem classes. The problem classes are used by
DuMu<sup>x</sup> during assembly of the system matrix and the residual (discrete equations).
For this to work, we have to tell the model which problem class to use. This is achieved through the
property system in `properties.hh` where the `Problem` property is specialized for the main models created in this example
(`NetworkTransportModel` and `TissueTransportModel`). These models are passed in the main function to the assembler.

Secondly, the spatial parameters (`spatialparams.hh`) are classes that specify (possibly) spatially varying parameter.
One such parameter is the radius field for the network. In the class `NetworkSpatialParams`, the radius field is
read from the grid file `network.dgf` (which is in the very simple, human-readable Dune Grid Format).
As for the problem, spatial parameters have to be added to the model by specializing the `SpatialParams` property
for the model in `properties.hh`.

Apart from problem and spatial params, the model (`properties.hh`) also has other configurable parameters.
(In fact most of the inner workings of the assembler can be configured like this.) One example is the grid type
used for each model. Dune provides specialized implementations for certain grid types behind a common interface.
In this exercise, we use a structured Cartesian grid (`YaspGrid`) for the tissue domain and an embedded network
grid (`FoamGrid`) for the network.
With the model configuration through the property system in mind, we can better understand the main program (`main.cc`)
and how the boundary conditions and parameter setting make their way into the assembler.

# Output

The simulation output VTK files at every time step that can be inspected in `ParaView`.
It also writes a simple text file `clearance_tracer_amounts.dat` containing the total tracer
amount at each time step. There is a small Python script `plot.py`, that visualizes this
data. (Execute `./plot.py` or `python3 plot.py`.)

## Folder layout and files

```
└── embedded_network_1d3d/
    ├── CMakeLists.txt          -> build system file
    ├── main.cc                 -> main program flow
    ├── params.input            -> runtime parameter configuration file
    ├── plot.py                 -> small Python script to plot the tracer amount in the tissue
    ├── problem.hh              -> boundary & initial conditions and source terms
    ├── properties.hh           -> compile time model configuration
    ├── solver.hh               -> specialized solver for linear multi-domain problems
    ├── spatialparams.hh        -> spatially varying model parameters
    └── tracerfluidsystem.hh    -> fluid system with a single tracer
```
