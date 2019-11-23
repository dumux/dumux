# Shallow water flow with bottom friction
This example shows how the shallow water flow model can be
applied to simulate steady subcritical flow including
bottom friction (bed shear stress).


## Shallow water model
The shallow water equations (SWEs) are given as:

$`
\frac{\partial \mathbf{U}}{\partial t} +
\frac{\partial \mathbf{F}}{\partial x} +
\frac{\partial \mathbf{G}}{\partial y} - \mathbf{S_b} - \mathbf{S_f} = 0
`$

where $U$, $F$ and $G$  defined as

$`
\mathbf{U} = \begin{bmatrix} h \\ uh \\ vh \end{bmatrix},
\mathbf{F} = \begin{bmatrix} hu \\ hu^2  + \frac{1}{2} gh^2 \\ huv \end{bmatrix},
\mathbf{G} = \begin{bmatrix} hv \\ huv \\ hv^2  + \frac{1}{2} gh^2 \end{bmatrix}
`$

Z is the bedSurface, h the water depth, u the velocity in
x-direction and v the velocity in y-direction, g is the constant of gravity.

The source terms for the bed friction $`S_b`$ and bed slope
$`S_f`$ are given as
$`
\mathbf{S_b} = \begin{bmatrix} 0 \\ -gh \frac{\partial z}{\partial x}
               \\ -gh \frac{\partial z}{\partial y}\end{bmatrix},
\mathbf{S_f} = \begin{bmatrix} 0 \\ -ghS_{fx} \\ -ghS_{fy}\end{bmatrix}.
`$

For this example, a cell-centered finite volume method (cctpfa) is applied to solve the SWEs
in combination with a fully-implicit time discretization. For cases where no sharp fronts or
traveling waves occur it is possible to apply time steps larger than CFL number = 1 to reduce
the computation time. Even if a steady state solution is considered, an implicit time stepping method
is applied.

## Problem set-up
The model domain is given by a rough channel with a slope of 0.001.
The domain is 500 meters long and 10 meters wide.
![Domain](img/domain.png).

Bottom friction is considered by applying
the friction law of Manning (Manning n = 0.025). At the lateral sides no friction is considered and  a
no-flow no slip boundary condition is applied. This is the default boundary condition for the shallow water model.


At the left border a discharge boundary condition
is applied as inflow boundary condition with q = -1.0 ($`m^2 s^{-1}`$). At the right border a water fixed depth boundary condition
is applied for the outflow. Normal flow is assumed, therefore the water depth at the right border is calculated after
the of Gaukler-Manning-Strickler equation:

 $` v_m = 1/n * R_{hy}^{2/3} * I_s^{1/2}`$

Where the mean velocity $`v_m`$ is given as

$`v_m = \frac{q}{h}`$

$`n`$ is the friction value after Manning. $`R_{hy}`$ the hydraulic radius, which is assumed to be equal to
the water depth. $`I_s`$ is the bed slope and $`q`$ the unity inflow discharge

The water depth h can be calculated as
$`h = \left(\frac{n*q}{\sqrt{I_s}} \right)^{3/5}`$

The formula of Gaukler Manning and Strickler is also used to calculate the analytic solution. All parameters
for the simulation are given in the file *params.input*.


## The file `spatialparams.hh`


We include the basic spatial parameters for finite volumes file from which we will inherit
```cpp
#include <dumux/material/spatialparams/fv.hh>
```
The parameters header is needed to retrieve run-time parameters.
```cpp
#include <dumux/common/parameters.hh>
```
We include all friction laws, between we can choose for the calculation of the bottom friction source.
```cpp
#include <dumux/material/fluidmatrixinteractions/frictionlaws/frictionlaw.hh>
#include <dumux/material/fluidmatrixinteractions/frictionlaws/manning.hh>
#include <dumux/material/fluidmatrixinteractions/frictionlaws/nikuradse.hh>
```
We enter the namespace Dumux. All Dumux functions and classes are in a namespace Dumux, to make sure they don`t clash with symbols from other libraries you may want to use in conjunction with Dumux.
```cpp
namespace Dumux {
```
In the RoughChannelSpatialParams class we define all functions needed to describe the spatial distributed parameters.
```cpp
template<class GridGeometry, class Scalar, class VolumeVariables>
class RoughChannelSpatialParams
: public FVSpatialParams<GridGeometry, Scalar,
                         RoughChannelSpatialParams<GridGeometry, Scalar, VolumeVariables>>
{
```
We introduce using declarations that are derived from the property system which we need in this class
```cpp
    using ThisType = RoughChannelSpatialParams<GridGeometry, Scalar, VolumeVariables>;
    using ParentType = FVSpatialParams<GridGeometry, Scalar, ThisType>;
    using GridView = typename GridGeometry::GridView;
    using FVElementGeometry = typename GridGeometry::LocalView;
    using SubControlVolume = typename FVElementGeometry::SubControlVolume;
    using Element = typename GridView::template Codim<0>::Entity;
    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;

public:
```
In the constructor be read some values from the `params.input` and initialize the friciton law.
```cpp
    RoughChannelSpatialParams(std::shared_ptr<const GridGeometry> gridGeometry)
    : ParentType(gridGeometry)
    {
        gravity_ = getParam<Scalar>("Problem.Gravity");
        bedSlope_ = getParam<Scalar>("Problem.BedSlope");
        frictionLawType_ = getParam<std::string>("Problem.FrictionLaw");
        initFrictionLaw();
    }
```
We initialize the friction law based on the law specified in `params.input`.
```cpp
    void initFrictionLaw()
    {
      if (frictionLawType_ == "Manning")
      {
          Scalar manningN = getParam<Scalar>("Problem.ManningN");
          frictionLaw_ = std::make_unique<FrictionLawManning<VolumeVariables>>(gravity_, manningN);
      }
      else if (frictionLawType_ == "Nikuradse")
      {
          Scalar ks = getParam<Scalar>("Problem.Ks");
          frictionLaw_ = std::make_unique<FrictionLawNikuradse<VolumeVariables>>(ks);
      }
      else
      {
          std::cout<<"The FrictionLaw in params.input is unknown. Valid entries are `Manning` and `Nikuradse`!"<<std::endl;
      }
    }
```
Use this function, if you want to vary the value for the gravity.
```cpp
    Scalar gravity(const GlobalPosition& globalPos) const
    {
        return gravity_;
    }
```
Use this function for a constant gravity.
```cpp
    Scalar gravity() const
    {
        return gravity_;
    }
```
This function returns an object of the friction law class, which is initialized with the appropriate friction values. If you want to use different friciton values or laws, you have to use a vector of unique_ptr for `frictionLaw_` and pick the right friction law instances via the `element` argument.
```cpp
    const FrictionLaw<VolumeVariables>& frictionLaw(const Element& element,
                                                    const SubControlVolume& scv) const
    {
        return *frictionLaw_;
    }
```
Define the bed surface based on the `bedSlope_`.
```cpp
    Scalar bedSurface(const Element& element,
                      const SubControlVolume& scv) const
    {
        return 10.0 - element.geometry().center()[0] * bedSlope_;
    }

private:
    Scalar gravity_;
    Scalar bedSlope_;
    std::string frictionLawType_;
    std::unique_ptr<FrictionLaw<VolumeVariables>> frictionLaw_;
};
```
end of namespace Dumux.
```cpp
}

```



## The file `problem.hh`


## Include files
We use the dune yasp grid.
```cpp
#include <dune/grid/yaspgrid.hh>
```
We include the cell centered, two-point-flux discretization scheme.
```cpp
#include <dumux/discretization/cctpfa.hh>
```
The parameters header is needed to retrieve run-time parameters.
```cpp
#include <dumux/common/parameters.hh>
```
We include the header which are needed for shallow water models.
```cpp
#include <dumux/freeflow/shallowwater/model.hh>
#include <dumux/freeflow/shallowwater/problem.hh>
#include <dumux/freeflow/shallowwater/boundaryfluxes.hh>
```
We include the header that specifies all spatially variable parameters.
```cpp
#include "spatialparams.hh"
```
## Define basic properties for our simulation
We enter the namespace Dumux. All Dumux functions and classes are in a namespace Dumux, to make sure they don't clash with symbols from other libraries you may want to use in conjunction with Dumux. One could use these functions and classes by prefixing every use of these names by ::, but that would quickly become cumbersome and annoying. Rather, we simply import the entire Dumux namespace for general use.
```cpp
namespace Dumux {
```
The problem class is forward declared.
```cpp
template <class TypeTag>
class RoughChannelProblem;
```
We enter the namespace Properties, which is a sub-namespace of the namespace Dumux.
```cpp
namespace Properties {
```
A TypeTag for our simulation is created which inherits from the shallow water model and the
cell centered, two-point-flux discretization scheme.
```cpp
namespace TTag {
struct RoughChannel { using InheritsFrom = std::tuple<ShallowWater, CCTpfaModel>; };
}
```
We define the grid of our simulation. We use a two-dimensional Yasp Grid.
```cpp
template<class TypeTag>
struct Grid<TypeTag, TTag::RoughChannel>
{ using type = Dune::YaspGrid<2, Dune::TensorProductCoordinates<GetPropType<TypeTag, Properties::Scalar>, 2> >; };
```
We set the problem. The problem class specifies initial and boundary conditions and is defined below.
```cpp
template<class TypeTag>
struct Problem<TypeTag, TTag::RoughChannel>
{ using type = Dumux::RoughChannelProblem<TypeTag>; };
```
We define the spatial parameters for our simulation. The values are specified in the corresponding spatialparameters header file, which is included above.
```cpp
template<class TypeTag>
struct SpatialParams<TypeTag, TTag::RoughChannel>
{
private:
```
We define convenient shortcuts to the properties GridGeometry, Scalar, ElementVolumeVariables and VolumeVariables:
```cpp
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using ElementVolumeVariables = typename GetPropType<TypeTag, Properties::GridVolumeVariables>::LocalView;
    using VolumeVariables = typename ElementVolumeVariables::VolumeVariables;
```
Finally we set the spatial parameters:
```cpp
public:
    using type = RoughChannelSpatialParams<GridGeometry, Scalar, VolumeVariables>;
};
```
We enable caching for the FV grid geometry and the grid volume variables. The cache
stores values that were already calculated for later usage. This makes the simulation faster.
```cpp
template<class TypeTag>
struct EnableGridGeometryCache<TypeTag, TTag::RoughChannel>
{ static constexpr bool value = true; };

template<class TypeTag>
struct EnableGridVolumeVariablesCache<TypeTag, TTag::RoughChannel>
{ static constexpr bool value = false; };
```
We leave the namespace Properties.
```cpp
}
```
## The problem class
We enter the problem class where all necessary boundary conditions and initial conditions are set for our simulation.
As this is a shallow water problem, we inherit from the basic ShallowWaterProblem.
```cpp
template <class TypeTag>
class RoughChannelProblem : public ShallowWaterProblem<TypeTag>
{
```
We use convenient declarations that we derive from the property system.
```cpp
    using ParentType = ShallowWaterProblem<TypeTag>;
    using PrimaryVariables = GetPropType<TypeTag, Properties::PrimaryVariables>;
    using BoundaryTypes = GetPropType<TypeTag, Properties::BoundaryTypes>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using Indices = typename GetPropType<TypeTag, Properties::ModelTraits>::Indices;
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using NeumannFluxes = GetPropType<TypeTag, Properties::NumEqVector>;
    using ElementVolumeVariables = typename GetPropType<TypeTag, Properties::GridVolumeVariables>::LocalView;
    using GridVariables = GetPropType<TypeTag, Properties::GridVariables>;
    using ElementFluxVariablesCache = typename GridVariables::GridFluxVariablesCache::LocalView;
    using VolumeVariables = typename ElementVolumeVariables::VolumeVariables;
    using FVElementGeometry = typename GetPropType<TypeTag, Properties::GridGeometry>::LocalView;
    using SubControlVolumeFace = typename FVElementGeometry::SubControlVolumeFace;
    using GridView = GetPropType<TypeTag, Properties::GridView>;
    using Element = typename GridView::template Codim<0>::Entity;
    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;
    using NumEqVector = GetPropType<TypeTag, Properties::NumEqVector>;
    using SubControlVolume = typename FVElementGeometry::SubControlVolume;

public:
```
This is the constructor of our problem class.
```cpp
    RoughChannelProblem(std::shared_ptr<const GridGeometry> gridGeometry)
    : ParentType(gridGeometry)
    {
```
We read the parameters from the params.input file.
```cpp
        name_ = getParam<std::string>("Problem.Name");
        constManningN_ = getParam<Scalar>("Problem.ManningN");
        bedSlope_ = getParam<Scalar>("Problem.BedSlope");
        discharge_ = getParam<Scalar>("Problem.Discharge");
```
We calculate the outflow boundary condition using the Gauckler-Manning-Strickler formula.
```cpp
        hBoundary_ = this->gauklerManningStrickler(discharge_,constManningN_,bedSlope_);
```
We initialize the analytic solution to a verctor of the appropriate size filled with zeros.
```cpp
        exactWaterDepth_.resize(gridGeometry->numDofs(), 0.0);
        exactVelocityX_.resize(gridGeometry->numDofs(), 0.0);
    }
```
Get the analytical water depth
```cpp
    const std::vector<Scalar>& getExactWaterDepth()
    {
        return exactWaterDepth_;
    }
```
Get the analytical velocity
```cpp
    const std::vector<Scalar>& getExactVelocityX()
    {
        return exactVelocityX_;
    }
```
Get the water depth with Gauckler-Manning-Strickler
```cpp
    Scalar gauklerManningStrickler(Scalar discharge, Scalar manningN, Scalar bedSlope)
    {
        using std::pow;
        using std::abs;
        using std::sqrt;

        return pow(abs(discharge)*manningN/sqrt(bedSlope), 0.6);
    }
```
Get the analytical solution
```cpp
    void analyticalSolution()
    {
        using std::abs;

        for (const auto& element : elements(this->gridGeometry().gridView()))
        {
            const Scalar h = this->gauklerManningStrickler(discharge_,constManningN_,bedSlope_);
            const Scalar u = abs(discharge_)/h;

            const auto eIdx = this->gridGeometry().elementMapper().index(element);
            exactWaterDepth_[eIdx] = h;
            exactVelocityX_[eIdx] = u;
        }
    }
```
Get the problem name. It is used as a prefix for files generated by the simulation.
```cpp
    const std::string& name() const
    {
        return name_;
    }
```
Get the source term.
```cpp
     NumEqVector source(const Element& element,
                        const FVElementGeometry& fvGeometry,
                        const ElementVolumeVariables& elemVolVars,
                        const SubControlVolume &scv) const
    {

        NumEqVector source (0.0);
```
In this model the bottom friction is the only source.
```cpp
        source += bottomFrictionSource(element, fvGeometry, elemVolVars, scv);

        return source;
    }
```
Get the source term due to bottom friction.
```cpp
     NumEqVector bottomFrictionSource(const Element& element,
                                      const FVElementGeometry& fvGeometry,
                                      const ElementVolumeVariables& elemVolVars,
                                      const SubControlVolume &scv) const
     {
        NumEqVector bottomFrictionSource(0.0);
        const auto& volVars = elemVolVars[scv];
```
For the calculation of the source term due to bottom friction the two-dimensional bottom shear stess vector is needed. This is the force per area, which works between the flow and the bed. It is calculated within the `FrictionLaw`, which is a spatialParameter. In this model the `FrictionLawManning` is used (see `params.input`).
```cpp
        Dune::FieldVector<Scalar, 2> bottomShearStress = this->spatialParams().frictionLaw(element, scv).shearStress(volVars);
```
The bottom shear stress causes a pure loss of momentum. Thus the first entry of the `bottomFrictionSource`, which is related to the mass balance equation is zero. The second entry of the `bottomFricitonSource` corresponds to the momentum equation in x-direction and is therefore equal to the first, the x-component, of the `bottomShearStress`. Accordingly the third entry of the `bottomFrictionSource` is equal to the second component of the `bottomShearStress`.
```cpp
        bottomFrictionSource[0] = 0.0;
        bottomFrictionSource[1] = bottomShearStress[0];
        bottomFrictionSource[2] = bottomShearStress[1];

        return bottomFrictionSource;
     }
```
We specify the boundary condition type.
```cpp
    BoundaryTypes boundaryTypesAtPos(const GlobalPosition &globalPos) const
    {
        BoundaryTypes bcTypes;
```
Since we use a weak imposition all boundary conditions are of Neumann type.
```cpp
        bcTypes.setAllNeumann();
        return bcTypes;
    }
```
We specify the neumann boundary. Due to the weak imposition we calculate the flux at the boundary, with a Rieman solver. For this the state of a virtual cell outside of the boundary is needed (`boundaryStateVariables`), wich is calculated with the Riemann invariants (see Yoon and Kang, Finite Volume Model for Two-Dimensional Shallow Water Flows on Unstructured Grids) . The calculation of the Riemann invariants differ depending on the type of the boundary (h, q or no-flow boundary).
```cpp
    NeumannFluxes neumann(const Element& element,
                          const FVElementGeometry& fvGeometry,
                          const ElementVolumeVariables& elemVolVars,
                          const ElementFluxVariablesCache& elemFluxVarsCache,
                          const SubControlVolumeFace& scvf) const
    {
        NeumannFluxes values(0.0);

        const auto& insideScv = fvGeometry.scv(scvf.insideScvIdx());
        const auto& insideVolVars = elemVolVars[insideScv];
        const auto& nxy = scvf.unitOuterNormal();
        const auto gravity = this->spatialParams().gravity(scvf.center());
        std::array<Scalar, 3> boundaryStateVariables;
```
Calculate the rieman invariants for imposed discharge at the left side.
```cpp
        if (scvf.center()[0] < 0.0 + eps_)
        {
            boundaryStateVariables = ShallowWater::fixedDischargeBoundary(discharge_,
                                                                          insideVolVars.waterDepth(),
                                                                          insideVolVars.velocity(0),
                                                                          insideVolVars.velocity(1),
                                                                          gravity,
                                                                          nxy);
        }
```
Calculate the rieman invariants for impose water depth at the right side.
```cpp
        else if (scvf.center()[0] > 100.0 - eps_)
        {
            boundaryStateVariables =  ShallowWater::fixedWaterDepthBoundary(hBoundary_,
                                                                            insideVolVars.waterDepth(),
                                                                            insideVolVars.velocity(0),
                                                                            insideVolVars.velocity(1),
                                                                            gravity,
                                                                            nxy);
        }
```
Calculate the rieman invarianty for the no-flow boundary.
```cpp
        else
        {
            boundaryStateVariables[0] = insideVolVars.waterDepth();
            boundaryStateVariables[1] = -insideVolVars.velocity(0);
            boundaryStateVariables[2] = -insideVolVars.velocity(1);
        }
```
We calculate the boundary fluxes based on a Riemann problem.
```cpp
        auto riemannFlux = ShallowWater::riemannProblem(insideVolVars.waterDepth(),
                                                        boundaryStateVariables[0],
                                                        insideVolVars.velocity(0),
                                                        boundaryStateVariables[1],
                                                        insideVolVars.velocity(1),
                                                        boundaryStateVariables[2],
                                                        insideVolVars.bedSurface(),
                                                        insideVolVars.bedSurface(),
                                                        gravity,
                                                        nxy);

        values[Indices::massBalanceIdx] = riemannFlux[0];
        values[Indices::velocityXIdx]   = riemannFlux[1];
        values[Indices::velocityYIdx]   = riemannFlux[2];

        return values;
    }
```
We set the initial conditions. In this example constant initial conditions are used. Therefore the argument `globalPos` is not needed. If you want to impose spatial variable initial conditions, you have to use the `globalPos`.
```cpp
    PrimaryVariables initialAtPos(const GlobalPosition &globalPos) const
    {
        PrimaryVariables values(0.0);
```
We set the initial water depth to one meter.
```cpp
        values[0] = 1.0;
```
We set the x-component of the initial velocity to zero.
```cpp
        values[1] = 0.0;
```
We set the y-component of the initial velocity to zero.
```cpp
        values[2] = 0.0;

        return values;
    };
```
\}
```cpp

private:
```
We declare the private variables of the problem. They are initialized in the problems constructor.
We declare the variable for the analytic solution.
```cpp
    std::vector<Scalar> exactWaterDepth_;
    std::vector<Scalar> exactVelocityX_;
```
constant friction value. An analytic solution is only available for const friction. If you want to run the simulation with a non constant friciton value (specified in the spatialParams) you have to remove the analytic solution.
```cpp
    Scalar constManningN_;
```
The constant bed slope.
```cpp
    Scalar bedSlope_;
```
The water depth at the outflow boundary.
```cpp
    Scalar hBoundary_;
```
The discharge at the inflow boundary.
```cpp
    Scalar discharge_;
```
eps is used as a small value for the definition of the boundry conditions
```cpp
    static constexpr Scalar eps_ = 1.0e-6;
    std::string name_;
};
```
We leave the namespace Dumux.
```cpp
}

```



## The file `main.cc`

## The main file
This is the main file for the shallow water example. Here we can see the programme sequence and how the system is solved using newton's method.
### Includes
```cpp
#include <config.h>
```
Standard header file for C++, to get time and date information.
```cpp
#include <ctime>
```
Standard header file for C++, for in- and output.
```cpp
#include <iostream>
```
Dumux is based on DUNE, the Distributed and Unified Numerics Environment, which provides several grid managers and linear solvers. So we need some includes from that.
```cpp
#include <dune/common/parallel/mpihelper.hh>
#include <dune/common/timer.hh>
#include <dune/grid/io/file/dgfparser/dgfexception.hh>
#include <dune/grid/io/file/vtk.hh>
```
We need the following class to simplify the writing of dumux simulation data to VTK format.
```cpp
#include <dumux/io/vtkoutputmodule.hh>
```
In Dumux a property system is used to specify the model. For this, different properties are defined containing type definitions, values and methods. All properties are declared in the file properties.hh.
```cpp
#include <dumux/common/properties.hh>
```
The following file contains the parameter class, which manages the definition of input parameters by a default value, the inputfile or the command line.
```cpp
#include <dumux/common/parameters.hh>
```
The file dumuxmessage.hh contains the class defining the start and end message of the simulation.
```cpp
#include <dumux/common/dumuxmessage.hh>
#include <dumux/common/defaultusagemessage.hh>
```
The gridmanager constructs a grid from the information in the input or grid file. There is a specification for the different supported grid managers.
```cpp
#include <dumux/io/grid/gridmanager.hh>
```
We include the linear solver to be used to solve the linear system
```cpp
#include <dumux/linear/amgbackend.hh>
```
We include the nonlinear newtons method
```cpp
#include <dumux/nonlinear/newtonsolver.hh>
```
Further we include assembler, which assembles the linear systems for finite volume schemes (box-scheme, tpfa-approximation, mpfa-approximation)
```cpp
#include <dumux/assembly/fvassembler.hh>
```
We include the problem file which defines initial and boundary conditions to describe our example problem
```cpp
#include "problem.hh"
```
### Beginning of the main function
```cpp
int main(int argc, char** argv) try
{
    using namespace Dumux;
```
We define the type tag for this problem
```cpp
    using TypeTag = Properties::TTag::RoughChannel;
```
We initialize MPI, finalize is done automatically on exit
```cpp
    const auto& mpiHelper = Dune::MPIHelper::instance(argc, argv);
```
We print dumux start message
```cpp
    if (mpiHelper.rank() == 0)
        DumuxMessage::print(/*firstCall=*/true);
```
We parse command line arguments and input file
```cpp
    Parameters::init(argc, argv);
```
### Create the grid
A gridmanager tries to create the grid either from a grid file or the input file.
```cpp
    GridManager<GetPropType<TypeTag, Properties::Grid>> gridManager;
    gridManager.init();
```
We compute on the leaf grid view
```cpp
    const auto& leafGridView = gridManager.grid().leafGridView();
```
### Setup and solving of the problem
#### Setup
We create and initialize the finite volume grid geometry, the problem, the linear system, including the jacobian matrix, the residual and the solution vector and the gridvariables.
We need the finite volume geometry to build up the subcontrolvolumes (scv) and subcontrolvolume faces (scvf) for each element of the grid partition.
```cpp
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    auto gridGeometry = std::make_shared<GridGeometry>(leafGridView);
    gridGeometry->update();
```
In the problem, we define the boundary and initial conditions.
```cpp
    using Problem = GetPropType<TypeTag, Properties::Problem>;
    auto problem = std::make_shared<Problem>(gridGeometry);
```
We initialize the solution vector
```cpp
    using SolutionVector = GetPropType<TypeTag, Properties::SolutionVector>;
    SolutionVector x(gridGeometry->numDofs());
    problem->applyInitialSolution(x);
    auto xOld = x;
```
And then use the solutionvector to intialize the gridVariables.
```cpp
    using GridVariables = GetPropType<TypeTag, Properties::GridVariables>;
    auto gridVariables = std::make_shared<GridVariables>(problem, gridGeometry);
    gridVariables->init(x);
```
We get some time loop parameters from the input file.
```cpp
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    const auto tEnd = getParam<Scalar>("TimeLoop.TEnd");
    const auto maxDt = getParam<Scalar>("TimeLoop.MaxTimeStepSize");
    auto dt = getParam<Scalar>("TimeLoop.DtInitial");
```
We intialize the vtk output module. Each model has a predefined model specific output with relevant parameters for that model.
```cpp
    using IOFields = GetPropType<TypeTag, Properties::IOFields>;
    VtkOutputModule<GridVariables, SolutionVector> vtkWriter(*gridVariables,x, problem->name());
```
We add the analytical solution ("exactWaterDepth" and "exactVelocityX") to the predefined specific output.
```cpp
    vtkWriter.addField(problem->getExactWaterDepth(), "exactWaterDepth");
    vtkWriter.addField(problem->getExactVelocityX(), "exactVelocityX");
```
We calculate the analytic solution.
```cpp
    problem->analyticalSolution();
    IOFields::initOutputModule(vtkWriter);
    vtkWriter.write(0.0);
```
We instantiate time loop.
```cpp
    auto timeLoop = std::make_shared<CheckPointTimeLoop<Scalar>>(0, dt, tEnd);
    timeLoop->setMaxTimeStepSize(maxDt);
```
we set the assembler with the time loop because we have an instationary problem.
```cpp
    using Assembler = FVAssembler<TypeTag, DiffMethod::numeric>;
    auto assembler = std::make_shared<Assembler>(problem, gridGeometry, gridVariables, timeLoop);
```
We set the linear solver.
```cpp
    using LinearSolver = Dumux::AMGBackend<TypeTag>;
    auto linearSolver = std::make_shared<LinearSolver>(leafGridView, gridGeometry->dofMapper());
```
Additionaly, we set the non-linear solver.
```cpp
    using NewtonSolver = Dumux::NewtonSolver<Assembler, LinearSolver>;
    NewtonSolver nonLinearSolver(assembler, linearSolver);
```
We set some check point at the end of the time loop. The check point is used to trigger the vtk output.
```cpp
    timeLoop->setCheckPoint(tEnd);
```
We start the time loop.
```cpp
    timeLoop->start(); do
    {
```
We start to calculate the new solution of that time step. First we define the old solution as the solution of the previous time step for storage evaluations.
```cpp
        assembler->setPreviousSolution(xOld);
```
We solve the non-linear system with time step control.
```cpp
        nonLinearSolver.solve(x,*timeLoop);
```
We make the new solution the old solution.
```cpp
        xOld = x;
        gridVariables->advanceTimeStep();
```
We advance to the time loop to the next step.
```cpp
        timeLoop->advanceTimeStep();
```
We write vtk output, if we reached the check point (end of time loop)
```cpp
        if (timeLoop->isCheckPoint())
            vtkWriter.write(timeLoop->time());
```
We report statistics of this time step.
```cpp
        timeLoop->reportTimeStep();
```
We set new dt as suggested by newton controller for the next time step.
```cpp
        timeLoop->setTimeStepSize(nonLinearSolver.suggestTimeStepSize(timeLoop->timeStepSize()));


    } while (!timeLoop->finished());

    timeLoop->finalize(leafGridView.comm());
```
### Final Output
We print dumux end message.
```cpp
    if (mpiHelper.rank() == 0)
    {
        Parameters::print();
        DumuxMessage::print(/*firstCall=*/false);
    }

    return 0;
} // end main

catch (const Dumux::ParameterException &e)
{
    std::cerr << std::endl << e << " ---> Abort!" << std::endl;
    return 1;
}
catch (const Dune::DGFException & e)
{
    std::cerr << "DGF exception thrown (" << e <<
                 "). Most likely, the DGF file name is wrong "
                 "or the DGF file is corrupted, "
                 "e.g. missing hash at end of file or wrong number (dimensions) of entries."
                 << " ---> Abort!" << std::endl;
    return 2;
}
catch (const Dune::Exception &e)
{
    std::cerr << "Dune reported error: " << e << " ---> Abort!" << std::endl;
    return 3;
}
catch (...)
{
    std::cerr << "Unknown exception thrown! ---> Abort!" << std::endl;
    return 4;
}
```

## Results
The solution and the analytical result for the problem is shown in ![Result Logo](img/result.png).
