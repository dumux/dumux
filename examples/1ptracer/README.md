This tutorial was copied from dumux/test/porousmediumflow/tracer/1ptracer.

# One-phase flow with random permeability distribution and a tracer model

## Problem set-up
This example contains a 2d simulation of a stationary groundwater flow. The permeability is distributed randomly. An contaminant initially concentrated at the domain bottom gets transported by the base groundwater flow.

The image below shows the simulation set-up. A pressure gradient between the top an the bottom boundary leads to a groundwater flux from the bottom to the top. Neumann no-flow boundaries are assigned to the left and right boundary. Initially, there is a contaminant concentration at the domain bottom.

 <img src="Plots/setup.png" width="500">

## Model description
Two different models are applied to simulate the system: In the first step, the groundwater velocity is evaluated under stationary conditions. Therefore the single phase model is applied. In a second step, the contaminant gets transported based on the volume fluxes of the single phase flow. It is assumed, that the dissolved contaminant does not affect density and viscosity of the groundwater. The tracer model is solved instationarily.

### 1p Model
The single phase model uses Darcy's law as the equation for the conservation of momentum:

$` \textbf v = - \frac{\textbf K}{\mu} \left(\textbf{grad}\, p - \varrho {\textbf g} \right) `$

With the darcy velocity $` \textbf v `$, the permeability $` \textbf K`$, the viscosity $` \mu`$, the pressure $`p`$, the density $`\rho`$ and the gravity $`\textbf g`$.
Darcy's law is inserted into the continuity equation:

$` \phi \frac{\partial \varrho}{\partial t} + \text{div} \textbf v = 0`$

with density $`\rho`$. It is solved for the pressure as primary variable. This equation is discretized using a cell-centered finite volume scheme as spatial discretization. For details on the discretization, we refer to the dumux handbook.

### Tracer Model
With the velocity field $`\textbf v`$ the transport of the contaminant component $`\kappa`$ is described by the following equation:

$` \phi \frac{ \partial X^\kappa}{\partial t} - \text{div} \left\lbrace X^\kappa {\textbf v}+ D^\kappa_\text{pm} \frac{M^\kappa}{M_\alpha} \textbf{grad} x^\kappa \right\rbrace = 0 `$

With the porosity $`\phi`$, the mass fraction of the component $`\kappa`$: $`X^\kappa`$, the binary diffusion coefficient in the porous medium $` D^\kappa_\text{pm} `$, the molar masses $` M `$ of the component $`\kappa`$ and the phase $`\alpha`$ and the mole fraction $`x`$.
The primary variable of this model is the mass fraction $`X^\kappa`$.

The porous medium diffusivity is yield out of the diffusion coefficient of the component, the porosity $`\phi `$ and the porous medium tortuosity $`\tau`$ by the following equation:

$`
D^\kappa_\text{pm}= \phi \tau D^\kappa
`$


## The file `spatialparams_1p.hh`


We want to generate a random permeability field. For this we use a random number generation of the C++ standard library.
```cpp
#include <random>
```
In the file properties.hh all properties are declared.
```cpp
#include <dumux/porousmediumflow/properties.hh>
```
We include the spatial parameters for single-phase, finite volumes from which we will inherit.
```cpp
#include <dumux/material/spatialparams/fv1p.hh>

namespace Dumux {
```
In the OnePTestSpatialParams class, we define all functions needed to describe the porous matrix, e.g. define porosity and permeability for the 1p_problem.
```cpp

template<class FVGridGeometry, class Scalar>
class OnePTestSpatialParams
: public FVSpatialParamsOneP<FVGridGeometry, Scalar,
                             OnePTestSpatialParams<FVGridGeometry, Scalar>>
{
```
We introduce using declarations that are derived from the property system which we need in this class.
```cpp
    using GridView = typename FVGridGeometry::GridView;
    using FVElementGeometry = typename FVGridGeometry::LocalView;
    using SubControlVolume = typename FVElementGeometry::SubControlVolume;
    using Element = typename GridView::template Codim<0>::Entity;
    using ParentType = FVSpatialParamsOneP<FVGridGeometry, Scalar,
                                           OnePTestSpatialParams<FVGridGeometry, Scalar>>;

    static constexpr int dimWorld = GridView::dimensionworld;
    using GlobalPosition = typename SubControlVolume::GlobalPosition;

public:
    using PermeabilityType = Scalar;
    OnePTestSpatialParams(std::shared_ptr<const FVGridGeometry> fvGridGeometry)
    : ParentType(fvGridGeometry), K_(fvGridGeometry->gridView().size(0), 0.0)
    {
```
### Generation of the random permeability field
We get the permeability of the domain and the lens from the params.input file.
```cpp
        permeability_ = getParam<Scalar>("SpatialParams.Permeability");
        permeabilityLens_ = getParam<Scalar>("SpatialParams.PermeabilityLens");
```
Further, we get the position of the lens, which is defined by the position of the lower left and the upper right corner.
```cpp
        lensLowerLeft_ = getParam<GlobalPosition>("SpatialParams.LensLowerLeft");
        lensUpperRight_ =getParam<GlobalPosition>("SpatialParams.LensUpperRight");
```
We generate random fields for the permeability using a lognormal distribution, with the permeability_ as mean value and 10 % of it as standard deviation. A seperate distribution is used for the lens using permeabilityLens_.
```cpp
        std::mt19937 rand(0);
        std::lognormal_distribution<Scalar> K(std::log(permeability_), std::log(permeability_)*0.1);
        std::lognormal_distribution<Scalar> KLens(std::log(permeabilityLens_), std::log(permeabilityLens_)*0.1);
        for (const auto& element : elements(fvGridGeometry->gridView()))
        {
            const auto eIdx = fvGridGeometry->elementMapper().index(element);
            const auto globalPos = element.geometry().center();
            K_[eIdx] = isInLens_(globalPos) ? KLens(rand) : K(rand);
        }
    }
```
### Properties of the porous matrix
We define the (intrinsic) permeability \f$[m^2]\f$ using the generated random permeability field. In this test, we use element-wise distributed permeabilities.
```cpp
    template<class ElementSolution>
    const PermeabilityType& permeability(const Element& element,
                                         const SubControlVolume& scv,
                                         const ElementSolution& elemSol) const
    {
        return K_[scv.dofIndex()];
    }

```
We set the porosity \f$[-]\f$ for the whole domain.
```cpp
    Scalar porosityAtPos(const GlobalPosition &globalPos) const
    { return 0.2; }
```
We reference to the permeability field. This is used in the main function to write an output for the permeability field.
```cpp
    const std::vector<Scalar>& getKField() const
    { return K_; }

private:
```
We have a convenience definition of the position of the lens.
```cpp
    bool isInLens_(const GlobalPosition &globalPos) const
    {
        for (int i = 0; i < dimWorld; ++i) {
            if (globalPos[i] < lensLowerLeft_[i] + eps_ || globalPos[i] > lensUpperRight_[i] - eps_)
                return false;
        }
        return true;
    }

    GlobalPosition lensLowerLeft_;
    GlobalPosition lensUpperRight_;

    Scalar permeability_, permeabilityLens_;

    std::vector<Scalar> K_;

    static constexpr Scalar eps_ = 1.5e-7;
};

} // end namespace Dumux

```



## The file `problem_1p.hh`


```cpp
#include <dune/grid/yaspgrid.hh>

#include <dumux/discretization/cctpfa.hh>
#include <dumux/porousmediumflow/1p/model.hh>
#include <dumux/porousmediumflow/problem.hh>
#include <dumux/material/components/simpleh2o.hh>
#include <dumux/material/fluidsystems/1pliquid.hh>
#include <dumux/porousmediumflow/1p/incompressiblelocalresidual.hh>

#include "spatialparams_1p.hh"

namespace Dumux {

/*!
 * \ingroup TracerTests
 * \brief The properties for the incompressible test
 */
```
forward declarations
```cpp
template<class TypeTag>
class OnePTestProblem;

namespace Properties {
```
Create new type tags
```cpp
namespace TTag {
struct IncompressibleTest { using InheritsFrom = std::tuple<OneP, CCTpfaModel>; };
} // end namespace TTag
```
Set the grid type
```cpp
template<class TypeTag>
struct Grid<TypeTag, TTag::IncompressibleTest> { using type = Dune::YaspGrid<2>; };
```
Set the problem type
```cpp
template<class TypeTag>
struct Problem<TypeTag, TTag::IncompressibleTest> { using type = OnePTestProblem<TypeTag>; };

template<class TypeTag>
struct SpatialParams<TypeTag, TTag::IncompressibleTest>
{
    using FVGridGeometry = GetPropType<TypeTag, Properties::FVGridGeometry>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using type = OnePTestSpatialParams<FVGridGeometry, Scalar>;
};

template<class TypeTag>
struct LocalResidual<TypeTag, TTag::IncompressibleTest> { using type = OnePIncompressibleLocalResidual<TypeTag>; };
```
the fluid system
```cpp
template<class TypeTag>
struct FluidSystem<TypeTag, TTag::IncompressibleTest>
{
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using type = FluidSystems::OnePLiquid<Scalar, Components::SimpleH2O<Scalar> >;
};
```
Enable caching
```cpp
template<class TypeTag>
struct EnableGridVolumeVariablesCache<TypeTag, TTag::IncompressibleTest> { static constexpr bool value = true; };
template<class TypeTag>
struct EnableGridFluxVariablesCache<TypeTag, TTag::IncompressibleTest> { static constexpr bool value = true; };
template<class TypeTag>
struct EnableFVGridGeometryCache<TypeTag, TTag::IncompressibleTest> { static constexpr bool value = true; };
} // end namespace Properties

template<class TypeTag>
class OnePTestProblem : public PorousMediumFlowProblem<TypeTag>
{
    using ParentType = PorousMediumFlowProblem<TypeTag>;
    using GridView = GetPropType<TypeTag, Properties::GridView>;
    using Element = typename GridView::template Codim<0>::Entity;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using PrimaryVariables = GetPropType<TypeTag, Properties::PrimaryVariables>;
    using FVElementGeometry = typename GetPropType<TypeTag, Properties::FVGridGeometry>::LocalView;
    using SubControlVolumeFace = typename FVElementGeometry::SubControlVolumeFace;
    using FVGridGeometry = GetPropType<TypeTag, Properties::FVGridGeometry>;
    using BoundaryTypes = GetPropType<TypeTag, Properties::BoundaryTypes>;

    static constexpr int dimWorld = GridView::dimensionworld;

public:
    OnePTestProblem(std::shared_ptr<const FVGridGeometry> fvGridGeometry)
    : ParentType(fvGridGeometry) {}

    /*!
     * \brief Specifies which kind of boundary condition should be
     *        used for which equation on a given boundary control volume.
     *
     * \param element The finite element
     * \param scvf The sub-control volume face
     */
    BoundaryTypes boundaryTypes(const Element &element,
                                const SubControlVolumeFace &scvf) const
    {
        BoundaryTypes values;
        const auto globalPos = scvf.ipGlobal();

        Scalar eps = 1.0e-6;
        if (globalPos[dimWorld-1] < eps || globalPos[dimWorld-1] > this->fvGridGeometry().bBoxMax()[dimWorld-1] - eps)
            values.setAllDirichlet();
        else
            values.setAllNeumann();

        return values;
    }

    /*!
     * \brief Evaluates the boundary conditions for a Dirichlet control volume.
     *
     * \param element The finite element
     * \param scvf The sub-control volume face
     *
     * For this method, the \a values parameter stores primary variables.
     */
    PrimaryVariables dirichlet(const Element &element,
                               const SubControlVolumeFace &scvf) const
    {
        const auto& pos = scvf.ipGlobal();
        PrimaryVariables values(0);
        values[0] = 1.0e+5*(1.1 - pos[dimWorld-1]*0.1);
        return values;
    }

    /*!
     * \brief Returns the temperature \f$\mathrm{[K]}\f$ for an isothermal problem.
     *
     * This is not specific to the discretization. By default it just
     * throws an exception so it must be overloaded by the problem if
     * no energy equation is used.
     */
    Scalar temperature() const
    {
        return 283.15; // 10°C
    }

};

} // end namespace Dumux

```



## The file `spatialparams_tracer.hh`


In the file properties.hh all properties are declared.
```cpp
#include <dumux/porousmediumflow/properties.hh>
```
As in the 1p spatialparams we inherit from the spatial parameters for single-phase, finite volumes, which we include here.
```cpp
#include <dumux/material/spatialparams/fv1p.hh>

namespace Dumux {
```
In the TracerTestSpatialParams class, we define all functions needed to describe spatially dependent parameters for the tracer_problem.
```cpp

template<class FVGridGeometry, class Scalar>
class TracerTestSpatialParams
: public FVSpatialParamsOneP<FVGridGeometry, Scalar,
                             TracerTestSpatialParams<FVGridGeometry, Scalar>>
{
    using GridView = typename FVGridGeometry::GridView;
    using FVElementGeometry = typename FVGridGeometry::LocalView;
    using SubControlVolume = typename FVElementGeometry::SubControlVolume;
    using SubControlVolumeFace = typename FVElementGeometry::SubControlVolumeFace;
    using Element = typename GridView::template Codim<0>::Entity;
    using ParentType = FVSpatialParamsOneP<FVGridGeometry, Scalar,
                                           TracerTestSpatialParams<FVGridGeometry, Scalar>>;

    static const int dimWorld = GridView::dimensionworld;
    using GlobalPosition = typename Dune::FieldVector<Scalar, dimWorld>;

public:

    TracerTestSpatialParams(std::shared_ptr<const FVGridGeometry> fvGridGeometry)
    : ParentType(fvGridGeometry) {}
```
### Properties of the porous matrix
We define the same porosity for the whole domain as in the 1p spatialparams.
```cpp
    Scalar porosityAtPos(const GlobalPosition& globalPos) const
    { return 0.2; }
```
We don't consider dispersivity for the tracer transport. So we set the dispersivity coefficient to zero.
```cpp
    template<class ElementSolution>
    Scalar dispersivity(const Element &element,
                        const SubControlVolume& scv,
                        const ElementSolution& elemSol) const
    { return 0; }
```
### Properties of the fluid system
In the following we define fluid properties that are spatial parameters in the tracer model. They can possible vary with space but are usually constants. Further spatially constant values of the fluid system are defined in the TracerFluidSystem class in problem.hh.
```cpp
```
We define the fluid density to a constant value of 1000 $`\frac{kg}{m^3}`$.
```cpp
    Scalar fluidDensity(const Element &element,
                        const SubControlVolume& scv) const
    { return 1000; }
```
We define the fluid molar mass.
```cpp
    Scalar fluidMolarMass(const Element &element,
                          const SubControlVolume& scv) const
    { return 18.0; }

    Scalar fluidMolarMass(const GlobalPosition &globalPos) const
    { return 18.0; }
```
### The volume fluxes
We define a function which returns the field of volume fluxes. This is e.g. used to calculate the transport of the tracer.
```cpp
    template<class ElementVolumeVariables>
    Scalar volumeFlux(const Element &element,
                      const FVElementGeometry& fvGeometry,
                      const ElementVolumeVariables& elemVolVars,
                      const SubControlVolumeFace& scvf) const
    {
        return volumeFlux_[scvf.index()];
    }
```
We define a function to set the volume flux. This is used in the main function to set the volume flux to the calculated value based on the solution of the 1p problem.
```cpp
    void setVolumeFlux(const std::vector<Scalar>& f)
    { volumeFlux_ = f; }

private:
    std::vector<Scalar> volumeFlux_;
};

} // end namespace Dumux

```



## The file `problem_tracer.hh`


```cpp
#include <dune/grid/yaspgrid.hh>

#include <dumux/discretization/cctpfa.hh>
#include <dumux/porousmediumflow/tracer/model.hh>
#include <dumux/porousmediumflow/problem.hh>
#include <dumux/material/fluidsystems/base.hh>

#include "spatialparams_tracer.hh"

namespace Dumux {
/**
 * \ingroup TracerTests
 * \brief Definition of a problem, for the tracer problem:
 * A rotating velocity field mixes a tracer band in a porous groundwater reservoir.
 */
template <class TypeTag>
class TracerTestProblem;

namespace Properties {
```
Create new type tags
```cpp
namespace TTag {
struct TracerTest { using InheritsFrom = std::tuple<Tracer>; };
struct TracerTestCC { using InheritsFrom = std::tuple<TracerTest, CCTpfaModel>; };
} // end namespace TTag
```
enable caching
```cpp
template<class TypeTag>
struct EnableGridVolumeVariablesCache<TypeTag, TTag::TracerTest> { static constexpr bool value = true; };
template<class TypeTag>
struct EnableGridFluxVariablesCache<TypeTag, TTag::TracerTest> { static constexpr bool value = true; };
template<class TypeTag>
struct EnableFVGridGeometryCache<TypeTag, TTag::TracerTest> { static constexpr bool value = true; };
```
Set the grid type
```cpp
template<class TypeTag>
struct Grid<TypeTag, TTag::TracerTest> { using type = Dune::YaspGrid<2>; };
```
Set the problem property
```cpp
template<class TypeTag>
struct Problem<TypeTag, TTag::TracerTest> { using type = TracerTestProblem<TypeTag>; };
```
Set the spatial parameters
```cpp
template<class TypeTag>
struct SpatialParams<TypeTag, TTag::TracerTest>
{
    using FVGridGeometry = GetPropType<TypeTag, Properties::FVGridGeometry>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using type = TracerTestSpatialParams<FVGridGeometry, Scalar>;
};
```
Define whether mole(true) or mass (false) fractions are used
```cpp
template<class TypeTag>
struct UseMoles<TypeTag, TTag::TracerTest> { static constexpr bool value = false; };
template<class TypeTag>
struct SolutionDependentMolecularDiffusion<TypeTag, TTag::TracerTestCC> { static constexpr bool value = false; };
```
! A simple fluid system with one tracer component
```cpp
template<class TypeTag>
class TracerFluidSystem : public FluidSystems::Base<GetPropType<TypeTag, Properties::Scalar>,
                                                               TracerFluidSystem<TypeTag>>
{
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using Problem = GetPropType<TypeTag, Properties::Problem>;
    using GridView = GetPropType<TypeTag, Properties::GridView>;
    using Element = typename GridView::template Codim<0>::Entity;
    using FVElementGeometry = typename GetPropType<TypeTag, Properties::FVGridGeometry>::LocalView;
    using SubControlVolume = typename FVElementGeometry::SubControlVolume;

public:
```
! If the fluid system only contains tracer components
```cpp
    static constexpr bool isTracerFluidSystem()
    { return true; }
```
! No component is the main component
```cpp
    static constexpr int getMainComponent(int phaseIdx)
    { return -1; }
```
! The number of components
```cpp
    static constexpr int numComponents = 1;
```
! Human readable component name (index compIdx) (for vtk output)
```cpp
    static std::string componentName(int compIdx)
    { return "tracer_" + std::to_string(compIdx); }
```
! Human readable phase name (index phaseIdx) (for velocity vtk output)
```cpp
    static std::string phaseName(int phaseIdx = 0)
    { return "Groundwater"; }
```
! Molar mass in kg/mol of the component with index compIdx
```cpp
    static Scalar molarMass(unsigned int compIdx)
    { return 0.300; }
```
! Binary diffusion coefficient
! (might depend on spatial parameters like pressure / temperature)
```cpp
    static Scalar binaryDiffusionCoefficient(unsigned int compIdx,
                                             const Problem& problem,
                                             const Element& element,
                                             const SubControlVolume& scv)
    { return 0.0; }
};

template<class TypeTag>
struct FluidSystem<TypeTag, TTag::TracerTest> { using type = TracerFluidSystem<TypeTag>; };

} // end namespace Properties


/*!
 * \ingroup TracerTests
 *
 * \brief Definition of a problem, for the tracer problem:
 * A lens of contaminant tracer is diluted by diffusion and a base groundwater flow
 *
 * This problem uses the \ref TracerModel model.
 *
 * To run the simulation execute the following line in shell:
 * <tt>./test_boxtracer -ParameterFile ./test_boxtracer.input</tt> or
 * <tt>./test_cctracer -ParameterFile ./test_cctracer.input</tt>
 */
template <class TypeTag>
class TracerTestProblem : public PorousMediumFlowProblem<TypeTag>
{
    using ParentType = PorousMediumFlowProblem<TypeTag>;

    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using Indices = typename GetPropType<TypeTag, Properties::ModelTraits>::Indices;
    using GridView = GetPropType<TypeTag, Properties::GridView>;
    using FVGridGeometry = GetPropType<TypeTag, Properties::FVGridGeometry>;
    using BoundaryTypes = GetPropType<TypeTag, Properties::BoundaryTypes>;
    using PrimaryVariables = GetPropType<TypeTag, Properties::PrimaryVariables>;
    using FluidSystem = GetPropType<TypeTag, Properties::FluidSystem>;
    using SpatialParams = GetPropType<TypeTag, Properties::SpatialParams>;
```
! property that defines whether mole or mass fractions are used
```cpp
    static constexpr bool useMoles = getPropValue<TypeTag, Properties::UseMoles>();

    using Element = typename FVGridGeometry::GridView::template Codim<0>::Entity;
    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;

public:
    TracerTestProblem(std::shared_ptr<const FVGridGeometry> fvGridGeom)
    : ParentType(fvGridGeom)
    {
```
stating in the console whether mole or mass fractions are used
```cpp
        if(useMoles)
            std::cout<<"problem uses mole fractions" << '\n';
        else
            std::cout<<"problem uses mass fractions" << '\n';
    }

    /*!
     * \name Boundary conditions
     */
```
\{
```cpp

    /*!
     * \brief Specifies which kind of boundary condition should be
     *        used for which equation on a given boundary segment.
     *
     * \param globalPos The position for which the bc type should be evaluated
     */
    BoundaryTypes boundaryTypesAtPos(const GlobalPosition &globalPos) const
    {
        BoundaryTypes values;
        values.setAllNeumann();
        return values;
    }
```
\}
```cpp

    /*!
     * \name Volume terms
     */
```
\{
```cpp

    /*!
     * \brief Evaluates the initial value for a control volume.
     *
     * \param globalPos The position for which the initial condition should be evaluated
     *
     * For this method, the \a values parameter stores primary
     * variables.
     */
    PrimaryVariables initialAtPos(const GlobalPosition &globalPos) const
    {
        PrimaryVariables initialValues(0.0);
        if (globalPos[1] < 0.1 + eps_)
        {
            if (useMoles)
                initialValues = 1e-9;
            else
                initialValues = 1e-9*FluidSystem::molarMass(0)/this->spatialParams().fluidMolarMass(globalPos);
        }
        return initialValues; }
```
\}
```cpp

private:
    static constexpr Scalar eps_ = 1e-6;
};

} // end namespace Dumux

```



## The file `main.cc`

## The main file
We look now at the main file for the tracer problem. We setup two problems in this file and solve them sequentially, first the 1p problem and afterwards the tracer problem. The result of the 1p problem is the pressure distribution in the problem domain. We use it to calculate the volume fluxes, which act as an input for the tracer problem. Based on this volume fluxes, we calculate the transport of a tracer in the following tracer problem.
```cpp
```
### Includes
```cpp
#include <config.h>
```
We include both problems in the main file, the problem_1p.hh and the problem_tracer.hh.
```cpp
#include "problem_1p.hh"
#include "problem_tracer.hh"
```
Further we include a standard header file for C++, to get time and date information
```cpp
#include <ctime>
```
and another one for in- and output.
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
```
The following file contains the class, which defines the sequential linear solver backends.
```cpp
#include <dumux/linear/seqsolverbackend.hh>
```
Further we include assembler, which assembles the linear systems for finite volume schemes (box-scheme, tpfa-approximation, mpfa-approximation).
```cpp
#include <dumux/assembly/fvassembler.hh>
```
The containing class in the following file defines the different differentiation methods used to compute the derivatives of the residual.
```cpp
#include <dumux/assembly/diffmethod.hh>
```
We need the following class to simplify the writing of dumux simulation data to VTK format.
```cpp
#include <dumux/io/vtkoutputmodule.hh>
```
The gridmanager constructs a grid from the information in the input or grid file. There is a specification for the different supported grid managers.
```cpp
#include <dumux/io/grid/gridmanager.hh>
```
### Beginning of the main function
```cpp
int main(int argc, char** argv) try
{
    using namespace Dumux;
```
We define the type tags for the two problems. They are created in the individual problem files.
```cpp
    using OnePTypeTag = Properties::TTag::IncompressibleTest;
    using TracerTypeTag = Properties::TTag::TracerTestCC;
```
We initialize MPI. Finalization is done automatically on exit.
```cpp
    const auto& mpiHelper = Dune::MPIHelper::instance(argc, argv);
```
We print the dumux start message.
```cpp
    if (mpiHelper.rank() == 0)
        DumuxMessage::print(/*firstCall=*/true);
```
We parse the command line arguments.
```cpp
    Parameters::init(argc, argv);
```
### Create the grid
```cpp
```
A gridmanager tries to create the grid either from a grid file or the input file. We solve both problems on the same grid. Hence, the grid is only created once using the type tag of the 1p problem.
```cpp
    GridManager<GetPropType<OnePTypeTag, Properties::Grid>> gridManager;
    gridManager.init();
```
We compute on the leaf grid view.
```cpp
    const auto& leafGridView = gridManager.grid().leafGridView();
```
### Setup and solving of the 1p problem
In the following section, we setup and solve the 1p problem. As the result of this problem, we obtain the pressure distribution in the problem domain.
```cpp
```
#### Setup
We create and initialize the finite volume grid geometry, the problem, the linear system, including the jacobian matrix, the residual and the solution vector and the gridvariables.
```cpp
```
We need the finite volume geometry to build up the subcontrolvolumes (scv) and subcontrolvolume faces (scvf) for each element of the grid partition.
```cpp
    using FVGridGeometry = GetPropType<OnePTypeTag, Properties::FVGridGeometry>;
    auto fvGridGeometry = std::make_shared<FVGridGeometry>(leafGridView);
    fvGridGeometry->update();
```
In the problem, we define the boundary and initial conditions.
```cpp
    using OnePProblem = GetPropType<OnePTypeTag, Properties::Problem>;
    auto problemOneP = std::make_shared<OnePProblem>(fvGridGeometry);
```
The jacobian matrix A, the solution vector p and the residual r are parts of the linear system.
```cpp
    using JacobianMatrix = GetPropType<OnePTypeTag, Properties::JacobianMatrix>;
    using SolutionVector = GetPropType<OnePTypeTag, Properties::SolutionVector>;
    SolutionVector p(leafGridView.size(0));

    auto A = std::make_shared<JacobianMatrix>();
    auto r = std::make_shared<SolutionVector>();
```
The grid variables store variables on scv and scvf (volume and flux variables).
```cpp
    using OnePGridVariables = GetPropType<OnePTypeTag, Properties::GridVariables>;
    auto onePGridVariables = std::make_shared<OnePGridVariables>(problemOneP, fvGridGeometry);
    onePGridVariables->init(p);
```
#### Assemling the linear system
We created and inizialize the assembler.
```cpp
    using OnePAssembler = FVAssembler<OnePTypeTag, DiffMethod::analytic>;
    auto assemblerOneP = std::make_shared<OnePAssembler>(problemOneP, fvGridGeometry, onePGridVariables);
    assemblerOneP->setLinearSystem(A, r);
```
We assemble the local jacobian and the residual and stop the time needed, which is displayed in the terminal output, using the assembly timer. Further, we start the timer to evaluate the total time of the assembly, solving and updating.
```cpp
    Dune::Timer timer;
    Dune::Timer assemblyTimer; std::cout << "Assembling linear system ..." << std::flush;
    assemblerOneP->assembleJacobianAndResidual(p);
    assemblyTimer.stop(); std::cout << " took " << assemblyTimer.elapsed() << " seconds." << std::endl;
```
We want to solve Ax = -r.
```cpp
    (*r) *= -1.0;
```
#### Solution
We set the linear solver "UMFPack" as the linear solver. Afterwards we solve the linear system. The time needed to solve the system is recorded by the solverTimer and displayed in the terminal output.
```cpp
    using LinearSolver = UMFPackBackend;
    Dune::Timer solverTimer; std::cout << "Solving linear system ..." << std::flush;
    auto linearSolver = std::make_shared<LinearSolver>();
    linearSolver->solve(*A, p, *r);
    solverTimer.stop(); std::cout << " took " << solverTimer.elapsed() << " seconds." << std::endl;
```
#### Update and output
We update the grid variables with the new solution.
```cpp
    Dune::Timer updateTimer; std::cout << "Updating variables ..." << std::flush;
    onePGridVariables->update(p);
    updateTimer.elapsed(); std::cout << " took " << updateTimer.elapsed() << std::endl;

```
We initialize the vtkoutput. Each model has a predefined model specific output with relevant parameters for that model. We add the pressure data from the solution vector p and the permeability field as output data.
```cpp
    using GridView = GetPropType<OnePTypeTag, Properties::GridView>;
    Dune::VTKWriter<GridView> onepWriter(leafGridView);
    onepWriter.addCellData(p, "p");
    const auto& k = problemOneP->spatialParams().getKField();
    onepWriter.addCellData(k, "permeability");
    onepWriter.write("1p");
```
We stop the timer and display the total time of the simulation as well as the cumulative CPU time.
```cpp
    timer.stop();

    const auto& comm = Dune::MPIHelper::getCollectiveCommunication();
    std::cout << "Simulation took " << timer.elapsed() << " seconds on "
              << comm.size() << " processes.\n"
              << "The cumulative CPU time was " << timer.elapsed()*comm.size() << " seconds.\n";

```
### Computation of the volume fluxes
We use the results of the 1p problem to calculate the the volume fluxes in the model domain.
```cpp

    using Scalar =  GetPropType<OnePTypeTag, Properties::Scalar>;
    std::vector<Scalar> volumeFlux(fvGridGeometry->numScvf(), 0.0);

    using FluxVariables =  GetPropType<OnePTypeTag, Properties::FluxVariables>;
    auto upwindTerm = [](const auto& volVars) { return volVars.mobility(0); };
```
We iterate over all elements.
```cpp
    for (const auto& element : elements(leafGridView))
    {
        auto fvGeometry = localView(*fvGridGeometry);
        fvGeometry.bind(element);

        auto elemVolVars = localView(onePGridVariables->curGridVolVars());
        elemVolVars.bind(element, fvGeometry, p);

        auto elemFluxVars = localView(onePGridVariables->gridFluxVarsCache());
        elemFluxVars.bind(element, fvGeometry, elemVolVars);
```
We calculate the volume flux for every subcontrolvolume face, which is not on a Neumann boundary (is not on the boundary or is on a Dirichlet boundary).
```cpp

        for (const auto& scvf : scvfs(fvGeometry))
        {
            const auto idx = scvf.index();

            if (!scvf.boundary())
            {
                FluxVariables fluxVars;
                fluxVars.init(*problemOneP, element, fvGeometry, elemVolVars, scvf, elemFluxVars);
                volumeFlux[idx] = fluxVars.advectiveFlux(0, upwindTerm);
            }
            else
            {
                const auto bcTypes = problemOneP->boundaryTypes(element, scvf);
                if (bcTypes.hasOnlyDirichlet())
                {
                    FluxVariables fluxVars;
                    fluxVars.init(*problemOneP, element, fvGeometry, elemVolVars, scvf, elemFluxVars);
                    volumeFlux[idx] = fluxVars.advectiveFlux(0, upwindTerm);
                }
            }
        }
    }

```
### Setup and solving of the tracer problem
```cpp
```
#### Setup
```cpp
```
Similar to the 1p problem, we first create and initialize the problem.
```cpp
    using TracerProblem = GetPropType<TracerTypeTag, Properties::Problem>;
    auto tracerProblem = std::make_shared<TracerProblem>(fvGridGeometry);
```
We use the the volume fluxes calculated in the previous section as input for the tracer model.
```cpp
    tracerProblem->spatialParams().setVolumeFlux(volumeFlux);
```
We create and initialize the solution vector. As the tracer problem is transient, the initial solution defined in the problem is applied to the solution vector.
```cpp
    SolutionVector x(leafGridView.size(0));
    tracerProblem->applyInitialSolution(x);
    auto xOld = x;
```
We create and initialize the grid variables.
```cpp
    using GridVariables = GetPropType<TracerTypeTag, Properties::GridVariables>;
    auto gridVariables = std::make_shared<GridVariables>(tracerProblem, fvGridGeometry);
    gridVariables->init(x);
```
We read in some time loop parameters from the input file. The parameter tEnd defines the duration of the simulation, dt the initial time step size and maxDt the maximal time step size.
```cpp
    const auto tEnd = getParam<Scalar>("TimeLoop.TEnd");
    auto dt = getParam<Scalar>("TimeLoop.DtInitial");
    const auto maxDt = getParam<Scalar>("TimeLoop.MaxTimeStepSize");
```
We instantiate the time loop.
```cpp
    auto timeLoop = std::make_shared<CheckPointTimeLoop<Scalar>>(0.0, dt, tEnd);
    timeLoop->setMaxTimeStepSize(maxDt);
```
We create and inizialize the assembler with time loop for instationary problem.
```cpp
    using TracerAssembler = FVAssembler<TracerTypeTag, DiffMethod::analytic, /*implicit=*/false>;
    auto assembler = std::make_shared<TracerAssembler>(tracerProblem, fvGridGeometry, gridVariables, timeLoop);
    assembler->setLinearSystem(A, r);
```
We initialize the vtk output module and add a velocity output.
```cpp
    VtkOutputModule<GridVariables, SolutionVector> vtkWriter(*gridVariables, x, tracerProblem->name());
    using IOFields = GetPropType<TracerTypeTag, Properties::IOFields>;
    IOFields::initOutputModule(vtkWriter);
    using VelocityOutput = GetPropType<TracerTypeTag, Properties::VelocityOutput>;
    vtkWriter.addVelocityOutput(std::make_shared<VelocityOutput>(*gridVariables));
    vtkWriter.write(0.0);

```
For the time loop we set a check point.
```cpp
    timeLoop->setPeriodicCheckPoint(tEnd/10.0);
```
We start the time loop and calculate a new time step as long as tEnd is reached. In every single time step the problem is assembled and solved.
```cpp
    timeLoop->start(); do
    {
```
First we define the old solution as the solution of the previous time step for storage evaluations.
```cpp
        assembler->setPreviousSolution(xOld);
```
Then the linear system is assembled.
```cpp
        Dune::Timer assembleTimer;
        assembler->assembleJacobianAndResidual(x);
        assembleTimer.stop();
```
We solve the linear system A(xOld-xNew) = r.
```cpp
        Dune::Timer solveTimer;
        SolutionVector xDelta(x);
        linearSolver->solve(*A, xDelta, *r);
        solveTimer.stop();
```
We calculate the actual solution and update it in the grid variables.
```cpp
        updateTimer.reset();
        x -= xDelta;
        gridVariables->update(x);
        updateTimer.stop();
```
We display the statistics of the actual time step.
```cpp
        const auto elapsedTot = assembleTimer.elapsed() + solveTimer.elapsed() + updateTimer.elapsed();
        std::cout << "Assemble/solve/update time: "
                  <<  assembleTimer.elapsed() << "(" << 100*assembleTimer.elapsed()/elapsedTot << "%)/"
                  <<  solveTimer.elapsed() << "(" << 100*solveTimer.elapsed()/elapsedTot << "%)/"
                  <<  updateTimer.elapsed() << "(" << 100*updateTimer.elapsed()/elapsedTot << "%)"
                  <<  std::endl;
```
The new solution is defined as the old solution.
```cpp
        xOld = x;
        gridVariables->advanceTimeStep();
```
We advance the time loop to the next time step.
```cpp
        timeLoop->advanceTimeStep();
```
We write the Vtk output on check points.
```cpp
        if (timeLoop->isCheckPoint())
            vtkWriter.write(timeLoop->time());
```
We report the statistics of this time step.
```cpp
        timeLoop->reportTimeStep();
```
We set the time step size dt of the next time step.
```cpp
        timeLoop->setTimeStepSize(dt);

    } while (!timeLoop->finished());

    timeLoop->finalize(leafGridView.comm());

```
### Final Output
```cpp

    if (mpiHelper.rank() == 0)
        DumuxMessage::print(/*firstCall=*/false);

    return 0;

}
catch (Dumux::ParameterException &e)
{
    std::cerr << std::endl << e << " ---> Abort!" << std::endl;
    return 1;
}
catch (Dune::DGFException & e)
{
    std::cerr << "DGF exception thrown (" << e <<
                 "). Most likely, the DGF file name is wrong "
                 "or the DGF file is corrupted, "
                 "e.g. missing hash at end of file or wrong number (dimensions) of entries."
                 << " ---> Abort!" << std::endl;
    return 2;
}
catch (Dune::Exception &e)
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

The 1-model calculated a stationary pressure distribution. It is shown in the following figure:
 <img src="Plots/pressure.png" width="300">


The random permeability distribution generates the velocity profile shown in the left plot of the next figure. The image in the middle illustrates the tracer distribution after 2500s and the image on the right after 5000s.

| ![](Plots/VelocityProfile.png)| ![](Plots/Tracer_2500.png) | ![](Plots/Tracer_5000.png)|
|:---:|:---:|:---:|
| velocity profile| tracer concentration after 2500s | tracer concentration after 5000s |
