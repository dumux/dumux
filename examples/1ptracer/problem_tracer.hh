// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   See the file COPYING for full copying permissions.                      *
 *                                                                           *
 *   This program is free software: you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation, either version 3 of the License, or       *
 *   (at your option) any later version.                                     *
 *                                                                           *
 *   This program is distributed in the hope that it will be useful,         *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of          *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the            *
 *   GNU General Public License for more details.                            *
 *                                                                           *
 *   You should have received a copy of the GNU General Public License       *
 *   along with this program.  If not, see <http://www.gnu.org/licenses/>.   *
 *****************************************************************************/

// ### Header guard
#ifndef DUMUX_TRACER_TEST_PROBLEM_HH
#define DUMUX_TRACER_TEST_PROBLEM_HH

//Before we enter the problem class containing initial and boundary conditions, we include necessary files and introduce properties.
// ### Include files
// Again, we have to include the dune grid interface:
#include <dune/grid/yaspgrid.hh>
// and the cell centered, two-point-flux discretization.
#include <dumux/discretization/cctpfa.hh>
// Here, we include the tracer model:
#include <dumux/porousmediumflow/tracer/model.hh>
// We include again the porous medium problem class that this class is derived from:
#include <dumux/porousmediumflow/problem.hh>
// and the base fluidsystem
#include <dumux/material/fluidsystems/base.hh>
// We include the header that specifies all spatially variable parameters for the tracer problem:
#include "spatialparams_tracer.hh"

// ### Define basic properties for our simulation
// We enter the namespace Dumux
namespace Dumux {

// The problem class is forward declared:
template <class TypeTag>
class TracerTestProblem;

// We enter the namespace Properties,
namespace Properties {
// A `TypeTag` for our simulation is created which inherits from the tracer model and the
// cell centered, two-point-flux discretization scheme.
namespace TTag {
struct TracerTest { using InheritsFrom = std::tuple<Tracer>; };
struct TracerTestCC { using InheritsFrom = std::tuple<TracerTest, CCTpfaModel>; };
}

// We enable caching for the grid volume variables, the flux variables and the FV grid geometry.
template<class TypeTag>
struct EnableGridVolumeVariablesCache<TypeTag, TTag::TracerTest> { static constexpr bool value = true; };
template<class TypeTag>
struct EnableGridFluxVariablesCache<TypeTag, TTag::TracerTest> { static constexpr bool value = true; };
template<class TypeTag>
struct EnableGridGeometryCache<TypeTag, TTag::TracerTest> { static constexpr bool value = true; };

// We use the same grid as in the stationary one-phase model, a structured 2D grid:
template<class TypeTag>
struct Grid<TypeTag, TTag::TracerTest> { using type = Dune::YaspGrid<2>; };

// The problem class specifies initial and boundary conditions:
template<class TypeTag>
struct Problem<TypeTag, TTag::TracerTest> { using type = TracerTestProblem<TypeTag>; };

// We define the spatial parameters for our tracer simulation:
template<class TypeTag>
struct SpatialParams<TypeTag, TTag::TracerTest>
{
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using type = TracerTestSpatialParams<GridGeometry, Scalar>;
};

// We define that mass fractions are used to define the concentrations
template<class TypeTag>
struct UseMoles<TypeTag, TTag::TracerTest> { static constexpr bool value = false; };
// We do not use a solution dependent molecular diffusion coefficient:
template<class TypeTag>
struct SolutionDependentMolecularDiffusion<TypeTag, TTag::TracerTestCC> { static constexpr bool value = false; };

// In the following, we create a new tracer fluid system and derive it from the base fluid system.
template<class TypeTag>
class TracerFluidSystem : public FluidSystems::Base<GetPropType<TypeTag, Properties::Scalar>,
                                                               TracerFluidSystem<TypeTag>>
{
    // We define convenient shortcuts to the properties `Scalar`, `Problem`, `GridView`,
    // `Element`, `FVElementGeometry` and `SubControlVolume`:
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using Problem = GetPropType<TypeTag, Properties::Problem>;
    using GridView = GetPropType<TypeTag, Properties::GridView>;
    using Element = typename GridView::template Codim<0>::Entity;
    using FVElementGeometry = typename GetPropType<TypeTag, Properties::GridGeometry>::LocalView;
    using SubControlVolume = typename FVElementGeometry::SubControlVolume;

public:
    // We specify, that the fluid system only contains tracer components,
    static constexpr bool isTracerFluidSystem()
    { return true; }

    // that no component is the main component
    static constexpr int getMainComponent(int phaseIdx)
    { return -1; }

    // and the number of components
    static constexpr int numComponents = 1;

    // We set the component name for the component index (`compIdx`) for the vtk output:
    static std::string componentName(int compIdx)
    { return "tracer_" + std::to_string(compIdx); }

    //We set the phase name for the phase index (`phaseIdx`) for velocity vtk output:
    static std::string phaseName(int phaseIdx = 0)
    { return "Groundwater"; }

    // We set the molar mass of the tracer component with index `compIdx`.
    static Scalar molarMass(unsigned int compIdx)
    { return 0.300; }

    // We set the value for the binary diffusion coefficient. This
    // might depend on spatial parameters like pressure / temperature. But for our case it is 0.0:
    static Scalar binaryDiffusionCoefficient(unsigned int compIdx,
                                             const Problem& problem,
                                             const Element& element,
                                             const SubControlVolume& scv)
    { return 0.0; }
};

// We set the above created tracer fluid system:
template<class TypeTag>
struct FluidSystem<TypeTag, TTag::TracerTest> { using type = TracerFluidSystem<TypeTag>; };

// We leave the namespace Properties.
}


// ### The problem class
// We enter the problem class where all necessary boundary conditions and initial conditions are set for our simulation.
// As this is a porous medium problem, we inherit from the basic `PorousMediumFlowProblem`.
template <class TypeTag>
class TracerTestProblem : public PorousMediumFlowProblem<TypeTag>
{
    // We use convenient declarations that we derive from the property system.
    using ParentType = PorousMediumFlowProblem<TypeTag>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using Indices = typename GetPropType<TypeTag, Properties::ModelTraits>::Indices;
    using GridView = GetPropType<TypeTag, Properties::GridView>;
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using BoundaryTypes = GetPropType<TypeTag, Properties::BoundaryTypes>;
    using PrimaryVariables = GetPropType<TypeTag, Properties::PrimaryVariables>;
    using FluidSystem = GetPropType<TypeTag, Properties::FluidSystem>;
    using SpatialParams = GetPropType<TypeTag, Properties::SpatialParams>;
    using Element = typename GridGeometry::GridView::template Codim<0>::Entity;
    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;
    using NumEqVector = GetPropType<TypeTag, Properties::NumEqVector>;
    using GridVariables = GetPropType<TypeTag, Properties::GridVariables>;
    using ElementVolumeVariables = typename GridVariables::GridVolumeVariables::LocalView;
    using ElementFluxVariablesCache = typename GridVariables::GridFluxVariablesCache::LocalView;
    using FVElementGeometry = typename GetPropType<TypeTag, Properties::GridGeometry>::LocalView;
    using SubControlVolumeFace = typename FVElementGeometry::SubControlVolumeFace;
    //We create a bool saying whether mole or mass fractions are used
    static constexpr bool useMoles = getPropValue<TypeTag, Properties::UseMoles>();
    //We create additional int to make dimWorld and numComponents available in the problem
    static constexpr int dimWorld = GridView::dimensionworld;
    static const int numComponents = FluidSystem::numComponents;

public:
    // This is the constructor of our problem class:
    TracerTestProblem(std::shared_ptr<const GridGeometry> gridGeometry)
    : ParentType(gridGeometry)
    {
        // We print out whether mole or mass fractions are used
        if(useMoles)
            std::cout<<"problem uses mole fractions" << '\n';
        else
            std::cout<<"problem uses mass fractions" << '\n';
    }

    // We define the type of boundary conditions depending on the location.
    // All boundaries are set to a neumann-type flow boundary condition.
    BoundaryTypes boundaryTypesAtPos(const GlobalPosition &globalPos) const
    {
        BoundaryTypes values;
        values.setAllNeumann();
        return values;
    }

    // We specify the initial conditions for the primary variable (tracer concentration) depending on the location.
    PrimaryVariables initialAtPos(const GlobalPosition &globalPos) const
    {
        PrimaryVariables initialValues(0.0);
        // The tracer concentration is located on the domain bottom:
        if (globalPos[1] < 0.1 + eps_)
        {
            // We assign different values, depending on wether mole concentrations or mass concentrations are used:
            if (useMoles)
                initialValues = 1e-9;
            else
                initialValues = 1e-9*FluidSystem::molarMass(0)/this->spatialParams().fluidMolarMass(globalPos);
        }
        return initialValues;
    }

        //We implement an outflow boundary on the top of the domain and prescribe zero-flux Neumann boundary conditions on all other boundaries.
        NumEqVector neumann(const Element& element,
                            const FVElementGeometry& fvGeometry,
                            const ElementVolumeVariables& elemVolVars,
                            const ElementFluxVariablesCache& elemFluxVarsCache,
                            const SubControlVolumeFace& scvf) const
        {
            NumEqVector values(0.0);
            const auto& volVars = elemVolVars[scvf.insideScvIdx()];
            const auto& globalPos = scvf.center();

            // This is the outflow boundary, where tracer is transported by advection
            // with the given flux field and diffusive flux is enforced to be zero
            if (globalPos[dimWorld-1] > this->gridGeometry().bBoxMax()[dimWorld-1] - eps_)
            {
                values = this->spatialParams().volumeFlux(element, fvGeometry, elemVolVars, scvf)
                         * volVars.massFraction(0, 0) * volVars.density(0)
                         / scvf.area();
                assert(values>=0.0 && "Volume flux at outflow boundary is expected to have a positive sign");
            }
            //prescribe zero-flux Neumann boundary conditions elsewhere
            else
                values = 0.0;

            return values;
        }

private:
// We assign a private global variable for the epsilon:
    static constexpr Scalar eps_ = 1e-6;

// This is everything the tracer problem class contains.
};

// We leave the namespace Dumux here.
} // end namespace Dumux

#endif
