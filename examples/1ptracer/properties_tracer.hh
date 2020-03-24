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
#ifndef DUMUX_TRACER_TEST_PROPERTIES_HH
#define DUMUX_TRACER_TEST_PROPERTIES_HH

// This file defines the `TypeTag` used for the tracer transport simulation, for
// which we then define the necessary properties.
//
// ### Include files
// As for the single-phase problem, a`TypeTag` is defined for this simulation.
// Here, we inherit all properties from the `Tracer` type tag, a convenience type tag
// that predefines most of the required properties for tracer transport flow simulations in DuMuX.
#include <dumux/porousmediumflow/tracer/model.hh>

// Again, we use YaspGrid, an implementation of the dune grid interface for structured grids:
#include <dune/grid/yaspgrid.hh>
// and the cell centered, two-point-flux discretization.
#include <dumux/discretization/cctpfa.hh>
// This includes the base class for fluid systems. We will define a custom fluid
// system that inherits from that class.
#include <dumux/material/fluidsystems/base.hh>

// We include the problem and spatial parameters headers used for this simulation.
#include "problem_tracer.hh"
#include "spatialparams_tracer.hh"

// ### Basic property definitions for the tracer transport problem
// We enter the namespace Dumux
namespace Dumux {

// In the following, we create a new tracer fluid system and derive from the base fluid system.
template<class TypeTag>
class TracerFluidSystem : public FluidSystems::Base<GetPropType<TypeTag, Properties::Scalar>,
                                                               TracerFluidSystem<TypeTag>>
{
    // We define some convenience aliases to be used inside this class.
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using Problem = GetPropType<TypeTag, Properties::Problem>;
    using GridView = GetPropType<TypeTag, Properties::GridView>;
    using Element = typename GridView::template Codim<0>::Entity;
    using FVElementGeometry = typename GetPropType<TypeTag, Properties::GridGeometry>::LocalView;
    using SubControlVolume = typename FVElementGeometry::SubControlVolume;

public:
    // We specify that the fluid system only contains tracer components,
    static constexpr bool isTracerFluidSystem()
    { return true; }

    // and that no component is the main component
    static constexpr int getMainComponent(int phaseIdx)
    { return -1; }

    // We define the number of components of this fluid system (one single tracer component)
    static constexpr int numComponents = 1;

    // This interface is designed to define the names of the components of the fluid system.
    // Here, we only have a single component, so `compIdx` should always be 0.
    // The component name is used for the vtk output.
    static std::string componentName(int compIdx = 0)
    { return "tracer_" + std::to_string(compIdx); }

    // We set the phase name for the phase index (`phaseIdx`) for velocity vtk output:
    // Here, we only have a single phase, so `phaseIdx` should always be zero.
    static std::string phaseName(int phaseIdx = 0)
    { return "Groundwater"; }

    // We set the molar mass of the tracer component with index `compIdx` (should again always be zero here).
    static Scalar molarMass(unsigned int compIdx = 0)
    { return 0.300; }

    // We set the value for the binary diffusion coefficient. This
    // might depend on spatial parameters like pressure / temperature.
    // But, in this case we neglect diffusion and return 0.0:
    static Scalar binaryDiffusionCoefficient(unsigned int compIdx,
                                             const Problem& problem,
                                             const Element& element,
                                             const SubControlVolume& scv)
    { return 0.0; }
};

// We enter the namespace Properties
namespace Properties {

// A `TypeTag` for our simulation is created which inherits from the tracer model and the
// cell centered discretization scheme using two-point flux approximation.
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

// The problem class that specifies initial and boundary conditions:
template<class TypeTag>
struct Problem<TypeTag, TTag::TracerTest> { using type = TracerTestProblem<TypeTag>; };

// We define the spatial parameters for our tracer simulation:
template<class TypeTag>
struct SpatialParams<TypeTag, TTag::TracerTest>
{
private:
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
public:
    using type = TracerTestSpatialParams<GridGeometry, Scalar>;
};

// One can choose between a formulation in terms of mass or mole fractions.
// Here, we are using mass fractions.
template<class TypeTag>
struct UseMoles<TypeTag, TTag::TracerTest> { static constexpr bool value = false; };

// We use solution-independent molecular diffusion coefficients. Per default, solution-dependent
// diffusion coefficients are assumed during the computation of the jacobian matrix entries. Specifying
// solution-independent diffusion coefficients can speed up computations:
template<class TypeTag>
struct SolutionDependentMolecularDiffusion<TypeTag, TTag::TracerTestCC>
{ static constexpr bool value = false; };

// We set the above created tracer fluid system:
template<class TypeTag>
struct FluidSystem<TypeTag, TTag::TracerTest> { using type = TracerFluidSystem<TypeTag>; };

// We leave the namespace Properties and Dumux.
} // end namespace Properties
} // end namespace Dumux

#endif
