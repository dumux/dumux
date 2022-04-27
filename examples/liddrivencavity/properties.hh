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

#ifndef DUMUX_LIDDRIVENCAVITY_EXAMPLE_PROPERTIES_HH
#define DUMUX_LIDDRIVENCAVITY_EXAMPLE_PROPERTIES_HH

// ## Compile-time settings (`properties.hh`)
//
// In this file, the type tag used for this simulation is defined,
// for which we then specialize properties (compile time options) to the needs of the desired setup.
//
// [[content]]
//
// ### Includes
// [[details]] includes
//
// The single-phase flow Navier-Stokes equations are solved by coupling a momentum balance model to a mass balance model.
#include <dumux/freeflow/navierstokes/momentum/model.hh>
#include <dumux/freeflow/navierstokes/mass/1p/model.hh>
#include <dumux/multidomain/traits.hh>
#include <dumux/multidomain/staggeredfreeflow/couplingmanager.hh>

// We want to use `YaspGrid`, an implementation of the dune grid interface for structured grids:
#include <dune/grid/yaspgrid.hh>

// In this example, we want to discretize the momentum and mass balance equations with the staggered-grid
// scheme which is so far the only available option for free-flow models in DuMux. Velocities are defined on the
// element edges while pressures are defined on the element centers.
#include <dumux/discretization/fcstaggered.hh>
#include <dumux/discretization/cctpfa.hh>

// The fluid properties are specified in the following headers (we use a liquid with constant properties as the fluid phase):
#include <dumux/material/components/constant.hh>
#include <dumux/material/fluidsystems/1pliquid.hh>

// We include the problem header used for this simulation.
#include "problem.hh"
// [[/details]]
//
// ### Type tag definition
//
// We define a type tag for our simulation with the name `LidDrivenCavityExample`. As we are dealing with a coupled model,
// we also define type tags for the momentum and mass model and inherit from the respective physical models (`NavierStokesMomentum` and `NavierStokesMassOneP`)
// as well as from the appropriate spatial discretization schemes (`FaceCenteredStaggeredModel` and `CCTpfaModel`).
// [[codeblock]]

namespace Dumux::Properties {

namespace TTag {
struct LidDrivenCavityExample {};
struct LidDrivenCavityExampleMomentum { using InheritsFrom = std::tuple<LidDrivenCavityExample, NavierStokesMomentum, FaceCenteredStaggeredModel>; };
struct LidDrivenCavityExampleMass { using InheritsFrom = std::tuple<LidDrivenCavityExample, NavierStokesMassOneP, CCTpfaModel>; };
} // end namespace TTag
// [[/codeblock]]

// ### Property specializations
//
// In the following piece of code, mandatory properties for which no meaningful
// default exist are specialized for our type tag `LidDrivenCavityExample`.
// [[codeblock]]
// This sets the fluid system to be used. Here, we use a liquid with constant properties as fluid phase.
template<class TypeTag>
struct FluidSystem<TypeTag, TTag::LidDrivenCavityExample>
{
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using type = FluidSystems::OnePLiquid<Scalar, Components::Constant<1, Scalar> >;
};

// We introduce the coupling manager to the properties system
template<class TypeTag>
struct CouplingManager<TypeTag, TTag::LidDrivenCavityExample>
{
private:
    using Traits = MultiDomainTraits<TTag::LidDrivenCavityExampleMomentum, TTag::LidDrivenCavityExampleMass>;
public:
    using type = StaggeredFreeFlowCouplingManager<Traits>;
};

// This sets the grid type used for the simulation. Here, we use a structured 2D grid.
template<class TypeTag>
struct Grid<TypeTag, TTag::LidDrivenCavityExample> { using type = Dune::YaspGrid<2>; };

// This sets our problem class (see problem.hh) containing initial and boundary conditions.
template<class TypeTag>
struct Problem<TypeTag, TTag::LidDrivenCavityExample> { using type = Dumux::LidDrivenCavityExampleProblem<TypeTag> ; };
// [[/codeblock]]

// We also set some properties related to memory management
// throughout the simulation.
// [[details]] caching properties
//
// In DuMux, one has the option to activate/deactivate the grid-wide caching of
// geometries and variables. If active, the CPU time can be significantly reduced
// as less dynamic memory allocation procedures are necessary. Per default, grid-wide
// caching is disabled to ensure minimal memory requirements, however, in this example we
// want to active all available caches, which significantly increases the memory
// demand but makes the simulation faster.
//
// [[codeblock]]
// This enables grid-wide caching of the volume variables.
template<class TypeTag>
struct EnableGridGeometryCache<TypeTag, TTag::LidDrivenCavityExample> { static constexpr bool value = true; };
//This enables grid wide caching for the flux variables.
template<class TypeTag>
struct EnableGridFluxVariablesCache<TypeTag, TTag::LidDrivenCavityExample> { static constexpr bool value = true; };
// This enables grid-wide caching for the finite volume grid geometry
template<class TypeTag>
struct EnableGridVolumeVariablesCache<TypeTag, TTag::LidDrivenCavityExample> { static constexpr bool value = true; };
} // end namespace Dumux::Properties
// [[/codeblock]]
// [[/details]]
// [[/content]]
#endif
