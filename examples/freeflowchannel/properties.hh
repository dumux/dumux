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

#ifndef DUMUX_EXAMPLES_FREEFLOW_CHANNEL_PROPERTIES_HH
#define DUMUX_EXAMPLES_FREEFLOW_CHANNEL_PROPERTIES_HH

// ## The file `properties.hh`
// In the following, we set the properties for our simulation
// (click [here](https://git.iws.uni-stuttgart.de/dumux-repositories/dumux-course/blob/master/slides/dumux-course-properties.pdf) for DuMux course slides on the property system).
// We start with includes
// <details><summary>Click to show the header includes</summary>
#include <dune/grid/yaspgrid.hh>

#include <dumux/common/properties.hh>
#include <dumux/discretization/staggered/freeflow/properties.hh>
#include <dumux/freeflow/navierstokes/model.hh>
#include <dumux/material/fluidsystems/1pliquid.hh>
#include <dumux/material/components/constant.hh>

#include "problem.hh"
// </details>
//
// Then we start setting the properties.
// 1. For every test problem, a new `TypeTag` has to be created, which is done within the namespace `TTag` (subnamespace of `Properties`). It inherits from the Navier-Stokes flow model and the staggered-grid discretization scheme.
// 2. The grid is chosen to be a two-dimensional Cartesian structured grid (`YaspGrid`).
// 3. We set the `FluidSystem` to be a one-phase liquid with a single component. The class `Component::Constant` refers to a component with constant fluid properties (density, viscosity, ...) that can be set via the input file in the group `[0.Component]` where the number is the identifier given as template argument to the class template `Component::Constant`.
// 4. The problem class `ChannelExampleProblem`, which is forward declared before we enter `namespace Dumux` and defined later in this file, is defined to be the problem used in this test problem (charaterized by the TypeTag `ChannelExample`). The fluid system, which contains information about the properties such as density, viscosity or diffusion coefficient of the fluid we're simulating, is set to a constant one phase liquid.
// 5. We enable caching for the following classes (which stores values that were already calculated for later usage and thus results in higher memory usage but improved CPU speed): the grid volume variables, the grid flux variables, the finite volume grid geometry.
// [[codeblock]]
namespace Dumux::Properties {

namespace TTag {
struct ChannelExample { using InheritsFrom = std::tuple<NavierStokes, StaggeredFreeFlowModel>; };
} // namespace TTag

template<class TypeTag>
struct Grid<TypeTag, TTag::ChannelExample> { using type = Dune::YaspGrid<2>; };

template<class TypeTag>
struct Problem<TypeTag, TTag::ChannelExample> { using type = Dumux::ChannelExampleProblem<TypeTag> ; };

template<class TypeTag>
struct FluidSystem<TypeTag, TTag::ChannelExample>
{
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using type = FluidSystems::OnePLiquid<Scalar, Components::Constant<1, Scalar> >;
};

template<class TypeTag>
struct EnableGridVolumeVariablesCache<TypeTag, TTag::ChannelExample> { static constexpr bool value = true; };

template<class TypeTag>
struct EnableGridFluxVariablesCache<TypeTag, TTag::ChannelExample> { static constexpr bool value = true; };

template<class TypeTag>
struct EnableGridGeometryCache<TypeTag, TTag::ChannelExample> { static constexpr bool value = true; };
} // end namespace Dumux::Properties
// [[/codeblock]]

#endif
