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
/*!
 * \file
 * \ingroup OnePTests
 * \brief The properties for the incompressible test
 */

#ifndef DUMUX_COMPRESSIBLE_ONEP_TEST_PROPERTIES_HH
#define DUMUX_COMPRESSIBLE_ONEP_TEST_PROPERTIES_HH

#include <dune/grid/yaspgrid.hh>

#include <dumux/material/fluidsystems/1pliquid.hh>

#include <dumux/discretization/cctpfa.hh>
#include <dumux/discretization/ccmpfa.hh>
#include <dumux/discretization/box.hh>
#include <dumux/discretization/fvgridvariables.hh>

#include <dumux/porousmediumflow/1p/model.hh>

#include "problem.hh"
#include "spatialparams.hh"

#ifndef EXPERIMENTAL
#define EXPERIMENTAL 0
#endif

namespace Dumux::Properties {
// create the type tag nodes
// Create new type tags
namespace TTag {
struct OnePCompressible { using InheritsFrom = std::tuple<OneP>; };
struct OnePCompressibleTpfa { using InheritsFrom = std::tuple<OnePCompressible, CCTpfaModel>; };
struct OnePCompressibleMpfa { using InheritsFrom = std::tuple<OnePCompressible, CCMpfaModel>; };
struct OnePCompressibleBox { using InheritsFrom = std::tuple<OnePCompressible, BoxModel>; };
} // end namespace TTag

// Set the grid type
template<class TypeTag>
struct Grid<TypeTag, TTag::OnePCompressible> { using type = Dune::YaspGrid<2>; };

// Set the problem type
template<class TypeTag>
struct Problem<TypeTag, TTag::OnePCompressible> { using type = OnePTestProblem<TypeTag>; };

// set the spatial params
template<class TypeTag>
struct SpatialParams<TypeTag, TTag::OnePCompressible>
{
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using type = OnePTestSpatialParams<GridGeometry, Scalar>;
};

// the fluid system
template<class TypeTag>
struct FluidSystem<TypeTag, TTag::OnePCompressible>
{
private:
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
public:
    using type = FluidSystems::OnePLiquid<Scalar, Components::TabulatedComponent<Components::H2O<Scalar>>>;
};

// for the experimental test, overwrite grid variables
// property in order to test grid variables-based assembly
#if EXPERIMENTAL
template<class TypeTag>
struct GridVariables<TypeTag, TTag::OnePCompressible>
{
private:
    using GVV = GetPropType<TypeTag, Properties::GridVolumeVariables>;
    using GFC = GetPropType<TypeTag, Properties::GridFluxVariablesCache>;
    using X = GetPropType<TypeTag, Properties::SolutionVector>;

public:
    using type = Dumux::Experimental::FVGridVariables<GVV, GFC, X>;
};
#endif

// Disable caching (for testing purposes)
template<class TypeTag>
struct EnableGridVolumeVariablesCache<TypeTag, TTag::OnePCompressible> { static constexpr bool value = false; };
template<class TypeTag>
struct EnableGridFluxVariablesCache<TypeTag, TTag::OnePCompressible> { static constexpr bool value = false; };
template<class TypeTag>
struct EnableGridGeometryCache<TypeTag, TTag::OnePCompressible> { static constexpr bool value = false; };

} // end namespace Dumux::Properties

#endif
