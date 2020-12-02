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
 * \ingroup SequentialTwoPTests
 * \brief test problem for the explicit transport model
 */
#ifndef DUMUX_TEST_TRANSPORT_PROBLEM_PROPERTIES_HH
#define DUMUX_TEST_TRANSPORT_PROBLEM_PROPERTIES_HH

#include <dune/grid/yaspgrid.hh>
#include <dune/grid/io/file/dgfparser/dgfyasp.hh>

#include <dumux/common/properties.hh>
#include <dumux/material/fluidsystems/1pliquid.hh>
#include <dumux/material/components/constant.hh>

#include <dumux/porousmediumflow/2p/sequential/transport/cellcentered/properties.hh>


#include "test_transportproblem.hh"
#include "test_transportspatialparams.hh"

namespace Dumux::Properties
{
// Create new type tags
namespace TTag {
struct TransportTest { using InheritsFrom = std::tuple<TestTransportSpatialParams, FVTransportTwoP>; };
} // end namespace TTag

// Set the grid type
template<class TypeTag>
struct Grid<TypeTag, TTag::TransportTest> { using type = Dune::YaspGrid<2>; };

// Set the problem property
template<class TypeTag>
struct Problem<TypeTag, TTag::TransportTest> { using type = TestTransportProblem<TypeTag>; };

// Set the fluid system
template<class TypeTag>
struct FluidSystem<TypeTag, TTag::TransportTest>
{
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using WettingPhase = FluidSystems::OnePLiquid<Scalar, Components::Constant<1, Scalar> >;
    using NonwettingPhase = FluidSystems::OnePLiquid<Scalar, Components::Constant<1, Scalar> >;
    using type = FluidSystems::TwoPImmiscible<Scalar, WettingPhase, NonwettingPhase>;
};

template<class TypeTag>
struct VelocityFormulation<TypeTag, TTag::TransportTest> { static constexpr int value = SequentialTwoPCommonIndices::velocityTotal; };
} // end namespace Properties

#endif
