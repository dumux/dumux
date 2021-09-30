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
 * \ingroup NavierStokesTests
 * \brief The properties of the channel flow test for the staggered grid (Navier-)Stokes model.
 */
#ifndef DUMUX_CHANNEL_TEST_PROPERTIES_HH
#define DUMUX_CHANNEL_TEST_PROPERTIES_HH

#include <dune/grid/yaspgrid.hh>

#include <dumux/discretization/fcstaggered.hh>
#include <dumux/discretization/cctpfa.hh>

#include <dumux/freeflow/navierstokes/momentum/model.hh>
#include <dumux/freeflow/navierstokes/mass/1p/model.hh>

#include <dumux/material/fluidsystems/1pliquid.hh>
#include <dumux/material/components/constant.hh>
#include <dumux/material/components/simpleh2o.hh>

#include <dumux/multidomain/staggeredfreeflow/couplingmanager.hh>

#include "problem.hh"

namespace Dumux::Properties {

// Create new type tags
namespace TTag {
struct ChannelTest {};
struct ChannelTestMomentum { using InheritsFrom = std::tuple<ChannelTest, NavierStokesMomentum, FaceCenteredStaggeredModel>; };
#if NONISOTHERMAL
struct ChannelTestMass { using InheritsFrom = std::tuple<ChannelTest, NavierStokesMassOnePNI, CCTpfaModel>; };
#else
struct ChannelTestMass { using InheritsFrom = std::tuple<ChannelTest, NavierStokesMassOneP, CCTpfaModel>; };
#endif
} // end namespace TTag

// Set the problem property
template<class TypeTag>
struct Problem<TypeTag, TTag::ChannelTest>
{ using type = ChannelTestProblem<TypeTag>; };

// the fluid system
template<class TypeTag>
struct FluidSystem<TypeTag, TTag::ChannelTest>
{
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
#if NONISOTHERMAL
    using type = FluidSystems::OnePLiquid<Scalar, Components::SimpleH2O<Scalar> >;
#else
    using type = FluidSystems::OnePLiquid<Scalar, Components::Constant<1, Scalar> >;
#endif
};

// Set the grid type
template<class TypeTag>
struct Grid<TypeTag, TTag::ChannelTest>
{ using type = Dune::YaspGrid<2>; };

template<class TypeTag>
struct EnableGridGeometryCache<TypeTag, TTag::ChannelTest>
{ static constexpr bool value = true; };

template<class TypeTag>
struct EnableGridFluxVariablesCache<TypeTag, TTag::ChannelTest>
{ static constexpr bool value = true; };

template<class TypeTag>
struct EnableGridVolumeVariablesCache<TypeTag, TTag::ChannelTest>
{ static constexpr bool value = true; };

// Set the problem property
template<class TypeTag>
struct CouplingManager<TypeTag, TTag::ChannelTest>
{
    using Traits = MultiDomainTraits<TTag::ChannelTestMomentum, TTag::ChannelTestMass>;
    using type = StaggeredFreeFlowCouplingManager<Traits>;
};

} // end namespace Dumux::Properties

#endif
