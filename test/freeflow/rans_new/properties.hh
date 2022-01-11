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
#include <dumux/freeflow/navierstokes/turbulence/twoeq/model.hh>

#include <dumux/material/fluidsystems/1pliquid.hh>
#include <dumux/material/components/constant.hh>

#include <dumux/multidomain/staggeredfreeflow/couplingmanager.hh>

#include "problem.hh"

namespace Dumux::Properties {

// Create new type tags
namespace TTag {
struct RANSTest {};
struct RANSTestMomentum { using InheritsFrom = std::tuple<RANSTest, NavierStokesMomentum, FaceCenteredStaggeredModel>; };
struct RANSTestMass { using InheritsFrom = std::tuple<RANSTest, RANSTwoEq, CCTpfaModel>; };
} // end namespace TTag

// Set the problem property
template<class TypeTag>
struct Problem<TypeTag, TTag::RANSTest>
{ using type = RansTestProblem<TypeTag>; };

// Set the grid type
template<class TypeTag>
struct Grid<TypeTag, TTag::RANSTest>
{ using type = Dune::YaspGrid<2, Dune::TensorProductCoordinates<GetPropType<TypeTag, Properties::Scalar>, 2> >; };

// the fluid system
template<class TypeTag>
struct FluidSystem<TypeTag, TTag::RANSTest>
{
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using type = FluidSystems::OnePLiquid<Scalar, Components::Constant<1, Scalar> >;
};

// Set the problem property
template<class TypeTag>
struct CouplingManager<TypeTag, TTag::RANSTest>
{
    using Traits = MultiDomainTraits<TTag::RANSTestMomentum, TTag::RANSTestMass>;
    using type = StaggeredFreeFlowCouplingManager<Traits>;
};

template<class TypeTag>
struct EnableGridGeometryCache<TypeTag, TTag::RANSTest> { static constexpr bool value = true; };
template<class TypeTag>
struct EnableGridFluxVariablesCache<TypeTag, TTag::RANSTest> { static constexpr bool value = true; };
template<class TypeTag>
struct EnableGridVolumeVariablesCache<TypeTag, TTag::RANSTest> { static constexpr bool value = true; };

} // end namespace Dumux::Properties

#endif
