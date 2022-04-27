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
 * \brief Stationary test for the staggered grid Navier-Stokes model with periodic BC
 */

#ifndef DUMUX_TEST_FREEFLOW_NAVIERSTOKES_PERIODIC_PROPERTIES_HH
#define DUMUX_TEST_FREEFLOW_NAVIERSTOKES_PERIODIC_PROPERTIES_HH


#include <dune/grid/spgrid.hh>

#include <dumux/material/fluidsystems/1pliquid.hh>
#include <dumux/material/components/constant.hh>

#include <dumux/freeflow/navierstokes/momentum/model.hh>
#include <dumux/freeflow/navierstokes/mass/1p/model.hh>

#include <dumux/multidomain/staggeredfreeflow/couplingmanager.hh>
#include <dumux/multidomain/traits.hh>

#include <dumux/discretization/fcstaggered.hh>
#include <dumux/discretization/cctpfa.hh>

#include "problem.hh"

namespace Dumux::Properties {

// Create new type tags
namespace TTag {
struct PeriodicTest {};
struct PeriodicTestMomentum { using InheritsFrom = std::tuple<PeriodicTest, NavierStokesMomentum, FaceCenteredStaggeredModel>; };
struct PeriodicTestMass { using InheritsFrom = std::tuple<PeriodicTest, NavierStokesMassOneP, CCTpfaModel>; };
} // end namespace TTag


// the fluid system
template<class TypeTag>
struct FluidSystem<TypeTag, TTag::PeriodicTest>
{
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using type = FluidSystems::OnePLiquid<Scalar, Components::Constant<1, Scalar> >;
};

// Set the grid type
template<class TypeTag>
struct Grid<TypeTag, TTag::PeriodicTest> { using type = Dune::SPGrid<double, 2>; };

// Set the problem property
template<class TypeTag>
struct Problem<TypeTag, TTag::PeriodicTest> { using type = PeriodicTestProblem<TypeTag> ; };

template<class TypeTag>
struct EnableGridGeometryCache<TypeTag, TTag::PeriodicTest> { static constexpr bool value = true; };
template<class TypeTag>
struct EnableGridFluxVariablesCache<TypeTag, TTag::PeriodicTest> { static constexpr bool value = true; };
template<class TypeTag>
struct EnableGridVolumeVariablesCache<TypeTag, TTag::PeriodicTest> { static constexpr bool value = true; };

template<class TypeTag>
struct CouplingManager<TypeTag, TTag::PeriodicTest>
{
    using Traits = MultiDomainTraits<TTag::PeriodicTestMomentum, TTag::PeriodicTestMass>;
    using type = StaggeredFreeFlowCouplingManager<Traits>;
};

} // end namespace Dumux::Properties

#endif
