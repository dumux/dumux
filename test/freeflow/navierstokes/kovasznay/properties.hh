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
 * \brief The properties of the test for the staggered grid Navier-Stokes model with analytical solution (Kovasznay 1948, \cite Kovasznay1948)
 */
#ifndef DUMUX_KOVASZNAY_TEST_PROPERTIES_HH
#define DUMUX_KOVASZNAY_TEST_PROPERTIES_HH

#ifndef UPWINDSCHEMEORDER
#define UPWINDSCHEMEORDER 0
#endif

#include <dune/grid/yaspgrid.hh>

#if HAVE_DUNE_SUBGRID
#include <dune/subgrid/subgrid.hh>
#endif

#include <dumux/discretization/fcstaggered.hh>
#include <dumux/discretization/cctpfa.hh>

#include <dumux/material/fluidsystems/1pliquid.hh>
#include <dumux/material/components/constant.hh>

#include <dumux/freeflow/navierstokes/momentum/model.hh>
#include <dumux/freeflow/navierstokes/mass/1p/model.hh>

#include <dumux/multidomain/staggeredfreeflow/couplingmanager.hh>
#include "problem.hh"

namespace Dumux::Properties {

// Create new type tags
namespace TTag {
struct KovasznayTest {};
struct KovasznayTestMomentum { using InheritsFrom = std::tuple<KovasznayTest, NavierStokesMomentum, FaceCenteredStaggeredModel>; };
struct KovasznayTestMass { using InheritsFrom = std::tuple<KovasznayTest, NavierStokesMassOneP, CCTpfaModel>; };
} // end namespace TTag

// the fluid system
template<class TypeTag>
struct FluidSystem<TypeTag, TTag::KovasznayTest>
{
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using type = FluidSystems::OnePLiquid<Scalar, Components::Constant<1, Scalar> >;
};

// Set the grid type
template<class TypeTag>
struct Grid<TypeTag, TTag::KovasznayTest>
{
    using HostGrid = Dune::YaspGrid<2, Dune::EquidistantOffsetCoordinates<GetPropType<TypeTag, Properties::Scalar>, 2> >;

#if HAVE_DUNE_SUBGRID
    using type = Dune::SubGrid<HostGrid::dimension, HostGrid>;
#else
    using type = HostGrid;
#endif
};

// set the coupling manager property
template<class TypeTag>
struct CouplingManager<TypeTag, TTag::KovasznayTest>
{
private:
    using Traits = MultiDomainTraits<TTag::KovasznayTestMomentum, TTag::KovasznayTestMass>;
public:
    using type = StaggeredFreeFlowCouplingManager<Traits>;
};

// Set the problem property
template<class TypeTag>
struct Problem<TypeTag, TTag::KovasznayTest> { using type = Dumux::KovasznayTestProblem<TypeTag> ; };

template<class TypeTag>
struct EnableGridGeometryCache<TypeTag, TTag::KovasznayTest> { static constexpr bool value = true; };
template<class TypeTag>
struct EnableGridFluxVariablesCache<TypeTag, TTag::KovasznayTest> { static constexpr bool value = true; };
template<class TypeTag>
struct EnableGridVolumeVariablesCache<TypeTag, TTag::KovasznayTest> { static constexpr bool value = true; };

template<class TypeTag>
struct UpwindSchemeOrder<TypeTag, TTag::KovasznayTest> { static constexpr int value = UPWINDSCHEMEORDER; };

} // end namespace Dumux::Properties

#endif
