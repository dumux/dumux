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
 * \brief The properties of the test for the staggered grid (Navier-)Stokes model with analytical solution (Donea 2003, \cite Donea2003).
 */
#ifndef DUMUX_DONEA_TEST_PROPERTIES_HH
#define DUMUX_DONEA_TEST_PROPERTIES_HH

#ifndef ENABLECACHING
#define ENABLECACHING 0
#endif

#ifndef NAVIER_STOKES_MOMENTUM_MODEL
#define NAVIER_STOKES_MOMENTUM_MODEL NavierStokesMomentum
#endif

#ifndef MOMENTUM_DISCRETIZATION_MODEL
#define MOMENTUM_DISCRETIZATION_MODEL FaceCenteredStaggeredModel
#endif

#ifndef MASS_DISCRETIZATION_MODEL
#define MASS_DISCRETIZATION_MODEL CCTpfaModel
#endif

#ifndef ALUGRID_CELL_TYPE
#define ALUGRID_CELL_TYPE cube
#endif

#if HAVE_DUNE_ALUGRID
#include <dune/alugrid/grid.hh>
#else
#include <dune/grid/yaspgrid.hh>
#endif

#include <dumux/discretization/fcstaggered.hh>
#include <dumux/discretization/cctpfa.hh>
#include <dumux/discretization/box.hh>

#include <dumux/freeflow/navierstokes/momentum/model.hh>
#include <dumux/freeflow/navierstokes/momentum/cvfe/model.hh>
#include <dumux/freeflow/navierstokes/mass/1p/model.hh>
#include <dumux/freeflow/navierstokes/momentum/problem.hh>
#include <dumux/freeflow/navierstokes/mass/problem.hh>
#include <dumux/multidomain/traits.hh>
#include <dumux/multidomain/freeflow/couplingmanager.hh>
#include <dumux/material/components/constant.hh>
#include <dumux/material/fluidsystems/1pliquid.hh>

#include "problem.hh"

namespace Dumux::Properties {

// Create new type tags
namespace TTag {
struct DoneaTest {};
struct DoneaTestMomentum { using InheritsFrom = std::tuple<DoneaTest, NAVIER_STOKES_MOMENTUM_MODEL, MOMENTUM_DISCRETIZATION_MODEL>; };
struct DoneaTestMass { using InheritsFrom = std::tuple<DoneaTest, NavierStokesMassOneP, MASS_DISCRETIZATION_MODEL>; };
} // end namespace TTag

template<class TypeTag>
struct Problem<TypeTag, TTag::DoneaTestMomentum>
{ using type = DoneaTestProblem<TypeTag, Dumux::NavierStokesMomentumProblem<TypeTag>>; };

template<class TypeTag>
struct Problem<TypeTag, TTag::DoneaTestMass>
{ using type = DoneaTestProblem<TypeTag, Dumux::NavierStokesMassProblem<TypeTag>>; };

template<class TypeTag>
struct FluidSystem<TypeTag, TTag::DoneaTest>
{
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using type = FluidSystems::OnePLiquid<Scalar, Components::Constant<1, Scalar> >;
};

template<class TypeTag>
struct Grid<TypeTag, TTag::DoneaTest>
#if HAVE_DUNE_ALUGRID
{ using type = Dune::ALUGrid<2, 2, Dune::ALUGRID_CELL_TYPE, Dune::nonconforming>; };
#else
{ using type = Dune::YaspGrid<2>; };
#endif

template<class TypeTag>
struct EnableGridGeometryCache<TypeTag, TTag::DoneaTest> { static constexpr bool value = ENABLECACHING; };
template<class TypeTag>
struct EnableGridFluxVariablesCache<TypeTag, TTag::DoneaTest> { static constexpr bool value = ENABLECACHING; };
template<class TypeTag>
struct EnableGridVolumeVariablesCache<TypeTag, TTag::DoneaTest> { static constexpr bool value = ENABLECACHING; };

// Set the problem property
template<class TypeTag>
struct CouplingManager<TypeTag, TTag::DoneaTest>
{
    using Traits = MultiDomainTraits<TTag::DoneaTestMomentum, TTag::DoneaTestMass>;
    using type = FreeFlowCouplingManager<Traits>;
};

} // end namespace Dumux::Properties

#endif
