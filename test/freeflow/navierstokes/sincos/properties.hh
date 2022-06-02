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
 * \ingroup NavierStokesNCTests
 * \brief The properties of the test for the staggered grid Navier-Stokes model with analytical solution.
 */
#ifndef DUMUX_SINCOS_TEST_PROPERTIES_HH
#define DUMUX_SINCOS_TEST_PROPERTIES_HH

#ifndef ENABLECACHING
#define ENABLECACHING true
#endif

#ifndef GRIDTYPE
#if HAVE_DUNE_ALUGRID
#define GRIDTYPE Dune::ALUGrid<2,2,Dune::cube,Dune::nonconforming>
#else
#define GRIDTYPE Dune::YaspGrid<2>
#endif
#endif

#ifndef NAVIER_STOKES_MODEL
#define NAVIER_STOKES_MODEL NavierStokesMomentum
#endif

#ifndef MOMENTUM_DISCRETIZATION_MODEL
#define MOMENTUM_DISCRETIZATION_MODEL FaceCenteredStaggeredModel
#endif

#ifndef MASS_DISCRETIZATION_MODEL
#define MASS_DISCRETIZATION_MODEL CCTpfaModel
#endif

#if HAVE_DUNE_ALUGRID
#include <dune/alugrid/grid.hh>
#else
#include <dune/grid/yaspgrid.hh>
#endif

#include <dumux/flux/fluxvariablescaching.hh>

#include <dumux/discretization/fcstaggered.hh>
#include <dumux/discretization/fcdiamond.hh>
#include <dumux/discretization/cctpfa.hh>
#include <dumux/discretization/box.hh>

#include <dumux/freeflow/navierstokes/momentum/model.hh>
#include <dumux/freeflow/navierstokes/momentum/diamond/model.hh>
#include <dumux/freeflow/navierstokes/mass/1p/model.hh>

#include <dumux/material/components/constant.hh>
#include <dumux/material/fluidsystems/1pliquid.hh>

#include <dumux/multidomain/staggeredfreeflow/couplingmanager.hh>

#include "problem.hh"

namespace Dumux::Properties {

// Create new type tags
namespace TTag {
struct SincosTest {};
struct SincosTestMomentum { using InheritsFrom = std::tuple<SincosTest, NAVIER_STOKES_MODEL, MOMENTUM_DISCRETIZATION_MODEL>; };
struct SincosTestMass { using InheritsFrom = std::tuple<SincosTest, NavierStokesMassOneP, MASS_DISCRETIZATION_MODEL>; };
} // end namespace TTag

template<class TypeTag>
struct Problem<TypeTag, TTag::SincosTest>
{ using type = SincosTestProblem<TypeTag>; };

template<class TypeTag>
struct FluidSystem<TypeTag, TTag::SincosTest>
{
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using type = FluidSystems::OnePLiquid<Scalar, Components::Constant<1, Scalar> >;
};

template<class TypeTag>
struct Grid<TypeTag, TTag::SincosTest>
{ using type = GRIDTYPE; };

template<class TypeTag>
struct EnableGridGeometryCache<TypeTag, TTag::SincosTest> { static constexpr bool value = ENABLECACHING; };
template<class TypeTag>
struct EnableGridFluxVariablesCache<TypeTag, TTag::SincosTest> { static constexpr bool value = ENABLECACHING; };
template<class TypeTag>
struct EnableGridVolumeVariablesCache<TypeTag, TTag::SincosTest> { static constexpr bool value = ENABLECACHING; };

template<class TypeTag>
struct FluxVariablesCache<TypeTag, TTag::SincosTestMomentum>
{
private:
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
public:
    using type = std::conditional_t<
        GridGeometry::discMethod == DiscretizationMethods::fcdiamond,
        FaceCenteredDiamondFluxVariablesCache<Scalar, GridGeometry>,
        FluxVariablesCaching::EmptyCache<Scalar>
    >;
};

// Set the problem property
template<class TypeTag>
struct CouplingManager<TypeTag, TTag::SincosTest>
{
    using Traits = MultiDomainTraits<TTag::SincosTestMomentum, TTag::SincosTestMass>;
    using type = StaggeredFreeFlowCouplingManager<Traits>;
};

} // end namespace Dumux::Properties

#endif
