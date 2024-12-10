// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup OnePTests
 * \ingroup MultiDomainTests
 * \brief The properties for the single-phase darcy-darcy mortar-coupling test
 */
#ifndef DUMUX_MORTAR_STOKES_DARCY_ONEP_TEST_PROPERTIES_HH
#define DUMUX_MORTAR_STOKES_DARCY_ONEP_TEST_PROPERTIES_HH

#include <dune/common/fvector.hh>
#include <dune/istl/bvector.hh>

#include <dune/grid/yaspgrid.hh>
#include <dune/foamgrid/foamgrid.hh>
#include <dumux/discretization/cctpfa.hh>
#include <dumux/discretization/box.hh>
#include <dumux/discretization/fcstaggered.hh>

#include <dumux/porousmediumflow/1p/model.hh>
#include <dumux/freeflow/navierstokes/momentum/model.hh>
#include <dumux/freeflow/navierstokes/mass/1p/model.hh>
#include <dumux/freeflow/navierstokes/momentum/problem.hh>
#include <dumux/freeflow/navierstokes/mass/problem.hh>

#include <dumux/material/components/constant.hh>
#include <dumux/material/fluidsystems/1pliquid.hh>

#include <dumux/multidomain/traits.hh>
#include <dumux/multidomain/freeflow/couplingmanager.hh>
#include <dumux/multidomain/mortar/properties.hh>

#include "problem_darcy.hh"
#include "problem_stokes.hh"
#include "spatialparams.hh"

namespace Dumux::Properties {

namespace TTag {
struct OnePDarcyMortar { using InheritsFrom = std::tuple<OneP>; };
struct OnePDarcyMortarTpfa { using InheritsFrom = std::tuple<OnePDarcyMortar, CCTpfaModel>; };
struct OnePStokesMortar {};
struct OnePStokesMortarMomentum { using InheritsFrom = std::tuple<OnePStokesMortar, NavierStokesMomentum, FaceCenteredStaggeredModel>; };
struct OnePStokesMortarMass { using InheritsFrom = std::tuple<OnePStokesMortar, NavierStokesMassOneP, CCTpfaModel>; };
} // namespace TTag

template<class TypeTag> struct MortarGrid<TypeTag, TTag::OnePDarcyMortar> { using type = Dune::FoamGrid<1, 2>; };
template<class TypeTag> struct MortarGrid<TypeTag, TTag::OnePStokesMortar> { using type = Dune::FoamGrid<1, 2>; };

template<class TypeTag> struct MortarSolutionVector<TypeTag, TTag::OnePDarcyMortar> { using type = Dune::BlockVector<Dune::FieldVector<double, 1>>; };
template<class TypeTag> struct MortarSolutionVector<TypeTag, TTag::OnePStokesMortar> { using type = Dune::BlockVector<Dune::FieldVector<double, 1>>; };

template<class TypeTag> struct Grid<TypeTag, TTag::OnePDarcyMortar> { using type = Dune::YaspGrid<2, Dune::EquidistantOffsetCoordinates<double, 2>>; };
template<class TypeTag> struct Grid<TypeTag, TTag::OnePStokesMortar> { using type = Dune::YaspGrid<2, Dune::EquidistantOffsetCoordinates<double, 2>>; };

template<class TypeTag> struct Problem<TypeTag, TTag::OnePDarcyMortar> { using type = DarcyProblem<TypeTag>; };
template<class TypeTag> struct Problem<TypeTag, TTag::OnePStokesMortarMomentum> { using type = StokesProblem<TypeTag, Dumux::NavierStokesMomentumProblem<TypeTag>>; };
template<class TypeTag> struct Problem<TypeTag, TTag::OnePStokesMortarMass> { using type = StokesProblem<TypeTag, Dumux::NavierStokesMassProblem<TypeTag>>; };

template<class TypeTag> struct EnableGridGeometryCache<TypeTag, TTag::OnePStokesMortar> { static constexpr bool value = true; };
template<class TypeTag> struct EnableGridFluxVariablesCache<TypeTag, TTag::OnePStokesMortar> { static constexpr bool value = true; };
template<class TypeTag> struct EnableGridVolumeVariablesCache<TypeTag, TTag::OnePStokesMortar> { static constexpr bool value = true; };

template<class TypeTag>
struct SpatialParams<TypeTag, TTag::OnePDarcyMortar>
{
    using type = DarcySpatialParams<
        GetPropType<TypeTag, Properties::GridGeometry>,
        GetPropType<TypeTag, Properties::Scalar>
    >;
};

template<class TypeTag>
struct FluidSystem<TypeTag, TTag::OnePDarcyMortar>
{
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using type = FluidSystems::OnePLiquid<Scalar, Components::Constant<0, Scalar> >;
};
template<class TypeTag>
struct FluidSystem<TypeTag, TTag::OnePStokesMortar>
{
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using type = FluidSystems::OnePLiquid<Scalar, Components::Constant<0, Scalar> >;
};

template<class TypeTag>
struct CouplingManager<TypeTag, TTag::OnePStokesMortar>
{
private:
    using Traits = MultiDomainTraits<TTag::OnePStokesMortarMomentum, TTag::OnePStokesMortarMass>;
public:
    using type = FreeFlowCouplingManager<Traits>;
};

}  // namespace Dumux::Properties

#endif
