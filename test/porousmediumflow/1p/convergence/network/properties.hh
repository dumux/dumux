// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup OnePTests
 * \brief The properties of the convergence test on lines and junctions
 */
#ifndef DUMUX_ONEP_NETWORK_CONVERGENCE_TEST_PROPERTIES_HH
#define DUMUX_ONEP_NETWORK_CONVERGENCE_TEST_PROPERTIES_HH

#include <dune/grid/onedgrid.hh>
#if HAVE_DUNE_FOAMGRID
#include <dune/foamgrid/foamgrid.hh>
#endif

#include <dumux/discretization/pq2.hh>
#include <dumux/discretization/pq3.hh>
#include <dumux/flux/cvfe/darcyslaw_.hh>
#include <dumux/porousmediumflow/1p/variables.hh>
#include <dumux/porousmediumflow/immiscible/localresidual_.hh>
#include <dumux/porousmediumflow/fvspatialparams1p.hh>

#include <dumux/material/components/constant.hh>
#include <dumux/material/fluidsystems/1pliquid.hh>

#include <dumux/porousmediumflow/1p/model.hh>

#include "problem.hh"

#ifndef GRIDTYPE
#define GRIDTYPE Dune::OneDGrid
#endif

namespace Dumux {

//! Unit permeability
template<class GridGeometry, class Scalar>
class OnePNetworkConvergenceSpatialParams
: public FVPorousMediumFlowSpatialParamsOneP<GridGeometry, Scalar,
                                             OnePNetworkConvergenceSpatialParams<GridGeometry, Scalar>>
{
    using ThisType = OnePNetworkConvergenceSpatialParams<GridGeometry, Scalar>;
    using ParentType = FVPorousMediumFlowSpatialParamsOneP<GridGeometry, Scalar, ThisType>;
    using GlobalPosition = typename GridGeometry::GridView::template Codim<0>::Geometry::GlobalCoordinate;
public:
    using PermeabilityType = Scalar;
    using ParentType::ParentType;

    PermeabilityType permeabilityAtPos(const GlobalPosition& globalPos) const
    { return 1.0; }

    Scalar porosityAtPos(const GlobalPosition& globalPos) const
    { return 0.4; }
};

} // end namespace Dumux

namespace Dumux::Properties {

namespace TTag {
struct OnePNetworkConvergence { using InheritsFrom = std::tuple<OneP>; };
struct OnePNetworkConvergencePQ2 { using InheritsFrom = std::tuple<OnePNetworkConvergence, PQ2HybridModel>; };
struct OnePNetworkConvergencePQ3 { using InheritsFrom = std::tuple<OnePNetworkConvergence, PQ3HybridModel>; };
} // end namespace TTag

template<class TypeTag>
struct Grid<TypeTag, TTag::OnePNetworkConvergence> { using type = GRIDTYPE; };

template<class TypeTag>
struct Problem<TypeTag, TTag::OnePNetworkConvergence>
{ using type = OnePNetworkConvergenceProblem<TypeTag>; };

template<class TypeTag>
struct FluidSystem<TypeTag, TTag::OnePNetworkConvergence>
{
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using type = FluidSystems::OnePLiquid<Scalar, Dumux::Components::Constant<1, Scalar>>;
};

template<class TypeTag>
struct SpatialParams<TypeTag, TTag::OnePNetworkConvergence>
{
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using type = OnePNetworkConvergenceSpatialParams<GridGeometry, Scalar>;
};

template<class TypeTag>
struct EnableGridVolumeVariablesCache<TypeTag, TTag::OnePNetworkConvergence> { static constexpr bool value = true; };
template<class TypeTag>
struct EnableGridFluxVariablesCache<TypeTag, TTag::OnePNetworkConvergence> { static constexpr bool value = true; };
template<class TypeTag>
struct EnableGridGeometryCache<TypeTag, TTag::OnePNetworkConvergence> { static constexpr bool value = true; };

template<class TypeTag>
struct GridVariables<TypeTag, TTag::OnePNetworkConvergence>
{
private:
    using GG = GetPropType<TypeTag, Properties::GridGeometry>;
    static constexpr bool enableCache = getPropValue<TypeTag, Properties::EnableGridVolumeVariablesCache>();
    using Problem = GetPropType<TypeTag, Properties::Problem>;

    using VariablesTraits = OnePVolumeVariablesTraits<
        GetPropType<TypeTag, Properties::PrimaryVariables>,
        GetPropType<TypeTag, Properties::FluidSystem>,
        GetPropType<TypeTag, Properties::FluidState>,
        GetPropType<TypeTag, Properties::SolidSystem>,
        GetPropType<TypeTag, Properties::SolidState>,
        typename GetPropType<TypeTag, Properties::SpatialParams>::PermeabilityType,
        GetPropType<TypeTag, Properties::ModelTraits>
    >;
    using Variables = Dumux::Experimental::OnePVariables<VariablesTraits>;

    using IPDataCache = Dumux::CVFE::LocalBasisInterpolationPointData<GG>;
    using Traits = Dumux::Experimental::CVFE::HybridCVFEDefaultGridVariablesCacheTraits<Problem, Variables, IPDataCache>;
    using GVC = Dumux::Experimental::CVFE::HybridCVFEGridVariablesCache<Traits, enableCache>;

public:
    using type = Dumux::Experimental::GridVariables<GG, GVC>;
};

template<class TypeTag>
struct AdvectionType<TypeTag, TTag::OnePNetworkConvergence>
{
    using type = Experimental::CVFEDarcysLawAtIp<GetPropType<TypeTag, Properties::Scalar>,
                                                 GetPropType<TypeTag, Properties::GridGeometry>>;
};

template<class TypeTag>
struct LocalResidual<TypeTag, TTag::OnePNetworkConvergence>
{ using type = Experimental::ImmiscibleLocalResidual<TypeTag>; };

} // end namespace Dumux::Properties

#endif
