// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup EmbeddedTests
 * \brief Properties for the parallel DGF-network 1d3d test (bulk reuses the tissue problem,
 *        the network uses an injected per-element radius migrated from a DGF file).
 */
#ifndef DUMUX_NETWORK_PARALLEL_PROPERTIES_HH
#define DUMUX_NETWORK_PARALLEL_PROPERTIES_HH

#ifndef BULKTYPETAG
#define BULKTYPETAG TissueCC
#endif
#ifndef LOWDIMTYPETAG
#define LOWDIMTYPETAG NetworkCC
#endif
#ifndef COUPLINGMODE
#define COUPLINGMODE Embedded1d3dCouplingMode::Average
#endif

#include <dune/grid/yaspgrid.hh>
#include <dune/foamgrid/foamgrid.hh>

#include <dumux/common/properties.hh>
#include <dumux/discretization/cctpfa.hh>

#include <dumux/porousmediumflow/1p/model.hh>
#include <dumux/material/components/constant.hh>
#include <dumux/material/fluidsystems/1pliquid.hh>
#include <dumux/porousmediumflow/1p/incompressiblelocalresidual.hh>

#include <dumux/multidomain/embedded/couplingmanager1d3d.hh>
#include <dumux/multidomain/traits.hh>

#include "problem_tissue.hh"
#include "spatialparams_tissue.hh"
#include "problem_network.hh"
#include "spatialparams_network.hh"

namespace Dumux::Properties {

// ---- bulk (tissue / soil): reuse the analytical tissue problem ----
namespace TTag {
struct Tissue { using InheritsFrom = std::tuple<OneP>; };
struct TissueCC { using InheritsFrom = std::tuple<Tissue, CCTpfaModel>; };
} // end namespace TTag

template<class TypeTag>
struct Grid<TypeTag, TTag::Tissue> { using type = Dune::YaspGrid<3, Dune::EquidistantOffsetCoordinates<GetPropType<TypeTag, Properties::Scalar>, 3> >; };

template<class TypeTag>
struct EnableGridGeometryCache<TypeTag, TTag::Tissue> { static constexpr bool value = true; };
template<class TypeTag>
struct EnableGridVolumeVariablesCache<TypeTag, TTag::Tissue> { static constexpr bool value = true; };
template<class TypeTag>
struct EnableGridFluxVariablesCache<TypeTag, TTag::Tissue> { static constexpr bool value = true; };

template<class TypeTag>
struct Problem<TypeTag, TTag::Tissue> { using type = TissueProblem<TypeTag>; };

template<class TypeTag>
struct LocalResidual<TypeTag, TTag::Tissue> { using type = OnePIncompressibleLocalResidual<TypeTag>; };

template<class TypeTag>
struct FluidSystem<TypeTag, TTag::Tissue>
{
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using type = FluidSystems::OnePLiquid<Scalar, Components::Constant<1, Scalar> >;
};

template<class TypeTag>
struct SpatialParams<TypeTag, TTag::Tissue>
{
    using type = TissueSpatialParams<GetPropType<TypeTag, Properties::GridGeometry>,
                                     GetPropType<TypeTag, Properties::Scalar>>;
};

// ---- low dim (network / root): injected per-element radius ----
namespace TTag {
struct Network { using InheritsFrom = std::tuple<OneP>; };
struct NetworkCC { using InheritsFrom = std::tuple<Network, CCTpfaModel>; };
} // end namespace TTag

template<class TypeTag>
struct Grid<TypeTag, TTag::Network> { using type = Dune::FoamGrid<1, 3>; };

template<class TypeTag>
struct EnableGridGeometryCache<TypeTag, TTag::Network> { static constexpr bool value = true; };
template<class TypeTag>
struct EnableGridVolumeVariablesCache<TypeTag, TTag::Network> { static constexpr bool value = true; };
template<class TypeTag>
struct EnableGridFluxVariablesCache<TypeTag, TTag::Network> { static constexpr bool value = true; };

template<class TypeTag>
struct Problem<TypeTag, TTag::Network> { using type = NetworkProblem<TypeTag>; };

template<class TypeTag>
struct LocalResidual<TypeTag, TTag::Network> { using type = OnePIncompressibleLocalResidual<TypeTag>; };

template<class TypeTag>
struct FluidSystem<TypeTag, TTag::Network>
{
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using type = FluidSystems::OnePLiquid<Scalar, Components::Constant<1, Scalar> >;
};

template<class TypeTag>
struct SpatialParams<TypeTag, TTag::Network>
{
    using type = NetworkSpatialParams<GetPropType<TypeTag, Properties::GridGeometry>,
                                      GetPropType<TypeTag, Properties::Scalar>>;
};

// ---- coupling manager wiring (same pattern as the analytical 1p_1p test) ----
template<class Traits>
using TheCouplingManager = Embedded1d3dCouplingManager<Traits, COUPLINGMODE>;

template<class TypeTag>
struct CouplingManager<TypeTag, TTag::BULKTYPETAG> { using type = TheCouplingManager<MultiDomainTraits<TypeTag, Properties::TTag::LOWDIMTYPETAG>>; };
template<class TypeTag>
struct PointSource<TypeTag, TTag::BULKTYPETAG> { using type = typename GetPropType<TypeTag, Properties::CouplingManager>::PointSourceTraits::template PointSource<0>; };
template<class TypeTag>
struct PointSourceHelper<TypeTag, TTag::BULKTYPETAG> { using type = typename GetPropType<TypeTag, Properties::CouplingManager>::PointSourceTraits::template PointSourceHelper<0>; };

template<class TypeTag>
struct CouplingManager<TypeTag, TTag::LOWDIMTYPETAG> { using type = TheCouplingManager<MultiDomainTraits<Properties::TTag::BULKTYPETAG, TypeTag>>; };
template<class TypeTag>
struct PointSource<TypeTag, TTag::LOWDIMTYPETAG> { using type = typename GetPropType<TypeTag, Properties::CouplingManager>::PointSourceTraits::template PointSource<1>; };
template<class TypeTag>
struct PointSourceHelper<TypeTag, TTag::LOWDIMTYPETAG> { using type = typename GetPropType<TypeTag, Properties::CouplingManager>::PointSourceTraits::template PointSourceHelper<1>; };

} // end namespace Dumux::Properties

#endif
