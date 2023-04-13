// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup FacetTests
 * \brief The properties for the single-phase facet coupling test.
 */

#ifndef DUMUX_TEST_FACETCOUPLING_THREEDOMAIN_PROPERTIES_HH
#define DUMUX_TEST_FACETCOUPLING_THREEDOMAIN_PROPERTIES_HH

#include <dune/alugrid/grid.hh>
#include <dune/foamgrid/foamgrid.hh>

#include <dumux/material/components/constant.hh>
#include <dumux/material/fluidsystems/1pliquid.hh>

#include <dumux/multidomain/facet/box/properties.hh>
#include <dumux/multidomain/facet/cellcentered/tpfa/properties.hh>
#include <dumux/multidomain/facet/cellcentered/mpfa/properties.hh>
#include <dumux/multidomain/facet/couplingmapper.hh>
#include <dumux/multidomain/facet/couplingmanager.hh>
#include <dumux/multidomain/traits.hh>
#include <dumux/porousmediumflow/1p/model.hh>

#include <dumux/discretization/box.hh>
#include <dumux/discretization/cctpfa.hh>
#include <dumux/discretization/ccmpfa.hh>

#include "spatialparams.hh"
#include "problem_bulk.hh"
#include "problem_facet.hh"
#include "problem_edge.hh"

namespace Dumux::Properties {

// create the type tag nodes
// Create new type tags
namespace TTag {
struct OnePBulk { using InheritsFrom = std::tuple<OneP>; };
struct OnePBulkTpfa { using InheritsFrom = std::tuple<CCTpfaFacetCouplingModel, OnePBulk>; };
struct OnePBulkMpfa { using InheritsFrom = std::tuple<CCMpfaFacetCouplingModel, OnePBulk>; };
struct OnePBulkBox { using InheritsFrom = std::tuple<BoxFacetCouplingModel, OnePBulk>; };
} // end namespace TTag

// Set the grid type
template<class TypeTag>
struct Grid<TypeTag, TTag::OnePBulk> { using type = Dune::ALUGrid<3, 3, Dune::simplex, Dune::nonconforming>; };
// Set the problem type
template<class TypeTag>
struct Problem<TypeTag, TTag::OnePBulk> { using type = OnePBulkProblem<TypeTag>; };
// set the spatial params
template<class TypeTag>
struct SpatialParams<TypeTag, TTag::OnePBulk>
{
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using type = OnePSpatialParams<GridGeometry, Scalar>;
};

// the fluid system
template<class TypeTag>
struct FluidSystem<TypeTag, TTag::OnePBulk>
{
private:
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
public:
    using type = FluidSystems::OnePLiquid< Scalar, Components::Constant<1, Scalar> >;
};

// Create new type tags
namespace TTag {
struct OnePFacet { using InheritsFrom = std::tuple<OneP>; };
struct OnePFacetTpfa { using InheritsFrom = std::tuple<CCTpfaFacetCouplingModel, OnePFacet>; };
struct OnePFacetMpfa { using InheritsFrom = std::tuple<CCMpfaFacetCouplingModel, OnePFacet>; };
struct OnePFacetBox { using InheritsFrom = std::tuple<BoxFacetCouplingModel, OnePFacet>; };
} // end namespace TTag

// Set the grid type
template<class TypeTag>
struct Grid<TypeTag, TTag::OnePFacet> { using type = Dune::FoamGrid<2, 3>; };
// Set the problem type
template<class TypeTag>
struct Problem<TypeTag, TTag::OnePFacet> { using type = OnePFacetProblem<TypeTag>; };
// set the spatial params
template<class TypeTag>
struct SpatialParams<TypeTag, TTag::OnePFacet>
{
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using type = OnePSpatialParams<GridGeometry, Scalar>;
};

// the fluid system
template<class TypeTag>
struct FluidSystem<TypeTag, TTag::OnePFacet>
{
private:
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
public:
    using type = FluidSystems::OnePLiquid< Scalar, Components::Constant<1, Scalar> >;
};

// Create new type tags
namespace TTag {
struct OnePEdge { using InheritsFrom = std::tuple<OneP>; };
struct OnePEdgeTpfa { using InheritsFrom = std::tuple<OnePEdge, CCTpfaModel>; };
struct OnePEdgeMpfa { using InheritsFrom = std::tuple<OnePEdge, CCTpfaModel>; };
struct OnePEdgeBox { using InheritsFrom = std::tuple<OnePEdge, BoxModel>; };
} // end namespace TTag

// Set the grid type
template<class TypeTag>
struct Grid<TypeTag, TTag::OnePEdge> { using type = Dune::FoamGrid<1, 3>; };
// Set the problem type
template<class TypeTag>
struct Problem<TypeTag, TTag::OnePEdge> { using type = OnePEdgeProblem<TypeTag>; };
// set the spatial params
template<class TypeTag>
struct SpatialParams<TypeTag, TTag::OnePEdge>
{
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using type = OnePSpatialParams<GridGeometry, Scalar>;
};

// the fluid system
template<class TypeTag>
struct FluidSystem<TypeTag, TTag::OnePEdge>
{
private:
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
public:
    using type = FluidSystems::OnePLiquid< Scalar, Components::Constant<1, Scalar> >;
};

template<class BulkTypeTag, class FacetTypeTag, class EdgeTypeTag>
class TestTraits
{
    using BulkFVG = Dumux::GetPropType<BulkTypeTag, Dumux::Properties::GridGeometry>;
    using FacetFVG = Dumux::GetPropType<FacetTypeTag, Dumux::Properties::GridGeometry>;
    using EdgeFVG = Dumux::GetPropType<EdgeTypeTag, Dumux::Properties::GridGeometry>;
public:
    using MDTraits = Dumux::MultiDomainTraits<BulkTypeTag, FacetTypeTag, EdgeTypeTag>;
    using CouplingMapper = Dumux::FacetCouplingThreeDomainMapper<BulkFVG, FacetFVG, EdgeFVG>;
    using CouplingManager = Dumux::FacetCouplingThreeDomainManager<MDTraits, CouplingMapper>;
};

using TpfaTraits = TestTraits<TTag::OnePBulkTpfa, TTag::OnePFacetTpfa, TTag::OnePEdgeTpfa>;
using MpfaTraits = TestTraits<TTag::OnePBulkMpfa, TTag::OnePFacetMpfa, TTag::OnePEdgeMpfa>;
using BoxTraits = TestTraits<TTag::OnePBulkBox, TTag::OnePFacetBox, TTag::OnePEdgeBox>;

template<class TypeTag>
struct CouplingManager<TypeTag, TTag::OnePBulkTpfa> { using type = typename TpfaTraits::CouplingManager; };
template<class TypeTag>
struct CouplingManager<TypeTag, TTag::OnePFacetTpfa> { using type = typename TpfaTraits::CouplingManager; };
template<class TypeTag>
struct CouplingManager<TypeTag, TTag::OnePEdgeTpfa> { using type = typename TpfaTraits::CouplingManager; };

template<class TypeTag>
struct CouplingManager<TypeTag, TTag::OnePBulkMpfa> { using type = typename MpfaTraits::CouplingManager; };
template<class TypeTag>
struct CouplingManager<TypeTag, TTag::OnePFacetMpfa> { using type = typename MpfaTraits::CouplingManager; };
template<class TypeTag>
struct CouplingManager<TypeTag, TTag::OnePEdgeMpfa> { using type = typename MpfaTraits::CouplingManager; };

template<class TypeTag>
struct CouplingManager<TypeTag, TTag::OnePBulkBox> { using type = typename BoxTraits::CouplingManager; };
template<class TypeTag>
struct CouplingManager<TypeTag, TTag::OnePFacetBox> { using type = typename BoxTraits::CouplingManager; };
template<class TypeTag>
struct CouplingManager<TypeTag, TTag::OnePEdgeBox> { using type = typename BoxTraits::CouplingManager; };

} // end namespace Dumux::Properties

#endif
