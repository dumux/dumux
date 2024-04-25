// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup FacetTests
 * \brief The properties for the box-facet coupling dof enrichment test.
 */
#ifndef DUMUX_TEST_FACET_BOX_DOF_ENRICHMENT_ONEP_PROPERTIES_HH
#define DUMUX_TEST_FACET_BOX_DOF_ENRICHMENT_ONEP_PROPERTIES_HH

#include <dune/alugrid/grid.hh>
#include <dune/foamgrid/foamgrid.hh>

#include <dumux/material/components/constant.hh>
#include <dumux/material/fluidsystems/1pliquid.hh>
#include <dumux/porousmediumflow/1p/model.hh>
#include <dumux/discretization/box.hh>
#include <dumux/discretization/cctpfa.hh>
#include <dumux/discretization/ccmpfa.hh>

#include <dumux/multidomain/facet/cellcentered/tpfa/properties.hh>
#include <dumux/multidomain/facet/cellcentered/mpfa/properties.hh>
#include <dumux/multidomain/facet/box/properties.hh>
#include <dumux/multidomain/facet/couplingmapper.hh>
#include <dumux/multidomain/facet/couplingmanager.hh>
#include <dumux/multidomain/traits.hh>

#include "spatialparams.hh"
#include "problem_bulk.hh"
#include "problem_lowdim.hh"

// default for the bulk grid type
#ifndef BULKGRIDTYPE
#define BULKGRIDTYPE Dune::ALUGrid<2, 2, Dune::simplex, Dune::conforming>
#endif

#ifndef LOWDIMGRIDTYPE
#define LOWDIMGRIDTYPE Dune::FoamGrid<1, 2>
#endif

#ifndef USEBULKGGCACHING
#define USEBULKGGCACHING 0
#endif

namespace Dumux::Properties {

// create the type tag nodes
namespace TTag {
struct OnePBulk { using InheritsFrom = std::tuple<OneP>; };
struct OnePBulkBox { using InheritsFrom = std::tuple<BoxFacetCouplingModel, OnePBulk>; };
struct OnePBulkTpfa { using InheritsFrom = std::tuple<CCTpfaFacetCouplingModel, OnePBulk>; };
struct OnePBulkMpfa { using InheritsFrom = std::tuple<CCMpfaFacetCouplingModel, OnePBulk>; };
} // end namespace TTag

// Set the grid type
template<class TypeTag>
struct Grid<TypeTag, TTag::OnePBulk> { using type = BULKGRIDTYPE; };
// Set the problem type
template<class TypeTag>
struct Problem<TypeTag, TTag::OnePBulk> { using type = OnePBulkProblem<TypeTag>; };
// set the spatial params
template<class TypeTag>
struct SpatialParams<TypeTag, TTag::OnePBulk>
{
    using type = OnePSpatialParams< GetPropType<TypeTag, Properties::GridGeometry>,
                                    GetPropType<TypeTag, Properties::Scalar> >;
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

// grid geometry caching
template<class TypeTag>
struct EnableGridGeometryCache<TypeTag, TTag::OnePBulk>
{
    static constexpr bool value = USEBULKGGCACHING;
};

// create the type tag nodes
namespace TTag {
struct OnePLowDim { using InheritsFrom = std::tuple<OneP>; };
struct OnePLowDimBox { using InheritsFrom = std::tuple<OnePLowDim, BoxModel>; };
struct OnePLowDimTpfa { using InheritsFrom = std::tuple<OnePLowDim, CCTpfaModel>; };
struct OnePLowDimMpfa { using InheritsFrom = std::tuple<OnePLowDim, CCTpfaModel>; };
} // end namespace TTag

// Set the grid type
template<class TypeTag>
struct Grid<TypeTag, TTag::OnePLowDim> { using type = LOWDIMGRIDTYPE; };
// Set the problem type
template<class TypeTag>
struct Problem<TypeTag, TTag::OnePLowDim> { using type = OnePLowDimProblem<TypeTag>; };
// set the spatial params
template<class TypeTag>
struct SpatialParams<TypeTag, TTag::OnePLowDim>
{
    using type = OnePSpatialParams< GetPropType<TypeTag, Properties::GridGeometry>,
                                    GetPropType<TypeTag, Properties::Scalar> >;
};

// the fluid system
template<class TypeTag>
struct FluidSystem<TypeTag, TTag::OnePLowDim>
{
private:
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
public:
    using type = FluidSystems::OnePLiquid< Scalar, Components::Constant<1, Scalar> >;
};

// obtain/define some types to be used below in the property definitions and in main
template< class BulkTypeTag, class LowDimTypeTag >
class TestTraits
{
    using BulkFVGridGeometry = GetPropType<BulkTypeTag, Properties::GridGeometry>;
    using LowDimFVGridGeometry = GetPropType<LowDimTypeTag, Properties::GridGeometry>;
public:
    using MDTraits = Dumux::MultiDomainTraits<BulkTypeTag, LowDimTypeTag>;
    using CouplingMapper = Dumux::FacetCouplingMapper<BulkFVGridGeometry, LowDimFVGridGeometry>;
    using CouplingManager = Dumux::FacetCouplingManager<MDTraits, CouplingMapper>;
};

// set cm property in the sub-problems
using BoxTraits = TestTraits<TTag::OnePBulkBox, TTag::OnePLowDimBox>;
template<class TypeTag> struct CouplingManager<TypeTag, TTag::OnePBulkBox> { using type = typename BoxTraits::CouplingManager; };
template<class TypeTag> struct CouplingManager<TypeTag, TTag::OnePLowDimBox> { using type = typename BoxTraits::CouplingManager; };

using TpfaTraits = TestTraits<TTag::OnePBulkTpfa, TTag::OnePLowDimTpfa>;
template<class TypeTag> struct CouplingManager<TypeTag, TTag::OnePBulkTpfa> { using type = typename TpfaTraits::CouplingManager; };
template<class TypeTag> struct CouplingManager<TypeTag, TTag::OnePLowDimTpfa> { using type = typename TpfaTraits::CouplingManager; };

using MpfaTraits = TestTraits<TTag::OnePBulkMpfa, TTag::OnePLowDimMpfa>;
template<class TypeTag> struct CouplingManager<TypeTag, TTag::OnePBulkMpfa> { using type = typename MpfaTraits::CouplingManager; };
template<class TypeTag> struct CouplingManager<TypeTag, TTag::OnePLowDimMpfa> { using type = typename MpfaTraits::CouplingManager; };
} // end namespace Dumux::Properties

#endif
