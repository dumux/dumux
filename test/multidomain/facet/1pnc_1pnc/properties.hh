// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup FacetTests
 * \brief Properties for for the bulk problem in the 1pnc facet coupling test.
 */
#ifndef DUMUX_TEST_FACETCOUPLING_ONEPNC_BULK_PROPERTIES_HH
#define DUMUX_TEST_FACETCOUPLING_ONEPNC_BULK_PROPERTIES_HH

#ifndef BULKTYPETAG
#define BULKTYPETAG OnePNCBulkTpfa
#endif
#ifndef FACETTYPETAG
#define FACETTYPETAG OnePNCFacetTpfa
#endif
#ifndef DIMWORLD
#define DIMWORLD 2
#endif

#include <dune/foamgrid/foamgrid.hh>
#include <dune/alugrid/grid.hh>

#include <dumux/material/fluidsystems/h2on2.hh>
#include <dumux/material/fluidsystems/1padapter.hh>

#include <dumux/discretization/box.hh>
#include <dumux/discretization/cctpfa.hh>
#include <dumux/multidomain/facet/box/properties.hh>
#include <dumux/multidomain/facet/cellcentered/tpfa/properties.hh>
#include <dumux/multidomain/facet/cellcentered/mpfa/properties.hh>
#include <dumux/multidomain/facet/couplingmapper.hh>
#include <dumux/multidomain/facet/couplingmanager.hh>
#include <dumux/multidomain/traits.hh>
#include <dumux/porousmediumflow/1pnc/model.hh>

#include "problem_facet.hh"
#include "problem_bulk.hh"
#include "spatialparams.hh"

namespace Dumux::Properties {

// create the type tag nodes
namespace TTag {
struct BaseBulk {};
struct OnePNCBulk { using InheritsFrom = std::tuple<OnePNC, BaseBulk>; };
struct OnePNCNIBulk { using InheritsFrom = std::tuple<OnePNCNI, BaseBulk>; };

struct OnePNCBulkTpfa { using InheritsFrom = std::tuple<CCTpfaFacetCouplingModel, OnePNCBulk>; };
struct OnePNCNIBulkTpfa { using InheritsFrom = std::tuple<CCTpfaFacetCouplingModel, OnePNCNIBulk>; };

struct OnePNCBulkMpfa { using InheritsFrom = std::tuple<CCMpfaFacetCouplingModel, OnePNCBulk>; };
struct OnePNCNIBulkMpfa { using InheritsFrom = std::tuple<CCMpfaFacetCouplingModel, OnePNCNIBulk>; };

struct OnePNCBulkBox { using InheritsFrom = std::tuple<BoxFacetCouplingModel, OnePNCBulk>; };
struct OnePNCNIBulkBox { using InheritsFrom = std::tuple<BoxFacetCouplingModel, OnePNCNIBulk>; };
} // end namespace TTag

// Set the grid type (DIMWORLD is defined in CMakeLists.txt)
template<class TypeTag>
struct Grid<TypeTag, TTag::BaseBulk> { using type = Dune::ALUGrid<2, DIMWORLD, Dune::cube, Dune::nonconforming>; };
// Set the problem type
template<class TypeTag>
struct Problem<TypeTag, TTag::BaseBulk> { using type = OnePNCBulkProblem<TypeTag>; };
// set the spatial params
template<class TypeTag>
struct SpatialParams<TypeTag, TTag::BaseBulk>
{
    using type = OnePSpatialParams< GetPropType<TypeTag, Properties::GridGeometry>,
                                    GetPropType<TypeTag, Properties::Scalar> >;
};

// Set fluid configuration
template<class TypeTag>
struct FluidSystem<TypeTag, TTag::BaseBulk>
{
private:
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using H2ON2 = FluidSystems::H2ON2<Scalar, FluidSystems::H2ON2DefaultPolicy</*simplified=*/true>>;
public:
    using type = FluidSystems::OnePAdapter<H2ON2, H2ON2::liquidPhaseIdx>;
};

// create the type tag nodes
namespace TTag {
struct BaseFacet {};
struct OnePNCFacet { using InheritsFrom = std::tuple<OnePNC, BaseFacet>; }; // Isothermal case
struct OnePNCNIFacet { using InheritsFrom = std::tuple<OnePNCNI, BaseFacet>; }; // Non-Isothermal case

struct OnePNCFacetTpfa { using InheritsFrom = std::tuple<OnePNCFacet, CCTpfaModel>; };
struct OnePNCNIFacetTpfa { using InheritsFrom = std::tuple<OnePNCNIFacet, CCTpfaModel>; };

struct OnePNCFacetBox { using InheritsFrom = std::tuple<OnePNCFacet, BoxModel>; };
struct OnePNCNIFacetBox { using InheritsFrom = std::tuple<OnePNCNIFacet, BoxModel>; };
} // end namespace TTag

// Set the grid type
template<class TypeTag>
struct Grid<TypeTag, TTag::BaseFacet> { using type = Dune::FoamGrid<1, DIMWORLD>; };
// Set the problem type
template<class TypeTag>
struct Problem<TypeTag, TTag::BaseFacet> { using type = OnePNCLowDimProblem<TypeTag>; };
// set the spatial params
template<class TypeTag>
struct SpatialParams<TypeTag, TTag::BaseFacet>
{
    using type = OnePSpatialParams< GetPropType<TypeTag, Properties::GridGeometry>,
                                    GetPropType<TypeTag, Properties::Scalar> >;
};

// the fluid system
template<class TypeTag>
struct FluidSystem<TypeTag, TTag::BaseFacet>
{
private:
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using H2ON2 = FluidSystems::H2ON2<Scalar, FluidSystems::H2ON2DefaultPolicy</*simplified=*/true>>;
public:
    using type = FluidSystems::OnePAdapter<H2ON2, H2ON2::liquidPhaseIdx>;
};

using BulkTypeTag = TTag::BULKTYPETAG;
using FacetTypeTag = TTag::FACETTYPETAG;

// obtain/define some types to be used below in the property definitions and in main
class TestTraits
{
    using BulkGridGeometry = GetPropType<BulkTypeTag, Properties::GridGeometry>;
    using FacetGridGeometry = GetPropType<FacetTypeTag, Properties::GridGeometry>;
public:
    using MDTraits = Dumux::MultiDomainTraits<BulkTypeTag, FacetTypeTag>;
    using CouplingMapper = Dumux::FacetCouplingMapper<BulkGridGeometry, FacetGridGeometry>;
    using CouplingManager = Dumux::FacetCouplingManager<MDTraits, CouplingMapper>;
};

// specify coupling manager property in sub-problems
template<class TypeTag> struct CouplingManager<TypeTag, BulkTypeTag> { using type = typename TestTraits::CouplingManager; };
template<class TypeTag> struct CouplingManager<TypeTag, FacetTypeTag> { using type = typename TestTraits::CouplingManager; };

} // end namespace Dumux::Properties

#endif
