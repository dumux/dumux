// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup EmbeddedTests
 * \brief The properties of the embedded 1d-3d test with the variables of the local dofs
 */
#ifndef DUMUX_EMBEDDED_1P1P_VARIABLES_PROPERTIES_HH
#define DUMUX_EMBEDDED_1P1P_VARIABLES_PROPERTIES_HH

#ifndef BULKTYPETAG
#define BULKTYPETAG TissueBox
#endif
#ifndef LOWDIMTYPETAG
#define LOWDIMTYPETAG BloodFlowBox
#endif
#ifndef COUPLINGMODE
#define COUPLINGMODE Embedded1d3dCouplingMode::Average
#endif

#include <dune/grid/yaspgrid.hh>
#include <dune/foamgrid/foamgrid.hh>

#include <dumux/common/properties.hh>
#include <dumux/discretization/box.hh>
#include <dumux/discretization/pq1.hh>
#include <dumux/discretization/pq1bubble.hh>
#include <dumux/discretization/pq2.hh>

#include <dumux/porousmediumflow/1p/model.hh>
#include <dumux/material/components/constant.hh>
#include <dumux/material/fluidsystems/1pliquid.hh>

#include <dumux/multidomain/embedded/couplingmanager1d3d.hh>
#include <dumux/multidomain/traits.hh>

#include "problem_tissue.hh"
#include "problem_bloodflow.hh"
#include "spatialparams.hh"

namespace Dumux::Properties {

namespace TTag {
struct Tissue { using InheritsFrom = std::tuple<OneP>; };
struct TissueBox { using InheritsFrom = std::tuple<Tissue, BoxModel>; };
struct TissuePQ1BubbleHybrid { using InheritsFrom = std::tuple<Tissue, PQ1BubbleHybridModel>; };
struct TissuePQ2Hybrid { using InheritsFrom = std::tuple<Tissue, PQ2HybridModel>; };
struct TissuePQ1FE { using InheritsFrom = std::tuple<Tissue, PQ1FEModel>; };
struct TissuePQ2FE { using InheritsFrom = std::tuple<Tissue, PQ2FEModel>; };

struct BloodFlow { using InheritsFrom = std::tuple<OneP>; };
struct BloodFlowBox { using InheritsFrom = std::tuple<BloodFlow, BoxModel>; };
struct BloodFlowPQ2Hybrid { using InheritsFrom = std::tuple<BloodFlow, PQ2HybridModel>; };
} // end namespace TTag

template<class TypeTag>
struct Grid<TypeTag, TTag::Tissue>
{ using type = Dune::YaspGrid<3, Dune::EquidistantOffsetCoordinates<GetPropType<TypeTag, Properties::Scalar>, 3> >; };

template<class TypeTag>
struct Grid<TypeTag, TTag::BloodFlow> { using type = Dune::FoamGrid<1, 3>; };

template<class TypeTag>
struct EnableGridGeometryCache<TypeTag, TTag::Tissue> { static constexpr bool value = true; };
template<class TypeTag>
struct EnableGridGeometryCache<TypeTag, TTag::BloodFlow> { static constexpr bool value = true; };

//! Assemble the schemes that can also be assembled with volume variables with the variables of the local dofs
template<class TypeTag>
struct Transition::UseVariables<TypeTag, TTag::TissueBox> { static constexpr bool value = true; };
template<class TypeTag>
struct Transition::UseVariables<TypeTag, TTag::BloodFlowBox> { static constexpr bool value = true; };

template<class TypeTag>
struct Problem<TypeTag, TTag::Tissue> { using type = TissueProblem<TypeTag>; };
template<class TypeTag>
struct Problem<TypeTag, TTag::BloodFlow> { using type = BloodFlowProblem<TypeTag>; };

template<class TypeTag>
struct FluidSystem<TypeTag, TTag::Tissue>
{
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using type = FluidSystems::OnePLiquid<Scalar, Components::Constant<1, Scalar> >;
};
template<class TypeTag>
struct FluidSystem<TypeTag, TTag::BloodFlow>
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
template<class TypeTag>
struct SpatialParams<TypeTag, TTag::BloodFlow>
{
    using type = BloodFlowSpatialParams<GetPropType<TypeTag, Properties::GridGeometry>,
                                        GetPropType<TypeTag, Properties::Scalar>>;
};

template<class Traits>
using TheCouplingManager = Embedded1d3dCouplingManager<Traits, COUPLINGMODE>;

template<class TypeTag>
struct CouplingManager<TypeTag, TTag::BULKTYPETAG> { using type = TheCouplingManager<MultiDomainTraits<TypeTag, Properties::TTag::LOWDIMTYPETAG>>; };
template<class TypeTag>
struct CouplingManager<TypeTag, TTag::LOWDIMTYPETAG> { using type = TheCouplingManager<MultiDomainTraits<Properties::TTag::BULKTYPETAG, TypeTag>>; };

} // end namespace Dumux::Properties

#endif
