// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup Discretization
 * \brief Defines a type tag and some properties for models using the diamond scheme.
          This scheme features degrees of freedom at the elements' centers and intersections (faces).
 */

#ifndef DUMUX_DISCRETIZATION_FACECENTERED_DIAMOND_HH
#define DUMUX_DISCRETIZATION_FACECENTERED_DIAMOND_HH

#include <concepts>
#include <type_traits>

#include <dumux/common/properties.hh>
#include <dumux/common/concepts/variables_.hh>
#include <dumux/common/typetraits/problem.hh>

#include <dumux/assembly/cvfelocalresidual.hh>
#include <dumux/assembly/cvfelocalresidual_.hh>
#include <dumux/discretization/method.hh>
#include <dumux/discretization/fvproperties.hh>
#include <dumux/discretization/defaultlocaloperator.hh>
#include <dumux/discretization/elementboundarytypes.hh>
#include <dumux/discretization/fvgridvariables.hh>
#include <dumux/flux/fluxvariablescaching.hh>

#include <dumux/discretization/facecentered/diamond/fvgridgeometry.hh>
#include <dumux/discretization/cvfe/gridvariablescache.hh>
#include <dumux/discretization/cvfe/variablesadapter.hh>
#include <dumux/discretization/cvfe/gridfluxvariablescache.hh>
#include <dumux/discretization/cvfe/fluxvariablescache.hh>
#include <dumux/discretization/cvfe/elementboundarytypes.hh>
#include <dumux/discretization/cvfe/elementsolution.hh>

namespace Dumux::Properties {

//! Type tag for the staggered scheme.
// Create new type tags
namespace TTag {
struct FaceCenteredDiamondModel { using InheritsFrom = std::tuple<FiniteVolumeModel>; };
} // end namespace TTag

//! Set the default for the grid geometry
template<class TypeTag>
struct GridGeometry<TypeTag, TTag::FaceCenteredDiamondModel>
{
private:
    static constexpr bool enableCache = getPropValue<TypeTag, Properties::EnableGridGeometryCache>();
    using GridView = typename GetPropType<TypeTag, Properties::Grid>::LeafGridView;
public:
    using type = FaceCenteredDiamondFVGridGeometry<GridView, enableCache>;
};

//! The grid volume variables vector class
template<class TypeTag>
struct GridVolumeVariables<TypeTag, TTag::FaceCenteredDiamondModel>
{
private:
    static constexpr bool enableCache = getPropValue<TypeTag, Properties::EnableGridVolumeVariablesCache>();
    using Problem = GetPropType<TypeTag, Properties::Problem>;
    using Variables = Dumux::Detail::CVFE::VariablesAdapter<GetPropType<TypeTag, Properties::VolumeVariables>>;
    using Traits = Dumux::Detail::CVFE::CVFEDefaultGridVariablesCacheTraits<Problem, Variables>;
public:
    using type = Dumux::Detail::CVFE::CVFEGridVariablesCache<Traits, enableCache>;
};

//! Set the global flux variables cache vector class
template<class TypeTag>
struct GridFluxVariablesCache<TypeTag, TTag::FaceCenteredDiamondModel>
{
private:
    static constexpr bool enableCache = getPropValue<TypeTag, Properties::EnableGridFluxVariablesCache>();
    using Problem = GetPropType<TypeTag, Properties::Problem>;

    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using FluxVariablesCache = GetPropTypeOr<TypeTag,
        Properties::FluxVariablesCache, FluxVariablesCaching::EmptyCache<Scalar>
    >;
public:
    using type = CVFEGridFluxVariablesCache<Problem, FluxVariablesCache, enableCache>;
};

//! Set the grid variables (volume, flux and face variables)
template<class TypeTag>
struct GridVariables<TypeTag, TTag::FaceCenteredDiamondModel>
{
private:
    using GG = GetPropType<TypeTag, Properties::GridGeometry>;
    using GVV = GetPropType<TypeTag, Properties::GridVolumeVariables>;
    using GFVC = GetPropType<TypeTag, Properties::GridFluxVariablesCache>;
public:
    using type = FVGridVariables<GG, GVV, GFVC>;
};

//! The flux variables cache type
template<class TypeTag>
struct FluxVariablesCache<TypeTag, TTag::FaceCenteredDiamondModel>
{
private:
    using S = GetPropType<TypeTag, Properties::Scalar>;
    using GG = GetPropType<TypeTag, Properties::GridGeometry>;
public:
    using type = CVFEFluxVariablesCache<S, GG>;
};

//! Set the default for the ElementBoundaryTypes
template<class TypeTag>
struct ElementBoundaryTypes<TypeTag, TTag::FaceCenteredDiamondModel>
{
private:
    using Problem = GetPropType<TypeTag, Properties::Problem>;
    using GG = std::decay_t<decltype(std::declval<Problem>().gridGeometry())>;
    using BoundaryTypes = typename ProblemTraits<Problem>::BoundaryTypes;
public:
    // Check if problem has new boundaryTypes interface
    // then use ElementIntersectionBoundaryTypes
    using type = std::conditional_t<
        Dumux::Detail::hasProblemBoundaryTypesForIntersectionFunction<Problem, typename GG::LocalView, typename GG::GridView::Intersection>(),
        Dumux::ElementIntersectionBoundaryTypes<BoundaryTypes>,
        Dumux::CVFEElementBoundaryTypes<BoundaryTypes>
    >;
};

} // namespace Dumux::Properties

namespace Dumux::Detail {

template<class Problem>
struct ProblemTraits<Problem, DiscretizationMethods::FCDiamond>
{
private:
    using GG = std::decay_t<decltype(std::declval<Problem>().gridGeometry())>;
public:
    using GridGeometry = GG;
    // Determine BoundaryTypes dependent on the used problem interface, either boundaryTypes(element, scv) or  boundaryTypes(element, intersection)
    using BoundaryTypes = Detail::BoundaryTypes<Problem, typename GG::LocalView, typename GG::GridView::Intersection>::type;
};

template<class TypeTag>
concept FaceCenteredDiamondModel = std::is_same_v<
    typename GetPropType<TypeTag, Properties::GridGeometry>::DiscretizationMethod,
    DiscretizationMethods::FCDiamond
>;

template<FaceCenteredDiamondModel TypeTag>
struct DiscretizationDefaultLocalOperator<TypeTag>
{
private:
    using GV = GetPropType<TypeTag, Properties::GridVariables>;
    static constexpr bool usesGeneralGridVariables =
        Dumux::Concept::GridVariables<GV> && !Dumux::Concept::FVGridVariables<GV>;
public:
    using type = std::conditional_t<usesGeneralGridVariables,
                                    Dumux::Experimental::CVFELocalResidual<TypeTag>,
                                    Dumux::CVFELocalResidual<TypeTag>>;
};

} // end namespace Dumux::Detail

#endif
