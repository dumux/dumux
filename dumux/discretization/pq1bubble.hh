// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup Discretization
 * \brief Defines a type tag and some properties for models using the pq1bubble scheme.
 */

#ifndef DUMUX_DISCRETIZTAION_PQ1BUBBLE_HH
#define DUMUX_DISCRETIZTAION_PQ1BUBBLE_HH

#include <concepts>
#include <type_traits>

#include <dune/common/fvector.hh>
#include <dune/geometry/multilineargeometry.hh>

#include <dumux/common/properties.hh>
#include <dumux/common/boundaryflag.hh>
#include <dumux/common/typetraits/problem.hh>

#include <dumux/assembly/cvfelocalresidual.hh>

#include <dumux/discretization/method.hh>
#include <dumux/discretization/fvproperties.hh>
#include <dumux/discretization/localdoftraits.hh>
#include <dumux/discretization/defaultlocaloperator.hh>

#include <dumux/discretization/cvfe/elementboundarytypes.hh>
#include <dumux/discretization/cvfe/gridfluxvariablescache.hh>
#include <dumux/discretization/cvfe/gridvariablescache.hh>
#include <dumux/discretization/cvfe/variablesadapter.hh>
#include <dumux/discretization/pq1bubble/fvgridgeometry.hh>
#include <dumux/discretization/cvfe/elementsolution.hh>
#include <dumux/discretization/cvfe/fluxvariablescache.hh>

#include <dumux/flux/fluxvariablescaching.hh>

namespace Dumux::Properties {

//! Type tag for the pq1bubble scheme.
// Create new type tags
namespace TTag {
struct PQ1BubbleModel { using InheritsFrom = std::tuple<FiniteVolumeModel>; };
} // end namespace TTag

//! Set the default for the grid geometry
template<class TypeTag>
struct GridGeometry<TypeTag, TTag::PQ1BubbleModel>
{
private:
    static constexpr bool enableCache = getPropValue<TypeTag, Properties::EnableGridGeometryCache>();
    using GridView = typename GetPropType<TypeTag, Properties::Grid>::LeafGridView;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
public:
    using type = PQ1BubbleFVGridGeometry<Scalar, GridView, enableCache>;
};

//! The grid volume variables vector class
template<class TypeTag>
struct GridVolumeVariables<TypeTag, TTag::PQ1BubbleModel>
{
private:
    static constexpr bool enableCache = getPropValue<TypeTag, Properties::EnableGridVolumeVariablesCache>();
    using Problem = GetPropType<TypeTag, Properties::Problem>;
    using Variables = Dumux::Detail::CVFE::VariablesAdapter<GetPropType<TypeTag, Properties::VolumeVariables>>;
    using Traits = Dumux::Detail::CVFE::CVFEDefaultGridVariablesCacheTraits<Problem, Variables>;
public:
    using type = Dumux::Detail::CVFE::CVFEGridVariablesCache<Traits, enableCache>;
};

//! The flux variables cache class
template<class TypeTag>
struct FluxVariablesCache<TypeTag, TTag::PQ1BubbleModel>
{
private:
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
public:
    using type = CVFEFluxVariablesCache<Scalar, GridGeometry>;
};

//! The grid flux variables cache vector class
template<class TypeTag>
struct GridFluxVariablesCache<TypeTag, TTag::PQ1BubbleModel>
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

//! Set the default for the ElementBoundaryTypes
template<class TypeTag>
struct ElementBoundaryTypes<TypeTag, TTag::PQ1BubbleModel>
{
private:
    using Problem = GetPropType<TypeTag, Properties::Problem>;
    using BoundaryTypes = typename ProblemTraits<Problem>::BoundaryTypes;
public:
    using type = CVFEElementBoundaryTypes<BoundaryTypes>;
};

} // namespace Dumux::Properties

namespace Dumux::Detail {

template<class Problem>
struct ProblemTraits<Problem, DiscretizationMethods::PQ1Bubble>
{
private:
    using GG = std::decay_t<decltype(std::declval<Problem>().gridGeometry())>;
    using Element = typename GG::GridView::template Codim<0>::Entity;
    using SubControlVolume = typename GG::SubControlVolume;
public:
    using GridGeometry = GG;
    // BoundaryTypes is whatever the problem returns from boundaryTypes(element, scv)
    using BoundaryTypes = std::decay_t<decltype(std::declval<Problem>().boundaryTypes(std::declval<Element>(), std::declval<SubControlVolume>()))>;
};

template<class GridView>
struct LocalDofTraits<GridView, DiscretizationMethods::PQ1Bubble>
{
    static constexpr int dim = GridView::dimension;
    // Dofs are located at the vertices and element
    static constexpr int numCubeElementDofs = (1<<dim) + 1;
};

template<class TypeTag>
concept PQ1BubbleModel = std::is_same_v<
    typename GetPropType<TypeTag, Properties::GridGeometry>::DiscretizationMethod,
    DiscretizationMethods::PQ1Bubble
>;

template<PQ1BubbleModel TypeTag>
struct DiscretizationDefaultLocalOperator<TypeTag>
{
    using type = CVFELocalResidual<TypeTag>;
};

} // end namespace Dumux::Detail

#endif
