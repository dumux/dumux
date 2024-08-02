// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup Discretization
 * \brief Defines a type tag and some properties for models using the pq2 scheme.
 */

#ifndef DUMUX_DISCRETIZTAION_PQ2_HH
#define DUMUX_DISCRETIZTAION_PQ2_HH

#include <concepts>
#include <type_traits>

#include <dune/common/fvector.hh>
#include <dune/geometry/multilineargeometry.hh>

#include <dumux/common/properties.hh>
#include <dumux/common/boundaryflag.hh>
#include <dumux/common/typetraits/problem.hh>
#include <dumux/common/typetraits/boundary_.hh>

#include <dumux/assembly/cvfelocalresidual.hh>

#include <dumux/discretization/method.hh>
#include <dumux/discretization/fvproperties.hh>
#include <dumux/discretization/localdoftraits.hh>
#include <dumux/discretization/defaultlocaloperator.hh>
#include <dumux/discretization/elementboundarytypes.hh>

#include <dumux/discretization/cvfe/elementboundarytypes.hh>
#include <dumux/discretization/cvfe/hybrid/gridfluxvariablescache.hh>
#include <dumux/discretization/cvfe/gridvariablescache.hh>
#include <dumux/discretization/cvfe/variablesadapter.hh>
#include <dumux/discretization/pq2/fvgridgeometry.hh>
#include <dumux/discretization/cvfe/elementsolution.hh>
#include <dumux/discretization/cvfe/hybrid/fluxvariablescache.hh>

#include <dumux/flux/fluxvariablescaching.hh>

namespace Dumux::Properties {

//! Type tag for the pq2 scheme.
// Create new type tags
namespace TTag {
struct PQ2HybridModel { using InheritsFrom = std::tuple<FiniteVolumeModel>; };
} // end namespace TTag

//! Set the default for the grid geometry
template<class TypeTag>
struct GridGeometry<TypeTag, TTag::PQ2HybridModel>
{
private:
    static constexpr bool enableCache = getPropValue<TypeTag, Properties::EnableGridGeometryCache>();
    using GridView = typename GetPropType<TypeTag, Properties::Grid>::LeafGridView;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
public:
    using type = PQ2FVGridGeometry<Scalar, GridView, enableCache>;
};

//! The grid volume variables vector class
template<class TypeTag>
struct GridVolumeVariables<TypeTag, TTag::PQ2HybridModel>
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
struct FluxVariablesCache<TypeTag, TTag::PQ2HybridModel>
{
private:
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
public:
    using type = HybridCVFEFluxVariablesCache<Scalar, GridGeometry>;
};

//! The grid flux variables cache vector class
template<class TypeTag>
struct GridFluxVariablesCache<TypeTag, TTag::PQ2HybridModel>
{
private:
    static constexpr bool enableCache = getPropValue<TypeTag, Properties::EnableGridFluxVariablesCache>();
    using Problem = GetPropType<TypeTag, Properties::Problem>;

    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using FluxVariablesCache = GetPropTypeOr<TypeTag,
        Properties::FluxVariablesCache, FluxVariablesCaching::EmptyCache<Scalar>
    >;
public:
    using type = HybridCVFEGridFluxVariablesCache<Problem, FluxVariablesCache, enableCache>;
};

//! Set the default for the ElementBoundaryTypes
template<class TypeTag>
struct ElementBoundaryTypes<TypeTag, TTag::PQ2HybridModel>
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
struct ProblemTraits<Problem, DiscretizationMethods::PQ2>
{
private:
    using GG = std::decay_t<decltype(std::declval<Problem>().gridGeometry())>;
public:
    using GridGeometry = GG;
    // Determine BoundaryTypes dependent on the used problem interface, either boundaryTypes(element, scv) or  boundaryTypes(element, intersection)
    using BoundaryTypes = Detail::BoundaryTypes<Problem, typename GG::LocalView, typename GG::GridView::Intersection>::type;
};

template<class GridView>
struct LocalDofTraits<GridView, DiscretizationMethods::PQ2>
{
    static constexpr int dim = GridView::dimension;
    // Dofs are located at the vertices and element
    static constexpr int numCubeElementDofs = (1<<dim)+ dim*(1<<(dim-1));
};

template<class TypeTag>
concept PQ2HybridModel = std::is_same_v<
    typename GetPropType<TypeTag, Properties::GridGeometry>::DiscretizationMethod,
    DiscretizationMethods::PQ2
>;

template<PQ2HybridModel TypeTag>
struct DiscretizationDefaultLocalOperator<TypeTag>
{
    using type = CVFELocalResidual<TypeTag>;
};

} // end namespace Dumux::Detail

#endif
