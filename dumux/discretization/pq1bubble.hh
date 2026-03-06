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
#include <dumux/common/concepts/variables_.hh>
#include <dumux/common/typetraits/problem.hh>
#include <dumux/common/typetraits/boundary_.hh>

#include <dumux/assembly/cvfelocalresidual.hh>
#include <dumux/assembly/cvfelocalresidual_.hh>

#include <dumux/discretization/method.hh>
#include <dumux/discretization/fvproperties.hh>
#include <dumux/discretization/defaultlocaloperator.hh>
#include <dumux/discretization/elementboundarytypes.hh>

#include <dumux/discretization/cvfe/elementboundarytypes.hh>
#include <dumux/discretization/cvfe/gridfluxvariablescache.hh>
#include <dumux/discretization/cvfe/hybrid/gridfluxvariablescache.hh>
#include <dumux/discretization/cvfe/gridvariablescache.hh>
#include <dumux/discretization/cvfe/variablesadapter.hh>
#include <dumux/discretization/pq1bubble/fvgridgeometry.hh>
#include <dumux/discretization/pq1bubble/fegridgeometry.hh>
#include <dumux/discretization/cvfe/elementsolution.hh>
#include <dumux/discretization/cvfe/fluxvariablescache.hh>
#include <dumux/discretization/cvfe/hybrid/fluxvariablescache.hh>

#include <dumux/assembly/localresidual.hh>
#include <dumux/discretization/fem/elementvariables.hh>
#include <dumux/discretization/fem/gridvariablescache.hh>
#include <dumux/discretization/cvfe/interpolationpointdata.hh>
#include <dumux/discretization/gridvariables.hh>

#include <dumux/flux/fluxvariablescaching.hh>

namespace Dumux::Properties {

//! Type tag for the pq1bubble scheme.
// Create new type tags
namespace TTag {
struct PQ1BubbleBase { using InheritsFrom = std::tuple<GridProperties>; };
struct PQ1BubbleFVBase { using InheritsFrom = std::tuple<FiniteVolumeModel,PQ1BubbleBase>; };
struct PQ1BubbleModel { using InheritsFrom = std::tuple<PQ1BubbleFVBase>; };
struct PQ1BubbleHybridModel { using InheritsFrom = std::tuple<PQ1BubbleFVBase>; };
struct PQ1BubbleFEModel { using InheritsFrom = std::tuple<PQ1BubbleBase>; };
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

//! Set the default for the grid geometry for hybrid model
template<class TypeTag>
struct GridGeometry<TypeTag, TTag::PQ1BubbleHybridModel>
{
private:
    static constexpr bool enableCache = getPropValue<TypeTag, Properties::EnableGridGeometryCache>();
    using GridView = typename GetPropType<TypeTag, Properties::Grid>::LeafGridView;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using QuadTraits = PQ1BubbleQuadratureTraits<GridView>;
    // For hybrid scheme we use two bubble functions on cubes
    // when only using one bubble function stability issues can arise, especially for coupled problems
    static constexpr std::size_t numCubeBubbles = 2;
    using Traits = HybridPQ1BubbleCVFEGridGeometryTraits<PQ1BubbleDefaultGridGeometryTraits<GridView, PQ1BubbleMapperTraits<GridView, numCubeBubbles>, QuadTraits>>;
public:
    using type = PQ1BubbleFVGridGeometry<Scalar, GridView, enableCache, Traits>;
};

//! Set the default for the grid geometry for fe model
template<class TypeTag>
struct GridGeometry<TypeTag, TTag::PQ1BubbleFEModel>
{
private:
    static constexpr bool enableCache = getPropValue<TypeTag, Properties::EnableGridGeometryCache>();
    using GridView = typename GetPropType<TypeTag, Properties::Grid>::LeafGridView;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
public:
    using type = PQ1BubbleFEGridGeometry<Scalar, GridView, enableCache>;
};

template<class TypeTag>
struct GridVariables<TypeTag, TTag::PQ1BubbleFEModel>
{
private:
    using GG = GetPropType<TypeTag, Properties::GridGeometry>;
    // ToDo: Do not determine enableCache by EnableGridVolumeVariablesCache
    static constexpr bool enableCache = getPropValue<TypeTag, Properties::EnableGridVolumeVariablesCache>();
    using Problem = GetPropType<TypeTag, Properties::Problem>;
    using Variables = Dumux::Detail::CVFE::VariablesAdapter<GetPropType<TypeTag, Properties::VolumeVariables>>;
    using IPDataCache = Dumux::CVFE::LocalBasisInterpolationPointData<GG>;
    using Traits = Dumux::Experimental::FE::FEDefaultGridVariablesCacheTraits<Problem, Variables, IPDataCache>;
    using GVC = Dumux::Experimental::FE::FEGridVariablesCache<Traits, enableCache>;
public:
    using type = Dumux::GridVariables<GG, GVC>;
};

//! TODO: Replace property
template<class TypeTag>
struct EnableGridVolumeVariablesCache<TypeTag, TTag::PQ1BubbleFEModel> { static constexpr bool value = false; };

//! TODO: Replace and move to LinearAlgebra traits
template<class TypeTag>
struct SolutionVector<TypeTag, TTag::PQ1BubbleFEModel> { using type = Dune::BlockVector<GetPropType<TypeTag, Properties::PrimaryVariables>>; };

//! TODO: Replace and move to LinearAlgebra traits
template<class TypeTag>
struct JacobianMatrix<TypeTag, TTag::PQ1BubbleFEModel>
{
private:
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    enum { numEq = GetPropType<TypeTag, Properties::ModelTraits>::numEq() };
    using MatrixBlock = typename Dune::FieldMatrix<Scalar, numEq, numEq>;
public:
    using type = typename Dune::BCRSMatrix<MatrixBlock>;
};

//! The grid volume variables vector class
template<class TypeTag>
struct GridVolumeVariables<TypeTag, TTag::PQ1BubbleFVBase>
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

//! The flux variables cache class
template<class TypeTag>
struct FluxVariablesCache<TypeTag, TTag::PQ1BubbleHybridModel>
{
private:
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
public:
    using type = HybridCVFEFluxVariablesCache<Scalar, GridGeometry>;
};

//! The grid flux variables cache vector class
template<class TypeTag>
struct GridFluxVariablesCache<TypeTag, TTag::PQ1BubbleHybridModel>
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
struct ElementBoundaryTypes<TypeTag, TTag::PQ1BubbleBase>
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
struct ProblemTraits<Problem, DiscretizationMethods::PQ1Bubble>
{
private:
    using GG = std::decay_t<decltype(std::declval<Problem>().gridGeometry())>;
public:
    using GridGeometry = GG;
    // Determine BoundaryTypes dependent on the used problem interface, either boundaryTypes(element, scv) or  boundaryTypes(element, intersection)
    using BoundaryTypes = Detail::BoundaryTypes<Problem, typename GG::LocalView, typename GG::GridView::Intersection>::type;
};

template<class TypeTag>
concept PQ1BubbleModel = std::is_same_v<
    typename GetPropType<TypeTag, Properties::GridGeometry>::DiscretizationMethod,
    DiscretizationMethods::PQ1Bubble
>;

template<PQ1BubbleModel TypeTag>
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

template<class T>
concept PQ1BubbleFEModel = PQ1BubbleModel<T> && Dumux::Properties::inheritsFrom<Properties::TTag::PQ1BubbleFEModel, T>();

template<PQ1BubbleFEModel TypeTag>
struct DiscretizationDefaultLocalOperator<TypeTag>
{ using type = LocalResidual<TypeTag>; };

} // end namespace Dumux::Detail

#endif
