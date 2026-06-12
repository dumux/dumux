// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup Discretization
 * \brief Defines a type tag and some properties for models using the pq3 scheme.
 */

#ifndef DUMUX_DISCRETIZATION_PQ3_HH
#define DUMUX_DISCRETIZATION_PQ3_HH

#include <concepts>
#include <type_traits>

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
#include <dumux/discretization/cvfe/hybrid/gridfluxvariablescache.hh>
#include <dumux/discretization/cvfe/gridvariablescache.hh>
#include <dumux/discretization/cvfe/variablesadapter.hh>
#include <dumux/discretization/pq3/fvgridgeometry.hh>
#include <dumux/discretization/pq3/fegriddiscretization.hh>
#include <dumux/discretization/cvfe/elementsolution.hh>
#include <dumux/discretization/cvfe/hybrid/fluxvariablescache.hh>

#include <dumux/assembly/localresidual.hh>
#include <dumux/discretization/fem/elementvariables.hh>
#include <dumux/discretization/fem/gridvariablescache.hh>
#include <dumux/discretization/cvfe/interpolationpointdata.hh>
#include <dumux/discretization/gridvariables.hh>

#include <dumux/flux/fluxvariablescaching.hh>

namespace Dumux::Properties {

namespace TTag {
struct PQ3Base { using InheritsFrom = std::tuple<GridProperties>; };
struct PQ3HybridModel { using InheritsFrom = std::tuple<FiniteVolumeModel, PQ3Base>; };
struct PQ3FEModel { using InheritsFrom = std::tuple<PQ3Base>; };
} // end namespace TTag

template<class TypeTag>
struct GridGeometry<TypeTag, TTag::PQ3HybridModel>
{
private:
    static constexpr bool enableCache = getPropValue<TypeTag, Properties::EnableGridGeometryCache>();
    using GridView = typename GetPropType<TypeTag, Properties::Grid>::LeafGridView;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
public:
    using type = PQ3FVGridGeometry<Scalar, GridView, enableCache>;
};

template<class TypeTag>
struct GridVolumeVariables<TypeTag, TTag::PQ3HybridModel>
{
private:
    static constexpr bool enableCache = getPropValue<TypeTag, Properties::EnableGridVolumeVariablesCache>();
    using Problem = GetPropType<TypeTag, Properties::Problem>;
    using Variables = Dumux::Detail::CVFE::VariablesAdapter<GetPropType<TypeTag, Properties::VolumeVariables>>;
    using Traits = Dumux::Detail::CVFE::CVFEDefaultGridVariablesCacheTraits<Problem, Variables>;
public:
    using type = Dumux::Detail::CVFE::CVFEGridVariablesCache<Traits, enableCache>;
};

template<class TypeTag>
struct FluxVariablesCache<TypeTag, TTag::PQ3HybridModel>
{
private:
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
public:
    using type = HybridCVFEFluxVariablesCache<Scalar, GridGeometry>;
};

template<class TypeTag>
struct GridFluxVariablesCache<TypeTag, TTag::PQ3HybridModel>
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

template<class TypeTag>
struct ElementBoundaryTypes<TypeTag, TTag::PQ3Base>
{
private:
    using Problem = GetPropType<TypeTag, Properties::Problem>;
    using GG = std::decay_t<decltype(std::declval<Problem>().gridDiscretization())>;
    using BoundaryTypes = typename ProblemTraits<Problem>::BoundaryTypes;
public:
    using type = std::conditional_t<
        Dumux::Detail::hasProblemBoundaryTypesForFaceFunction<Problem, typename GG::LocalView>(),
        Dumux::ElementIntersectionBoundaryTypes<BoundaryTypes>,
        Dumux::CVFEElementBoundaryTypes<BoundaryTypes>
    >;
};

//! Set the default FE grid discretization
template<class TypeTag>
struct GridGeometry<TypeTag, TTag::PQ3FEModel>
{
private:
    static constexpr bool enableCache = getPropValue<TypeTag, Properties::EnableGridGeometryCache>();
    using GridView = typename GetPropType<TypeTag, Properties::Grid>::LeafGridView;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
public:
    using type = Dumux::Experimental::PQ3FEGridDiscretization<Scalar, GridView, enableCache>;
};

template<class TypeTag>
struct GridVariables<TypeTag, TTag::PQ3FEModel>
{
private:
    using GG = GetPropType<TypeTag, Properties::GridGeometry>;
    static constexpr bool enableCache = getPropValue<TypeTag, Properties::EnableGridVolumeVariablesCache>();
    using Problem = GetPropType<TypeTag, Properties::Problem>;
    using Variables = Dumux::Detail::CVFE::VariablesAdapter<GetPropType<TypeTag, Properties::VolumeVariables>>;
    using IPDataCache = Dumux::CVFE::LocalBasisInterpolationPointData<GG>;
    using Traits = Dumux::Experimental::FEDefaultGridVariablesCacheTraits<Problem, Variables, IPDataCache>;
    using GVC = Dumux::Experimental::FEGridVariablesCache<Traits, enableCache>;
public:
    using type = Dumux::Experimental::GridVariables<GG, GVC>;
};

//! TODO: Replace property
template<class TypeTag>
struct EnableGridVolumeVariablesCache<TypeTag, TTag::PQ3FEModel> { static constexpr bool value = false; };

//! TODO: Replace and move to LinearAlgebra traits
template<class TypeTag>
struct SolutionVector<TypeTag, TTag::PQ3FEModel> { using type = Dune::BlockVector<GetPropType<TypeTag, Properties::PrimaryVariables>>; };

//! TODO: Replace and move to LinearAlgebra traits
template<class TypeTag>
struct JacobianMatrix<TypeTag, TTag::PQ3FEModel>
{
private:
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    enum { numEq = GetPropType<TypeTag, Properties::ModelTraits>::numEq() };
    using MatrixBlock = typename Dune::FieldMatrix<Scalar, numEq, numEq>;
public:
    using type = typename Dune::BCRSMatrix<MatrixBlock>;
};

} // namespace Dumux::Properties

namespace Dumux::Detail {

template<class Problem>
struct ProblemTraits<Problem, DiscretizationMethods::PQ3>
{
private:
    using GG = std::decay_t<decltype(std::declval<Problem>().gridDiscretization())>;
public:
    using GridGeometry = GG;
    using BoundaryTypes = Detail::BoundaryTypes<Problem, typename GG::LocalView>::type;
};

template<class TypeTag>
concept PQ3Model = std::is_same_v<
    typename GetPropType<TypeTag, Properties::GridGeometry>::DiscretizationMethod,
    DiscretizationMethods::PQ3
>;

template<class T>
concept PQ3HybridModel = PQ3Model<T> && Dumux::Properties::inheritsFrom<Properties::TTag::PQ3HybridModel, T>();

template<PQ3HybridModel TypeTag>
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
concept PQ3FEModel = PQ3Model<T> && Dumux::Properties::inheritsFrom<Properties::TTag::PQ3FEModel, T>();

template<PQ3FEModel TypeTag>
struct DiscretizationDefaultLocalOperator<TypeTag>
{ using type = Dumux::Experimental::LocalResidual<TypeTag>; };

} // end namespace Dumux::Detail

#endif
