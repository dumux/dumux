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

#ifndef DUMUX_DISCRETIZATION_PQ2_HH
#define DUMUX_DISCRETIZATION_PQ2_HH

#include <concepts>
#include <type_traits>
#include <utility>

#include <dumux/common/properties.hh>
#include <dumux/common/boundaryflag.hh>
#include <dumux/common/typetraits/problem.hh>
#include <dumux/common/typetraits/boundary_.hh>

#include <dumux/assembly/cvfelocalresidual_.hh>

#include <dumux/discretization/method.hh>
#include <dumux/discretization/fvproperties.hh>
#include <dumux/discretization/defaultlocaloperator.hh>
#include <dumux/discretization/elementboundarytypes.hh>

#include <dumux/discretization/cvfe/elementboundarytypes.hh>
#include <dumux/discretization/cvfe/hybrid/gridvariablescache.hh>
#include <dumux/discretization/cvfe/variablesadapter.hh>
#include <dumux/discretization/pq2/fvgridgeometry.hh>
#include <dumux/discretization/pq2/fegriddiscretization.hh>
#include <dumux/discretization/cvfe/elementsolution.hh>

#include <dumux/assembly/localresidual.hh>
#include <dumux/discretization/fem/elementvariables.hh>
#include <dumux/discretization/fem/gridvariablescache.hh>
#include <dumux/discretization/cvfe/interpolationpointdata.hh>
#include <dumux/discretization/gridvariables.hh>


namespace Dumux::Properties {

//! Type tag for the pq2 scheme.
// Create new type tags
namespace TTag {
struct PQ2Base { using InheritsFrom = std::tuple<GridProperties>; };
struct PQ2HybridModel { using InheritsFrom = std::tuple<PQ2Base>; };
struct PQ2FEModel { using InheritsFrom = std::tuple<PQ2Base>; };
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

//! The grid variables for the hybrid model
template<class TypeTag>
struct GridVariables<TypeTag, TTag::PQ2HybridModel>
{
private:
    using GG = GetPropType<TypeTag, Properties::GridGeometry>;
    // ToDo: Do not determine enableCache by EnableGridVolumeVariablesCache
    static constexpr bool enableCache = getPropValue<TypeTag, Properties::EnableGridVolumeVariablesCache>();
    using Problem = GetPropType<TypeTag, Properties::Problem>;
    using Variables = Dumux::Detail::CVFE::VariablesAdapter<GetPropType<TypeTag, Properties::VolumeVariables>>;
    using IPDataCache = Dumux::CVFE::LocalBasisInterpolationPointData<GG>;
    using Traits = Dumux::Experimental::CVFE::HybridCVFEDefaultGridVariablesCacheTraits<Problem, Variables, IPDataCache>;
    using GVC = Dumux::Experimental::CVFE::HybridCVFEGridVariablesCache<Traits, enableCache>;
public:
    using type = Dumux::Experimental::GridVariables<GG, GVC>;
};

//! TODO: Replace property
template<class TypeTag>
struct EnableGridGeometryCache<TypeTag, TTag::PQ2HybridModel> { static constexpr bool value = false; };

//! TODO: Replace property
template<class TypeTag>
struct EnableGridVolumeVariablesCache<TypeTag, TTag::PQ2HybridModel> { static constexpr bool value = false; };

//! TODO: Replace and move to LinearAlgebra traits
template<class TypeTag>
struct SolutionVector<TypeTag, TTag::PQ2HybridModel> { using type = Dune::BlockVector<GetPropType<TypeTag, Properties::PrimaryVariables>>; };

//! TODO: Replace and move to LinearAlgebra traits
template<class TypeTag>
struct JacobianMatrix<TypeTag, TTag::PQ2HybridModel>
{
private:
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    enum { numEq = GetPropType<TypeTag, Properties::ModelTraits>::numEq() };
    using MatrixBlock = typename Dune::FieldMatrix<Scalar, numEq, numEq>;
public:
    using type = typename Dune::BCRSMatrix<MatrixBlock>;
};

//! Set the default for the ElementBoundaryTypes
template<class TypeTag>
struct ElementBoundaryTypes<TypeTag, TTag::PQ2Base>
{
private:
    using Problem = GetPropType<TypeTag, Properties::Problem>;
    using GD = typename ProblemTraits<Problem>::GridDiscretization;
    using BoundaryTypes = typename ProblemTraits<Problem>::BoundaryTypes;
public:
    // Check if problem has new boundaryTypes interface
    // then use ElementIntersectionBoundaryTypes
    using type = std::conditional_t<
        Dumux::Detail::hasProblemBoundaryTypesForFaceFunction<Problem, typename GD::LocalView>(),
        Dumux::ElementIntersectionBoundaryTypes<BoundaryTypes>,
        Dumux::CVFEElementBoundaryTypes<BoundaryTypes>
    >;
};

//! Set the default FE grid discretization
template<class TypeTag>
struct GridGeometry<TypeTag, TTag::PQ2FEModel>
{
private:
    static constexpr bool enableCache = getPropValue<TypeTag, Properties::EnableGridGeometryCache>();
    using GridView = typename GetPropType<TypeTag, Properties::Grid>::LeafGridView;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
public:
    using type = Dumux::Experimental::PQ2FEGridDiscretization<Scalar, GridView, enableCache>;
};

template<class TypeTag>
struct GridVariables<TypeTag, TTag::PQ2FEModel>
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
struct EnableGridVolumeVariablesCache<TypeTag, TTag::PQ2FEModel> { static constexpr bool value = false; };

//! TODO: Replace and move to LinearAlgebra traits
template<class TypeTag>
struct SolutionVector<TypeTag, TTag::PQ2FEModel> { using type = Dune::BlockVector<GetPropType<TypeTag, Properties::PrimaryVariables>>; };

//! TODO: Replace and move to LinearAlgebra traits
template<class TypeTag>
struct JacobianMatrix<TypeTag, TTag::PQ2FEModel>
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
struct ProblemTraits<Problem, DiscretizationMethods::PQ2>
{
private:
    using GD = std::decay_t<decltype(std::declval<Problem>().gridDiscretization())>;
public:
    using GridDiscretization = GD;
    // Determine BoundaryTypes dependent on the used problem interface, either boundaryTypes(element, scv) or  boundaryTypes(element, boundaryFace)
    using BoundaryTypes = Detail::BoundaryTypes<Problem, typename GD::LocalView>::type;
};

template<class TypeTag>
concept PQ2Model = std::is_same_v<
    typename GetPropType<TypeTag, Properties::GridGeometry>::DiscretizationMethod,
    DiscretizationMethods::PQ2
>;

template<class T>
concept PQ2HybridModel = PQ2Model<T> && Dumux::Properties::inheritsFrom<Properties::TTag::PQ2HybridModel, T>();

template<PQ2HybridModel TypeTag>
struct DiscretizationDefaultLocalOperator<TypeTag>
{ using type = Dumux::Experimental::CVFELocalResidual<TypeTag>; };

template<class T>
concept PQ2FEModel = PQ2Model<T> && Dumux::Properties::inheritsFrom<Properties::TTag::PQ2FEModel, T>();

template<PQ2FEModel TypeTag>
struct DiscretizationDefaultLocalOperator<TypeTag>
{ using type = Dumux::Experimental::LocalResidual<TypeTag>; };

} // end namespace Dumux::Detail

#endif
