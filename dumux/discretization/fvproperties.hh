// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup Discretization
 * \brief Declares properties required for finite-volume models models.
 */

#ifndef DUMUX_FV_PROPERTIES_HH
#define DUMUX_FV_PROPERTIES_HH

#include <dune/istl/bvector.hh>
#include <dune/istl/bcrsmatrix.hh>

#include <dumux/common/properties.hh>
#include <dumux/common/properties/grid.hh>
#include <dumux/common/boundarytypes.hh>

#include <dumux/discretization/fvgridvariables.hh>

namespace Dumux {
namespace Properties {

//! Type tag for finite-volume schemes.
// Create new type tags
namespace TTag {
struct FiniteVolumeModel { using InheritsFrom = std::tuple<GridProperties>; };
} // end namespace TTag

//! The grid variables
template<class TypeTag>
struct GridVariables<TypeTag, TTag::FiniteVolumeModel>
{
private:
    using GG = GetPropType<TypeTag, Properties::GridGeometry>;
    using GVV = GetPropType<TypeTag, Properties::GridVolumeVariables>;
    using GFVC = GetPropType<TypeTag, Properties::GridFluxVariablesCache>;
public:
    using type = FVGridVariables<GG, GVV, GFVC>;
};

//! We do not store the FVGeometry by default
template<class TypeTag>
struct EnableGridGeometryCache<TypeTag, TTag::FiniteVolumeModel> { static constexpr bool value = false; };

//! We do not store the volume variables by default
template<class TypeTag>
struct EnableGridVolumeVariablesCache<TypeTag, TTag::FiniteVolumeModel> { static constexpr bool value = false; };

//! disable flux variables data caching by default
template<class TypeTag>
struct EnableGridFluxVariablesCache<TypeTag, TTag::FiniteVolumeModel> { static constexpr bool value = false; };

// TODO: bundle SolutionVector, JacobianMatrix
//       in LinearAlgebra traits

//! The type of a solution for the whole grid at a fixed time TODO: move to LinearAlgebra traits
template<class TypeTag>
struct SolutionVector<TypeTag, TTag::FiniteVolumeModel> { using type = Dune::BlockVector<GetPropType<TypeTag, Properties::PrimaryVariables>>; };

//! Set the type of a global jacobian matrix from the solution types TODO: move to LinearAlgebra traits
template<class TypeTag>
struct JacobianMatrix<TypeTag, TTag::FiniteVolumeModel>
{
private:
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    enum { numEq = GetPropType<TypeTag, Properties::ModelTraits>::numEq() };
    using MatrixBlock = typename Dune::FieldMatrix<Scalar, numEq, numEq>;
public:
    using type = typename Dune::BCRSMatrix<MatrixBlock>;
};

} // namespace Properties
} // namespace Dumux

 #endif
