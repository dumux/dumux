// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   See the file COPYING for full copying permissions.                      *
 *                                                                           *
 *   This program is free software: you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation, either version 3 of the License, or       *
 *   (at your option) any later version.                                     *
 *                                                                           *
 *   This program is distributed in the hope that it will be useful,         *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of          *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the            *
 *   GNU General Public License for more details.                            *
 *                                                                           *
 *   You should have received a copy of the GNU General Public License       *
 *   along with this program.  If not, see <http://www.gnu.org/licenses/>.   *
 *****************************************************************************/
/*!
 * \file
 * \ingroup Discretization
 * \brief Declares properties required for finite-volume models models.
 */

#ifndef DUMUX_FV_PROPERTIES_HH
#define DUMUX_FV_PROPERTIES_HH

#include <dune/istl/bvector.hh>
#include <dune/istl/bcrsmatrix.hh>

#include <dumux/common/deprecated.hh>
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

DUNE_NO_DEPRECATED_BEGIN
//! Boundary types at a single degree of freedom
template<class TypeTag>
struct BoundaryTypes<TypeTag, TTag::FiniteVolumeModel>
{ using type = Dumux::BoundaryTypes<GetPropType<TypeTag, Properties::ModelTraits>::numEq()>; };
DUNE_NO_DEPRECATED_END
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
