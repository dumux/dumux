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
 * \brief Defines a type tag and some properties for models using the staggered scheme.
          This scheme features degrees of freedom at the elements' centers and intersections (faces).
 * TODO: detailed documentation and figures
 */

#ifndef DUMUX_DISCRETIZATION_STAGGERD_HH
#define DUMUX_DISCRETIZATION_STAGGERD_HH

#include <dumux/common/properties.hh>
#include <dumux/common/typetraits/problem.hh>

#include <dumux/discretization/method.hh>
#include <dumux/discretization/fvproperties.hh>
#include <dumux/flux/fluxvariablescaching.hh>

#include <dumux/discretization/cellcentered/elementboundarytypes.hh>
#include <dumux/assembly/staggeredlocalresidual.hh>

#include <dumux/discretization/cellcentered/subcontrolvolume.hh>
#include <dumux/discretization/staggered/gridvariables.hh>
#include <dumux/discretization/staggered/gridfluxvariablescache.hh>
#include <dumux/discretization/staggered/fvgridgeometry.hh>
#include <dumux/discretization/staggered/gridfacevariables.hh>
#include <dumux/discretization/staggered/facesolution.hh>
#include <dumux/discretization/staggered/subcontrolvolumeface.hh>

#include <dune/istl/multitypeblockvector.hh>
#include <dune/istl/multitypeblockmatrix.hh>
#include <dumux/common/deprecated.hh>

namespace Dumux {

// forward declarations
class CCElementBoundaryTypes;

namespace Properties
{
//! Type tag for the staggered scheme.
// Create new type tags
namespace TTag {
struct StaggeredModel { using InheritsFrom = std::tuple<FiniteVolumeModel>; };
} // end namespace TTag

//! Set the default global face variables cache vector class
template<class TypeTag>
struct GridFaceVariables<TypeTag, TTag::StaggeredModel>
{
private:
    using Problem = GetPropType<TypeTag, Properties::Problem>;
    using FaceVariables = GetPropType<TypeTag, Properties::FaceVariables>;
    static constexpr auto enableCache = getPropValue<TypeTag, Properties::EnableGridFaceVariablesCache>();
public:
    using type = StaggeredGridFaceVariables<Problem, FaceVariables, enableCache>;
};

//! Cache the face variables per default
template<class TypeTag>
struct EnableGridFaceVariablesCache<TypeTag, TTag::StaggeredModel> { static constexpr bool value = true; };

//! Set the global flux variables cache vector class
template<class TypeTag>
struct GridFluxVariablesCache<TypeTag, TTag::StaggeredModel>
{
private:
    using Problem = GetPropType<TypeTag, Properties::Problem>;
    using FluxVariablesCache = GetPropType<TypeTag, Properties::FluxVariablesCache>;
    using FluxVariablesCacheFiller =  GetPropType<TypeTag, Properties::FluxVariablesCacheFiller>;
    static constexpr auto enableCache = getPropValue<TypeTag, Properties::EnableGridFluxVariablesCache>();
    static constexpr auto upwindSchemeOrder = getPropValue<TypeTag, Properties::UpwindSchemeOrder>();
public:
    using type = StaggeredGridFluxVariablesCache<Problem, FluxVariablesCache, FluxVariablesCacheFiller, enableCache, upwindSchemeOrder>;
};

//! Set the face solution type
template<class TypeTag>
struct StaggeredFaceSolution<TypeTag, TTag::StaggeredModel>
{
private:
    using FaceSolutionVector = GetPropType<TypeTag, Properties::FaceSolutionVector>;
public:
    using type = Dumux::StaggeredFaceSolution<FaceSolutionVector>;
};

//! Set the grid variables (volume, flux and face variables)
template<class TypeTag>
struct GridVariables<TypeTag, TTag::StaggeredModel>
{
private:
    using GG = GetPropType<TypeTag, Properties::GridGeometry>;
    using GVV = GetPropType<TypeTag, Properties::GridVolumeVariables>;
    using GFVC = GetPropType<TypeTag, Properties::GridFluxVariablesCache>;
    using GFV = GetPropType<TypeTag, Properties::GridFaceVariables>;
public:
    using type = StaggeredGridVariables<GG, GVV, GFVC, GFV>;
};

//! Use the cell center element boundary types per default
template<class TypeTag>
struct ElementBoundaryTypes<TypeTag, TTag::StaggeredModel> { using type = CCElementBoundaryTypes; };

//! Set the BaseLocalResidual to StaggeredLocalResidual
template<class TypeTag>
struct BaseLocalResidual<TypeTag, TTag::StaggeredModel> { using type = StaggeredLocalResidual<TypeTag>; };

//! The cell center primary variables
template<class TypeTag>
struct CellCenterPrimaryVariables<TypeTag, TTag::StaggeredModel>
{
    using type = Dune::FieldVector<GetPropType<TypeTag, Properties::Scalar>,
                                   getPropValue<TypeTag, Properties::NumEqCellCenter>()>;
};

//! The face primary variables
template<class TypeTag>
struct FacePrimaryVariables<TypeTag, TTag::StaggeredModel>
{
    using type = Dune::FieldVector<GetPropType<TypeTag, Properties::Scalar>,
                                   getPropValue<TypeTag, Properties::NumEqFace>()>;
};

DUNE_NO_DEPRECATED_BEGIN
//! Boundary types at a single degree of freedom
template<class TypeTag>
struct BoundaryTypes<TypeTag, TTag::StaggeredModel> { using type = Dumux::BoundaryTypes<GetPropType<TypeTag, Properties::ModelTraits>::numEq()>; };
DUNE_NO_DEPRECATED_END
// TODO: bundle SolutionVector, JacobianMatrix
//       in LinearAlgebra traits

//! The type of a solution for the whole grid at a fixed time TODO: move to LinearAlgebra traits
template<class TypeTag>
struct CellCenterSolutionVector<TypeTag, TTag::StaggeredModel>
{ using type = Dune::BlockVector<GetPropType<TypeTag, Properties::CellCenterPrimaryVariables>>; };

//! The type of a solution for the whole grid at a fixed time TODO: move to LinearAlgebra traits
template<class TypeTag>
struct FaceSolutionVector<TypeTag, TTag::StaggeredModel>
{ using type = Dune::BlockVector<GetPropType<TypeTag, Properties::FacePrimaryVariables>>; };

//! default property value for the solution vector only used for monolithic solver TODO: move to LinearAlgebra traits
template<class TypeTag>
struct SolutionVector<TypeTag, TTag::StaggeredModel>
{
private:
    using FaceSolutionVector = GetPropType<TypeTag, Properties::FaceSolutionVector>;
    using CellCenterSolutionVector = GetPropType<TypeTag, Properties::CellCenterSolutionVector>;
public:
    using type = Dune::MultiTypeBlockVector<FaceSolutionVector, CellCenterSolutionVector>;
};

//! Set the type of a global jacobian matrix from the solution types TODO: move to LinearAlgebra traits
template<class TypeTag>
struct JacobianMatrix<TypeTag, TTag::StaggeredModel>
{
private:
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;

    static constexpr auto numEqCellCenter = getPropValue<TypeTag, Properties::NumEqCellCenter>();
    static constexpr auto numEqFace = getPropValue<TypeTag, Properties::NumEqFace>();

public:
    // the sub-blocks
    using MatrixLittleBlockCCToCC = typename Dune::FieldMatrix<Scalar, numEqCellCenter, numEqCellCenter>;
    using MatrixLittleBlockCCToFace = typename Dune::FieldMatrix<Scalar, numEqCellCenter, numEqFace>;

    using MatrixLittleBlockFaceToFace = typename Dune::FieldMatrix<Scalar, numEqFace, numEqFace>;
    using MatrixLittleBlockFaceToCC = typename Dune::FieldMatrix<Scalar, numEqFace, numEqCellCenter>;

    // the BCRS matrices of the subproblems as big blocks
    using MatrixBlockCCToCC = typename Dune::BCRSMatrix<MatrixLittleBlockCCToCC>;
    using MatrixBlockCCToFace = typename Dune::BCRSMatrix<MatrixLittleBlockCCToFace>;

    using MatrixBlockFaceToFace = typename Dune::BCRSMatrix<MatrixLittleBlockFaceToFace>;
    using MatrixBlockFaceToCC = typename Dune::BCRSMatrix<MatrixLittleBlockFaceToCC>;

    // the row types
    using RowFace = typename Dune::MultiTypeBlockVector<MatrixBlockFaceToFace, MatrixBlockFaceToCC>;
    using RowCellCenter = typename Dune::MultiTypeBlockVector<MatrixBlockCCToFace, MatrixBlockCCToCC>;

    // the jacobian matrix
    using type = typename Dune::MultiTypeBlockMatrix<RowFace, RowCellCenter>;
};

} // namespace Properties

namespace Detail {

template<class Problem>
struct ProblemTraits<Problem, DiscretizationMethod::staggered>
{
private:
    using GG = std::decay_t<decltype(std::declval<Problem>().gridGeometry())>;
    using Element = typename GG::GridView::template Codim<0>::Entity;
    using SubControlVolumeFace = typename GG::SubControlVolumeFace;
public:
    using GridGeometry = GG;
    // BoundaryTypes is whatever the problem returns from boundaryTypes(element, scvf)
    using BoundaryTypes = std::decay_t<decltype(std::declval<Problem>().boundaryTypes(std::declval<Element>(), std::declval<SubControlVolumeFace>()))>;
};

} // end namespace Detail

} // namespace Dumux

#endif
