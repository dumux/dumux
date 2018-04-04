// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   See the file COPYING for full copying permissions.                      *
 *                                                                           *
 *   This program is free software: you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation, either version 2 of the License, or       *
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
 * \ingroup StaggeredDiscretization
 * \file
 *
 * \brief Defines a type tag and some properties for models using the staggered scheme.
          This scheme features degrees of freedom at the elements' centers and intersections (faces).
 * TODO: detailed documentation and figures
 */

#ifndef DUMUX_STAGGERD_PROPERTIES_HH
#define DUMUX_STAGGERD_PROPERTIES_HH

#include <dumux/common/properties.hh>

#include <dumux/discretization/methods.hh>
#include <dumux/discretization/fvproperties.hh>

#include <dumux/discretization/cellcentered/elementboundarytypes.hh>
#include <dumux/assembly/staggeredlocalresidual.hh>

#include <dumux/discretization/cellcentered/subcontrolvolume.hh>
#include <dumux/discretization/staggered/gridvariables.hh>
#include <dumux/discretization/staggered/gridfluxvariablescache.hh>
#include <dumux/discretization/staggered/gridvolumevariables.hh>
#include <dumux/discretization/staggered/fvgridgeometry.hh>
#include <dumux/discretization/staggered/gridfacevariables.hh>
#include <dumux/discretization/staggered/facesolution.hh>
#include <dumux/discretization/staggered/subcontrolvolumeface.hh>

#include <dune/istl/multitypeblockvector.hh>
#include <dune/istl/multitypeblockmatrix.hh>

#include <dumux/discretization/staggered/gridvariablestraits.hh>

namespace Dumux {

// forward declarations
class CCElementBoundaryTypes;

namespace Properties
{
//! Type tag for the staggered scheme.
NEW_TYPE_TAG(StaggeredModel, INHERITS_FROM(FiniteVolumeModel));

//! Set the default global face variables cache vector class
SET_PROP(StaggeredModel, GridFaceVariables)
{
private:
    using FaceVariables = typename GET_PROP_TYPE(TypeTag, FaceVariables);
    using FVGridGeometry = typename GET_PROP_TYPE(TypeTag, FVGridGeometry);
    using Problem = typename GET_PROP_TYPE(TypeTag, Problem);
    static constexpr auto enableCache = GET_PROP_VALUE(TypeTag, EnableGridFaceVariablesCache);

    using Traits = StaggeredGridFaceVariablesTraits<FaceVariables, Problem>;

public:
    using type = StaggeredGridFaceVariables<FVGridGeometry, Traits, enableCache>;
};

//! Cache the face variables per default
SET_BOOL_PROP(StaggeredModel, EnableGridFaceVariablesCache, true);

//! Set the default global volume variables cache vector class
SET_PROP(StaggeredModel, GridVolumeVariables)
{
private:
    using CellCenterPrimaryVariables = typename GET_PROP_TYPE(TypeTag, CellCenterPrimaryVariables);
    using VolumeVariables = typename GET_PROP_TYPE(TypeTag, VolumeVariables);
    using FVGridGeometry = typename GET_PROP_TYPE(TypeTag, FVGridGeometry);
    using Indices = typename GET_PROP_TYPE(TypeTag, Indices); // TODO extract indices from volumevariables
    using Problem = typename GET_PROP_TYPE(TypeTag, Problem);

    static constexpr auto enableCache = GET_PROP_VALUE(TypeTag, EnableGridVolumeVariablesCache);

    using Traits = StaggeredGridVolumeVariablesTraits<Problem, VolumeVariables, Indices>; // TODO extract indices from volumevariables

public:
    using type = StaggeredGridVolumeVariables<Traits, enableCache>;
};

//! Set the global flux variables cache vector class
SET_PROP(StaggeredModel, GridFluxVariablesCache)
{
private:
    static constexpr auto enableCache = GET_PROP_VALUE(TypeTag, EnableGridFluxVariablesCache);
    using Problem = typename GET_PROP_TYPE(TypeTag, Problem);
    using FluxVariablesCache = typename GET_PROP_TYPE(TypeTag, FluxVariablesCache);
public:
    using type = StaggeredGridFluxVariablesCache<Problem, FluxVariablesCache, enableCache>;
};

//! Set the face solution type
SET_PROP(StaggeredModel, StaggeredFaceSolution)
{
private:
    using FaceSolutionVector = typename GET_PROP_TYPE(TypeTag, FaceSolutionVector);
public:
    using type = StaggeredFaceSolution<FaceSolutionVector>;
};

//! Set the grid variables (volume, flux and face variables)
SET_PROP(StaggeredModel, GridVariables)
{
private:
    using GG = typename GET_PROP_TYPE(TypeTag, FVGridGeometry);
    using GVV = typename GET_PROP_TYPE(TypeTag, GridVolumeVariables);
    using GFVC = typename GET_PROP_TYPE(TypeTag, GridFluxVariablesCache);
    using GFV = typename GET_PROP_TYPE(TypeTag, GridFaceVariables);
public:
    using type = StaggeredGridVariables<GG, GVV, GFVC, GFV>;
};

//! Use the cell center element boundary types per default
SET_TYPE_PROP(StaggeredModel, ElementBoundaryTypes, CCElementBoundaryTypes);

//! Set the BaseLocalResidual to StaggeredLocalResidual
SET_TYPE_PROP(StaggeredModel, BaseLocalResidual, StaggeredLocalResidual<TypeTag>);

//! Definition of the indices for cell center and face dofs in the global solution vector
SET_PROP(StaggeredModel, DofTypeIndices)
{
    using CellCenterIdx = Dune::index_constant<0>;
    using FaceIdx = Dune::index_constant<1>;
};

//! The cell center primary variables
SET_TYPE_PROP(StaggeredModel,
              CellCenterPrimaryVariables,
              Dune::FieldVector<typename GET_PROP_TYPE(TypeTag, Scalar),
                                GET_PROP_VALUE(TypeTag, NumEqCellCenter)>);
//! The face primary variables
SET_TYPE_PROP(StaggeredModel,
              FacePrimaryVariables,
              Dune::FieldVector<typename GET_PROP_TYPE(TypeTag, Scalar),
                                GET_PROP_VALUE(TypeTag, NumEqFace)>);

//! Boundary types at a single degree of freedom
SET_PROP(StaggeredModel, BoundaryTypes)
{
private:
    static constexpr auto numEqCellCenter = GET_PROP_VALUE(TypeTag, NumEqCellCenter);
    static constexpr auto numEqFace = GET_PROP_VALUE(TypeTag, NumEqFace);
public:
    using type = BoundaryTypes<numEqCellCenter + numEqFace>;
};

//! Set one or different base epsilons for the calculations of the localJacobian's derivatives
SET_PROP(StaggeredModel, BaseEpsilon)
{
private:
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    static constexpr Scalar dCCdCC = 1e-8;
    static constexpr Scalar dCCdFace = 1e-8;
    static constexpr Scalar dFacedCC = 1e-8;
    static constexpr Scalar dFacedFace = 1e-8;

public:
    static constexpr auto getEps()
    {
        return std::array<std::array<Scalar, 2>, 2>{{{dCCdCC, dCCdFace},
                                                     {dFacedCC, dFacedFace}}};
    }
};

// TODO: bundle SolutionVector, JacobianMatrix
//       in LinearAlgebra traits

//! The type of a solution for the whole grid at a fixed time TODO: move to LinearAlgebra traits
SET_TYPE_PROP(StaggeredModel,
              CellCenterSolutionVector,
              Dune::BlockVector<typename GET_PROP_TYPE(TypeTag, CellCenterPrimaryVariables)>);

//! The type of a solution for the whole grid at a fixed time TODO: move to LinearAlgebra traits
SET_TYPE_PROP(StaggeredModel,
              FaceSolutionVector,
              Dune::BlockVector<typename GET_PROP_TYPE(TypeTag, FacePrimaryVariables)>);

//! default property value for the solution vector only used for monolithic solver TODO: move to LinearAlgebra traits
SET_PROP(StaggeredModel, SolutionVector)
{
private:
    using CellCenterSolutionVector = typename GET_PROP_TYPE(TypeTag, CellCenterSolutionVector);
    using FaceSolutionVector = typename GET_PROP_TYPE(TypeTag, FaceSolutionVector);
public:
    using type = Dune::MultiTypeBlockVector<CellCenterSolutionVector, FaceSolutionVector>;
};

//! Set the type of a global jacobian matrix from the solution types TODO: move to LinearAlgebra traits
SET_PROP(StaggeredModel, JacobianMatrix)
{
private:
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);

    static constexpr auto numEqCellCenter = GET_PROP_VALUE(TypeTag, NumEqCellCenter);
    static constexpr auto numEqFace = GET_PROP_VALUE(TypeTag, NumEqFace);

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
    using RowCellCenter = typename Dune::MultiTypeBlockVector<MatrixBlockCCToCC, MatrixBlockCCToFace>;
    using RowFace = typename Dune::MultiTypeBlockVector<MatrixBlockFaceToCC, MatrixBlockFaceToFace>;

    // the jacobian matrix
    using type = typename Dune::MultiTypeBlockMatrix<RowCellCenter, RowFace>;
};

} // namespace Properties
} // namespace Dumux

#endif
