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
 * \file
 * \ingroup MixedDimension
 * \brief Base properties for problems of mixed dimension
 * two sub problems
 */

#ifndef DUMUX_MIXEDDIMENSION_PROPERTIES_HH
#define DUMUX_MIXEDDIMENSION_PROPERTIES_HH

#include <dune/common/fvector.hh>
#include <dune/common/fmatrix.hh>
#include <dune/common/indices.hh>
#include <dune/istl/bcrsmatrix.hh>
#include <dune/istl/multitypeblockvector.hh>
#include <dune/istl/multitypeblockmatrix.hh>

#include <dumux/common/basicproperties.hh>
#include <dumux/nonlinear/newtonmethod.hh>
#include <dumux/common/timemanager.hh>

#include <dumux/mixeddimension/subproblemlocaljacobian.hh>

namespace Dumux
{

// forward declarations
template <class TypeTag> class MixedDimensionModel;
// template <class TypeTag> class BulkLocalJacobian;
// template <class TypeTag> class LowDimLocalJacobian;
template <class TypeTag> class MixedDimensionNewtonController;
template <class TypeTag> class MixedDimensionAssembler;

namespace Properties
{
// NumericModel provides Scalar, GridCreator, ParameterTree
NEW_TYPE_TAG(MixedDimension, INHERITS_FROM(NewtonMethod, LinearSolverTypeTag, NumericModel));

NEW_PROP_TAG(Model); //!< The type of the base class of the model
NEW_PROP_TAG(BulkLocalJacobian); //!< The type of the bulk local jacobian operator
NEW_PROP_TAG(LowDimLocalJacobian); //!< The type of the low dim local jacobian operator

NEW_PROP_TAG(SolutionVector); //!< Vector containing all primary variable vector of the grid

NEW_PROP_TAG(JacobianAssembler); //!< Assembles the global jacobian matrix
NEW_PROP_TAG(JacobianMatrix); //!< Type of the global jacobian matrix

// high level simulation control
NEW_PROP_TAG(TimeManager);  //!< Manages the simulation time

/*!
 * \brief Specify which kind of method should be used to numerically
 * calculate the partial derivatives of the residual.
 *
 * -1 means backward differences, 0 means central differences, 1 means
 * forward differences. By default we use central differences.
 */
NEW_PROP_TAG(ImplicitNumericDifferenceMethod);

//! the maximum allowed number of timestep divisions for the
//! Newton solver
NEW_PROP_TAG(ImplicitMaxTimeStepDivisions);

// property tags that will be set in the problem at hand
NEW_PROP_TAG(BulkProblemTypeTag);
NEW_PROP_TAG(LowDimProblemTypeTag);
NEW_PROP_TAG(Problem);
NEW_PROP_TAG(CouplingManager);
NEW_PROP_TAG(MixedDimensionUseIterativeSolver);
NEW_PROP_TAG(SubProblemBlockIndices);

// forward declarations
NEW_PROP_TAG(NumEq);

//////////////////////////////////////////////////////////////////
// Properties
//////////////////////////////////////////////////////////////////
//! use the plain newton method by default
SET_TYPE_PROP(MixedDimension, NewtonMethod, NewtonMethod<TypeTag>);

//! set default values
SET_TYPE_PROP(MixedDimension, NewtonController, MixedDimensionNewtonController<TypeTag>);

//! Set the assembler
SET_TYPE_PROP(MixedDimension, JacobianAssembler, MixedDimensionAssembler<TypeTag>);

//! Set the BaseModel to MixedDimensionModel
SET_TYPE_PROP(MixedDimension, Model, MixedDimensionModel<TypeTag>);

//! The local jacobian operators
SET_TYPE_PROP(MixedDimension, BulkLocalJacobian, BulkLocalJacobian<TypeTag>);
SET_TYPE_PROP(MixedDimension, LowDimLocalJacobian, LowDimLocalJacobian<TypeTag>);

//! use forward differences to calculate the jacobian by default
SET_INT_PROP(MixedDimension, ImplicitNumericDifferenceMethod, +1);

//! default property value for the time manager
SET_TYPE_PROP(MixedDimension, TimeManager, TimeManager<TypeTag>);

//! default property is monolithic solver
SET_BOOL_PROP(MixedDimension, MixedDimensionUseIterativeSolver, false);

//! default property value for the solution vector only used for monolithic solver
SET_PROP(MixedDimension, SolutionVector)
{
private:
    using BulkProblemTypeTag = typename GET_PROP_TYPE(TypeTag, BulkProblemTypeTag);
    using LowDimProblemTypeTag = typename GET_PROP_TYPE(TypeTag, LowDimProblemTypeTag);
public:
    using SolutionVectorBulk = typename GET_PROP_TYPE(BulkProblemTypeTag, SolutionVector);
    using SolutionVectorLowDim = typename GET_PROP_TYPE(LowDimProblemTypeTag, SolutionVector);
    using type = Dune::MultiTypeBlockVector<SolutionVectorBulk, SolutionVectorLowDim>;
};

//! Set the type of a global jacobian matrix from the solution types
SET_PROP(MixedDimension, JacobianMatrix)
{
private:
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using BulkProblemTypeTag = typename GET_PROP_TYPE(TypeTag, BulkProblemTypeTag);
    using LowDimProblemTypeTag = typename GET_PROP_TYPE(TypeTag, LowDimProblemTypeTag);
    enum {
        numEqBulk = GET_PROP_VALUE(BulkProblemTypeTag, NumEq),
        numEqLowDim = GET_PROP_VALUE(LowDimProblemTypeTag, NumEq)
    };

public:
    // the sub-blocks
    using MatrixLittleBlockBulk = Dune::FieldMatrix<Scalar, numEqBulk, numEqBulk>;
    using MatrixLittleBlockBulkCoupling = Dune::FieldMatrix<Scalar, numEqBulk, numEqLowDim>;
    using MatrixLittleBlockLowDim = Dune::FieldMatrix<Scalar, numEqLowDim, numEqLowDim>;
    using MatrixLittleBlockLowDimCoupling = Dune::FieldMatrix<Scalar, numEqLowDim, numEqBulk>;

    // the BCRS matrices of the subproblems as big blocks
    using MatrixBlockBulk = Dune::BCRSMatrix<MatrixLittleBlockBulk>;
    using MatrixBlockBulkCoupling = Dune::BCRSMatrix<MatrixLittleBlockBulkCoupling>;
    using MatrixBlockLowDim = Dune::BCRSMatrix<MatrixLittleBlockLowDim>;
    using MatrixBlockLowDimCoupling = Dune::BCRSMatrix<MatrixLittleBlockLowDimCoupling>;

    // the row types
    using RowBulk = Dune::MultiTypeBlockVector<MatrixBlockBulk, MatrixBlockBulkCoupling>;
    using RowLowDim = Dune::MultiTypeBlockVector<MatrixBlockLowDimCoupling, MatrixBlockLowDim>;

    // the jacobian matrix
    using type = Dune::MultiTypeBlockMatrix<RowBulk, RowLowDim>;
};

//! Definition of the indices of the subproblems in the global solution vector
SET_PROP(MixedDimension, SubProblemBlockIndices)
{
    using BulkIdx = Dune::index_constant<0>;
    using LowDimIdx = Dune::index_constant<1>;
};

//! set default solver
SET_TYPE_PROP(MixedDimension, LinearSolver, GSBiCGSTABBackend);

//! if the deflection of the newton method is large, we do not
//! need to solve the linear approximation accurately. Assuming
//! that the initial value for the delta vector u is quite
//! close to the final value, a reduction of 6 orders of
//! magnitude in the defect should be sufficient...
SET_SCALAR_PROP(MixedDimension, LinearSolverResidualReduction, 1e-6);

//! set the default number of maximum iterations for the linear solver
SET_INT_PROP(MixedDimension, LinearSolverMaxIterations, 250);

//! set number of equations of the mathematical model as default
SET_INT_PROP(MixedDimension, LinearSolverBlockSize, 1);

//! set number of maximum timestep divisions to 10
SET_INT_PROP(MixedDimension, ImplicitMaxTimeStepDivisions, 10);

}//end namespace Properties

}//end namespace Dumux

#endif
