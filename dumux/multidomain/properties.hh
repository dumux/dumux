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
 * \ingroup MultiDomain
 * \brief Base properties for problems of mixed dimension
 * two sub problems
 */

#ifndef DUMUX_MULTIDOMAIN_PROPERTIES_HH
#define DUMUX_MULTIDOMAIN_PROPERTIES_HH

#include <dune/common/fvector.hh>
#include <dune/common/fmatrix.hh>
#include <dune/common/indices.hh>
#include <dune/istl/bcrsmatrix.hh>
#include <dune/istl/multitypeblockvector.hh>
#include <dune/istl/multitypeblockmatrix.hh>

#include <dumux/common/basicproperties.hh>
#include <dumux/linear/linearsolverproperties.hh>
#include <dumux/nonlinear/newtonmethod.hh>
#include <dumux/common/timemanager.hh>

#include <dumux/multidimension/subproblemlocaljacobian.hh>

namespace Dumux
{

// forward declarations
template <class TypeTag> class MultiDomainModel;
// template <class TypeTag> class StokesLocalJacobian;
// template <class TypeTag> class DarcyLocalJacobian;
template <class TypeTag> class MultiDomainNewtonController;
template <class TypeTag> class MultiDomainAssembler;

namespace Properties
{
// NumericModel provides Scalar, GridCreator, ParameterTree
NEW_TYPE_TAG(MultiDomain, INHERITS_FROM(NewtonMethod, LinearSolverTypeTag, NumericModel)); // TODO

NEW_PROP_TAG(Model); //!< The type of the base class of the model
NEW_PROP_TAG(StokesLocalJacobian); //!< The type of the Stokes local jacobian operator
NEW_PROP_TAG(DarcyLocalJacobian); //!< The type of the Darcy local jacobian operator

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
NEW_PROP_TAG(StokesProblemTypeTag);
NEW_PROP_TAG(DarcyProblemTypeTag);
NEW_PROP_TAG(Problem);
NEW_PROP_TAG(CouplingManager);
NEW_PROP_TAG(MultiDomainUseIterativeSolver);
NEW_PROP_TAG(SubProblemBlockIndices);

// forward declarations
NEW_PROP_TAG(NumEq);

//////////////////////////////////////////////////////////////////
// Properties
//////////////////////////////////////////////////////////////////
//! use the plain newton method by default
SET_TYPE_PROP(MultiDomain, NewtonMethod, NewtonMethod<TypeTag>);

//! set default values
SET_TYPE_PROP(MultiDomain, NewtonController, MultiDomainNewtonController<TypeTag>);

//! Set the assembler
SET_TYPE_PROP(MultiDomain, JacobianAssembler, MultiDomainAssembler<TypeTag>);

//! Set the BaseModel to MultiDomainModel
SET_TYPE_PROP(MultiDomain, Model, MultiDomainModel<TypeTag>);

//! The local jacobian operators
SET_TYPE_PROP(MultiDomain, StokesLocalJacobian, StokesLocalJacobian<TypeTag>);
SET_TYPE_PROP(MultiDomain, DarcyLocalJacobian, DarcyLocalJacobian<TypeTag>);

//! use forward differences to calculate the jacobian by default
SET_INT_PROP(MultiDomain, ImplicitNumericDifferenceMethod, +1);

//! default property value for the time manager
SET_TYPE_PROP(MultiDomain, TimeManager, TimeManager<TypeTag>);

//! default property is monolithic solver
SET_BOOL_PROP(MultiDomain, MultiDomainUseIterativeSolver, false);

//! default property value for the solution vector only used for monolithic solver
SET_PROP(MultiDomain, SolutionVector)
{
private:
    typedef typename GET_PROP_TYPE(TypeTag, StokesProblemTypeTag) StokesProblemTypeTag;
    typedef typename GET_PROP_TYPE(TypeTag, DarcyProblemTypeTag) DarcyProblemTypeTag;
public:
    typedef typename GET_PROP_TYPE(StokesProblemTypeTag, SolutionVector) SolutionVectorStokes;
    typedef typename GET_PROP_TYPE(DarcyProblemTypeTag, SolutionVector) SolutionVectorDarcy;
    typedef typename Dune::MultiTypeBlockVector<SolutionVectorStokes, SolutionVectorDarcy> type;
};

//! Set the type of a global jacobian matrix from the solution types
SET_PROP(MultiDomain, JacobianMatrix)
{
private:
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, StokesProblemTypeTag) StokesProblemTypeTag;
    typedef typename GET_PROP_TYPE(TypeTag, DarcyProblemTypeTag) DarcyProblemTypeTag;
    enum {
        numEqStokes = GET_PROP_VALUE(StokesProblemTypeTag, NumEq),
        numEqDarcy = GET_PROP_VALUE(DarcyProblemTypeTag, NumEq)
    };

public:
    // the sub-blocks
    typedef typename Dune::FieldMatrix<Scalar, numEqStokes, numEqStokes> MatrixLittleBlockStokes;
    typedef typename Dune::FieldMatrix<Scalar, numEqStokes, numEqDarcy> MatrixLittleBlockStokesCoupling;
    typedef typename Dune::FieldMatrix<Scalar, numEqDarcy, numEqDarcy> MatrixLittleBlockDarcy;
    typedef typename Dune::FieldMatrix<Scalar, numEqDarcy, numEqStokes> MatrixLittleBlockDarcyCoupling;

    // the BCRS matrices of the subproblems as big blocks
    typedef typename Dune::BCRSMatrix<MatrixLittleBlockStokes> MatrixBlockStokes;
    typedef typename Dune::BCRSMatrix<MatrixLittleBlockStokesCoupling> MatrixBlockStokesCoupling;
    typedef typename Dune::BCRSMatrix<MatrixLittleBlockDarcy> MatrixBlockDarcy;
    typedef typename Dune::BCRSMatrix<MatrixLittleBlockDarcyCoupling> MatrixBlockDarcyCoupling;

    // the row types
    typedef typename Dune::MultiTypeBlockVector<MatrixBlockStokes, MatrixBlockStokesCoupling> RowStokes;
    typedef typename Dune::MultiTypeBlockVector<MatrixBlockDarcyCoupling, MatrixBlockDarcy> RowDarcy;

    // the jacobian matrix
    typedef typename Dune::MultiTypeBlockMatrix<RowStokes, RowDarcy> type;
};

//! Definition of the indices of the subproblems in the global solution vector
SET_PROP(MultiDomain, SubProblemBlockIndices)
{
    using StokesIdx = Dune::index_constant<0>;
    using DarcyIdx = Dune::index_constant<1>;
};

//! set default solver
SET_TYPE_PROP(MultiDomain, LinearSolver, GSBiCGSTABBackend<TypeTag>);

//! set the block level to 2, suitable for e.g. the Dune::MultiTypeBlockMatrix
SET_INT_PROP(MultiDomain, LinearSolverPreconditionerBlockLevel, 2);

//! if the deflection of the newton method is large, we do not
//! need to solve the linear approximation accurately. Assuming
//! that the initial value for the delta vector u is quite
//! close to the final value, a reduction of 6 orders of
//! magnitude in the defect should be sufficient...
SET_SCALAR_PROP(MultiDomain, LinearSolverResidualReduction, 1e-6);

//! set the default number of maximum iterations for the linear solver
SET_INT_PROP(MultiDomain, LinearSolverMaxIterations, 250);

//! set number of equations of the mathematical model as default
SET_INT_PROP(MultiDomain, LinearSolverBlockSize, 1);

//! set number of maximum timestep divisions to 10
SET_INT_PROP(MultiDomain, ImplicitMaxTimeStepDivisions, 10);

}//end namespace Properties

}//end namespace Dumux

#endif
