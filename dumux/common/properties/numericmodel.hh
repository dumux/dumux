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
 * \ingroup Properties
 * \file
 *
 * \brief Defines a type tags and some fundamental properties for the numeric model
 */
#ifndef DUMUX_NUMERIC_MODEL_PROPERTIES_HH
#define DUMUX_NUMERIC_MODEL_PROPERTIES_HH

#include <dune/common/fvector.hh>
#include <dune/istl/bcrsmatrix.hh>
#include <dumux/common/properties.hh>
#include <dumux/common/properties/basic.hh>
#include <dumux/common/boundarytypes.hh>


namespace Dumux
{

namespace Properties
{
//! Type tag for numeric models.
NEW_TYPE_TAG(NumericModel, INHERITS_FROM(BasicProperties));

//! Set the default vector with size number of equations to a field vector
SET_TYPE_PROP(NumericModel, NumEqVector, Dune::FieldVector<typename GET_PROP_TYPE(TypeTag, Scalar), GET_PROP_VALUE(TypeTag, NumEq)>);

//! Set the default primary variable vector to a vector of size of number of equations
SET_TYPE_PROP(NumericModel, PrimaryVariables, typename GET_PROP_TYPE(TypeTag, NumEqVector));

//! The type of a solution for the whole grid at a fixed time
SET_TYPE_PROP(NumericModel, SolutionVector, Dune::BlockVector<typename GET_PROP_TYPE(TypeTag, PrimaryVariables)>);

//! Set the type of a global jacobian matrix from the solution types
SET_PROP(NumericModel, JacobianMatrix)
{
private:
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    enum { numEq = GET_PROP_VALUE(TypeTag, NumEq) };
    using MatrixBlock = typename Dune::FieldMatrix<Scalar, numEq, numEq>;
public:
    using type = typename Dune::BCRSMatrix<MatrixBlock>;
};

//! set the block level to 1, suitable for e.g. a simple Dune::BCRSMatrix.
// Set this to more than one if the matrix to solve is nested multiple times
// e.g. for Dune::MultiTypeBlockMatrix'es.
SET_INT_PROP(NumericModel, LinearSolverPreconditionerBlockLevel, 1);

//! set the block size to number of equations as default
/*!
 * The number of different types of equations which build the system of equations to solve
 * can differ from the number of equations given by the mathematical/physical model (e.g. IMPES).
 * Thus, the block size does not have to be equal to NumEq.
 * (Especially important for the SuperLU solver!)
 */
SET_INT_PROP(NumericModel, LinearSolverBlockSize, GET_PROP_VALUE(TypeTag, NumEq));

//! Boundary types at a single degree of freedom
SET_TYPE_PROP(NumericModel, BoundaryTypes, BoundaryTypes<GET_PROP_VALUE(TypeTag, NumEq)>);

} // namespace Properties
} // namespace Dumux

#endif
