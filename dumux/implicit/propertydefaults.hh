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
 * \ingroup ImplicitProperties
 * \file
 * \brief Default properties for implicit models
 */
#ifndef DUMUX_IMPLICIT_PROPERTY_DEFAULTS_HH
#define DUMUX_IMPLICIT_PROPERTY_DEFAULTS_HH

#include <dune/common/fvector.hh>
#include <dune/common/fmatrix.hh>
#include <dune/grid/common/mcmgmapper.hh>
#include <dune/istl/bcrsmatrix.hh>

#include <dumux/nonlinear/newtonmethod.hh>
#include <dumux/nonlinear/newtoncontroller.hh>
#include <dumux/common/boundarytypes.hh>
#include <dumux/common/timemanager.hh>
#include <dumux/linear/amgbackend.hh>

#include "properties.hh"
#include "model.hh"
#include "localjacobian.hh"
#include "volumevariables.hh"

namespace Dumux {

// forward declarations
template <class TypeTag> class NewtonController;
template <class TypeTag> class AMGBackend;

namespace Properties {
//////////////////////////////////////////////////////////////////
// Some defaults for very fundamental properties
//////////////////////////////////////////////////////////////////

//! Set the default type for the time manager
SET_TYPE_PROP(ImplicitBase, TimeManager, TimeManager<TypeTag>);

//////////////////////////////////////////////////////////////////
// Properties
//////////////////////////////////////////////////////////////////

//! Use the leaf grid view if not defined otherwise
SET_TYPE_PROP(ImplicitBase,
              GridView,
              typename GET_PROP_TYPE(TypeTag, Grid)::LeafGridView);

//! use the plain newton method by default
SET_TYPE_PROP(ImplicitBase, NewtonMethod, NewtonMethod<TypeTag>);

//! use the plain newton controller by default
SET_TYPE_PROP(ImplicitBase, NewtonController, NewtonController<TypeTag>);

//! Mapper for the grid view's vertices.
SET_TYPE_PROP(ImplicitBase,
              VertexMapper,
              Dune::MultipleCodimMultipleGeomTypeMapper<typename GET_PROP_TYPE(TypeTag, GridView),
                                                        Dune::MCMGVertexLayout>);

//! Mapper for the grid view's elements.
SET_TYPE_PROP(ImplicitBase,
              ElementMapper,
              Dune::MultipleCodimMultipleGeomTypeMapper<typename GET_PROP_TYPE(TypeTag, GridView),
                                                        Dune::MCMGElementLayout>);

//! Set the BaseModel to ImplicitModel
SET_TYPE_PROP(ImplicitBase, BaseModel, ImplicitModel<TypeTag>);

//! The volume variable class, to be overloaded by the model
SET_TYPE_PROP(ImplicitBase, VolumeVariables, ImplicitVolumeVariables<TypeTag>);

//! The local jacobian operator
SET_TYPE_PROP(ImplicitBase, LocalJacobian, ImplicitLocalJacobian<TypeTag>);

//! The type of a solution for the whole grid at a fixed time
SET_TYPE_PROP(ImplicitBase,
              SolutionVector,
              Dune::BlockVector<typename GET_PROP_TYPE(TypeTag, PrimaryVariables)>);

//! The type of a solution for a whole element
SET_TYPE_PROP(ImplicitBase,
              ElementSolutionVector,
              Dune::BlockVector<typename GET_PROP_TYPE(TypeTag, PrimaryVariables)>);

//! A vector of primary variables
SET_TYPE_PROP(ImplicitBase,
              PrimaryVariables,
              Dune::FieldVector<typename GET_PROP_TYPE(TypeTag, Scalar),
                                GET_PROP_VALUE(TypeTag, NumEq)>);

//! Boundary types at a single degree of freedom
SET_TYPE_PROP(ImplicitBase,
              BoundaryTypes,
              BoundaryTypes<GET_PROP_VALUE(TypeTag, NumEq)>);

//! use forward differences to calculate the jacobian by default
SET_INT_PROP(ImplicitBase, ImplicitNumericDifferenceMethod, +1);

//! do not use hints by default
SET_BOOL_PROP(ImplicitBase, ImplicitEnableHints, false);

//! disable jacobian matrix recycling by default
SET_BOOL_PROP(ImplicitBase, ImplicitEnableJacobianRecycling, false);

//! disable partial reassembling by default
SET_BOOL_PROP(ImplicitBase, ImplicitEnablePartialReassemble, false);

//! Set the type of a global jacobian matrix from the solution types
SET_PROP(ImplicitBase, JacobianMatrix)
{
private:
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    enum { numEq = GET_PROP_VALUE(TypeTag, NumEq) };
    typedef typename Dune::FieldMatrix<Scalar, numEq, numEq> MatrixBlock;
public:
    typedef typename Dune::BCRSMatrix<MatrixBlock> type;
};

//! use the AMGBackend solver by default
SET_TYPE_PROP(ImplicitBase, LinearSolver, Dumux::AMGBackend<TypeTag> );

// if the deflection of the newton method is large, we do not
// need to solve the linear approximation accurately. Assuming
// that the initial value for the delta vector u is quite
// close to the final value, a reduction of 6 orders of
// magnitude in the defect should be sufficient...
SET_SCALAR_PROP(ImplicitBase, LinearSolverResidualReduction, 1e-6);

//! set the default number of maximum iterations for the linear solver
SET_INT_PROP(ImplicitBase, LinearSolverMaxIterations, 250);

//! set number of equations of the mathematical model as default
SET_INT_PROP(ImplicitBase, LinearSolverBlockSize, GET_PROP_VALUE(TypeTag, NumEq));

//! set number of maximum timestep divisions to 10
SET_INT_PROP(ImplicitBase, ImplicitMaxTimeStepDivisions, 10);

} // namespace Properties
} // namespace Dumux

#endif
