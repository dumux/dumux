// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   Copyright (C) 2008-2010 by Andreas Lauser                               *
 *   Copyright (C) 2008-2010 by Bernd Flemisch                               *
 *   Institute for Modelling Hydraulic and Environmental Systems             *
 *   University of Stuttgart, Germany                                        *
 *   email: <givenname>.<name>@iws.uni-stuttgart.de                          *
 *                                                                           *
 *   This program is free software: you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation, either version 2 of the License, or       *
 *   (at your option) any later version.                                     *
 *                                                                           *
 *   This program is distributed in the hope that it will be useful,         *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of          *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the           *
 *   GNU General Public License for more details.                            *
 *                                                                           *
 *   You should have received a copy of the GNU General Public License       *
 *   along with this program.  If not, see <http://www.gnu.org/licenses/>.   *
 *****************************************************************************/
/*!
 * \ingroup Properties
 * \ingroup CCProperties
 * \ingroup CCModel
 * \file
 *
 * \brief Default properties for box models
 */
#ifndef DUMUX_CC_PROPERTY_DEFAULTS_HH
#define DUMUX_CC_PROPERTY_DEFAULTS_HH

#include <dumux/nonlinear/newtonmethod.hh>
#include <dumux/nonlinear/newtoncontroller.hh>

#include "ccassembler.hh"
#include "ccmodel.hh"
#include "ccfvelementgeometry.hh"
#include "ccelementboundarytypes.hh"
#include "cclocaljacobian.hh"
#include "cclocalresidual.hh"
#include "ccelementvolumevariables.hh"
#include <dumux/boxmodels/common/boxvolumevariables.hh>

#include <dumux/common/boundarytypes.hh>
#include <dumux/common/timemanager.hh>

#include "ccproperties.hh"

#include <limits>

namespace Dumux {

// forward declaration
template<class TypeTag>
class CCModel;

namespace Properties {
//////////////////////////////////////////////////////////////////
// Some defaults for very fundamental properties
//////////////////////////////////////////////////////////////////

//! Set the default type for the time manager
SET_TYPE_PROP(CCModel, TimeManager, Dumux::TimeManager<TypeTag>);

//////////////////////////////////////////////////////////////////
// Properties
//////////////////////////////////////////////////////////////////

//! Use the leaf grid view if not defined otherwise
SET_TYPE_PROP(CCModel,
              GridView,
              typename GET_PROP_TYPE(TypeTag, Grid)::LeafGridView);

//! Set the default for the FVElementGeometry
SET_TYPE_PROP(CCModel, FVElementGeometry, Dumux::CCFVElementGeometry<TypeTag>);

//! Set the default for the ElementBoundaryTypes
SET_TYPE_PROP(CCModel, ElementBoundaryTypes, Dumux::CCElementBoundaryTypes<TypeTag>);

//! use the plain newton method for the box scheme by default
SET_TYPE_PROP(CCModel, NewtonMethod, Dumux::NewtonMethod<TypeTag>);

//! use the plain newton controller for the box scheme by default
SET_TYPE_PROP(CCModel, NewtonController, Dumux::NewtonController<TypeTag>);

//! Mapper for the grid view's vertices.
SET_TYPE_PROP(CCModel,
              VertexMapper,
              Dune::MultipleCodimMultipleGeomTypeMapper<typename GET_PROP_TYPE(TypeTag, GridView),
                                                        Dune::MCMGVertexLayout>);

//! Mapper for the grid view's elements.
SET_TYPE_PROP(CCModel,
              ElementMapper,
              Dune::MultipleCodimMultipleGeomTypeMapper<typename GET_PROP_TYPE(TypeTag, GridView),
                                                        Dune::MCMGElementLayout>);

//! Mapper for the degrees of freedoms.
SET_TYPE_PROP(CCModel, DofMapper, typename GET_PROP_TYPE(TypeTag, ElementMapper));

//! Set the BaseLocalResidual to CCLocalResidual
SET_TYPE_PROP(CCModel, BaseLocalResidual, Dumux::CCLocalResidual<TypeTag>);

//! Set the BaseModel to CCModel
SET_TYPE_PROP(CCModel, BaseModel, Dumux::CCModel<TypeTag>);

//! The local jacobian operator for the box scheme
SET_TYPE_PROP(CCModel, LocalJacobian, Dumux::CCLocalJacobian<TypeTag>);

/*!
 * \brief The type of a solution for the whole grid at a fixed time.
 */
SET_TYPE_PROP(CCModel,
              SolutionVector,
              Dune::BlockVector<typename GET_PROP_TYPE(TypeTag, PrimaryVariables)>);

/*!
 * \brief The type of a solution for a whole element.
 */
SET_TYPE_PROP(CCModel,
              ElementSolutionVector,
              Dune::BlockVector<typename GET_PROP_TYPE(TypeTag, PrimaryVariables)>);

/*!
 * \brief A vector of primary variables.
 */
SET_TYPE_PROP(CCModel,
              PrimaryVariables,
              Dune::FieldVector<typename GET_PROP_TYPE(TypeTag, Scalar),
                                GET_PROP_VALUE(TypeTag, NumEq)>);

/*!
 * \brief The volume variable class.
 *
 * This should almost certainly be overloaded by the model...
 */
SET_TYPE_PROP(CCModel, VolumeVariables, Dumux::BoxVolumeVariables<TypeTag>);

/*!
 * \brief An array of secondary variable containers.
 */
SET_TYPE_PROP(CCModel, ElementVolumeVariables, Dumux::CCElementVolumeVariables<TypeTag>);

/*!
 * \brief Boundary types at a single degree of freedom.
 */
SET_TYPE_PROP(CCModel,
              BoundaryTypes,
              Dumux::BoundaryTypes<GET_PROP_VALUE(TypeTag, NumEq)>);

/*!
 * \brief Assembler for the global jacobian matrix.
 */
SET_TYPE_PROP(CCModel, JacobianAssembler, Dumux::CCAssembler<TypeTag>);

//! use an unlimited time step size by default
#if 0
// requires GCC 4.6 and above to call the constexpr function of
// numeric_limits
SET_SCALAR_PROP(CCModel, TimeManagerMaxTimeStepSize, std::numeric_limits<Scalar>::infinity());
#else
SET_SCALAR_PROP(CCModel, TimeManagerMaxTimeStepSize, 1e100);
#endif

//! use forward differences to calculate the jacobian by default
SET_INT_PROP(CCModel, ImplicitNumericDifferenceMethod, +1);

//! do not use hints by default
SET_BOOL_PROP(CCModel, ImplicitEnableHints, false);

// disable jacobian matrix recycling by default
SET_BOOL_PROP(CCModel, ImplicitEnableJacobianRecycling, false);

// disable partial reassembling by default
SET_BOOL_PROP(CCModel, ImplicitEnablePartialReassemble, false);

//! Set the type of a global jacobian matrix from the solution types
SET_PROP(CCModel, JacobianMatrix)
{
private:
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    enum { numEq = GET_PROP_VALUE(TypeTag, NumEq) };
    typedef typename Dune::FieldMatrix<Scalar, numEq, numEq> MatrixBlock;
public:
    typedef typename Dune::BCRSMatrix<MatrixBlock> type;
};

// use the stabilized BiCG solver preconditioned by the ILU-0 by default
SET_TYPE_PROP(CCModel, LinearSolver, Dumux::BoxBiCGStabILU0Solver<TypeTag> );

// if the deflection of the newton method is large, we do not
// need to solve the linear approximation accurately. Assuming
// that the initial value for the delta vector u is quite
// close to the final value, a reduction of 6 orders of
// magnitude in the defect should be sufficient...
SET_SCALAR_PROP(CCModel, LinearSolverResidualReduction, 1e-6);

//! set the default number of maximum iterations for the linear solver
SET_INT_PROP(CCModel, LinearSolverMaxIterations, 250);

//! set number of equations of the mathematical model as default
SET_INT_PROP(CCModel, LinearSolverBlockSize, GET_PROP_VALUE(TypeTag, NumEq));

} // namespace Properties
} // namespace Dumux

#endif
