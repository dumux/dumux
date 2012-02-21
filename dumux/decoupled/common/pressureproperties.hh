// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   Copyright (C) 2011 by Markus Wolff                                      *
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
#ifndef DUMUX_PRESSURE_PROPERTIES_HH
#define DUMUX_PRESSURE_PROPERTIES_HH

//Dune-includes
#include <dune/istl/bcrsmatrix.hh>
#include <dune/istl/bvector.hh>

#include "decoupledproperties.hh"
#include <dumux/linear/linearsolverproperties.hh>
#include <dumux/linear/seqsolverbackend.hh>


/*!
 * \ingroup Sequential
 * \ingroup IMPETProperties
 */
/*!
 * \file
 * \brief Base file for properties related to sequential IMPET algorithms
 */
namespace Dumux
{

template<class TypeTag>
class IMPET;

namespace Properties
{
/*!
 *
 * \brief General properties for sequential IMPET algorithms
 *
 * This class holds properties necessary for the sequential IMPET solution.
 */

//////////////////////////////////////////////////////////////////
// Type tags tags
//////////////////////////////////////////////////////////////////

//! The type tag for models based on the diffusion-scheme
NEW_TYPE_TAG(Pressure, INHERITS_FROM(LinearSolverTypeTag, DecoupledModel));

//////////////////////////////////////////////////////////////////
// Property tags
//////////////////////////////////////////////////////////////////
//Properties for linear solvers
NEW_PROP_TAG(PressureCoefficientMatrix);//!< Type of the coefficient matrix given to the linear solver
NEW_PROP_TAG(PressureRHSVector);//!< Type of the right hand side vector given to the linear solver
NEW_PROP_TAG( VisitFacesOnlyOnce); //!< Indicates if faces are only regarded from one side
NEW_PROP_TAG(Velocity);
}
}

#include <dumux/decoupled/common/fv/velocitydefault.hh>

namespace Dumux
{
namespace Properties
{
//! Faces are only regarded from one side and not from both cells
SET_BOOL_PROP(Pressure, VisitFacesOnlyOnce, false);

//Set defaults
SET_PROP(Pressure, PressureCoefficientMatrix)
{
private:
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef Dune::FieldMatrix<Scalar, 1, 1> MB;

public:
    typedef Dune::BCRSMatrix<MB> type;
};
SET_PROP(Pressure, PressureRHSVector)
{
private:
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;

public:
    typedef Dune::BlockVector<Dune::FieldVector<Scalar, 1> > type;
};

// use the stabilized BiCG solver preconditioned by the ILU-0 by default
SET_TYPE_PROP(Pressure, LinearSolver, Dumux::ILU0BiCGSTABBackend<TypeTag> );

//! set the default for the reduction of the initial residual
SET_SCALAR_PROP(Pressure, LinearSolverResidualReduction, 1e-13);

//! set the default number of maximum iterations for the linear solver
SET_INT_PROP(Pressure, LinearSolverMaxIterations, 500);

//! set the default number of maximum iterations for the linear solver
SET_INT_PROP(Pressure, LinearSolverBlockSize, 1);

SET_TYPE_PROP( Pressure, Velocity, FVVelocityDefault<TypeTag>);

}
}

#endif
