/*****************************************************************************
 *   Copyright (C) 2010 by Markus Wolff                                      *
 *   Copyright (C) 2008-2009 by Andreas Lauser                               *
 *   Copyright (C) 2008 by Bernd Flemisch                                    *
 *   Institute of Hydraulic Engineering                                      *
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
 * \ingroup OnePhase
 * \file
 *
 * \brief Defines the properties required for the single phase sequential model.
 */

#ifndef DUMUX_1PPROPERTIES_HH
#define DUMUX_1PPROPERTIES_HH

#define OnePModel

//Dune-includes
#include <dune/istl/operators.hh>
#include <dune/istl/solvers.hh>
#include <dune/istl/preconditioners.hh>

//Dumux-includes
#include <dumux/decoupled/common/decoupledproperties.hh>
#include <dumux/linear/seqsolverbackend.hh>

namespace Dumux
{

////////////////////////////////
// forward declarations
////////////////////////////////

template<class TypeTag>
class VariableClass;

////////////////////////////////
// properties
////////////////////////////////
namespace Properties
{
//////////////////////////////////////////////////////////////////
// Type tags
//////////////////////////////////////////////////////////////////

//! The type tag for the single-phase problem
NEW_TYPE_TAG(DecoupledOneP, INHERITS_FROM(DecoupledModel));

//////////////////////////////////////////////////////////////////
// Property tags
//////////////////////////////////////////////////////////////////

NEW_PROP_TAG( SpatialParameters )
; //!< The type of the spatial parameters object
NEW_PROP_TAG( EnableGravity)
; //!< Returns whether gravity is considered in the problem
NEW_PROP_TAG( Fluid )
; //!< The fluid for one-phase models

//Properties for linear solvers
NEW_PROP_TAG(PressureCoefficientMatrix);//!< Type of the coefficient matrix given to the linear solver
NEW_PROP_TAG(PressureRHSVector);//!< Type of the right hand side vector given to the linear solver
NEW_PROP_TAG( LinearSolver );//!< Type of linear solver

//////////////////////////////////////////////////////////////////
// Properties
//////////////////////////////////////////////////////////////////

SET_INT_PROP(DecoupledOneP, NumEq, 1);

SET_INT_PROP(DecoupledOneP, NumPhases, 1)//!< Single phase system
;
SET_INT_PROP(DecoupledOneP, NumComponents, 1); //!< Each phase consists of 1 pure component

SET_TYPE_PROP(DecoupledOneP, Variables, VariableClass<TypeTag>);

SET_PROP(DecoupledOneP, PressureCoefficientMatrix)
{
private:
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(Scalar)) Scalar;
    typedef Dune::FieldMatrix<Scalar, 1, 1> MB;

public:
    typedef Dune::BCRSMatrix<MB> type;
};
SET_PROP(DecoupledOneP, PressureRHSVector)
{
private:
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(Scalar)) Scalar;

public:
    typedef Dune::BlockVector<Dune::FieldVector<Scalar, 1> > type;
};

// use the stabilized BiCG solver preconditioned by the ILU-0 by default
SET_TYPE_PROP(DecoupledOneP, LinearSolver, Dumux::ILU0BiCGSTABBackend<TypeTag> );
}
}
#endif
