// $Id$
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

//Dune-includes
#include <dune/istl/operators.hh>
#include <dune/istl/solvers.hh>
#include <dune/istl/preconditioners.hh>

//Dumux-includes
#include <dumux/decoupled/common/decoupledproperties.hh>

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
NEW_PROP_TAG(PressureCoefficientMatrix)
;//!< Type of the coefficient matrix given to the linear solver
NEW_PROP_TAG(PressureRHSVector)
;//!< Type of the right hand side vector given to the linear solver
NEW_PROP_TAG( PressurePreconditioner )
;//!< Type of the preconditioner
NEW_PROP_TAG( PressureSolver )
;//!< Type of the solver
NEW_PROP_TAG( SolverParameters )
;//!< Container for the solver parameters
NEW_PROP_TAG(ReductionSolver)
;//!< Reduction value for linear solver
NEW_PROP_TAG(MaxIterationNumberSolver)
;//!< Maximum iteration number for linear solver
NEW_PROP_TAG(IterationNumberPreconditioner)
;//!< Iteration number for preconditioner
NEW_PROP_TAG(VerboseLevelSolver)
;//!<Verbose level of linear solver and preconditioner
NEW_PROP_TAG(RelaxationPreconditioner)
;//!< Relaxation factor for preconditioner

//////////////////////////////////////////////////////////////////
// Properties
//////////////////////////////////////////////////////////////////

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

SET_PROP(DecoupledOneP, PressurePreconditioner)
{
private:
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(PressureCoefficientMatrix)) Matrix;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(PressureRHSVector)) Vector;
public:
    typedef Dune::SeqILUn<Matrix, Vector, Vector> type;
};
SET_PROP(DecoupledOneP, PressureSolver)
{
private:
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(PressureCoefficientMatrix)) Matrix;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(PressureRHSVector)) Vector;
public:
    typedef Dune::BiCGSTABSolver<Vector> type;
};

SET_SCALAR_PROP(DecoupledOneP, ReductionSolver, 1E-12);
SET_INT_PROP(DecoupledOneP, MaxIterationNumberSolver, 10000);
SET_INT_PROP(DecoupledOneP, IterationNumberPreconditioner, 0);
SET_INT_PROP(DecoupledOneP, VerboseLevelSolver, 1);
SET_SCALAR_PROP(DecoupledOneP, RelaxationPreconditioner, 1.0);

SET_PROP(DecoupledOneP, SolverParameters)
{
public:
    //solver parameters
    static const double reductionSolver = GET_PROP_VALUE(TypeTag, PTAG(ReductionSolver));
    static const int maxIterationNumberSolver = GET_PROP_VALUE(TypeTag, PTAG(MaxIterationNumberSolver));
    static const int iterationNumberPreconditioner = GET_PROP_VALUE(TypeTag, PTAG(IterationNumberPreconditioner));
    static const int verboseLevelSolver = GET_PROP_VALUE(TypeTag, PTAG(VerboseLevelSolver));
    static const double relaxationPreconditioner = GET_PROP_VALUE(TypeTag, PTAG(RelaxationPreconditioner));
};
}
}
#endif
