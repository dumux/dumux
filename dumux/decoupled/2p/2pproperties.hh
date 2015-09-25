// $Id$
/*****************************************************************************
 *   Copyright (C) 2008-2009 by Markus Wolff                                 *
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
 * \ingroup IMPES
 * \ingroup Properties
 */
/*!
 * \file
 *
 * \brief Defines the properties required for (immiscible) twophase sequential models.
 */

#ifndef DUMUX_2PPROPERTIES_HH
#define DUMUX_2PPROPERTIES_HH

//Dune-includes
#include <dune/istl/operators.hh>
#include <dune/istl/solvers.hh>
#include <dune/istl/preconditioners.hh>

//Dumux-includes
#include <dumux/decoupled/common/impetproperties.hh>
#include <dumux/decoupled/2p/transport/transportproperties.hh>

namespace Dumux
{

////////////////////////////////
// forward declarations
////////////////////////////////

template<class TypeTag>
class VariableClass2P;

template<class TypeTag>
class FluidSystem2P;

template<class TypeTag>
class TwoPFluidState;

template<class TypeTag>
struct TwoPCommonIndices;

////////////////////////////////
// properties
////////////////////////////////
namespace Properties
{
//////////////////////////////////////////////////////////////////
// Type tags
//////////////////////////////////////////////////////////////////

//! The type tag for the two-phase problems
NEW_TYPE_TAG(DecoupledTwoP, INHERITS_FROM(IMPET, Transport))
;

//////////////////////////////////////////////////////////////////
// Property tags
//////////////////////////////////////////////////////////////////

NEW_PROP_TAG ( TwoPIndices )
;
NEW_PROP_TAG( SpatialParameters )
; //!< The type of the spatial parameters object
NEW_PROP_TAG( EnableGravity)
; //!< Returns whether gravity is considered in the problem
NEW_PROP_TAG( PressureFormulation)
; //!< The formulation of the model
NEW_PROP_TAG( SaturationFormulation)
; //!< The formulation of the model
NEW_PROP_TAG( VelocityFormulation)
; //!< The formulation of the model
NEW_PROP_TAG( EnableCompressibility)
;// !< Returns whether compressibility is allowed
NEW_PROP_TAG( WettingPhase)
; //!< The wetting phase for two-phase models
NEW_PROP_TAG( NonwettingPhase)
; //!< The non-wetting phase for two-phase models
NEW_PROP_TAG( FluidSystem )//!< Defines the fluid system
;
NEW_PROP_TAG( FluidState )//!< Defines the fluid state
;

//Properties for linear solvers
NEW_PROP_TAG(PressureCoefficientMatrix);//!< Type of the coefficient matrix given to the linear solver
NEW_PROP_TAG(PressureRHSVector);//!< Type of the right hand side vector given to the linear solver
NEW_PROP_TAG( PressurePreconditioner );//!< Type of the pressure preconditioner
NEW_PROP_TAG( PressureSolver );//!< Type of the pressure solver
NEW_PROP_TAG( SolverParameters );//!< Container for the solver parameters
NEW_PROP_TAG(ReductionSolver);//!< Reduction value for linear solver
NEW_PROP_TAG(MaxIterationNumberSolver);//!< Maximum iteration number for linear solver
NEW_PROP_TAG(IterationNumberPreconditioner);//!< Iteration number for preconditioner
NEW_PROP_TAG(VerboseLevelSolver);//!<Verbose level of linear solver and preconditioner
NEW_PROP_TAG(RelaxationPreconditioner);//!< Relaxation factor for preconditioner

//////////////////////////////////////////////////////////////////
// Properties
//////////////////////////////////////////////////////////////////

SET_PROP(DecoupledTwoP, NumPhases)
//!< The number of phases in the 2p model is 2
{
private:
typedef typename GET_PROP_TYPE(TypeTag, PTAG(FluidSystem)) FluidSystem;

public:
static const int value = FluidSystem::numPhases;
static_assert(value == 2,
        "Only fluid systems with 2 phases are supported by the 2p model!");
};
SET_INT_PROP(DecoupledTwoP, NumComponents, 1); //!< Each phase consists of 1 pure component


SET_PROP(DecoupledTwoP, TwoPIndices)
{
typedef TwoPCommonIndices<TypeTag> type;
};

//! Set the default formulation
SET_INT_PROP(DecoupledTwoP,
    PressureFormulation,
    TwoPCommonIndices<TypeTag>::pressureW);

SET_INT_PROP(DecoupledTwoP,
    SaturationFormulation,
    TwoPCommonIndices<TypeTag>::saturationW);

SET_INT_PROP(DecoupledTwoP,
    VelocityFormulation,
    TwoPCommonIndices<TypeTag>::velocityTotal);

SET_BOOL_PROP(DecoupledTwoP, EnableCompressibility, false);

SET_TYPE_PROP(DecoupledTwoP, Variables, VariableClass2P<TypeTag>);

SET_TYPE_PROP(DecoupledTwoP, FluidSystem, FluidSystem2P<TypeTag>);

SET_TYPE_PROP(DecoupledTwoP, FluidState, TwoPFluidState<TypeTag>);

SET_PROP(DecoupledTwoP, PressureCoefficientMatrix)
{
private:
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(Scalar)) Scalar;
    typedef Dune::FieldMatrix<Scalar, 1, 1> MB;

public:
    typedef Dune::BCRSMatrix<MB> type;
};
SET_PROP(DecoupledTwoP, PressureRHSVector)
{
private:
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(Scalar)) Scalar;

public:
    typedef Dune::BlockVector<Dune::FieldVector<Scalar, 1> > type;
};

SET_PROP(DecoupledTwoP, PressurePreconditioner)
{
private:
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(PressureCoefficientMatrix)) Matrix;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(PressureRHSVector)) Vector;
public:
    typedef Dune::SeqILUn<Matrix, Vector, Vector> type;
};
SET_PROP(DecoupledTwoP, PressureSolver)
{
private:
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(PressureCoefficientMatrix)) Matrix;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(PressureRHSVector)) Vector;
public:
    typedef Dune::BiCGSTABSolver<Vector> type;
//    typedef Dune::CGSolver<Vector> type;
};
SET_PROP_DEFAULT(LocalStiffness)
{
private:
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(GridView)) GridView;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(Scalar)) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(GridView)) Variables;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(Problem)) Problem;

public:
    typedef MimeticGroundwaterEquationLocalStiffness<GridView,Scalar,Variables, Problem> type;
};
//SET_INT_PROP(DecoupledTwoP, PressurePreconditioner, SolverIndices::seqILU0);
//SET_INT_PROP(DecoupledTwoP, PressureSolver, SolverIndices::biCGSTAB);
SET_SCALAR_PROP(DecoupledTwoP, ReductionSolver, 1E-12);
SET_INT_PROP(DecoupledTwoP, MaxIterationNumberSolver, 10000);
SET_INT_PROP(DecoupledTwoP, IterationNumberPreconditioner, 0);
SET_INT_PROP(DecoupledTwoP, VerboseLevelSolver, 1);
SET_SCALAR_PROP(DecoupledTwoP, RelaxationPreconditioner, 1.0);
SET_PROP(DecoupledTwoP, SolverParameters)
{
public:
    //solver parameters
    static const double reductionSolver = GET_PROP_VALUE(TypeTag, PTAG(ReductionSolver));
    static const int maxIterationNumberSolver = GET_PROP_VALUE(TypeTag, PTAG(MaxIterationNumberSolver));
    static const int iterationNumberPreconditioner = GET_PROP_VALUE(TypeTag, PTAG(IterationNumberPreconditioner));
    static const int verboseLevelSolver = GET_PROP_VALUE(TypeTag, PTAG(VerboseLevelSolver));
    static const double relaxationPreconditioner = GET_PROP_VALUE(TypeTag, PTAG(RelaxationPreconditioner));
};

// \}
}

/*!
 * \brief The common indices for the two-phase model.
 */
template <class TypeTag>
struct TwoPCommonIndices
{
private:
typedef typename GET_PROP_TYPE(TypeTag, PTAG(FluidSystem)) FluidSystem;

public:
// Formulations
//saturation flags
static const int saturationW = 0;
static const int saturationNW = 1;
//pressure flags
static const int pressureW = 0;
static const int pressureNW = 1;
static const int pressureGlobal = 2;
//velocity flags
static const int velocityW = 0;
static const int velocityNW = 1;
static const int velocityTotal = 2;

// Phase indices
static const int wPhaseIdx = FluidSystem::wPhaseIdx; //!< Index of the wetting phase in a phase vector
static const int nPhaseIdx = FluidSystem::nPhaseIdx; //!< Index of the non-wetting phase in a phase vector
};


}

#endif
