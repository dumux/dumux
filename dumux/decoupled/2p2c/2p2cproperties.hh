// $Id:
/*****************************************************************************
 *   Copyright (C) 20010 by Benjamin Faigle                                  *
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
 * \ingroup IMPEC Properties
 * \file
 *
 * \brief Defines the properties required for the decoupled 2p2c models.
 */
#ifndef DUMUX_2P2CPROPERTIES_HH
#define DUMUX_2P2CPROPERTIES_HH

//Dune-includes
#include <dune/istl/operators.hh>
#include <dune/istl/solvers.hh>
#include <dune/istl/preconditioners.hh>

//DUMUX includes
#include <dumux/decoupled/common/impetproperties.hh>

namespace Dumux
{

////////////////////////////////
// forward declarations
////////////////////////////////


template<class TypeTag>
class VariableClass2P2C;

template<class TypeTag>
class DecTwoPTwoCFluidState;

template <class TypeTag>
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
NEW_TYPE_TAG(DecoupledTwoPTwoC, INHERITS_FROM(IMPET));

//////////////////////////////////////////////////////////////////
// Property tags
//////////////////////////////////////////////////////////////////

NEW_PROP_TAG ( TwoPIndices );
NEW_PROP_TAG( SpatialParameters ); //!< The type of the soil properties object
NEW_PROP_TAG( EnableGravity); //!< Returns whether gravity is considered in the problem
NEW_PROP_TAG( PressureFormulation); //!< The formulation of the model
NEW_PROP_TAG( SaturationFormulation); //!< The formulation of the model
NEW_PROP_TAG( VelocityFormulation); //!< The formulation of the model
NEW_PROP_TAG( EnableCompressibility);// ! Returns whether compressibility is allowed
NEW_PROP_TAG(EnableCapillarity); //!< Returns whether capillarity is regarded
NEW_PROP_TAG( BoundaryMobility );
NEW_PROP_TAG( NumDensityTransport );
NEW_PROP_TAG( FluidSystem );
NEW_PROP_TAG( FluidState );

//Properties for linear solvers
NEW_PROP_TAG(PressureCoefficientMatrix);
NEW_PROP_TAG(PressureRHSVector);
NEW_PROP_TAG( PressurePreconditioner );
NEW_PROP_TAG( PressureSolver );
NEW_PROP_TAG( SolverParameters );
NEW_PROP_TAG(ReductionSolver);
NEW_PROP_TAG(MaxIterationNumberSolver);
NEW_PROP_TAG(IterationNumberPreconditioner);
NEW_PROP_TAG(VerboseLevelSolver);
NEW_PROP_TAG(RelaxationPreconditioner);

//////////////////////////////////////////////////////////////////
// Properties
//////////////////////////////////////////////////////////////////
SET_PROP(DecoupledTwoPTwoC, TwoPIndices)
{
  typedef TwoPCommonIndices<TypeTag> type;
};

// set fluid/component information
SET_PROP(DecoupledTwoPTwoC, NumPhases) //!< The number of phases in the 2p model is 2
{
    // the property is created in decoupledproperties.hh
private:
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(FluidSystem)) FluidSystem;

public:
    static const int value = FluidSystem::numPhases;
    static_assert(value == 2,
                  "Only fluid systems with 2 phases are supported by the 2p2c model!");
};

SET_PROP(DecoupledTwoPTwoC, NumComponents) //!< The number of components in the 2p2c model is 2
{
private:
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(FluidSystem)) FluidSystem;

public:
    static const int value = FluidSystem::numComponents;
    static_assert(value == 2,
                  "Only fluid systems with 2 components are supported by the 2p2c model!");
};

//! Set the default formulation
SET_INT_PROP(DecoupledTwoPTwoC,
        PressureFormulation,
        TwoPCommonIndices<TypeTag>::pressureW);

SET_INT_PROP(DecoupledTwoPTwoC,
        SaturationFormulation,
        TwoPCommonIndices<TypeTag>::saturationW);

SET_INT_PROP(DecoupledTwoPTwoC,
        VelocityFormulation,
        TwoPCommonIndices<TypeTag>::velocityW);

SET_PROP(DecoupledTwoPTwoC, TransportSolutionType)
{
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(Scalar)) Scalar;
    typedef typename Dune::BlockVector<Dune::BlockVector<Dune::FieldVector<Scalar,1> > > type;//!<type for vector of vector (of scalars)

};

SET_BOOL_PROP(DecoupledTwoPTwoC, EnableCompressibility, true);

SET_BOOL_PROP(DecoupledTwoPTwoC, EnableCapillarity, false);


SET_PROP_DEFAULT(BoundaryMobility)
{
    static const int value = TwoPCommonIndices<TypeTag>::satDependent;
};
SET_PROP_DEFAULT(NumDensityTransport)
{
    static const bool value = false;
};

SET_TYPE_PROP(DecoupledTwoPTwoC, Variables, VariableClass2P2C<TypeTag>);

SET_TYPE_PROP(DecoupledTwoPTwoC, FluidState, DecTwoPTwoCFluidState<TypeTag>);

// solver stuff
SET_PROP(DecoupledTwoPTwoC, PressureCoefficientMatrix)
{
private:
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(Scalar)) Scalar;
    typedef Dune::FieldMatrix<Scalar, 1, 1> MB;

public:
    typedef Dune::BCRSMatrix<MB> type;
};
SET_PROP(DecoupledTwoPTwoC, PressureRHSVector)
{
private:
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(Scalar)) Scalar;

public:
    typedef Dune::BlockVector<Dune::FieldVector<Scalar, 1> > type;
};

SET_PROP(DecoupledTwoPTwoC, PressurePreconditioner)
{
private:
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(PressureCoefficientMatrix)) Matrix;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(PressureRHSVector)) Vector;
public:
    typedef Dune::SeqILUn<Matrix, Vector, Vector> type;
};
SET_PROP(DecoupledTwoPTwoC, PressureSolver)
{
private:
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(PressureCoefficientMatrix)) Matrix;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(PressureRHSVector)) Vector;
public:
    typedef Dune::BiCGSTABSolver<Vector> type;
};
//SET_INT_PROP(DecoupledTwoPTwoC, PressurePreconditioner, SolverIndices::seqILU0);
//SET_INT_PROP(DecoupledTwoPTwoC, PressureSolver, SolverIndices::biCGSTAB);
SET_SCALAR_PROP(DecoupledTwoPTwoC, ReductionSolver, 1E-12);
SET_INT_PROP(DecoupledTwoPTwoC, MaxIterationNumberSolver, 10000);
SET_INT_PROP(DecoupledTwoPTwoC, IterationNumberPreconditioner, 0);
SET_INT_PROP(DecoupledTwoPTwoC, VerboseLevelSolver, 1);
SET_SCALAR_PROP(DecoupledTwoPTwoC, RelaxationPreconditioner, 1.0);
SET_PROP(DecoupledTwoPTwoC, SolverParameters)
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

/*!
 * \brief The common indices for the 2p2c are the same as for the two-phase model.
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

    // BoundaryCondition flags
    static const int satDependent = 0;
    static const int permDependent = 1;
};

// \}

}

#endif
