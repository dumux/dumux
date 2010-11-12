// $Id: test_2p_injectionproblem.hh 3456 2010-04-09 12:11:51Z mwolff $
/*****************************************************************************
 *   Copyright (C) 2007-2008 by Klaus Mosthaf                                *
 *   Copyright (C) 2007-2008 by Bernd Flemisch                               *
 *   Copyright (C) 2008-2009 by Andreas Lauser                               *
 *   Institute of Hydraulic Engineering                                      *
 *   University of Stuttgart, Germany                                        *
 *   email: <givenname>.<name>@iws.uni-stuttgart.de                          *
 *                                                                           *
 *   This program is free software; you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation; either version 2 of the License, or       *
 *   (at your option) any later version, as long as this copyright notice    *
 *   is included in its original form.                                       *
 *                                                                           *
 *   This program is distributed WITHOUT ANY WARRANTY.                       *
 *****************************************************************************/
/*!
 * \file
 *
 * \brief test problem for the multiphysics 2p2c model
 */
#ifndef DUMUX_TEST_2P2C_PROBLEM_HH
#define DUMUX_TEST_2P2C_PROBLEM_HH

#if HAVE_UG
#include <dune/grid/uggrid.hh>
#endif

#include <dune/grid/yaspgrid.hh>
#include <dune/grid/sgrid.hh>

// fluid properties
//#include <dumux/material/fluidsystems/simple_h2o_n2_system.hh>
#include <dumux/material/fluidsystems/h2o_n2_system.hh>

#include <dumux/decoupled/2p2c/2p2cproblem.hh>
#include <dumux/decoupled/2p2c/fvpressure2p2cmultiphysics.hh>
#include <dumux/decoupled/2p2c/fvtransport2p2cmultiphysics.hh>

#include "test_dec2p2c_spatialparams.hh"

namespace Dumux
{

template<class TypeTag>
class TestTwoPTwoCProblem;

//////////
// Specify the properties for the lens problem
//////////
namespace Properties
{
NEW_TYPE_TAG(TestTwoPTwoCProblem, INHERITS_FROM(DecoupledTwoPTwoC/*, Transport*/));

// Set the grid type
SET_PROP(TestTwoPTwoCProblem, Grid)
{
    //    typedef Dune::YaspGrid<2> type;
    typedef Dune::SGrid<3, 3> type;
};

// Set the problem property
SET_PROP(TestTwoPTwoCProblem, Problem)
{
    typedef Dumux::TestTwoPTwoCProblem<TTAG(TestTwoPTwoCProblem)> type;
};

// Set the model properties
SET_PROP(TestTwoPTwoCProblem, TransportModel)
{
    typedef Dumux::FVTransport2P2CMultiPhysics<TTAG(TestTwoPTwoCProblem)> type;
};

SET_PROP(TestTwoPTwoCProblem, PressureModel)
{
    typedef Dumux::FVPressure2P2CMultiPhysics<TTAG(TestTwoPTwoCProblem)> type;
};

SET_INT_PROP(TestTwoPTwoCProblem, VelocityFormulation,
        GET_PROP_TYPE(TypeTag, PTAG(TwoPIndices))::velocityW);

SET_INT_PROP(TestTwoPTwoCProblem, PressureFormulation,
        GET_PROP_TYPE(TypeTag, PTAG(TwoPIndices))::pressureW);


// Select fluid system
SET_PROP(TestTwoPTwoCProblem, FluidSystem)
{
    typedef Dumux::H2O_N2_System<TypeTag> type;
};

// Select water formulation
SET_PROP(TestTwoPTwoCProblem, Components) : public GET_PROP(TypeTag, PTAG(DefaultComponents))
{
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(Scalar)) Scalar;
//    typedef Dumux::TabulatedComponent<Scalar, typename Dumux::H2O<Scalar> > type;
        typedef Dumux::H2O<Scalar> H2O;
//        static void init(){};
};

// Set the soil properties
SET_PROP(TestTwoPTwoCProblem, SpatialParameters)
{
private:
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(Grid)) Grid;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(Scalar)) Scalar;

public:
    typedef Dumux::Test2P2CSpatialParams<TypeTag> type;
};

// Enable gravity
SET_BOOL_PROP(TestTwoPTwoCProblem, EnableGravity, true);
SET_INT_PROP(DecoupledTwoPTwoC,
        BoundaryMobility,
        TwoPCommonIndices<TypeTag>::satDependent);
SET_SCALAR_PROP(TestTwoPTwoCProblem, CFLFactor, 0.8);
}

/*!
 * \ingroup DecoupledProblems
 */
template<class TypeTag = TTAG(TestTwoPTwoCProblem)>
class TestTwoPTwoCProblem: public IMPETProblem2P2C<TypeTag, TestTwoPTwoCProblem<TypeTag> >
{
typedef TestTwoPTwoCProblem<TypeTag> ThisType;
typedef IMPETProblem2P2C<TypeTag, ThisType> ParentType;
typedef typename GET_PROP_TYPE(TypeTag, PTAG(GridView)) GridView;

typedef typename GET_PROP_TYPE(TypeTag, PTAG(TwoPIndices)) Indices;

typedef typename GET_PROP_TYPE(TypeTag, PTAG(FluidSystem)) FluidSystem;
typedef typename GET_PROP_TYPE(TypeTag, PTAG(FluidState)) FluidState;

enum
{
    dim = GridView::dimension, dimWorld = GridView::dimensionworld
};

enum
{
    wPhaseIdx = Indices::wPhaseIdx, nPhaseIdx = Indices::nPhaseIdx
};

typedef typename GET_PROP_TYPE(TypeTag, PTAG(Scalar)) Scalar;

typedef typename GridView::Traits::template Codim<0>::Entity Element;
typedef typename GridView::Intersection Intersection;
typedef Dune::FieldVector<Scalar, dimWorld> GlobalPosition;
typedef Dune::FieldVector<Scalar, dim> LocalPosition;

public:
TestTwoPTwoCProblem(const GridView &gridView, const GlobalPosition lowerLeft = 0, const GlobalPosition upperRight = 0) :
ParentType(gridView), lowerLeft_(lowerLeft), upperRight_(upperRight)
{
    // initialize the tables of the fluid system
//    WaterFormulation::init(273.15, 623.15, 100,
//                            -10,   20e6, 200);
    FluidSystem::init();
}

/*!
 * \name Problem parameters
 */
// \{

/*!
 * \brief The problem name.
 *
 * This is used as a prefix for files generated by the simulation.
 */
const char *name() const
{
    return "test_multiphysics2p2c";
}

bool shouldWriteRestartFile() const
{
    return false;
}

/*!
 * \brief Returns the temperature within the domain.
 *
 * This problem assumes a temperature of 10 degrees Celsius.
 */
Scalar temperature(const GlobalPosition& globalPos, const Element& element) const
{
    return 273.15 + 10; // -> 10°C
}

// \}

Scalar referencePressure(const GlobalPosition& globalPos, const Element& element) const
{
    return 1e6; // TODO: proper description!!
}

typename BoundaryConditions::Flags bcTypePress(const GlobalPosition& globalPos, const Intersection& intersection) const
{
    if (globalPos[0] > 10-1E-6 || globalPos[0] < 1e-6)
        return BoundaryConditions::dirichlet;
    // all other boundaries
    return BoundaryConditions::neumann;
}

typename BoundaryConditions::Flags bcTypeTransport(const GlobalPosition& globalPos, const Intersection& intersection) const
{
    return bcTypePress(globalPos, intersection);
}

BoundaryConditions2p2c::Flags bcFormulation(const GlobalPosition& globalPos, const Intersection& intersection) const
{
	return BoundaryConditions2p2c::concentration;
}

Scalar dirichletPress(const GlobalPosition& globalPos, const Intersection& intersection) const
{
    const Element& element = *(intersection.inside());

    Scalar pRef = referencePressure(globalPos, element);
    Scalar temp = temperature(globalPos, element);

    // all other boundaries
    return (globalPos[0] < 1e-6) ? (2.5e5 - FluidSystem::H2O::liquidDensity(temp, pRef) * this->gravity()[dim-1])
            : (2e5 - FluidSystem::H2O::liquidDensity(temp, pRef) * this->gravity()[dim-1]);
}

Scalar dirichletTransport(const GlobalPosition& globalPos, const Intersection& intersection) const
{
    return 1.;
}

Dune::FieldVector<Scalar,2> neumann(const GlobalPosition& globalPos, const Intersection& intersection) const
{
	Dune::FieldVector<Scalar,2> neumannFlux(0.0);
//    if (globalPos[1] < 15 && globalPos[1]> 5)
//    {
//        neumannFlux[nPhaseIdx] = -0.015;
//    }
    return neumannFlux;
}

Dune::FieldVector<Scalar,2> source(const GlobalPosition& globalPos, const Element& element)
{
    Dune::FieldVector<Scalar,2> q_(0);
        if (fabs(globalPos[0] - 4.5) < 1 && fabs(globalPos[1] - 4.5) < 1) q_[1] = 0.0001;
    return q_;
}

const BoundaryConditions2p2c::Flags initFormulation (const GlobalPosition& globalPos, const Element& element) const
{
    return BoundaryConditions2p2c::concentration;
}

Scalar initSat(const GlobalPosition& globalPos, const Element& element) const
{
    return 1;
}
Scalar initConcentration(const GlobalPosition& globalPos, const Element& element) const
{
    return 1;
}

private:
GlobalPosition lowerLeft_;
GlobalPosition upperRight_;

static const Scalar eps_ = 1e-6;
static const Scalar depthBOR_ = 1000;
};
} //end namespace

#endif
