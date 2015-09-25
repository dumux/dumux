// $Id: test_2p_injectionproblem.hh 3456 2010-04-09 12:11:51Z mwolff $
/*****************************************************************************
 *   Copyright (C) 2007-2008 by Klaus Mosthaf                                *
 *   Copyright (C) 2007-2008 by Bernd Flemisch                               *
 *   Copyright (C) 2008-2009 by Andreas Lauser                               *
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
//#include <dumux/material/fluidsystems/brine_co2_system.hh>
#include <dumux/material/fluidsystems/h2o_n2_system.hh>

#include <dumux/decoupled/2p2c/2p2cproblem.hh>
#include <dumux/decoupled/2p2c/fvpressure2p2cmultiphysics.hh>
#include <dumux/decoupled/2p2c/fvtransport2p2cmultiphysics.hh>

#include "test_dec2p2c_spatialparams.hh"

namespace Dumux
{

template<class TypeTag>
class TestMultTwoPTwoCProblem;

// Specify the properties
namespace Properties
{
NEW_TYPE_TAG(TestMultTwoPTwoCProblem, INHERITS_FROM(DecoupledTwoPTwoC/*, Transport*/));

// Set the grid type
SET_PROP(TestMultTwoPTwoCProblem, Grid)
{
    //    typedef Dune::YaspGrid<2> type;
    typedef Dune::SGrid<3, 3> type;
};

// Set the problem property
SET_PROP(TestMultTwoPTwoCProblem, Problem)
{
    typedef Dumux::TestMultTwoPTwoCProblem<TTAG(TestMultTwoPTwoCProblem)> type;
};

// Set the model properties
SET_PROP(TestMultTwoPTwoCProblem, TransportModel)
{
    typedef Dumux::FVTransport2P2CMultiPhysics<TTAG(TestMultTwoPTwoCProblem)> type;
};

SET_PROP(TestMultTwoPTwoCProblem, PressureModel)
{
    typedef Dumux::FVPressure2P2CMultiPhysics<TTAG(TestMultTwoPTwoCProblem)> type;
};

SET_INT_PROP(TestMultTwoPTwoCProblem, VelocityFormulation,
        GET_PROP_TYPE(TypeTag, PTAG(TwoPIndices))::velocityW);

SET_INT_PROP(TestMultTwoPTwoCProblem, PressureFormulation,
        GET_PROP_TYPE(TypeTag, PTAG(TwoPIndices))::pressureW);

//// Select fluid system
//SET_PROP(TestMultTwoPTwoCProblem, FluidSystem)
//{
//    typedef Dumux::Brine_CO2_System<TypeTag, Dumux::Benchmark3::CO2Tables> type;
//};
// Select fluid system
SET_PROP(TestMultTwoPTwoCProblem, FluidSystem)
{
    typedef Dumux::H2O_N2_System<TypeTag> type;
};

// Select water formulation
SET_PROP(TestMultTwoPTwoCProblem, Components) : public GET_PROP(TypeTag, PTAG(DefaultComponents))
{
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(Scalar)) Scalar;
//    typedef Dumux::TabulatedComponent<Scalar, typename Dumux::H2O<Scalar> > H20;
        typedef Dumux::H2O<Scalar> H2O;
};

// Set the soil properties
SET_PROP(TestMultTwoPTwoCProblem, SpatialParameters)
{
private:
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(Grid)) Grid;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(Scalar)) Scalar;

public:
    typedef Dumux::Test2P2CSpatialParams<TypeTag> type;
};

// Enable gravity
SET_BOOL_PROP(TestMultTwoPTwoCProblem, EnableGravity, true);
SET_INT_PROP(DecoupledTwoPTwoC,
        BoundaryMobility,
        TwoPCommonIndices<TypeTag>::satDependent);
SET_SCALAR_PROP(TestMultTwoPTwoCProblem, CFLFactor, 0.8);
}

/*!
 * \ingroup IMPETtests
 *
 * \brief test problem for the sequential 2p2c model with multiphysics
 *
 * The domain is box shaped (3D). All sides are closed (Neumann 0 boundary)
 * except the top and bottom boundaries (Dirichlet). A Gas (Nitrogen)
 * is injected over a vertical well in the center of the domain.
 *
 * A multiphysics approach is used to adapt model complexity (see
 * description in the pressure module)
 *
 * To run the simulation execute the following line in shell:
 * <tt>./test_dec2p2c</tt>
 * Optionally, simulation endtime and first timestep size can be
 * specified by programm arguments.
 */
template<class TypeTag = TTAG(TestMultTwoPTwoCProblem)>
class TestMultTwoPTwoCProblem: public IMPETProblem2P2C<TypeTag, TestMultTwoPTwoCProblem<TypeTag> >
{
typedef TestMultTwoPTwoCProblem<TypeTag> ThisType;
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
TestMultTwoPTwoCProblem(const GridView &gridView, const GlobalPosition lowerLeft = 0, const GlobalPosition upperRight = 0) :
ParentType(gridView), lowerLeft_(lowerLeft), upperRight_(upperRight)
{
    // initialize the tables of the fluid system
//    FluidSystem::init();
}

/*!
 * \name Problem parameters
 */
// \{

//! The problem name.
/*! This is used as a prefix for files generated by the simulation.
*/
const char *name() const
{
    return "test_multiphysics2p2c";
}
//!  Returns true if a restart file should be written.
/* The default behaviour is to write no restart file.
 */
bool shouldWriteRestartFile() const
{
    return false;
}

//! Returns the temperature within the domain.
/*! This problem assumes a temperature of 10 degrees Celsius.
 */
Scalar temperature(const GlobalPosition& globalPos, const Element& element) const
{
    return 273.15 + 10; // -> 10°C
}

// \}
/*!
 * \copydoc Dumux::TestDecTwoPTwoCProblem::referencePressure()
 */
Scalar referencePressure(const GlobalPosition& globalPos, const Element& element) const
{
    return 1e6;
}
/*!
 * \copydoc Dumux::TestDecTwoPTwoCProblem::bcTypePress()
 */
typename BoundaryConditions::Flags bcTypePress(const GlobalPosition& globalPos, const Intersection& intersection) const
{
    if (globalPos[0] > 10-1E-6 || globalPos[0] < 1e-6)
        return BoundaryConditions::dirichlet;
    // all other boundaries
    return BoundaryConditions::neumann;
}
/*!
 * \copydoc Dumux::TestDecTwoPTwoCProblem::bcTypeTransport()
 */
typename BoundaryConditions::Flags bcTypeTransport(const GlobalPosition& globalPos, const Intersection& intersection) const
{
    return bcTypePress(globalPos, intersection);
}
/*!
 * \copydoc Dumux::TestDecTwoPTwoCProblem::bcFormulation()
 */
BoundaryConditions2p2c::Flags bcFormulation(const GlobalPosition& globalPos, const Intersection& intersection) const
{
    return BoundaryConditions2p2c::concentration;
}
/*!
 * \copydoc Dumux::TestDecTwoPTwoCProblem::dirichletPress()
 */
Scalar dirichletPress(const GlobalPosition& globalPos, const Intersection& intersection) const
{
    const Element& element = *(intersection.inside());

    Scalar pRef = referencePressure(globalPos, element);
    Scalar temp = temperature(globalPos, element);

    // all other boundaries
    return (globalPos[0] < 1e-6) ? (2.5e5 - FluidSystem::H2O::liquidDensity(temp, pRef) * this->gravity()[dim-1])
            : (2e5 - FluidSystem::H2O::liquidDensity(temp, pRef) * this->gravity()[dim-1]);
}
/*!
 * \copydoc Dumux::TestDecTwoPTwoCProblem::dirichletTransport()
 */
Scalar dirichletTransport(const GlobalPosition& globalPos, const Intersection& intersection) const
{
    return 1.;
}
/*!
 * \copydoc Dumux::TestDecTwoPTwoCProblem::neumann()
 */
Dune::FieldVector<Scalar,2> neumann(const GlobalPosition& globalPos, const Intersection& intersection) const
{
    Dune::FieldVector<Scalar,2> neumannFlux(0.0);
//    if (globalPos[1] < 15 && globalPos[1]> 5)
//    {
//        neumannFlux[nPhaseIdx] = -0.015;
//    }
    return neumannFlux;
}
/*!
 * \copydoc Dumux::TestDecTwoPTwoCProblem::source()
 */
Dune::FieldVector<Scalar,2> source(const GlobalPosition& globalPos, const Element& element)
{
    Dune::FieldVector<Scalar,2> q_(0);
        if (fabs(globalPos[0] - 4.5) < 1 && fabs(globalPos[1] - 4.5) < 1) q_[1] = 0.0001;
    return q_;
}
/*!
 * \copydoc Dumux::TestDecTwoPTwoCProblem::initFormulation()
 */
const BoundaryConditions2p2c::Flags initFormulation (const GlobalPosition& globalPos, const Element& element) const
{
    return BoundaryConditions2p2c::concentration;
}
/*!
 * \copydoc Dumux::TestDecTwoPTwoCProblem::initSat()
 */
Scalar initSat(const GlobalPosition& globalPos, const Element& element) const
{
    return 1;
}
/*!
 * \copydoc Dumux::TestDecTwoPTwoCProblem::initConcentration()
 */
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
