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
 * \brief test problem for the sequential 2p2c model
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
#include <dumux/decoupled/2p2c/fvpressure2p2c.hh>
#include <dumux/decoupled/2p2c/fvtransport2p2c.hh>

#include "test_dec2p2c_spatialparams.hh"

namespace Dumux
{

template<class TypeTag>
class TestDecTwoPTwoCProblem;

// Specify the properties
namespace Properties
{
NEW_TYPE_TAG(TestDecTwoPTwoCProblem, INHERITS_FROM(DecoupledTwoPTwoC/*, Transport*/));

// Set the grid type
SET_PROP(TestDecTwoPTwoCProblem, Grid)
{
    //    typedef Dune::YaspGrid<2> type;
    typedef Dune::SGrid<3, 3> type;
};

// Set the problem property
SET_PROP(TestDecTwoPTwoCProblem, Problem)
{
    typedef Dumux::TestDecTwoPTwoCProblem<TTAG(TestDecTwoPTwoCProblem)> type;
};

// Set the model properties
SET_PROP(TestDecTwoPTwoCProblem, TransportModel)
{
    typedef Dumux::FVTransport2P2C<TTAG(TestDecTwoPTwoCProblem)> type;
};

SET_PROP(TestDecTwoPTwoCProblem, PressureModel)
{
    typedef Dumux::FVPressure2P2C<TTAG(TestDecTwoPTwoCProblem)> type;
};

SET_INT_PROP(TestDecTwoPTwoCProblem, VelocityFormulation,
        GET_PROP_TYPE(TypeTag, PTAG(TwoPIndices))::velocityW);

SET_INT_PROP(TestDecTwoPTwoCProblem, PressureFormulation,
        GET_PROP_TYPE(TypeTag, PTAG(TwoPIndices))::pressureW);


// Select fluid system
SET_PROP(TestDecTwoPTwoCProblem, FluidSystem)
{
    typedef Dumux::H2O_N2_System<TypeTag> type;
};

// Select water formulation
SET_PROP(TestDecTwoPTwoCProblem, Components) : public GET_PROP(TypeTag, PTAG(DefaultComponents))
{
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(Scalar)) Scalar;
//    typedef Dumux::TabulatedComponent<Scalar, typename Dumux::H2O<Scalar> > H20;
        typedef Dumux::SimpleH2O<Scalar> H2O;
};

// Set the soil properties
SET_PROP(TestDecTwoPTwoCProblem, SpatialParameters)
{
private:
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(Grid)) Grid;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(Scalar)) Scalar;

public:
    typedef Dumux::Test2P2CSpatialParams<TypeTag> type;
};

// Enable gravity
SET_BOOL_PROP(TestDecTwoPTwoCProblem, EnableGravity, true);
SET_INT_PROP(DecoupledTwoPTwoC,
        BoundaryMobility,
        TwoPCommonIndices<TypeTag>::satDependent);
SET_SCALAR_PROP(TestDecTwoPTwoCProblem, CFLFactor, 0.8);
}

/*!
 * \ingroup IMPETtests
 *
 * \brief test problem for the sequential 2p2c model
 *
 * The domain is box shaped (3D). All sides are closed (Neumann 0 boundary)
 * except the top and bottom boundaries (Dirichlet). A Gas (Nitrogen)
 * is injected over a vertical well in the center of the domain.
 *
 * To run the simulation execute the following line in shell:
 * <tt>./test_dec2p2c</tt>
 * Optionally, simulation endtime and first timestep size can be
 * specified by programm arguments.
 */
template<class TypeTag = TTAG(TestDecTwoPTwoCProblem)>
class TestDecTwoPTwoCProblem: public IMPETProblem2P2C<TypeTag, TestDecTwoPTwoCProblem<TypeTag> >
{
typedef TestDecTwoPTwoCProblem<TypeTag> ThisType;
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
TestDecTwoPTwoCProblem(const GridView &gridView, const GlobalPosition lowerLeft = 0, const GlobalPosition upperRight = 0) :
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
    return "test_dec2p2c";
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
//! Returns the reference pressure.
 /*This pressure is used in order to calculate the material properties
 * at the beginning of the initialization routine. It should lie within
 * a reasonable pressure range for the current problem.
 */
Scalar referencePressure(const GlobalPosition& globalPos, const Element& element) const
{
    return 1e6;
}
//! Type of pressure boundary condition.
/*! Defines the type the boundary condition for the pressure equation,
 *  either pressure (dirichlet) or flux (neumann).
 */
typename BoundaryConditions::Flags bcTypePress(const GlobalPosition& globalPos, const Intersection& intersection) const
{
    if (globalPos[0] > 10-1E-6 || globalPos[0] < 1e-6)
        return BoundaryConditions::dirichlet;
    // all other boundaries
    return BoundaryConditions::neumann;
}
//! Type of Transport boundary condition.
/*! Defines the type the boundary condition for the transport equation,
 *  either saturation (dirichlet) or flux (neumann).
 */
typename BoundaryConditions::Flags bcTypeTransport(const GlobalPosition& globalPos, const Intersection& intersection) const
{
    return bcTypePress(globalPos, intersection);
}
//! Flag for the type of Dirichlet conditions
/*! The Dirichlet BCs can be specified by a given concentration (mass of
 * a component per total mass inside the control volume) or by means
 * of a saturation.
 */
BoundaryConditions2p2c::Flags bcFormulation(const GlobalPosition& globalPos, const Intersection& intersection) const
{
    return BoundaryConditions2p2c::concentration;
}
//! Value for dirichlet pressure boundary condition \f$ [Pa] \f$.
/*! In case of a dirichlet BC for the pressure equation, the pressure
 *  have to be defined on boundaries.
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
//! Value for transport dirichlet boundary condition (dimensionless).
/*! In case of a dirichlet BC for the transport equation,
 *  has to be defined on boundaries.
 */
Scalar dirichletTransport(const GlobalPosition& globalPos, const Intersection& intersection) const
{
    return 1.;
}
//! Value for pressure neumann boundary condition \f$ [\frac{kg}{m^3 \cdot s}] \f$.
/*! In case of a neumann boundary condition, the flux of matter
 *  is returned as a vector.
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
//! Source of mass \f$ [\frac{kg}{m^3 \cdot s}] \f$
/*! Evaluate the source term for all phases within a given
 *  volume. The method returns the mass generated (positive) or
 *  annihilated (negative) per volume unit.
 */
Dune::FieldVector<Scalar,2> source(const GlobalPosition& globalPos, const Element& element)
{
    Dune::FieldVector<Scalar,2> q_(0);
        if (fabs(globalPos[0] - 4.5) < 1 && fabs(globalPos[1] - 4.5) < 1) q_[1] = 0.0001;
    return q_;
}
//! Flag for the type of initial conditions
/*! The problem can be initialized by a given concentration (mass of
 * a component per total mass inside the control volume) or by means
 * of a saturation.
 */
const BoundaryConditions2p2c::Flags initFormulation (const GlobalPosition& globalPos, const Element& element) const
{
    return BoundaryConditions2p2c::concentration;
}
//! Saturation initial condition (dimensionless)
/*! The problem is initialized with the following saturation. Both
 * phases are assumed to contain an equilibrium concentration of the
 * correspondingly other component.
 */
Scalar initSat(const GlobalPosition& globalPos, const Element& element) const
{
    return 1;
}
//! Concentration initial condition (dimensionless)
/*! The problem is initialized with the following concentration.
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
