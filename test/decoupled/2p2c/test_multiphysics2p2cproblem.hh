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
#include <dumux/material/old_fluidsystems/h2o_n2_system.hh>

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
NEW_TYPE_TAG(TestMultTwoPTwoCProblem, INHERITS_FROM(DecoupledTwoPTwoC));

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

SET_INT_PROP(TestMultTwoPTwoCProblem, PressureFormulation,
        GET_PROP_TYPE(TypeTag, PTAG(TwoPTwoCIndices))::pressureNW);

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
SET_BOOL_PROP(TestMultTwoPTwoCProblem, EnableCapillarity, true);
SET_INT_PROP(TestMultTwoPTwoCProblem,
             BoundaryMobility,
             GET_PROP_TYPE(TypeTag, PTAG(TwoPTwoCIndices))::satDependent);
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
 * <tt>./test_multiphyiscs2p2c</tt>
 * Optionally, simulation endtime and first timestep size can be
 * specified by programm arguments.
 */
template<class TypeTag = TTAG(TestMultTwoPTwoCProblem)>
class TestMultTwoPTwoCProblem: public IMPETProblem2P2C<TypeTag>
{
typedef IMPETProblem2P2C<TypeTag> ParentType;
typedef typename GET_PROP_TYPE(TypeTag, PTAG(GridView)) GridView;
typedef typename GET_PROP_TYPE(TypeTag, PTAG(TimeManager)) TimeManager;
typedef typename GET_PROP_TYPE(TypeTag, PTAG(TwoPTwoCIndices)) Indices;

typedef typename GET_PROP_TYPE(TypeTag, PTAG(FluidSystem)) FluidSystem;
typedef typename GET_PROP_TYPE(TypeTag, PTAG(FluidState)) FluidState;

// boundary typedefs
typedef typename GET_PROP_TYPE(TypeTag, PTAG(BoundaryTypes)) BoundaryTypes;
typedef typename GET_PROP_TYPE(TypeTag, PTAG(PrimaryVariables)) PrimaryVariables;

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
TestMultTwoPTwoCProblem(TimeManager &timeManager, const GridView &gridView, const GlobalPosition lowerLeft = 0, const GlobalPosition upperRight = 0) :
ParentType(timeManager, gridView), lowerLeft_(lowerLeft), upperRight_(upperRight)
{
    // Specifies how many time-steps are done before output will be written.
    this->setOutputInterval(10);

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
Scalar temperatureAtPos(const GlobalPosition& globalPos) const
{
    return 273.15 + 10; // -> 10°C
}

// \}
/*!
 * \copydoc Dumux::TestDecTwoPTwoCProblem::referencePressureAtPos()
 */
Scalar referencePressureAtPos(const GlobalPosition& globalPos) const
{
    return 1e6;
}
/*!
 * \copydoc Dumux::TestDecTwoPTwoCProblem::boundaryTypesAtPos()
 */
void boundaryTypesAtPos(BoundaryTypes &bcTypes, const GlobalPosition& globalPos) const
{
    if (globalPos[0] > 10-1E-6 || globalPos[0] < 1e-6)
        bcTypes.setAllDirichlet();
    else
        // all other boundaries
        bcTypes.setAllNeumann();
}

/*!
 * \copydoc Dumux::TestDecTwoPTwoCProblem::boundaryFormulation()
 */
const void boundaryFormulation(typename Indices::BoundaryFormulation &bcFormulation, const Intersection& intersection) const
{
    bcFormulation = Indices::concentration;
}
/*!
 * \copydoc Dumux::TestDecTwoPTwoCProblem::dirichletAtPos()
 */
void dirichletAtPos(PrimaryVariables &bcValues ,const GlobalPosition& globalPos) const
{
    Scalar pRef = referencePressureAtPos(globalPos);
    Scalar temp = temperatureAtPos(globalPos);

    // Dirichlet for pressure equation
    bcValues[Indices::pressureEqIdx] = (globalPos[0] < 1e-6) ? (2.5e5 - FluidSystem::H2O::liquidDensity(temp, pRef) * this->gravity()[dim-1])
            : (2e5 - FluidSystem::H2O::liquidDensity(temp, pRef) * this->gravity()[dim-1]);

    // Dirichlet values for transport equations
    bcValues[Indices::contiWEqIdx] = 1.;
    bcValues[Indices::contiNEqIdx] = 1.- bcValues[Indices::contiWEqIdx];

}
/*!
 * \copydoc Dumux::TestDecTwoPTwoCProblem::neumannAtPos()
 */
void neumannAtPos(PrimaryVariables &neumannValues, const GlobalPosition& globalPos) const
{
    neumannValues[Indices::contiNEqIdx] = 0.;
    neumannValues[Indices::contiWEqIdx] = 0.;
//    if (globalPos[1] < 15 && globalPos[1]> 5)
//    {
//        neumannValues[Indices::contiNEqIdx] = -0.015;
//    }
}
/*!
 * \copydoc Dumux::TestDecTwoPTwoCProblem::sourceAtPos()
 */
void sourceAtPos(PrimaryVariables &sourceValues, const GlobalPosition& globalPos) const
{
    setZero(sourceValues);
    if (fabs(globalPos[0] - 4.5) < 1 && fabs(globalPos[1] - 4.5) < 1)
        sourceValues[Indices::contiNEqIdx] = 0.0001;
}
/*!
 * \copydoc Dumux::TestDecTwoPTwoCProblem::initialFormulation()
 */
const void initialFormulation(typename Indices::BoundaryFormulation &initialFormulation, const Element& element) const
{
    initialFormulation = Indices::concentration;
}
/*!
 * \copydoc Dumux::TestDecTwoPTwoCProblem::initConcentrationAtPos()
 */
Scalar initConcentrationAtPos(const GlobalPosition& globalPos) const
{
    return 1;
}

private:
GlobalPosition lowerLeft_;
GlobalPosition upperRight_;

static constexpr Scalar eps_ = 1e-6;
static constexpr Scalar depthBOR_ = 1000;
};
} //end namespace

#endif
