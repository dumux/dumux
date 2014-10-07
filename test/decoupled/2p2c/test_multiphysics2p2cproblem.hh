// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   See the file COPYING for full copying permissions.                      *
 *                                                                           *
 *   This program is free software: you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation, either version 2 of the License, or       *
 *   (at your option) any later version.                                     *
 *                                                                           *
 *   This program is distributed in the hope that it will be useful,         *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of          *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the            *
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

#include <dune/grid/yaspgrid.hh>

#include <dumux/decoupled/2p2c/2p2cproblem.hh>
#include <dumux/decoupled/2p2c/fvpressure2p2cmultiphysics.hh>
#include <dumux/decoupled/2p2c/fvtransport2p2cmultiphysics.hh>
#include <dumux/decoupled/2p2c/celldata2p2cmultiphysics.hh>
// fluid properties
//#include <dumux/material/fluidsystems/brine_co2_system.hh>
#include <dumux/material/fluidsystems/h2on2fluidsystem.hh>
#include <dumux/material/fluidsystems/h2oairfluidsystem.hh>

#include "test_dec2p2c_spatialparams.hh"

namespace Dumux
{

template<class TypeTag>
class TestMultTwoPTwoCProblem;

// Specify the properties
namespace Properties
{
NEW_TYPE_TAG(TestMultTwoPTwoCProblem, INHERITS_FROM(DecoupledTwoPTwoC, Test2P2CSpatialParams));

SET_TYPE_PROP(TestMultTwoPTwoCProblem, CellData, Dumux::CellData2P2Cmultiphysics<TypeTag>);

// Set the grid type
SET_PROP(TestMultTwoPTwoCProblem, Grid)
{
    typedef Dune::YaspGrid<3> type;
};

// Set the problem property
SET_PROP(TestMultTwoPTwoCProblem, Problem)
{
    typedef Dumux::TestMultTwoPTwoCProblem<TypeTag> type;
};

// Set the model properties
SET_PROP(TestMultTwoPTwoCProblem, TransportModel)
{
    typedef Dumux::FVTransport2P2CMultiPhysics<TypeTag> type;
};

SET_PROP(TestMultTwoPTwoCProblem, PressureModel)
{
    typedef Dumux::FVPressure2P2CMultiPhysics<TypeTag> type;
};

SET_INT_PROP(TestMultTwoPTwoCProblem, PressureFormulation,
        GET_PROP_TYPE(TypeTag, Indices)::pressureNw);


//// Select fluid system
//SET_PROP(TestMultTwoPTwoCProblem, FluidSystem)
//{
//    typedef Dumux::Brine_CO2_System<TypeTag, Dumux::Benchmark3::CO2Tables> type;
//};
// Select fluid system
SET_PROP(TestMultTwoPTwoCProblem, FluidSystem)
{
    typedef Dumux::H2OAirFluidSystem<TypeTag> type;
};
// Select fluid system
SET_BOOL_PROP(TestMultTwoPTwoCProblem, EnableComplicatedFluidSystem, true);

// Select water formulation
SET_PROP(TestMultTwoPTwoCProblem, Components) : public GET_PROP(TypeTag, DefaultComponents)
{
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
//    typedef Dumux::TabulatedComponent<Scalar, typename Dumux::H2O<Scalar> > H20;
        typedef Dumux::H2O<Scalar> H2O;
};

// Enable gravity
SET_BOOL_PROP(TestMultTwoPTwoCProblem, ProblemEnableGravity, true);
SET_BOOL_PROP(TestMultTwoPTwoCProblem, EnableCapillarity, true);
SET_INT_PROP(TestMultTwoPTwoCProblem,
             BoundaryMobility,
             GET_PROP_TYPE(TypeTag, Indices)::satDependent);
SET_SCALAR_PROP(TestMultTwoPTwoCProblem, ImpetCFLFactor, 0.8);
//SET_SCALAR_PROP(TestMultTwoPTwoCProblem, ImpetSubCFLFactor, 0.8);//can be defined to use sub-time-stepping for the transport
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
template<class TypeTag>
class TestMultTwoPTwoCProblem: public IMPETProblem2P2C<TypeTag>
{
typedef IMPETProblem2P2C<TypeTag> ParentType;
typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
typedef typename GET_PROP_TYPE(TypeTag, TimeManager) TimeManager;
typedef typename GET_PROP_TYPE(TypeTag, Indices) Indices;

typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;

// boundary typedefs
typedef typename GET_PROP_TYPE(TypeTag, BoundaryTypes) BoundaryTypes;
typedef typename GET_PROP_TYPE(TypeTag, PrimaryVariables) PrimaryVariables;

enum
{
    dim = GridView::dimension, dimWorld = GridView::dimensionworld
};

typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;

typedef typename GridView::Traits::template Codim<0>::Entity Element;
typedef typename GridView::Intersection Intersection;
typedef Dune::FieldVector<Scalar, dimWorld> GlobalPosition;

public:
TestMultTwoPTwoCProblem(TimeManager &timeManager, const GridView &gridView, const GlobalPosition upperRight = 0) :
ParentType(timeManager, gridView), lowerLeft_(0), upperRight_(upperRight), eps_(1e-6), depthBOR_(1000.0)
{
    // Specifies how many time-steps are done before output will be written.
//    this->setOutputInterval(20);

//    // initialize the tables of the fluid system
//    FluidSystem::init(/*tempMin=*/280,
//            /*tempMax=*/290,
//            /*numTemp=*/10,
//            /*pMin=*/190000,
//            /*pMax=*/280000,
//            /*numP=*/400);
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
 * \param globalPos The global Position
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
    this->setZero(sourceValues);
    if (fabs(globalPos[0] - 4.8) < 0.5 && fabs(globalPos[1] - 4.8) < 0.5)
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

const Scalar eps_;
const Scalar depthBOR_;
};
} //end namespace

#endif
