// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   See the file COPYING for full copying permissions.                      *
 *                                                                           *
 *   This program is free software: you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation, either version 3 of the License, or       *
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
 * \ingroup SequentialTwoPTwoCTests
 * \brief test problem for the multiphysics 2p2c model
 */
#ifndef DUMUX_TEST_2P2C_PROBLEM_HH
#define DUMUX_TEST_2P2C_PROBLEM_HH

#include <dune/grid/yaspgrid.hh>

#include <dumux/porousmediumflow/2p2c/sequential/problem.hh>
#include <dumux/porousmediumflow/2p2c/sequential/fvpressuremultiphysics.hh>
#include <dumux/porousmediumflow/2p2c/sequential/fvtransportmultiphysics.hh>
#include <dumux/porousmediumflow/2p2c/sequential/celldatamultiphysics.hh>

#include <dumux/material/fluidsystems/h2on2.hh>
#include <dumux/material/fluidsystems/h2oair.hh>

#include "test_dec2p2c_spatialparams.hh"

namespace Dumux
{
/*!
 * \ingroup SequentialTwoPTwoCTests
 * \brief test problem for the multiphysics 2p2c model
 */
template<class TypeTag>
class TestMultTwoPTwoCProblem;

// Specify the properties
namespace Properties
{
// Create new type tags
namespace TTag {
struct TestMultTwoPTwoC { using InheritsFrom = std::tuple<Test2P2CSpatialParams, SequentialTwoPTwoC>; };
} // end namespace TTag

template<class TypeTag>
struct CellData<TypeTag, TTag::TestMultTwoPTwoC> { using type = CellData2P2CMultiPhysics<TypeTag>; };

// Set the grid type
template<class TypeTag>
struct Grid<TypeTag, TTag::TestMultTwoPTwoC> { using type = Dune::YaspGrid<3>; };

// Set the problem property
template<class TypeTag>
struct Problem<TypeTag, TTag::TestMultTwoPTwoC> { using type = TestMultTwoPTwoCProblem<TypeTag>; };

// Set the model properties
template<class TypeTag>
struct TransportModel<TypeTag, TTag::TestMultTwoPTwoC> { using type = FVTransport2P2CMultiPhysics<TypeTag>; };

template<class TypeTag>
struct PressureModel<TypeTag, TTag::TestMultTwoPTwoC> { using type = FVPressure2P2CMultiPhysics<TypeTag>; };

template<class TypeTag>
struct PressureFormulation<TypeTag, TTag::TestMultTwoPTwoC> { static constexpr int value = GetPropType<TypeTag, Properties::Indices>::pressureNw; };

// Select fluid system
template<class TypeTag>
struct FluidSystem<TypeTag, TTag::TestMultTwoPTwoC>
{
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using type = FluidSystems::H2OAir<Scalar, Components::H2O<Scalar>>;
};

template<class TypeTag>
struct EnableCapillarity<TypeTag, TTag::TestMultTwoPTwoC> { static constexpr bool value = true; };
template<class TypeTag>
struct BoundaryMobility<TypeTag, TTag::TestMultTwoPTwoC> { static constexpr int value = GetPropType<TypeTag, Properties::Indices>::satDependent; };
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
using ParentType = IMPETProblem2P2C<TypeTag>;
using GridView = typename GetPropType<TypeTag, Properties::GridGeometry>::GridView;
using Grid = typename GridView::Grid;
using TimeManager = GetPropType<TypeTag, Properties::TimeManager>;
using Indices = typename GetPropType<TypeTag, Properties::ModelTraits>::Indices;

using FluidSystem = GetPropType<TypeTag, Properties::FluidSystem>;

using BoundaryTypes = GetPropType<TypeTag, Properties::SequentialBoundaryTypes>;
using PrimaryVariables = GetPropType<TypeTag, Properties::PrimaryVariables>;

enum
{
    dim = GridView::dimension, dimWorld = GridView::dimensionworld
};

using Scalar = GetPropType<TypeTag, Properties::Scalar>;

using Element = typename GridView::Traits::template Codim<0>::Entity;
using Intersection = typename GridView::Intersection;
using GlobalPosition = typename Element::Geometry::GlobalCoordinate;

public:
TestMultTwoPTwoCProblem(TimeManager &timeManager, Grid &grid, const GlobalPosition& upperRight = 0) :
ParentType(timeManager, grid), lowerLeft_(0), upperRight_(upperRight), depthBOR_(1000.0)
{}

/*!
 * \name Problem parameters
 */
// \{

//! The problem name.
/*! This is used as a prefix for files generated by the simulation.
*/
std::string name() const
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
 * \copydoc TestDecTwoPTwoCProblem::referencePressureAtPos()
 */
Scalar referencePressureAtPos(const GlobalPosition& globalPos) const
{
    return 1e6;
}
/*!
 * \copydoc TestDecTwoPTwoCProblem::boundaryTypesAtPos()
 */
void boundaryTypesAtPos(BoundaryTypes &bcTypes, const GlobalPosition& globalPos) const
{
    if (globalPos[0] > 10 - eps_ || globalPos[0] < eps_)
        bcTypes.setAllDirichlet();
    else
        // all other boundaries
        bcTypes.setAllNeumann();
}

/*!
 * \copydoc TestDecTwoPTwoCProblem::boundaryFormulation()
 */
void boundaryFormulation(typename Indices::BoundaryFormulation &bcFormulation, const Intersection& intersection) const
{
    bcFormulation = Indices::concentration;
}
/*!
 * \copydoc TestDecTwoPTwoCProblem::dirichletAtPos()
 */
void dirichletAtPos(PrimaryVariables &bcValues ,const GlobalPosition& globalPos) const
{
    Scalar pRef = referencePressureAtPos(globalPos);
    Scalar temp = temperatureAtPos(globalPos);

    // Dirichlet for pressure equation
    bcValues[Indices::pressureEqIdx] = (globalPos[0] < eps_) ? (2.5e5 - FluidSystem::H2O::liquidDensity(temp, pRef) * this->gravity()[dim-1])
            : (2e5 - FluidSystem::H2O::liquidDensity(temp, pRef) * this->gravity()[dim-1]);

    // Dirichlet values for transport equations
    bcValues[Indices::contiWEqIdx] = 1.;
    bcValues[Indices::contiNEqIdx] = 1.- bcValues[Indices::contiWEqIdx];

}
/*!
 * \copydoc TestDecTwoPTwoCProblem::neumannAtPos()
 */
void neumannAtPos(PrimaryVariables &neumannValues, const GlobalPosition& globalPos) const
{
    neumannValues[Indices::contiNEqIdx] = 0.;
    neumannValues[Indices::contiWEqIdx] = 0.;
}
/*!
 * \copydoc TestDecTwoPTwoCProblem::sourceAtPos()
 */
void sourceAtPos(PrimaryVariables &sourceValues, const GlobalPosition& globalPos) const
{
    this->setZero(sourceValues);
    using std::abs;
    if (abs(globalPos[0] - 4.8) < 0.5 + eps_ && abs(globalPos[1] - 4.8) < 0.5 + eps_)
        sourceValues[Indices::contiNEqIdx] = 0.0001;
}
/*!
 * \copydoc TestDecTwoPTwoCProblem::initialFormulation()
 */
void initialFormulation(typename Indices::BoundaryFormulation &initialFormulation, const Element& element) const
{
    initialFormulation = Indices::concentration;
}
/*!
 * \copydoc TestDecTwoPTwoCProblem::initConcentrationAtPos()
 */
Scalar initConcentrationAtPos(const GlobalPosition& globalPos) const
{
    return 1;
}

private:
GlobalPosition lowerLeft_;
GlobalPosition upperRight_;

static constexpr Scalar eps_ = 1e-6;
const Scalar depthBOR_;
};
} //end namespace

#endif
