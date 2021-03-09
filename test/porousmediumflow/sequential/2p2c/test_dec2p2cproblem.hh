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
 * \brief test problem for the sequential 2p2c model
 */
#ifndef DUMUX_TEST_2P2C_PROBLEM_HH
#define DUMUX_TEST_2P2C_PROBLEM_HH

#include <dune/grid/yaspgrid.hh>

#include <dumux/porousmediumflow/2p2c/sequential/problem.hh>
#include <dumux/porousmediumflow/2p2c/sequential/fvpressure.hh>
#include <dumux/porousmediumflow/2p2c/sequential/fvtransport.hh>

// fluid properties
#include <dumux/material/fluidsystems/h2oair.hh>

#include "test_dec2p2c_spatialparams.hh"

namespace Dumux
{
/*!
 * \ingroup SequentialTwoPTwoCTests
 * \brief test problem for the sequential 2p2c model
 */
template<class TypeTag>
class TestDecTwoPTwoCProblem;

// Specify the properties
namespace Properties
{
// Create new type tags
namespace TTag {
struct TestDecTwoPTwoC { using InheritsFrom = std::tuple<Test2P2CSpatialParams, SequentialTwoPTwoC>; };
} // end namespace TTag

// Set the grid type
template<class TypeTag>
struct Grid<TypeTag, TTag::TestDecTwoPTwoC> { using type = Dune::YaspGrid<3>; };

// Set the problem property
template<class TypeTag>
struct Problem<TypeTag, TTag::TestDecTwoPTwoC> { using type = TestDecTwoPTwoCProblem<TypeTag>; };

// Set the model properties
template<class TypeTag>
struct TransportModel<TypeTag, TTag::TestDecTwoPTwoC> { using type = FVTransport2P2C<TypeTag>; };

template<class TypeTag>
struct PressureModel<TypeTag, TTag::TestDecTwoPTwoC> { using type = FVPressure2P2C<TypeTag>; };


template<class TypeTag>
struct PressureFormulation<TypeTag, TTag::TestDecTwoPTwoC> { static constexpr int value = GetPropType<TypeTag, Properties::Indices>::pressureN; };

// Select fluid system
template<class TypeTag>
struct FluidSystem<TypeTag, TTag::TestDecTwoPTwoC>
{
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using type = FluidSystems::H2OAir<Scalar, Components::H2O<Scalar>>;
};

template<class TypeTag>
struct EnableCapillarity<TypeTag, TTag::TestDecTwoPTwoC> { static constexpr bool value = true; };
template<class TypeTag>
struct BoundaryMobility<TypeTag, TTag::TestDecTwoPTwoC> { static constexpr int value = GetPropType<TypeTag, Properties::Indices>::satDependent; };
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
 * <tt>./test_dec2p2c -parameterFile ./test_dec2p2c.input</tt>
 */
template<class TypeTag>
class TestDecTwoPTwoCProblem: public IMPETProblem2P2C<TypeTag>
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
TestDecTwoPTwoCProblem(TimeManager& timeManager, Grid& grid) :
ParentType(timeManager, grid), depthBOR_(1000.0)
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
 * \param globalPos The global Position
 */
Scalar temperatureAtPos(const GlobalPosition& globalPos) const
{
    return 273.15 + 10; // -> 10°C
}

// \}
//! Returns the reference pressure.
 /*This pressure is used in order to calculate the material properties
 * at the beginning of the initialization routine. It should lie within
 * a reasonable pressure range for the current problem.
 * \param globalPos The global Position
 */
Scalar referencePressureAtPos(const GlobalPosition& globalPos) const
{
    return 1e6;
}
//! Type of boundary condition.
/*! Defines the type the boundary condition for the conservation equation,
 *  e.g. Dirichlet or Neumann. The Pressure equation is acessed via
 *  Indices::pressureEqIdx, while the transport (or mass conservation -)
 *  equations are reached with Indices::contiWEqIdx and Indices::contiNEqIdx
 *
 * \param bcTypes The boundary types for the conservation equations
 * \param globalPos The global Position
 */
void boundaryTypesAtPos(BoundaryTypes &bcTypes, const GlobalPosition& globalPos) const
{
    if (globalPos[0] > this->bBoxMax()[0] - eps_ || globalPos[0] < eps_)
        bcTypes.setAllDirichlet();
    else
        // all other boundaries
        bcTypes.setAllNeumann();
}

//! Flag for the type of Dirichlet conditions
/*! The Dirichlet BCs can be specified by a given concentration (mass of
 * a component per total mass inside the control volume) or by means
 * of a saturation.
 *
 * \param bcFormulation The boundary formulation for the conservation equations.
 * \param intersection The intersection on the boundary.
 */
void boundaryFormulation(typename Indices::BoundaryFormulation &bcFormulation, const Intersection& intersection) const
{
    bcFormulation = Indices::concentration;
}
//! Values for dirichlet boundary condition \f$ [Pa] \f$ for pressure and \f$ \frac{mass}{totalmass} \f$ or \f$ S_{\alpha} \f$ for transport.
/*! In case of a dirichlet BC, values for all primary variables have to be set. In the sequential 2p2c model, a pressure
 * is required for the pressure equation and component fractions for the transport equations. Although one BC for the two
 * transport equations can be deduced by the other, it is seperately defined for consistency reasons with other models.
 * Depending on the boundary Formulation, either saturation or total mass fractions can be defined.
 *
 * \param bcValues Vector holding values for all equations (pressure and transport).
 * \param globalPos The global Position
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

//! Value for neumann boundary condition \f$ [\frac{kg}{m^3 \cdot s}] \f$.
/*! In case of a neumann boundary condition, the flux of matter
 *  is returned. Both pressure and transport module regard the neumann-bc values
 *  with the tranport (mass conservation -) equation indices. So the first entry
 *  of the vector is superflous.
 *  An influx into the domain has negative sign.
 *
 * \param neumannValues Vector holding values for all equations (pressure and transport).
 * \param globalPos The global Position
 */
void neumannAtPos(PrimaryVariables &neumannValues, const GlobalPosition& globalPos) const
{
    this->setZero(neumannValues);
}
//! Source of mass \f$ [\frac{kg}{m^3 \cdot s}] \f$
/*! Evaluate the source term for all phases within a given
 *  volume. The method returns the mass generated (positive) or
 *  annihilated (negative) per volume unit.
 *  Both pressure and transport module regard the neumann-bc values
 *  with the tranport (mass conservation -) equation indices. So the first entry
 *  of the vector is superflous.
 *
 * \param sourceValues Vector holding values for all equations (pressure and transport).
 * \param globalPos The global Position
 */
void sourceAtPos(PrimaryVariables &sourceValues, const GlobalPosition& globalPos) const
{
    this->setZero(sourceValues);
    using std::abs;
    if (abs(globalPos[0] - 4.8) < 0.5 + eps_ && abs(globalPos[1] - 4.8) < 0.5 + eps_)
        sourceValues[Indices::contiNEqIdx] = 0.0001;
}
//! Flag for the type of initial conditions
/*! The problem can be initialized by a given concentration (mass of
 * a component per total mass inside the control volume) or by means
 * of a saturation.
 */
void initialFormulation(typename Indices::BoundaryFormulation &initialFormulation, const Element& element) const
{
    initialFormulation = Indices::concentration;
}
//! Concentration initial condition (dimensionless)
/*! The problem is initialized with the following concentration.
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
