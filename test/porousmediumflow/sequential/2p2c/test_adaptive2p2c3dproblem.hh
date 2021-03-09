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
 * \brief test problem for the grid-adaptive 3d 2p2c model
 */
#ifndef DUMUX_TEST_ADAPTIVE3D_2P2C_PROBLEM_HH
#define DUMUX_TEST_ADAPTIVE3D_2P2C_PROBLEM_HH

#if HAVE_UG
#include <dune/grid/uggrid.hh>
#endif
#include <dune/grid/yaspgrid.hh>

#include <dumux/common/math.hh>
#include <dumux/porousmediumflow/2p2c/sequential/adaptiveproperties.hh>
#include <dumux/porousmediumflow/2p2c/sequential/problem.hh>
#include <dumux/porousmediumflow/2p/sequential/impes/gridadaptionindicator.hh>
#include <dumux/material/fluidmatrixinteractions/2p/linearmaterial.hh>
#include <dumux/material/fluidsystems/h2oair.hh>

#include "test_dec2p2c_spatialparams.hh"

namespace Dumux
{
/*!
 * \ingroup SequentialTwoPTwoCTests
 * \brief test problem for the grid-adaptive 3d 2p2c model
 */
template<class TypeTag>
class Adaptive2p2c3d;

namespace Properties
{
// Create new type tags
namespace TTag {
struct Adaptive2p2c3d { using InheritsFrom = std::tuple<Test2P2CSpatialParams, MPFAProperties, SequentialTwoPTwoCAdaptive>; };
} // end namespace TTag

// Set the grid type
#if HAVE_UG
template<class TypeTag>
struct Grid<TypeTag, TTag::Adaptive2p2c3d> { using type = Dune::UGGrid<3>; };
#else
template<class TypeTag>
struct Grid<TypeTag, TTag::Adaptive2p2c3d> { using type = Dune::YaspGrid<3>; };
#endif

// Set the problem property
template<class TypeTag>
struct Problem<TypeTag, TTag::Adaptive2p2c3d> { using type = Adaptive2p2c3d<Properties::TTag::Adaptive2p2c3d>; };

// Set the model properties
template<class TypeTag>
struct TransportModel<TypeTag, TTag::Adaptive2p2c3d> { using type = FV3dTransport2P2CAdaptive<Properties::TTag::Adaptive2p2c3d>; };

template<class TypeTag>
struct PressureModel<TypeTag, TTag::Adaptive2p2c3d> { using type = FV3dPressure2P2CAdaptive<Properties::TTag::Adaptive2p2c3d>; };

// Select fluid system
template<class TypeTag>
struct FluidSystem<TypeTag, TTag::Adaptive2p2c3d>
{
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using type = FluidSystems::H2OAir<Scalar, Components::H2O<Scalar>, FluidSystems::H2OAirDefaultPolicy</*fastButSimplifiedRelations=*/true>>;
};

// Specify indicator
template<class TypeTag>
struct AdaptionIndicator<TypeTag, TTag::Adaptive2p2c3d> { using type = GridAdaptionIndicator2P<TypeTag>; };

template<class TypeTag>
struct EnableCapillarity<TypeTag, TTag::Adaptive2p2c3d> { static constexpr bool value = true; };
template<class TypeTag>
struct PressureFormulation<TypeTag, TTag::Adaptive2p2c3d> { static constexpr int value = GetPropType<TypeTag, Properties::Indices>::pressureN; };

}

/*!
 * \ingroup Adaptive2p2c
 * \ingroup IMPETtests
 *
 * \brief test problem for the grid-adaptive sequential 2p2c model
 *
 * The domain is box shaped (2D). All sides are closed (Neumann 0 boundary)
 * except the top and bottom boundaries (Dirichlet). A Gas (Air)
 * is injected over a vertical well in the center of the domain.
 *
 * It is the same test as for the non-adaptive compositional models.
 * The input file can alter the grid adaptation in group [GridAdapt]:
 * MinLevel and MaxLevel define the refinement levels of the grid.
 * Afterwards, the refinement tolerances can be specified. On hanging
 * nodes, a MPFA can be used (EnableMultiPointFluxApproximation = 1), at
 * best taken both interaction regions into account (MaxInteractionVolumes = 2).
 *
 * To run the simulation execute the following line in shell:
 * <tt>./test_adaptive2p2c -parameterFile ./test_adaptive2p2c.input</tt>
 */
template<class TypeTag = Properties::TTag::Adaptive2p2c3d>
class Adaptive2p2c3d: public IMPETProblem2P2C<TypeTag>
{
using ParentType = IMPETProblem2P2C<TypeTag>;
using GridView = typename GetPropType<TypeTag, Properties::GridGeometry>::GridView;
using Grid = GetPropType<TypeTag, Properties::Grid>;
using TimeManager = GetPropType<TypeTag, Properties::TimeManager>;

using Indices = typename GetPropType<TypeTag, Properties::ModelTraits>::Indices;
using FluidSystem = GetPropType<TypeTag, Properties::FluidSystem>;
using SpatialParams = GetPropType<TypeTag, Properties::SpatialParams>;

using BoundaryTypes = GetPropType<TypeTag, Properties::SequentialBoundaryTypes>;
using PrimaryVariables = GetPropType<TypeTag, Properties::PrimaryVariables>;
enum
{
    dim = GridView::dimension, dimWorld = GridView::dimensionworld
};

enum
{
    wPhaseIdx = Indices::wPhaseIdx, nPhaseIdx = Indices::nPhaseIdx
};

using Scalar = GetPropType<TypeTag, Properties::Scalar>;

using Element = typename GridView::Traits::template Codim<0>::Entity;
using Intersection = typename GridView::Intersection;
using GlobalPosition = typename Element::Geometry::GlobalCoordinate;

public:
Adaptive2p2c3d(TimeManager& timeManager, Grid& grid) :
    ParentType(timeManager, grid),
            debugWriter_(grid.leafGridView(), "gridAfterAdapt")
{
    grid.globalRefine(getParam<int>("GridAdapt.MaxLevel"));

    //Process parameter file
    //Simulation Control
    const int outputInterval = getParam<int>("Problem.OutputInterval");
    this->setOutputInterval(outputInterval);

    injectionrate_ = getParam<Scalar>("BoundaryConditions.Injectionrate");
}

/*!
 * \name Problem parameters
 */
// \{
//! @copydoc TestDecTwoPTwoCProblem::name()
std::string name() const
{
    return getParam<std::string>("Problem.Name");
}

//! @copydoc TestDecTwoPTwoCProblem::shouldWriteRestartFile()
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
    bcFormulation = Indices::BoundaryFormulation::concentration;
}

/*!
 * \copydoc TestDecTwoPTwoCProblem::dirichletAtPos()
 */
void dirichletAtPos(PrimaryVariables &bcValues, const GlobalPosition& globalPos) const
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
    this->setZero(neumannValues);
}

/*!
 * \copydoc TestDecTwoPTwoCProblem::sourceAtPos()
 */
void sourceAtPos(PrimaryVariables &sourceValues, const GlobalPosition& globalPos) const
{
    this->setZero(sourceValues);
    using std::abs;
    if (abs(globalPos[0] - 4.8) < 0.5 + eps_ && abs(globalPos[1] - 4.8) < 0.5 + eps_)
        sourceValues[Indices::contiNEqIdx] = injectionrate_;
}

/*!
 * \copydoc TestDecTwoPTwoCProblem::initialFormulation()
 */
void initialFormulation(typename Indices::BoundaryFormulation &initialFormulation, const Element& element) const
{
    initialFormulation = Indices::BoundaryFormulation::concentration;
}

/*!
 * \copydoc TestDecTwoPTwoCProblem::initConcentrationAtPos()
 */
Scalar initConcentrationAtPos(const GlobalPosition& globalPos) const
{
    return 1.0;
}

private:
VtkMultiWriter<GridView> debugWriter_;
Scalar injectionrate_;
static constexpr Scalar eps_ = 1e-6;
};
} //end namespace

#endif
