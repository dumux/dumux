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
 * \brief test problem for the grid-adaptive 3d 2p2c model
 */
#ifndef DUMUX_TEST_ADAPTIVE3D_2P2C_PROBLEM_HH
#define DUMUX_TEST_ADAPTIVE3D_2P2C_PROBLEM_HH

#if HAVE_ALUGRID
#include <dune/grid/alugrid.hh>
#elif HAVE_DUNE_ALUGRID
#include <dune/alugrid/grid.hh>
#else
#include <dune/grid/yaspgrid.hh>
#endif

#include <dumux/common/math.hh>
#include <dumux/decoupled/2p2c/2p2cadaptiveproperties.hh>
#include <dumux/decoupled/2p2c/2p2cproblem.hh>
#include <dumux/decoupled/2p/impes/gridadaptionindicator2p.hh>
#include <dumux/io/cubegridcreator.hh>
#include <dumux/material/fluidmatrixinteractions/2p/linearmaterial.hh>
#include <dumux/material/fluidsystems/h2oairfluidsystem.hh>

#include "test_dec2p2c_spatialparams.hh"

namespace Dumux
{

template<class TypeTag>
class Adaptive2p2c3d;

namespace Properties
{
NEW_TYPE_TAG(Adaptive2p2c3d, INHERITS_FROM(DecoupledTwoPTwoCAdaptive,Test2P2CSpatialParams, MPFAProperties));

// Set the grid type
#if HAVE_ALUGRID || HAVE_DUNE_ALUGRID
SET_TYPE_PROP(Adaptive2p2c3d, Grid, Dune::ALUGrid<3, 3, Dune::cube, Dune::nonconforming>);
#else
SET_TYPE_PROP(Adaptive2p2c3d, Grid, Dune::YaspGrid<3>);
#endif

// set the GridCreator property
SET_TYPE_PROP(Adaptive2p2c3d, GridCreator, CubeGridCreator<TypeTag>);

// Set the problem property
SET_TYPE_PROP(Adaptive2p2c3d, Problem, Dumux::Adaptive2p2c3d<TTAG(Adaptive2p2c3d)>);

// Set the model properties
SET_TYPE_PROP(Adaptive2p2c3d, TransportModel, Dumux::FV3dTransport2P2CAdaptive<TTAG(Adaptive2p2c3d)>);

SET_TYPE_PROP(Adaptive2p2c3d, PressureModel, Dumux::FV3dPressure2P2CAdaptive<TTAG(Adaptive2p2c3d)>);

// Select fluid system
SET_TYPE_PROP(Adaptive2p2c3d, FluidSystem, Dumux::H2OAirFluidSystem<TypeTag>);

SET_BOOL_PROP(Adaptive2p2c3d, EnableComplicatedFluidSystem, false);

// Select water formulation
SET_PROP(Adaptive2p2c3d, Components) : public GET_PROP(TypeTag, DefaultComponents)
{
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef Dumux::H2O<Scalar> H2O;
};

// Specify indicator
SET_TYPE_PROP(Adaptive2p2c3d, AdaptionIndicator, GridAdaptionIndicator2P<TypeTag>);

SET_BOOL_PROP(Adaptive2p2c3d, EnableCapillarity, true);
SET_BOOL_PROP(Adaptive2p2c3d, AdaptiveGrid, true);
SET_INT_PROP(Adaptive2p2c3d, PressureFormulation, GET_PROP_TYPE(TypeTag, Indices)::pressureN);

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
template<class TypeTag = TTAG(Adaptive2p2c3d)>
class Adaptive2p2c3d: public IMPETProblem2P2C<TypeTag>
{
typedef IMPETProblem2P2C<TypeTag> ParentType;
typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
typedef typename GET_PROP_TYPE(TypeTag, Grid) Grid;
typedef typename GET_PROP_TYPE(TypeTag, GridCreator) GridCreator;
typedef typename GET_PROP_TYPE(TypeTag, TimeManager) TimeManager;

typedef typename GET_PROP_TYPE(TypeTag, Indices) Indices;
typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;
typedef typename GET_PROP_TYPE(TypeTag, SpatialParams)    SpatialParams;

// boundary typedefs
typedef typename GET_PROP_TYPE(TypeTag, BoundaryTypes) BoundaryTypes;
typedef typename GET_PROP_TYPE(TypeTag, PrimaryVariables) PrimaryVariables;
enum
{
    dim = GridView::dimension, dimWorld = GridView::dimensionworld
};

enum
{
    wPhaseIdx = Indices::wPhaseIdx, nPhaseIdx = Indices::nPhaseIdx
};

typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;

typedef typename GridView::Traits::template Codim<0>::Entity Element;
typedef typename GridView::Intersection Intersection;
typedef Dune::FieldVector<Scalar, dimWorld> GlobalPosition;

public:
Adaptive2p2c3d(TimeManager &timeManager, const GridView& gridView) :
    ParentType(timeManager, gridView),
            debugWriter_(gridView, "gridAfterAdapt")
{
    GridCreator::grid().globalRefine(GET_PARAM_FROM_GROUP(TypeTag, int, GridAdapt, MaxLevel));
    this->setGrid(GridCreator::grid());

    //Process parameter file 
    //Simulation Control
    const int outputInterval = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, int, Problem, OutputInterval);
    this->setOutputInterval(outputInterval);

    injectionrate_ = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, double, BoundaryConditions, Injectionrate);
}

/*!
 * \name Problem parameters
 */
// \{
//! @copydoc Dumux::TestDecTwoPTwoCProblem::name()
const char *name() const
{
    return GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, std::string, Problem, Name).c_str();
}

//! @copydoc Dumux::TestDecTwoPTwoCProblem::shouldWriteRestartFile()
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
    bcFormulation = Indices::BoundaryFormulation::concentration;
}

/*!
 * \copydoc Dumux::TestDecTwoPTwoCProblem::dirichletAtPos()
 */
void dirichletAtPos(PrimaryVariables &bcValues, const GlobalPosition& globalPos) const
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
    this->setZero(neumannValues);
}

/*!
 * \copydoc Dumux::TestDecTwoPTwoCProblem::sourceAtPos()
 */
void sourceAtPos(PrimaryVariables &sourceValues, const GlobalPosition& globalPos) const
{
    this->setZero(sourceValues);
    if (fabs(globalPos[0] - 4.8) < 0.5 && fabs(globalPos[1] - 4.8) < 0.5)
        sourceValues[Indices::contiNEqIdx] = injectionrate_;
}

/*!
 * \copydoc Dumux::TestDecTwoPTwoCProblem::initialFormulation()
 */
const void initialFormulation(typename Indices::BoundaryFormulation &initialFormulation, const Element& element) const
{
    initialFormulation = Indices::BoundaryFormulation::concentration;
}

/*!
 * \copydoc Dumux::TestDecTwoPTwoCProblem::initConcentrationAtPos()
 */
Scalar initConcentrationAtPos(const GlobalPosition& globalPos) const
{
    return 1.0;
}

private:
Grid grid_;
Dumux::VtkMultiWriter<GridView> debugWriter_;
Scalar injectionrate_;
};
} //end namespace

#endif // DUMUX_TEST_ADAPTIVE3D_2P2C_PROBLEM_HH
