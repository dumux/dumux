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
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the           *
 *   GNU General Public License for more details.                            *
 *                                                                           *
 *   You should have received a copy of the GNU General Public License       *
 *   along with this program.  If not, see <http://www.gnu.org/licenses/>.   *
 *****************************************************************************/
/*!
 * \file
 *
 * \brief Test problem for the adaptive sequential 2p2c model
 */
#ifndef DUMUX_TEST_ADAPTIVE_2P2C_PROBLEM_HH
#define DUMUX_TEST_ADAPTIVE_2P2C_PROBLEM_HH

#include <dune/grid/alugrid.hh>
// #include <dune/grid/uggrid.hh>

#include <dumux/io/cubegridcreator.hh>

#include <dumux/common/math.hh>
#include <dumux/decoupled/2p2c/2p2cadaptiveproperties.hh>
#include <dumux/decoupled/2p2c/2p2cproblem.hh>
#include <dumux/decoupled/2p/impes/gridadaptionindicator2p.hh>
#include <dumux/material/fluidmatrixinteractions/2p/linearmaterial.hh>

#include <dumux/material/fluidsystems/h2oairfluidsystem.hh>
#include <dumux/material/fluidsystems/h2on2fluidsystem.hh>

#include "test_dec2p2c_spatialparams.hh"

namespace Dumux
{

template<class TypeTag>
class Adaptive2p2c;

namespace Properties
{
NEW_TYPE_TAG(Adaptive2p2c, INHERITS_FROM(DecoupledTwoPTwoCAdaptive,Test2P2CSpatialParams));

// Set the grid type
SET_PROP(Adaptive2p2c, Grid)
{
#if HAVE_ALUGRID
    typedef Dune::ALUGrid<2, 2, Dune::cube, Dune::nonconforming> type;
// #elif HAVE_UG
//     typedef Dune::UGGrid<2> type;
#endif
    
};
// set the GridCreator property
SET_TYPE_PROP(Adaptive2p2c, GridCreator, CubeGridCreator<TypeTag>);

// Set the problem property
SET_PROP(Adaptive2p2c, Problem)
{
    typedef Dumux::Adaptive2p2c<TTAG(Adaptive2p2c)> type;
};

// Select fluid system
SET_PROP(Adaptive2p2c, FluidSystem)
{
    typedef Dumux::H2OAirFluidSystem<TypeTag> type;
//    typedef Dumux::H2ON2FluidSystem<TypeTag> type;
//    typedef Dumux::Brine_CO2_System<TypeTag, Dumux::Benchmark3::CO2Tables> type;
};

SET_BOOL_PROP(Adaptive2p2c, EnableComplicatedFluidSystem, true);

// Select water formulation
SET_PROP(Adaptive2p2c, Components) : public GET_PROP(TypeTag, DefaultComponents)
{
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
//    typedef Dumux::TabulatedComponent<Scalar, typename Dumux::H2O<Scalar> > H20;
        typedef Dumux::H2O<Scalar> H2O;
};
// Specify indicator
SET_TYPE_PROP(Adaptive2p2c, AdaptionIndicator, GridAdaptionIndicator2P<TypeTag>);

// Enable gravity
SET_BOOL_PROP(Adaptive2p2c, ProblemEnableGravity, true);
SET_INT_PROP(Adaptive2p2c,
        BoundaryMobility,
        GET_PROP_TYPE(TypeTag, Indices)::permDependent);
SET_BOOL_PROP(Adaptive2p2c, EnableCapillarity, true);
SET_INT_PROP(Adaptive2p2c, PressureFormulation,
        GET_PROP_TYPE(TypeTag, Indices)::pressureNw);

}

/*!
 * \ingroup Adaptive2p2c
 * \ingroup IMPETtests
 * 
 * \brief test problem for the grid-adaptive sequential 2p2c model
 *
 * The domain is box shaped (2D), as the gird-adaptive model can only run in
 * two dimensions!. All sides are closed (Neumann 0 boundary)
 * except the top and bottom boundaries (Dirichlet). A Gas (Air)
 * is injected over a vertical well in the center of the domain.
 *
 * It is the same test as for the non-adaptive compositional models.
 * The input file can alter the grid adaptation in group [GridAdapt]:
 * MinLevel and MaxLevel define the refinement levels of the grid.
 * Afterwards, the refinement tolerances can be specified. On hanging
 * nodes, an mpfa can be used (EnableMultiPointFluxApproximation = 1), at
 * best taken both half-edges into account (EnableSecondHalfEdge = 1).
 *
 * To run the simulation execute the following line in shell:
 * <tt>./test_adaptive2p2c -parameterFile ./test_adaptive2p2c.input</tt>
 */
template<class TypeTag = TTAG(Adaptive2p2c)>
class Adaptive2p2c: public IMPETProblem2P2C<TypeTag>
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
typedef typename Grid::Traits::template Codim<0>::EntityPointer ElementPointer;
typedef typename GridView::Intersection Intersection;
typedef Dune::FieldVector<Scalar, dimWorld> GlobalPosition;

public:
Adaptive2p2c(TimeManager &timeManager, const GridView& gridView) :
    ParentType(timeManager, gridView),
            debugWriter_(gridView, "gridAfterAdapt")
{
    this->setGrid(GridCreator::grid());
    std::string s = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, std::string, Problem, SimulationName);
    this->setName(s.c_str());
    this->setOutputInterval(GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, int, Problem, OutputInterval));
    // initialize the tables of the fluid system
//    WaterFormulation::init(273.15, 623.15, 100,
//                            -10,   20e6, 200);
//    FluidSystem::init();
}

//void preTimeStep()
//{
//    ParentType::preTimeStep();
//            // use second writer
//            debugWriter_.gridChanged();
//            // write
//            debugWriter_.beginWrite(this->timeManager().time());
//            //write stuff out
//            typedef typename GET_PROP(TypeTag, PTAG(SolutionTypes))::ScalarSolution ScalarSolutionType;
//            typedef typename GET_PROP_TYPE(TypeTag, PTAG(CellData)) CellData;
//            int size = this->gridView().size(0);
//            ScalarSolutionType *pressureW = debugWriter_.allocateManagedBuffer (size);
//            ScalarSolutionType *pressureN = debugWriter_.allocateManagedBuffer (size);
//            ScalarSolutionType *totalConcentration1 = debugWriter_.allocateManagedBuffer (size);
//            ScalarSolutionType *totalConcentration2 = debugWriter_.allocateManagedBuffer (size);
//            for (int i = 0; i < size; i++)
//            {
//                CellData& cellData = this->variables().cellData(i);
//                (*pressureW)[i] = cellData.pressure(wPhaseIdx);
//                (*pressureN)[i] = cellData.pressure(nPhaseIdx);
//                (*totalConcentration1)[i] = cellData.massConcentration(wPhaseIdx);
//                (*totalConcentration2)[i] = cellData.massConcentration(nPhaseIdx);
//            }
//            debugWriter_.attachCellData(*pressureW, "wetting pressure");
//            debugWriter_.attachCellData(*pressureN, "nonwetting pressure");
//            debugWriter_.attachCellData(*totalConcentration1, "C^w from cellData");
//            debugWriter_.attachCellData(*totalConcentration2, "C^n from cellData");
//            debugWriter_.endWrite();
//            return;
//}

/*!
 * \name Problem parameters
 */
// \{

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
    if (globalPos[0] > this->bBoxMax()[0]-1E-6 || globalPos[0] < 1e-6)
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
 * \copydoc Dumux::IMPETProblem::source()
 */
void source(PrimaryVariables &values, const Element &element)
{
    this->setZero(values);
    ElementPointer father(element);
    // access level 1 entity
    while (father->level() != this->gridAdapt().getMinLevel())
    {
        father = father->father();
    }
    GlobalPosition globalPos = father->geometry().center();
    if (fabs(globalPos[0] - 4.8) < 0.5 && fabs(globalPos[1] - 4.8) < 0.5)
        values[Indices::contiNEqIdx] = 0.0001;
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
    return 1.0;
}

private:
Grid grid_;
Dumux::VtkMultiWriter<GridView> debugWriter_;
};
} //end namespace

#endif
