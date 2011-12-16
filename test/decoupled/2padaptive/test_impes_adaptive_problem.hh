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
 * \brief test problem for the sequential 2p model
 */
#ifndef DUMUX_TEST_IMPES_ADAPTIVE_PROBLEM_HH
#define DUMUX_TEST_IMPES_ADAPTIVE_PROBLEM_HH

#include <dune/grid/uggrid.hh>

//#include <dune/grid/yaspgrid.hh>
//#include <dune/grid/sgrid.hh>
#include <dune/grid/io/file/dgfparser/dgfug.hh>

#include <dumux/material/fluidsystems/liquidphase.hh>
#include <dumux/material/components/simpleh2o.hh>
#include <dumux/material/components/oil.hh>

#include <dumux/decoupled/2p/impes/impesproblem2p.hh>
#include <dumux/decoupled/2p/diffusion/fv/fvvelocity2padaptive.hh>
#include <dumux/decoupled/2p/transport/fv/fvsaturation2p.hh>
#include <dumux/decoupled/2p/transport/fv/capillarydiffusion.hh>
#include <dumux/decoupled/2p/transport/fv/gravitypart.hh>
#include <dumux/decoupled/2p/variableclass2p_gridadapt.hh>

#include "test_impes_adaptive_spatialparams.hh"

#include<dumux/decoupled/2p/transport/fv/evalcflflux_coats.hh>

namespace Dumux
{

template<class TypeTag>
class TestIMPESAdaptiveProblem;

//////////
// Specify the properties
//////////
namespace Properties
{
NEW_TYPE_TAG(TestIMPESAdaptiveProblem, INHERITS_FROM(DecoupledTwoP, Transport, TestIMPESAdaptiveSpatialParams));

// Set the grid type
SET_PROP(TestIMPESAdaptiveProblem, Grid)
{
    typedef Dune::UGGrid<2> type;
};

// Set the problem property
SET_TYPE_PROP(TestIMPESAdaptiveProblem, Problem, Dumux::TestIMPESAdaptiveProblem<TTAG(TestIMPESAdaptiveProblem)>);

// Set the model properties
SET_TYPE_PROP(TestIMPESAdaptiveProblem, TransportModel, Dumux::FVSaturation2P<TTAG(TestIMPESAdaptiveProblem)>);

SET_TYPE_PROP(TestIMPESAdaptiveProblem, DiffusivePart, Dumux::CapillaryDiffusion<TypeTag>);
SET_TYPE_PROP(TestIMPESAdaptiveProblem, ConvectivePart, Dumux::GravityPart<TypeTag>);

SET_PROP(TestIMPESAdaptiveProblem, PressureModel)
{
    typedef Dumux::FVVelocity2Padaptive<TTAG(TestIMPESAdaptiveProblem)> type;
};

//SET_INT_PROP(TestIMPESAdaptiveProblem, Formulation,
//        DecoupledTwoPCommonIndices::pnSn);


// Set the wetting phase
SET_PROP(TestIMPESAdaptiveProblem, WettingPhase)
{
private:
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(Scalar)) Scalar;
public:
    typedef Dumux::LiquidPhase<Scalar, Dumux::SimpleH2O<Scalar> > type;
};

// Set the non-wetting phase
SET_PROP(TestIMPESAdaptiveProblem, NonWettingPhase)
{
private:
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(Scalar)) Scalar;
public:
    typedef Dumux::LiquidPhase<Scalar, Dumux::SimpleH2O<Scalar> > type;
};

SET_PROP(TestIMPESAdaptiveProblem, Variables)
{
    typedef Dumux::VariableClass2PGridAdapt<TypeTag> type;
};

// Enable gravity
SET_BOOL_PROP(TestIMPESAdaptiveProblem, EnableGravity, false);
SET_BOOL_PROP(TestIMPESAdaptiveProblem, AdaptiveGrid, true);

SET_TYPE_PROP(TestIMPESAdaptiveProblem, EvalCflFluxFunction, Dumux::EvalCflFluxCoats<TypeTag>);

SET_SCALAR_PROP(TestIMPESAdaptiveProblem, CFLFactor, 0.95);
}

/*!
 * \ingroup IMPETtests
 *
 * \brief test problem for the sequential 2p model
 *
 * Water is injected from the left side into a rectangular 2D domain also
 * filled with water. Upper and lower boundary is closed (Neumann = 0),
 * and there is free outflow on the right side.
 *
 * To run the simulation execute the following line in shell:
 * <tt>./test_impes 1e8</tt>,
 * where the argument defines the simulation endtime.
 */
template<class TypeTag = TTAG(TestIMPESAdaptiveProblem)>
class TestIMPESAdaptiveProblem: public IMPESProblem2P<TypeTag>//IMPESProblem2Padaptive<TypeTag>
{
//typedef IMPESProblem2Padaptive<TypeTag> ParentType;
typedef IMPESProblem2P<TypeTag> ParentType;
typedef typename GET_PROP_TYPE(TypeTag, PTAG(Grid)) Grid;
typedef typename GET_PROP_TYPE(TypeTag, PTAG(GridView)) GridView;

typedef typename GET_PROP_TYPE(TypeTag, PTAG(TwoPIndices)) Indices;

typedef typename GET_PROP_TYPE(TypeTag, PTAG(FluidSystem)) FluidSystem;
typedef typename GET_PROP_TYPE(TypeTag, PTAG(FluidState)) FluidState;

typedef typename GET_PROP_TYPE(TypeTag, PTAG(TimeManager)) TimeManager;

enum
{
    dim = GridView::dimension, dimWorld = GridView::dimensionworld
};

enum
{
    wPhaseIdx = Indices::wPhaseIdx, nPhaseIdx = Indices::nPhaseIdx,
    pWIdx = Indices::pwIdx,
    SwIdx = Indices::SwIdx,
    eqIdxPress = Indices::pressEqIdx,
    eqIdxSat = Indices::satEqIdx
};

typedef typename GET_PROP_TYPE(TypeTag, PTAG(Scalar)) Scalar;

typedef typename GridView::Traits::template Codim<0>::Entity Element;
typedef typename GridView::Intersection Intersection;
typedef Dune::FieldVector<Scalar, dimWorld> GlobalPosition;

typedef typename GET_PROP_TYPE(TypeTag, PTAG(BoundaryTypes)) BoundaryTypes;
typedef typename GET_PROP(TypeTag, PTAG(SolutionTypes)) SolutionTypes;
    typedef typename SolutionTypes::PrimaryVariables PrimaryVariables;

public:
TestIMPESAdaptiveProblem(TimeManager &timeManager, const GridView &gridView) :
ParentType(timeManager, gridView)
{
this->setOutputInterval(10);
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
    return "output2padaptive";
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
Scalar temperatureAtPos(const GlobalPosition& globalPos) const
{
    return 273.15 + 10; // -> 10°C
}

// \}

//! Returns the reference pressure for evaluation of constitutive relations
Scalar referencePressureAtPos(const GlobalPosition& globalPos) const
{
    return 1e5; // -> 10°C
}

void source(PrimaryVariables &values,const Element& element) const
{
    values = 0;
}

/*!
* \brief Returns the type of boundary condition.
*
* BC for pressure equation can be dirichlet (pressure) or neumann (flux).
*
* BC for saturation equation can be dirichlet (saturation), neumann (flux), or outflow.
*/
void boundaryTypesAtPos(BoundaryTypes &bcTypes, const GlobalPosition& globalPos) const
{
        if (globalPos[0] < eps_)
        {
            bcTypes.setAllDirichlet();
        }
        else if (globalPos[0] > this->bboxMax()[0] - eps_)
        {
            bcTypes.setNeumann(eqIdxPress);
            bcTypes.setOutflow(eqIdxSat);
        }
        // all other boundaries
        else
        {
            bcTypes.setAllNeumann();
        }
}

//! set dirichlet condition  (pressure [Pa], saturation [-])
void dirichletAtPos(PrimaryVariables &values, const GlobalPosition& globalPos) const
{
    values = 0;
    if (globalPos[0] < eps_)
    {
        if (GET_PARAM(TypeTag, bool, EnableGravity))
        {
            Scalar pRef = referencePressureAtPos(globalPos);
            Scalar temp = temperatureAtPos(globalPos);
            Scalar sat = 1;

            FluidState fluidState;
            fluidState.update(sat, pRef, pRef, temp);
            values[pWIdx] = (2e5 + (this->bboxMax()[dim-1] - globalPos[dim-1]) * FluidSystem::phaseDensity(wPhaseIdx, temp, pRef, fluidState) * this->gravity().two_norm());
        }
        else
        {
            values[pWIdx] = 2e5;
        }
        values[SwIdx] = 0.8;
    }
    else
    {
        values[pWIdx] = 2e5;
        values[SwIdx] = 0.2;
    }
}

//! set neumann condition for phases (flux, [kg/(m^2 s)])
void neumannAtPos(PrimaryVariables &values, const GlobalPosition& globalPos) const
{
    values = 0;
    if (globalPos[0] > this->bboxMax()[0] - eps_)
    {
        values[nPhaseIdx] = 3e-4;
    }
}
//! return initial solution -> only saturation values have to be given!
void initial(PrimaryVariables &values,
        const Element& element) const
{
    values[pWIdx] = 0;
    values[SwIdx] = 0.2;
}

private:

static constexpr Scalar eps_ = 1e-6;
};
} //end namespace

#endif
