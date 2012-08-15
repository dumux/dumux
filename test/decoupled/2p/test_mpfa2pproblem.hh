/*****************************************************************************
 *   Copyright (C) 2012 by Markus Wolff                                      *
 *   Institute for Modelling Hydraulic and Environmental Systems             *
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
#ifndef DUMUX_TEST_MPFA2P_PROBLEM_HH
#define DUMUX_TEST_MPFA2P_PROBLEM_HH

#include <dune/grid/alugrid/2d/alugrid.hh>

#include <dumux/common/cubegridcreator.hh>

#include <dumux/material/fluidsystems/liquidphase.hh>
#include <dumux/material/components/simpleh2o.hh>
#include <dumux/material/components/simplednapl.hh>

#include <dumux/decoupled/2p/diffusion/fv/fvpressureproperties2p.hh>
#include <dumux/decoupled/2p/diffusion/fv/fvpressureproperties2padaptive.hh>
#include <dumux/decoupled/2p/diffusion/fvmpfa/omethod/fvmpfaopressureproperties2p.hh>
#include <dumux/decoupled/2p/diffusion/fvmpfa/lmethod/fvmpfalpressureproperties2p.hh>
#include <dumux/decoupled/2p/diffusion/fvmpfa/lmethod/fvmpfalpressureproperties2padaptive.hh>
#include <dumux/decoupled/2p/transport/fv/fvtransportproperties2p.hh>
#include <dumux/decoupled/2p/impes/impesproblem2p.hh>


#include<dumux/decoupled/2p/transport/fv/evalcflfluxcoats.hh>
#include<dumux/decoupled/2p/impes/gridadaptionindicator2plocal.hh>

#include "test_mpfa2pspatialparams.hh"

namespace Dumux
{

template<class TypeTag>
class MPFATwoPTestProblem;

//////////
// Specify the properties
//////////
namespace Properties
{

NEW_TYPE_TAG(MPFATwoPTestProblem, INHERITS_FROM(Test2PSpatialParams));

// Set the grid type
SET_PROP(MPFATwoPTestProblem, Grid)
{
    typedef Dune::ALUGrid<2, 2, Dune::cube, Dune::nonconforming> type;
};

// set the GridCreator property
SET_TYPE_PROP(MPFATwoPTestProblem, GridCreator, CubeGridCreator<TypeTag>);


// Set the problem property
SET_PROP(MPFATwoPTestProblem, Problem)
{
public:
    typedef Dumux::MPFATwoPTestProblem<TypeTag> type;
};

// Set the wetting phase
SET_PROP(MPFATwoPTestProblem, WettingPhase)
{
private:
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
public:
    typedef Dumux::LiquidPhase<Scalar, Dumux::SimpleH2O<Scalar> > type;
};

// Set the non-wetting phase
SET_PROP(MPFATwoPTestProblem, NonwettingPhase)
{
private:
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
public:
    typedef Dumux::LiquidPhase<Scalar, Dumux::SimpleDNAPL<Scalar> > type;
};

// Enable gravity
SET_BOOL_PROP(MPFATwoPTestProblem, ProblemEnableGravity, true);

SET_TYPE_PROP(MPFATwoPTestProblem, EvalCflFluxFunction, Dumux::EvalCflFluxCoats<TypeTag>);
SET_SCALAR_PROP(MPFATwoPTestProblem, ImpetCFLFactor, 1.0);
SET_TYPE_PROP(MPFATwoPTestProblem, AdaptionIndicator, Dumux::GridAdaptionIndicator2PLocal<TypeTag>);

NEW_TYPE_TAG(FVTwoPTestProblem, INHERITS_FROM(FVPressureTwoP, FVTransportTwoP, IMPESTwoP, MPFATwoPTestProblem));
NEW_TYPE_TAG(FVAdaptiveTwoPTestProblem, INHERITS_FROM(FVPressureTwoPAdaptive, FVTransportTwoP, IMPESTwoPAdaptive, MPFATwoPTestProblem));
NEW_TYPE_TAG(MPFAOTwoPTestProblem, INHERITS_FROM(FVMPFAOPressureTwoP, FVTransportTwoP, IMPESTwoP, MPFATwoPTestProblem));
NEW_TYPE_TAG(MPFALTwoPTestProblem, INHERITS_FROM(FVMPFALPressureTwoP, FVTransportTwoP, IMPESTwoP, MPFATwoPTestProblem));
NEW_TYPE_TAG(MPFALAdaptiveTwoPTestProblem, INHERITS_FROM(FVMPFALPressureTwoPAdaptive, FVTransportTwoP, IMPESTwoPAdaptive, MPFATwoPTestProblem));

}

/*!
 * \ingroup DecoupledProblems
 */
template<class TypeTag>
class MPFATwoPTestProblem: public IMPESProblem2P<TypeTag>
{
typedef IMPESProblem2P<TypeTag> ParentType;
typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;

typedef typename GET_PROP_TYPE(TypeTag, TwoPIndices) Indices;

typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;
typedef typename GET_PROP_TYPE(TypeTag, FluidState) FluidState;
typedef typename GET_PROP_TYPE(TypeTag, WettingPhase) WettingPhase;

typedef typename GET_PROP_TYPE(TypeTag, TimeManager) TimeManager;

typedef typename GET_PROP_TYPE(TypeTag, BoundaryTypes) BoundaryTypes;
typedef typename GET_PROP(TypeTag, SolutionTypes)::PrimaryVariables PrimaryVariables;
typedef typename GET_PROP(TypeTag, SolutionTypes)::ScalarSolution ScalarSolutionType;

typedef typename GET_PROP_TYPE(TypeTag, GridCreator) GridCreator;
typedef typename GET_PROP(TypeTag, ParameterTree) ParameterTree;

enum
{
    dim = GridView::dimension, dimWorld = GridView::dimensionworld
};

enum
{
    wPhaseIdx = Indices::wPhaseIdx,
    nPhaseIdx = Indices::nPhaseIdx,
    pWIdx = Indices::pwIdx,
    SwIdx = Indices::SwIdx,
    eqIdxPress = Indices::pressEqIdx,
    eqIdxSat = Indices::satEqIdx
};

typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;

typedef typename GridView::Traits::template Codim<0>::Entity Element;
typedef typename GridView::template Codim<0>::Iterator ElementIterator;
typedef typename GridView::Intersection Intersection;
typedef Dune::FieldVector<Scalar, dimWorld> GlobalPosition;
typedef Dune::FieldVector<Scalar, dim> LocalPosition;

public:
MPFATwoPTestProblem(TimeManager &timeManager,const GridView &gridView) :
ParentType(timeManager, gridView)
{
    this->setGrid(GridCreator::grid());

    Scalar inletWidth = 1.0;
    if (ParameterTree::tree().hasKey("Problem.InletWidth"))
    {
        inletWidth = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, Scalar, Problem, InletWidth);
    }
    GlobalPosition inletCenter = this->bboxMax();
    inletCenter[0] *= 0.5;

    inletLeftCoord_ = inletCenter;
    inletLeftCoord_[0] -=0.5*inletWidth;
    inletRightCoord_ = inletCenter;
    inletRightCoord_[0] +=0.5*inletWidth;

    inFlux_ = 1e-4;
    if (ParameterTree::tree().hasKey("Problem.InjectionFlux"))
    {
        inFlux_ = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, Scalar, Problem, InjectionFlux);
    }

    int outputInterval = 0;
    if (ParameterTree::tree().hasKey("Problem.OutputInterval"))
    {
        outputInterval = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, int, Problem, OutputInterval);
    }
    this->setOutputInterval(outputInterval);

    Scalar outputTimeInterval = 1e6;
    if (ParameterTree::tree().hasKey("Problem.OutputTimeInterval"))
    {
        outputTimeInterval = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, Scalar, Problem, OutputTimeInterval);
    }
    this->setOutputTimeInterval(outputTimeInterval);
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
    if (ParameterTree::tree().hasKey("Problem.OutputfileName"))
    {
        std::string fileName(GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, std::string, Problem, OutputfileName));
        return fileName.c_str();
    }
    else
    {
        return "testmpfa2p";
    }
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


Scalar referencePressureAtPos(const GlobalPosition& globalPos) const
{
    return 1e5; // -> 10°C
}

void source(PrimaryVariables &values,const Element& element) const
{
    values = 0;
}

void boundaryTypesAtPos(BoundaryTypes &bcTypes, const GlobalPosition& globalPos) const
{
    if (isInlet(globalPos))
    {
        bcTypes.setNeumann(eqIdxPress);
        bcTypes.setDirichlet(eqIdxSat);
    }
    else if (isBottom(globalPos) || isTop(globalPos))
    {
        bcTypes.setAllNeumann();
    }
    else
    {
        bcTypes.setDirichlet(eqIdxPress);
        bcTypes.setOutflow(eqIdxSat);
    }
}

void dirichletAtPos(PrimaryVariables &values, const GlobalPosition& globalPos) const
{
    Scalar pRef = referencePressureAtPos(globalPos);
    Scalar temp = temperatureAtPos(globalPos);

    values[pWIdx] = (1e5 - (this->bboxMax()- globalPos) * this->gravity() * WettingPhase::density(temp, pRef));
    values[SwIdx] = 1.0;

    if (isInlet(globalPos))
    {
        values[SwIdx] = 0.0;
    }
}

void neumannAtPos(PrimaryVariables &values, const GlobalPosition& globalPos) const
{
    values = 0;
    if (isInlet(globalPos))
    {
        values[nPhaseIdx] = -inFlux_;
    }

}

void initialAtPos(PrimaryVariables &values,
        const GlobalPosition& globalPos) const
{
    values[pWIdx] = 0;
    values[SwIdx] = 1.0;
}

private:

bool isInlet(const GlobalPosition& globalPos) const
{
        if (!isTop(globalPos))
            return false;

    for (int i = 0; i < dim; i++)
    {
        if (globalPos[i] < inletLeftCoord_[i] - eps_)
            return false;
        if (globalPos[i] > inletRightCoord_[i] + eps_)
            return false;
    }
    return true;
}

bool isTop(const GlobalPosition& globalPos) const
{
    if (dim == 2)
    {
        if (globalPos[1] > this->bboxMax()[1] - eps_)
            return true;
    }
    if (dim == 3)
    {
        if (globalPos[2] > this->bboxMax()[2] - eps_)
            return true;
    }
    return false;
}

bool isBottom(const GlobalPosition& globalPos) const
{
    if (dim == 2)
    {
        if (globalPos[1] < eps_)
            return true;
    }
    if (dim == 3)
    {
        if (globalPos[2] < eps_)
            return true;
    }
    return false;
}

Scalar inFlux_;
GlobalPosition inletLeftCoord_;
GlobalPosition inletRightCoord_;
static const Scalar eps_ = 1e-6;
};
} //end namespace

#endif
