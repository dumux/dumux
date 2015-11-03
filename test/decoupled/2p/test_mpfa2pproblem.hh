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
 * \brief test problem for sequential 2p models
 */

#ifndef DUMUX_TEST_MPFA2P_PROBLEM_HH
#define DUMUX_TEST_MPFA2P_PROBLEM_HH

#if HAVE_DUNE_ALUGRID
#include <dune/alugrid/grid.hh>
#endif

#include <dumux/io/cubegridcreator.hh>

#include <dumux/material/fluidsystems/liquidphase.hh>
#include <dumux/material/components/simpleh2o.hh>
#include <dumux/material/components/dnapl.hh>

#include <dumux/decoupled/2p/diffusion/fv/fvpressureproperties2p.hh>
#include <dumux/decoupled/2p/diffusion/fv/fvpressureproperties2padaptive.hh>
#include <dumux/decoupled/2p/diffusion/fvmpfa/omethod/fvmpfao2dpressureproperties2p.hh>
#include <dumux/decoupled/2p/diffusion/fvmpfa/lmethod/fvmpfal2dpressureproperties2p.hh>
#include <dumux/decoupled/2p/diffusion/fvmpfa/lmethod/fvmpfal2dpressureproperties2padaptive.hh>
#include <dumux/decoupled/2p/transport/fv/fvtransportproperties2p.hh>
#include <dumux/decoupled/2p/impes/impesproblem2p.hh>


#include <dumux/decoupled/2p/transport/fv/evalcflfluxcoats.hh>
#include <dumux/decoupled/2p/impes/gridadaptionindicator2plocal.hh>

#include "test_mpfa2pspatialparams.hh"
#include "buckleyleverettanalyticsolution.hh"
#include "mcwhorteranalyticsolution.hh"

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
#if HAVE_DUNE_ALUGRID
SET_TYPE_PROP(MPFATwoPTestProblem, Grid, Dune::ALUGrid<2, 2, Dune::cube, Dune::nonconforming>);
#endif

#if PROBLEM == 2
//set the GridCreator property
SET_TYPE_PROP(MPFATwoPTestProblem, GridCreator, CubeGridCreator<TypeTag>);
#endif


// Set the problem property
SET_TYPE_PROP(MPFATwoPTestProblem, Problem, Dumux::MPFATwoPTestProblem<TypeTag>);

// Set the wetting phase
SET_PROP(MPFATwoPTestProblem, WettingPhase)
{
private:
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
public:
    typedef Dumux::LiquidPhase<Scalar, Dumux::SimpleH2O<Scalar> > type;
};

#if PROBLEM == 2
// Set the non-wetting phase
SET_PROP(MPFATwoPTestProblem, NonwettingPhase)
{
private:
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
public:
    typedef Dumux::LiquidPhase<Scalar, Dumux::DNAPL<Scalar> > type;
};
#else
// Set the non-wetting phase
SET_PROP(MPFATwoPTestProblem, NonwettingPhase)
{
private:
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
public:
    typedef Dumux::LiquidPhase<Scalar, Dumux::SimpleH2O<Scalar> > type;
};
#endif

#if PROBLEM == 1
SET_INT_PROP(MPFATwoPTestProblem, Formulation, DecoupledTwoPCommonIndices::pnsw);
#endif

#if PROBLEM == 2
// Enable gravity
SET_BOOL_PROP(MPFATwoPTestProblem, ProblemEnableGravity, true);
#else
SET_BOOL_PROP(MPFATwoPTestProblem, ProblemEnableGravity, false);
#endif

SET_TYPE_PROP(MPFATwoPTestProblem, EvalCflFluxFunction, Dumux::EvalCflFluxCoats<TypeTag>);
SET_SCALAR_PROP(MPFATwoPTestProblem, ImpetCFLFactor, 1.0);
SET_TYPE_PROP(MPFATwoPTestProblem, AdaptionIndicator, Dumux::GridAdaptionIndicator2PLocal<TypeTag>);

NEW_TYPE_TAG(FVTwoPTestProblem, INHERITS_FROM(FVPressureTwoP, FVTransportTwoP, IMPESTwoP, MPFATwoPTestProblem));
NEW_TYPE_TAG(FVAdaptiveTwoPTestProblem, INHERITS_FROM(FVPressureTwoPAdaptive, FVTransportTwoP, IMPESTwoPAdaptive, MPFATwoPTestProblem));
NEW_TYPE_TAG(MPFAOTwoPTestProblem, INHERITS_FROM(FvMpfaO2dPressureTwoP, FVTransportTwoP, IMPESTwoP, MPFATwoPTestProblem));
NEW_TYPE_TAG(MPFALTwoPTestProblem, INHERITS_FROM(FvMpfaL2dPressureTwoP, FVTransportTwoP, IMPESTwoP, MPFATwoPTestProblem));
NEW_TYPE_TAG(MPFALAdaptiveTwoPTestProblem, INHERITS_FROM(FvMpfaL2dPressureTwoPAdaptive, FVTransportTwoP, IMPESTwoPAdaptive, MPFATwoPTestProblem));

}

/*!
 * \ingroup IMPETtests
 *
 * \brief test problem for sequential 2p models
 *
 * A DNAPL is injected from the top into a rectangular 2D domain saturated by water.
 * The remaining upper and the lower boundary is closed (Neumann = 0). At the sides a hydrostatic pressure condition
 * and free outflow for saturation are set. The domain is heterogeneous with a backround material and three lenses.
 *
 * To run the simulation execute the following line in shell:
 *
 * <tt>./test_mpfa2p</tt>,
 *
 * Additionally, the numerical model can be switched by executing with the parameter "ModelType":
 *
 * <tt>./test_mpfa2p --ModelType=...</tt>,
 *
 * where ModelType can be:
 * - FV (standard finite volume)
 * - FVAdaptive (adaptive finite volume)
 * - MPFAO (MPFA o-method)
 * - MPFAL (MPFA l-method)
 * - MPFALAdaptive (adaptive MPFA l-method)
 */
template<class TypeTag>
class MPFATwoPTestProblem: public IMPESProblem2P<TypeTag>
{
typedef IMPESProblem2P<TypeTag> ParentType;
typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;

typedef typename GET_PROP_TYPE(TypeTag, Indices) Indices;

typedef typename GET_PROP_TYPE(TypeTag, WettingPhase) WettingPhase;

typedef typename GET_PROP_TYPE(TypeTag, TimeManager) TimeManager;

typedef typename GET_PROP_TYPE(TypeTag, BoundaryTypes) BoundaryTypes;
typedef typename GET_PROP(TypeTag, SolutionTypes) SolutionTypes;
typedef typename SolutionTypes::PrimaryVariables PrimaryVariables;

typedef typename GET_PROP_TYPE(TypeTag, GridCreator) GridCreator;
typedef typename GET_PROP(TypeTag, ParameterTree) ParameterTree;

enum
{
    dim = GridView::dimension, dimWorld = GridView::dimensionworld
};

enum
{
    nPhaseIdx = Indices::nPhaseIdx,
#if PROBLEM == 1
    pnIdx = Indices::pnIdx,
#else
    pwIdx = Indices::pwIdx,
#endif
    swIdx = Indices::swIdx,
    eqIdxPress = Indices::pressureEqIdx,
    eqIdxSat = Indices::satEqIdx
};

typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;

typedef typename GridView::Traits::template Codim<0>::Entity Element;
typedef Dune::FieldVector<Scalar, dimWorld> GlobalPosition;

public:
MPFATwoPTestProblem(TimeManager &timeManager,const GridView &gridView) :
ParentType(timeManager, gridView)
#if PROBLEM != 2
, analyticSolution_(*this)
#endif
{
    this->setGrid(GridCreator::grid());

    int refinementFactor = 0;
    if (ParameterTree::tree().hasKey("Grid.RefinementFactor"))
    {
        refinementFactor = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, Scalar, Grid, RefinementFactor);
    }

    this->grid().globalRefine(refinementFactor);

    Scalar inletWidth = 1.0;
    if (ParameterTree::tree().hasKey("Problem.InletWidth"))
    {
        inletWidth = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, Scalar, Problem, InletWidth);
    }
    GlobalPosition inletCenter = this->bBoxMax();
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

#if PROBLEM != 2
void init()
{
    ParentType::init();

#if PROBLEM == 0
    Scalar vTot = 3e-6;
    analyticSolution_.initialize(vTot);
#endif
#if PROBLEM == 1
    analyticSolution_.initialize();
#endif
}

void addOutputVtkFields()
{
    analyticSolution_.calculateAnalyticSolution();

    ParentType::addOutputVtkFields();
    analyticSolution_.addOutputVtkFields(this->resultWriter());
}
#endif

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
    if (ParameterTree::tree().hasKey("Problem.OutputFileName"))
    {
        std::string fileName(GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, std::string, Problem, OutputFileName));
        return fileName.c_str();
    }
    else
    {
        return "test_mpfa2p";
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
#if PROBLEM == 2
    if (isInlet(globalPos))
    {
        bcTypes.setNeumann(eqIdxPress);
        bcTypes.setDirichlet(swIdx);
    }
    else if (isBottom(globalPos) || isTop(globalPos))
    {
        bcTypes.setAllNeumann();
    }
    else
    {
        bcTypes.setDirichlet(pwIdx);
        bcTypes.setOutflow(eqIdxSat);
    }
#elif  PROBLEM == 0
    if (globalPos[0] < eps_)
    {
        bcTypes.setAllDirichlet();
    }
    else if (globalPos[0] > this->bBoxMax()[0] - eps_)
    {
        bcTypes.setNeumann(eqIdxPress);
        bcTypes.setOutflow(eqIdxSat);
    }
    else
    {
        bcTypes.setAllNeumann();
    }
#else
    if (globalPos[0] < eps_)
    {
        bcTypes.setAllDirichlet();
    }
    else
    {
        bcTypes.setAllNeumann();
    }
#endif
}

//! set dirichlet condition  (pressure [Pa], saturation [-])
void dirichletAtPos(PrimaryVariables &values, const GlobalPosition& globalPos) const
{
    values = 0;
#if PROBLEM == 2
    Scalar pRef = referencePressureAtPos(globalPos);
    Scalar temp = temperatureAtPos(globalPos);

    values[pwIdx] = (1e5 - (this->bBoxMax()- globalPos) * this->gravity() * WettingPhase::density(temp, pRef));
    values[swIdx] = 1.0;

    if (isInlet(globalPos))
    {
        values[swIdx] = 0.0;
    }
#elif PROBLEM == 0
    if (globalPos[0] < eps_)
    {
        values[swIdx] = 0.8;
        values[pwIdx] = 1;
    }
#else
    if (globalPos[0] < eps_)
    {
        values[swIdx] = 1.0;
        values[pnIdx] = 1e5;
    }
#endif
}

//! set neumann condition for phases (flux, [kg/(m^2 s)])
void neumannAtPos(PrimaryVariables &values, const GlobalPosition& globalPos) const
{
    values = 0;
#if PROBLEM == 2
    if (isInlet(globalPos))
    {
        values[nPhaseIdx] = -inFlux_;
    }
#elif PROBLEM == 0
    if (globalPos[0] > this->bBoxMax()[0] - eps_)
    {
        values[nPhaseIdx] = 3e-3;
    }
#endif
}

//! return initial solution -> only saturation values have to be given!
void initialAtPos(PrimaryVariables &values,
        const GlobalPosition& globalPos) const
{
#if PROBLEM == 2
    values[pwIdx] = 0;
    values[swIdx] = 1.0;
#elif PROBLEM == 0
    values[pwIdx] = 0;
    values[swIdx] = 0.2;
#else
    values[pnIdx] = 0;
    values[swIdx] = 0.0;
#endif
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
        if (globalPos[1] > this->bBoxMax()[1] - eps_)
            return true;
    }
    if (dim == 3)
    {
        if (globalPos[2] > this->bBoxMax()[2] - eps_)
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
static constexpr Scalar eps_ = 1e-6;
#if PROBLEM == 0
BuckleyLeverettAnalytic<TypeTag> analyticSolution_;
#endif
#if PROBLEM == 1
McWhorterAnalytic<TypeTag> analyticSolution_;
#endif
};
} //end namespace

#endif // DUMUX_TEST_MPFA2P_PROBLEM_HH
