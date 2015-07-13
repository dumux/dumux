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
 * \brief test problem for sequential 2p models in 3d
 */

#ifndef DUMUX_TEST_3D2P_PROBLEM_HH
#define DUMUX_TEST_3D2P_PROBLEM_HH

#if HAVE_ALUGRID
#include <dune/grid/alugrid/3d/alugrid.hh>
#elif HAVE_DUNE_ALUGRID
#include <dune/alugrid/grid.hh>
#endif

#include <dumux/io/cubegridcreator.hh>

#include <dumux/material/fluidsystems/liquidphase.hh>
#include <dumux/material/fluidsystems/gasphase.hh>
#include <dumux/material/components/simpleh2o.hh>

#include <dumux/decoupled/2p/diffusion/fvmpfa/lmethod/fvmpfal3dpressureproperties2p.hh>
#include <dumux/decoupled/2p/diffusion/fvmpfa/lmethod/fvmpfal3dpressureproperties2padaptive.hh>
#include <dumux/decoupled/2p/diffusion/fv/fvpressureproperties2p.hh>
#include <dumux/decoupled/2p/diffusion/fv/fvpressureproperties2padaptive.hh>
#include <dumux/decoupled/2p/diffusion/mimetic/mimeticpressureproperties2p.hh>
#include <dumux/decoupled/2p/diffusion/mimetic/mimeticpressureproperties2padaptive.hh>
#include <dumux/decoupled/2p/impes/impesproblem2p.hh>
#include <dumux/decoupled/2p/transport/fv/fvtransportproperties2p.hh>

#include "test_3d2pspatialparams.hh"

#include <dumux/decoupled/2p/transport/fv/evalcflfluxcoats.hh>
#include <dumux/decoupled/2p/impes/gridadaptionindicator2plocal.hh>


namespace Dumux
{
template<class TypeTag>
class Test3D2PProblem;

//////////
// Specify the properties
//////////
namespace Properties
{
NEW_TYPE_TAG(ThreeDTwoPTestProblem, INHERITS_FROM(Test3d2pSpatialParams));

// Set the grid type
#if HAVE_ALUGRID || HAVE_DUNE_ALUGRID
SET_TYPE_PROP(ThreeDTwoPTestProblem, Grid, Dune::ALUGrid<3, 3, Dune::cube, Dune::nonconforming>);
#endif

// Set the problem property
SET_TYPE_PROP(ThreeDTwoPTestProblem, Problem, Dumux::Test3D2PProblem<TypeTag>);

// Set the wetting phase
SET_PROP(ThreeDTwoPTestProblem, WettingPhase)
{
private:
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
public:
    typedef Dumux::LiquidPhase<Scalar, Dumux::SimpleH2O<Scalar> > type;
};

// Set the non-wetting phase
SET_PROP(ThreeDTwoPTestProblem, NonwettingPhase)
{
private:
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
public:
    typedef Dumux::LiquidPhase<Scalar, Dumux::SimpleH2O<Scalar> > type;
};

#if PROBLEM == 1
SET_INT_PROP(ThreeDTwoPTestProblem, Formulation, DecoupledTwoPCommonIndices::pnSw);
#endif

// Set the spatial parameters
SET_PROP(ThreeDTwoPTestProblem, SpatialParams)
{
private:
    typedef typename GET_PROP_TYPE(TypeTag, Grid) Grid;
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;

public:
    typedef Dumux::Test3d2pSpatialParams<TypeTag> type;
};

#if PROBLEM == 1
SET_TYPE_PROP(ThreeDTwoPTestProblem, EvalCflFluxFunction, Dumux::EvalCflFluxCoats<TypeTag>);
SET_SCALAR_PROP(ThreeDTwoPTestProblem, ImpetCFLFactor, 1.0);
#else
SET_SCALAR_PROP(ThreeDTwoPTestProblem, ImpetCFLFactor, 0.95);
#endif

SET_TYPE_PROP(ThreeDTwoPTestProblem, AdaptionIndicator, Dumux::GridAdaptionIndicator2PLocal<TypeTag>);

NEW_TYPE_TAG(FVTwoPTestProblem, INHERITS_FROM(FVPressureTwoP, FVTransportTwoP, IMPESTwoP, ThreeDTwoPTestProblem));
NEW_TYPE_TAG(FVAdaptiveTwoPTestProblem, INHERITS_FROM(FVPressureTwoPAdaptive, FVTransportTwoP, IMPESTwoPAdaptive, ThreeDTwoPTestProblem));
NEW_TYPE_TAG(MPFALTwoPTestProblem, INHERITS_FROM(FvMpfaL3dPressureTwoP, FVTransportTwoP, IMPESTwoP, ThreeDTwoPTestProblem));
NEW_TYPE_TAG(MPFALAdaptiveTwoPTestProblem, INHERITS_FROM(FvMpfaL3dPressureTwoPAdaptive, FVTransportTwoP, IMPESTwoPAdaptive, ThreeDTwoPTestProblem));
NEW_TYPE_TAG(MimeticTwoPTestProblem, INHERITS_FROM(MimeticPressureTwoP, FVTransportTwoP, IMPESTwoP, ThreeDTwoPTestProblem));
NEW_TYPE_TAG(MimeticAdaptiveTwoPTestProblem, INHERITS_FROM(MimeticPressureTwoPAdaptive, FVTransportTwoP, IMPESTwoPAdaptive, ThreeDTwoPTestProblem));
}

/*!
 * \ingroup IMPETtests
 *
 * \brief test problem for the sequential 2p models in 3d
 *
 * Standard setting: one eighth of a nine-spot waterflood problem. The domain of size 1 x 1 x 1 m is initially saturated by a non-wetting fluid.
 * Water is injected in the corner of the origin (0,0,0) and non-wetting fluid produced in the upper corner at (1,1,1).
 *
 * To run the simulation execute the following line in shell:
 * <tt>./test_3d2p</tt>.
 */
template<class TypeTag>
class Test3D2PProblem: public IMPESProblem2P<TypeTag>
{
typedef Test3D2PProblem<TypeTag> ThisType;
typedef IMPESProblem2P<TypeTag> ParentType;
typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;

typedef typename GET_PROP_TYPE(TypeTag, Indices) Indices;

typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;
typedef typename GET_PROP_TYPE(TypeTag, FluidState) FluidState;

typedef typename GET_PROP_TYPE(TypeTag, TimeManager) TimeManager;

typedef typename GET_PROP_TYPE(TypeTag, BoundaryTypes) BoundaryTypes;
typedef typename GET_PROP(TypeTag, SolutionTypes)::PrimaryVariables PrimaryVariables;

enum
{
    dim = GridView::dimension, dimWorld = GridView::dimensionworld
};

enum
{
    wPhaseIdx = Indices::wPhaseIdx,
    nPhaseIdx = Indices::nPhaseIdx,
#if PROBLEM == 1
    pNIdx = Indices::pnIdx,
#else
    pWIdx = Indices::pwIdx,
#endif
    swIdx = Indices::swIdx,
    pressureEqIdx = Indices::pressureEqIdx,
    satEqIdx = Indices::satEqIdx
};

typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;

typedef typename GridView::Traits::template Codim<0>::Entity Element;
typedef typename GridView::Traits::template Codim<0>::EntityPointer ElementPointer;
typedef typename GridView::Traits::template Codim<dim>::EntityPointer VertexPointer;
typedef typename GridView::Intersection Intersection;
typedef Dune::FieldVector<Scalar, dimWorld> GlobalPosition;
typedef Dune::FieldVector<Scalar, dim> LocalPosition;
typedef typename GridView::template Codim<dim>::Iterator VertexIterator;

typedef typename GET_PROP_TYPE(TypeTag, GridCreator) GridCreator;
typedef typename GET_PROP(TypeTag, ParameterTree) ParameterTree;

public:

Test3D2PProblem(TimeManager &timeManager,const GridView &gridView) :
ParentType(timeManager, gridView), inflowEdge_(0), outflowEdge_(0)
{
    this->setGrid(GridCreator::grid());

    int refinementFactor = 0;
    if (ParameterTree::tree().hasKey("Grid.RefinementFactor") && !GET_PROP_VALUE(TypeTag, AdaptiveGrid))
    {
        refinementFactor = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, Scalar, Grid, RefinementFactor);
        this->grid().globalRefine(refinementFactor);
    }

    Scalar minDist = this->bBoxMax().two_norm();
    Scalar maxDist = this->bBoxMin().two_norm();

    // calculate the bounding box of the grid view
    VertexIterator vIt = gridView.template begin<dim>();
    const VertexIterator vEndIt = gridView.template end<dim>();
    for (; vIt!=vEndIt; ++vIt) {
        GlobalPosition vertexCoord(vIt->geometry().center());

        Scalar dist = vertexCoord.two_norm();
        if (dist > maxDist)
        {
            maxDist = dist;
            outflowEdge_ = vIt->geometry().center();
        }
        if (dist < minDist)
        {
            minDist = dist;
            inflowEdge_ = vIt->geometry().center();
        }
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
    if (ParameterTree::tree().hasKey("Problem.OutputName"))
    {
        std::string name(GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, std::string, Problem, OutputName));
        return name.c_str();
    }
    else
        return "test_3d2p";
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

void sourceAtPos(PrimaryVariables &values,const GlobalPosition& globalPos) const
{
    values = 0;

#if PROBLEM == 2 //Nine-Spot
        if (globalPos[0] < 0.1 * outflowEdge_[0] + eps_ && globalPos[1] < 0.1 * outflowEdge_[1] && globalPos[2] < 0.1 * outflowEdge_[2])
        {
            values[wPhaseIdx] = 1;
        }
#endif
}

void source(PrimaryVariables &values,const Element& element) const
{
    values = 0;

#if PROBLEM == 2 //Nine-Spot
    int numVertices = element.geometry().corners();
    for(int i = 0; i < numVertices; i++)
    {
        GlobalPosition globalPos(element.template subEntity<dim>(i)->geometry().center());

        if (globalPos[0] < inflowEdge_[0] + eps_ && globalPos[1] < inflowEdge_[1] + eps_ && globalPos[2] < inflowEdge_[2] + eps_)
        {
            values[wPhaseIdx] = 1;
            break;
        }
    }
#endif
}

void boundaryTypesAtPos(BoundaryTypes &bcTypes, const GlobalPosition& globalPos) const
{
#if PROBLEM == 2 //Nine-Spot
        if (globalPos[0] > 0.9 * outflowEdge_[0]  && globalPos[1] > 0.9 * outflowEdge_[1] && globalPos[2] > 0.9 * outflowEdge_[2])
        {
            bcTypes.setAllDirichlet();
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
        else if (globalPos[0] > this->bBoxMax()[0] - eps_)
        {
            bcTypes.setNeumann(pressureEqIdx);
            bcTypes.setOutflow(satEqIdx);
        }
        // all other boundaries
        else
        {
            bcTypes.setAllNeumann();
        }
#endif
}

void boundaryTypes(BoundaryTypes &bcTypes, const Intersection& intersection) const
{
#if PROBLEM == 2 //Nine-Spot
    ElementPointer element = intersection.inside();
    int numVertices = element->geometry().corners();
    for(int i = 0; i < numVertices; i++)
    {
        GlobalPosition globalPos(element->template subEntity<dim>(i)->geometry().center());

        if (globalPos[0] > outflowEdge_[0] - eps_ && globalPos[1] > outflowEdge_[1] - eps_ && globalPos[2] > outflowEdge_[2] - eps_)
        {
            bcTypes.setAllDirichlet();
            break;
        }
        else
        {
            bcTypes.setAllNeumann();
        }
    }
#else
    GlobalPosition globalPos(intersection.geometry().center());

    boundaryTypesAtPos(bcTypes,globalPos);
#endif
}

void dirichletAtPos(PrimaryVariables &values, const GlobalPosition& globalPos) const
{
    values = 0;

#if PROBLEM == 2 //Nine-Spot
            values[swIdx] = 0.8;
            values[pWIdx] = 1;
#elif PROBLEM == 0
    if (globalPos[0] < eps_)
    {
        values[swIdx] = 0.8;
        typedef typename  GET_PROP_TYPE(TypeTag, NonwettingPhase) NonwettingPhase;
        values[pWIdx] = 1;
    }
#elif PROBLEM == 1
    if (globalPos[0] < eps_)
    {
        values[swIdx] = 1.0;
        values[pNIdx] = 1e5;
    }
#endif
}

void neumannAtPos(PrimaryVariables &values, const GlobalPosition& globalPos) const
{
    values = 0;

#if PROBLEM == 0 //Buckley-Leverett
    if (globalPos[0] > this->bBoxMax()[0] - eps_)
    {
        values[nPhaseIdx] = 3e-3;
    }
#endif
}

void initialAtPos(PrimaryVariables &values,
        const GlobalPosition& globalPos) const
{
#if PROBLEM == 1
    values[pNIdx] = 0;
    values[swIdx] = 0.0;
#else
    values[pWIdx] = 0;
    values[swIdx] = 0.2;
#endif
}

private:
static constexpr Scalar eps_ = 1e-6;
GlobalPosition inflowEdge_;
GlobalPosition outflowEdge_;
};
} //end namespace

#endif
