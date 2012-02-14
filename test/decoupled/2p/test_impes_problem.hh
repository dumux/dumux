// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   Copyright (C) 2007-2008 by Klaus Mosthaf                                *
 *   Copyright (C) 2007-2008 by Bernd Flemisch                               *
 *   Copyright (C) 2008-2009 by Andreas Lauser                               *
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
/*!
 * \file
 *
 * \brief test problem for the sequential 2p model
 */
#ifndef DUMUX_TEST_IMPES_PROBLEM_HH
#define DUMUX_TEST_IMPES_PROBLEM_HH

#if HAVE_UG
#include <dune/grid/uggrid.hh>
#endif

#include <dune/grid/yaspgrid.hh>
#include <dune/grid/sgrid.hh>
#include <dumux/common/structuredcubegridcreator.hh>

#include <dumux/material/fluidsystems/liquidphase.hh>
#include <dumux/material/components/simpleh2o.hh>
#include <dumux/material/components/simplednapl.hh>
#include <dumux/material/components/oil.hh>

#include <dumux/decoupled/2p/diffusion/fv/fvpressureproperties2p.hh>
#include <dumux/decoupled/2p/transport/fv/fvtransportproperties2p.hh>
#include <dumux/decoupled/2p/impes/impesproblem2p.hh>

//following includes are only needed if a global pressure formulation is chosen! Then only a total velocity can be reconstructed for the transport step
#include <dumux/decoupled/2p/transport/fv/capillarydiffusion.hh>
#include <dumux/decoupled/2p/transport/fv/gravitypart.hh>

#include "test_impes_spatialparams.hh"

#include<dumux/decoupled/2p/transport/fv/evalcflflux_coats.hh>

namespace Dumux
{

template<class TypeTag>
class IMPESTestProblem;

//////////
// Specify the properties
//////////
namespace Properties
{
NEW_TYPE_TAG(IMPESTestProblem, INHERITS_FROM(FVPressureTwoP, FVTransportTwoP, IMPESTwoP, TestIMPESSpatialParams));

// set the GridCreator property
SET_TYPE_PROP(IMPESTestProblem, GridCreator, StructuredCubeGridCreator<TypeTag>);

// Set the grid type
SET_PROP(IMPESTestProblem, Grid)
{
#if HAVE_UG
    typedef Dune::UGGrid<2> type;
#else
    typedef Dune::YaspGrid<2> type;
    //typedef Dune::SGrid<2, 2> type;
#endif
};

// Set the problem property
SET_TYPE_PROP(IMPESTestProblem, Problem, Dumux::IMPESTestProblem<TypeTag>);

////////////////////////////////////////////////////////////////////////
//Switch to a p_n-S_w formulation
//
//SET_INT_PROP(IMPESTestProblem, Formulation,
//        DecoupledTwoPCommonIndices::pnSn);
//
////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////
//Switch to a p_global-S_w formulation
//
//SET_INT_PROP(IMPESTestProblem, Formulation,
//        DecoupledTwoPCommonIndices::pGlobalSw);
//
//Define the capillary pressure term in the transport equation -> only needed in case of a p_global-S_w formulation!
//SET_TYPE_PROP(IMPESTestProblem, CapillaryFlux, CapillaryDiffusion<TypeTag>);
//
//Define the gravity term in the transport equation -> only needed in case of a p_global-S_w formulation!
//SET_TYPE_PROP(IMPESTestProblem, GravityFlux, GravityPart<TypeTag>);
//
////////////////////////////////////////////////////////////////////////

// Set the wetting phase
SET_PROP(IMPESTestProblem, WettingPhase)
{
private:
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
public:
    typedef Dumux::LiquidPhase<Scalar, Dumux::SimpleH2O<Scalar> > type;
};

// Set the non-wetting phase
SET_PROP(IMPESTestProblem, NonwettingPhase)
{
private:
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
public:
    typedef Dumux::LiquidPhase<Scalar, Dumux::SimpleH2O<Scalar> > type;
};

// Enable gravity
SET_BOOL_PROP(IMPESTestProblem, EnableGravity, false);

SET_TYPE_PROP(IMPESTestProblem, EvalCflFluxFunction, Dumux::EvalCflFluxCoats<TypeTag>);

SET_SCALAR_PROP(IMPESTestProblem, CFLFactor, 0.95);
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
 * <tt>./test_impes -parameterFile ./test_impes.input</tt>,
 * where the arguments define the parameter file..
 */
template<class TypeTag>
class IMPESTestProblem: public IMPESProblem2P<TypeTag>
{
typedef IMPESProblem2P<TypeTag> ParentType;
typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;

typedef typename GET_PROP_TYPE(TypeTag, Indices) Indices;

typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;
typedef typename GET_PROP_TYPE(TypeTag, FluidState) FluidState;

typedef typename GET_PROP_TYPE(TypeTag, WettingPhase) WettingPhase;

typedef typename GET_PROP_TYPE(TypeTag, TimeManager) TimeManager;

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

typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;

typedef typename GridView::Traits::template Codim<0>::Entity Element;
typedef typename GridView::Intersection Intersection;
typedef Dune::FieldVector<Scalar, dimWorld> GlobalPosition;

typedef typename GET_PROP_TYPE(TypeTag, BoundaryTypes) BoundaryTypes;
typedef typename GET_PROP(TypeTag, SolutionTypes) SolutionTypes;
    typedef typename SolutionTypes::PrimaryVariables PrimaryVariables;

public:
IMPESTestProblem(TimeManager &timeManager, const GridView &gridView) :
ParentType(timeManager, gridView), eps_(1e-6)
{}

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
    return "test_impes";
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

            values[pWIdx] = (2e5 + (this->bboxMax()[dim-1] - globalPos[dim-1]) * WettingPhase::density(temp, pRef) * this->gravity().two_norm());
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

const Scalar eps_;
};
} //end namespace

#endif
