// $Id$
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
#ifndef DUMUX_TEST_IMPES_PROBLEM_HH
#define DUMUX_TEST_IMPES_PROBLEM_HH

#if HAVE_UG
#include <dune/grid/uggrid.hh>
#endif

#include <dune/grid/yaspgrid.hh>
#include <dune/grid/sgrid.hh>

#include <dumux/material/fluidsystems/liquidphase.hh>
#include <dumux/material/components/simpleh2o.hh>
#include <dumux/material/components/oil.hh>

#include <dumux/decoupled/2p/impes/impesproblem2p.hh>
#include <dumux/decoupled/2p/diffusion/fv/fvvelocity2p.hh>
#include <dumux/decoupled/2p/diffusion/fvmpfa/fvmpfaovelocity2p.hh>
#include <dumux/decoupled/2p/transport/fv/fvsaturation2p.hh>
#include <dumux/decoupled/2p/transport/fv/capillarydiffusion.hh>
#include <dumux/decoupled/2p/transport/fv/gravitypart.hh>

#include "test_impes_spatialparams.hh"

//#include<dumux/decoupled/2p/transport/fv/evalcflflux_coats.hh>

namespace Dumux
{

template<class TypeTag>
class TestIMPESProblem;

//////////
// Specify the properties
//////////
namespace Properties
{
NEW_TYPE_TAG(IMPESTestProblem, INHERITS_FROM(DecoupledTwoP, MPFAProperties, Transport));

// Set the grid type
SET_PROP(IMPESTestProblem, Grid)
{
    //    typedef Dune::YaspGrid<2> type;
    typedef Dune::SGrid<2, 2> type;
};

// Set the problem property
SET_PROP(IMPESTestProblem, Problem)
{
public:
    typedef Dumux::TestIMPESProblem<TTAG(IMPESTestProblem)> type;
};

// Set the model properties
SET_PROP(IMPESTestProblem, TransportModel)
{
    typedef Dumux::FVSaturation2P<TTAG(IMPESTestProblem)> type;
};
SET_TYPE_PROP(IMPESTestProblem, DiffusivePart, Dumux::CapillaryDiffusion<TypeTag>);
SET_TYPE_PROP(IMPESTestProblem, ConvectivePart, Dumux::GravityPart<TypeTag>);

SET_PROP(IMPESTestProblem, PressureModel)
{
    typedef Dumux::FVVelocity2P<TTAG(IMPESTestProblem)> type;
//    typedef Dumux::FVMPFAOVelocity2P<TTAG(IMPESTestProblem)> type;
};

//SET_INT_PROP(IMPESTestProblem, VelocityFormulation,
//        GET_PROP_TYPE(TypeTag, PTAG(TwoPIndices))::velocityW);

//SET_INT_PROP(IMPESTestProblem, PressureFormulation,
//        GET_PROP_TYPE(TypeTag, PTAG(TwoPIndices))::pressureGlobal);

// Set the wetting phase
SET_PROP(IMPESTestProblem, WettingPhase)
{
private:
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(Scalar)) Scalar;
public:
    typedef Dumux::LiquidPhase<Scalar, Dumux::SimpleH2O<Scalar> > type;
};

// Set the non-wetting phase
SET_PROP(IMPESTestProblem, NonwettingPhase)
{
private:
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(Scalar)) Scalar;
public:
    typedef Dumux::LiquidPhase<Scalar, Dumux::SimpleH2O<Scalar> > type;
};

// Set the spatial parameters
SET_PROP(IMPESTestProblem, SpatialParameters)
{
private:
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(Grid)) Grid;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(Scalar)) Scalar;

public:
    typedef Dumux::TestIMPESSpatialParams<TypeTag> type;
};

// Enable gravity
SET_BOOL_PROP(IMPESTestProblem, EnableGravity, false);

//SET_TYPE_PROP(IMPESTestProblem, EvalCflFluxFunction, Dumux::EvalCflFluxCoats<TypeTag>);

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
 * <tt>./test_impes 1e8</tt>,
 * where the argument defines the simulation endtime.
 */
template<class TypeTag = TTAG(IMPESTestProblem)>
class TestIMPESProblem: public IMPESProblem2P<TypeTag, TestIMPESProblem<TypeTag> >
{
typedef TestIMPESProblem<TypeTag> ThisType;
typedef IMPESProblem2P<TypeTag, ThisType> ParentType;
typedef typename GET_PROP_TYPE(TypeTag, PTAG(GridView)) GridView;

typedef typename GET_PROP_TYPE(TypeTag, PTAG(TwoPIndices)) Indices;

typedef typename GET_PROP_TYPE(TypeTag, PTAG(FluidSystem)) FluidSystem;
typedef typename GET_PROP_TYPE(TypeTag, PTAG(FluidState)) FluidState;

enum
{
    dim = GridView::dimension, dimWorld = GridView::dimensionworld
};

enum
{
    wPhaseIdx = Indices::wPhaseIdx, nPhaseIdx = Indices::nPhaseIdx
};

typedef typename GET_PROP_TYPE(TypeTag, PTAG(Scalar)) Scalar;

typedef typename GridView::Traits::template Codim<0>::Entity Element;
typedef typename GridView::Intersection Intersection;
typedef Dune::FieldVector<Scalar, dimWorld> GlobalPosition;
typedef Dune::FieldVector<Scalar, dim> LocalPosition;

public:
TestIMPESProblem(const GridView &gridView, const GlobalPosition lowerLeft = 0, const GlobalPosition upperRight = 0) :
ParentType(gridView), lowerLeft_(lowerLeft), upperRight_(upperRight)
{
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
    return "test2p";
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
Scalar temperature(const GlobalPosition& globalPos, const Element& element) const
{
    return 273.15 + 10; // -> 10°C
}

// \}

//! Returns the reference pressure for evaluation of constitutive relations
Scalar referencePressure(const GlobalPosition& globalPos, const Element& element) const
{
    return 1e5; // -> 10°C
}

//!source term [kg/(m^3 s)]
std::vector<Scalar> source(const GlobalPosition& globalPos, const Element& element)
{
    return std::vector<Scalar>(2, 0.0);
}

/*!
* \brief Returns the type of boundary condition for the pressure equation.
*
* BC can be dirichlet (pressure) or neumann (flux).
*/
typename BoundaryConditions::Flags bctypePress(const GlobalPosition& globalPos, const Intersection& intersection) const
{
    if ((globalPos[0] < eps_))
    return BoundaryConditions::dirichlet;
    // all other boundaries
    return BoundaryConditions::neumann;
}

/*!
* \brief Returns the type of boundary condition for the saturation equation.
*
* BC can be dirichlet (saturation), neumann (flux), or outflow.
*/
BoundaryConditions::Flags bctypeSat(const GlobalPosition& globalPos, const Intersection& intersection) const
{
    if (globalPos[0] < eps_)
    return Dumux::BoundaryConditions::dirichlet;
    else if (globalPos[0] > upperRight_[0] - eps_)
        return Dumux::BoundaryConditions::outflow;
    else
    return Dumux::BoundaryConditions::neumann;
}

//! return dirichlet condition  (pressure, [Pa])
Scalar dirichletPress(const GlobalPosition& globalPos, const Intersection& intersection) const
{
    if (globalPos[0] < eps_)
    {
        if (GET_PROP_VALUE(TypeTag, PTAG(EnableGravity)))
        {
            const Element& element = *(intersection.inside());

            Scalar pRef = referencePressure(globalPos, element);
            Scalar temp = temperature(globalPos, element);
            Scalar sat = this->variables().satElement(element);

            FluidState fluidState;
            fluidState.update(sat, pRef, pRef, temp);
            return (2e5 + (upperRight_[dim-1] - globalPos[dim-1]) * FluidSystem::phaseDensity(wPhaseIdx, temp, pRef, fluidState) * this->gravity().two_norm());
        }
        else
        return 2e5;
    }
    // all other boundaries
    return 2e5;
}

//! return dirichlet condition  (saturation, [-])
Scalar dirichletSat(const GlobalPosition& globalPos, const Intersection& intersection) const
{
    if (globalPos[0] < eps_)
    return 0.8;
    // all other boundaries
    return 0.2;
}

//! return neumann condition  (flux, [kg/(m^2 s)])
std::vector<Scalar> neumann(const GlobalPosition& globalPos, const Intersection& intersection) const
{
    std::vector<Scalar> neumannFlux(2, 0.0);
    if (globalPos[0] > upperRight_[0] - eps_)
    {
        neumannFlux[nPhaseIdx] = 3e-4;
    }
    return neumannFlux;
}

//! initial condition for saturation
Scalar initSat(const GlobalPosition& globalPos, const Element& element) const
{
    return 0.2;
}

private:
GlobalPosition lowerLeft_;
GlobalPosition upperRight_;

static const Scalar eps_ = 1e-6;
};
} //end namespace

#endif
