// $Id: test_2p_problem.hh 3783 2010-06-24 11:33:53Z bernd $
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
 * \brief test problem for the explicit transport model
 */
#ifndef DUMUX_TEST_TRANSPORT_PROBLEM_HH
#define DUMUX_TEST_TRANSPORT_PROBLEM_HH

#include <dune/grid/io/file/dgfparser/dgfyasp.hh>

#include <dumux/material/fluidsystems/liquidphase.hh>
#include <dumux/material/components/unit.hh>

#include <dumux/decoupled/2p/transport/transportproblem2p.hh>
#include <dumux/decoupled/2p/diffusion/fv/fvvelocity2p.hh>
#include <dumux/decoupled/2p/transport/fv/fvsaturation2p.hh>

#include "test_transport_spatialparams.hh"

namespace Dumux
{

template<class TypeTag>
class TestTransportProblem;

//////////
// Specify the properties
//////////
namespace Properties
{
NEW_TYPE_TAG(TransportTestProblem, INHERITS_FROM(DecoupledModel, Transport));


// Set the grid type
SET_PROP(TransportTestProblem, Grid)
{
    typedef Dune::YaspGrid<2> type;
};

// Set the problem property
SET_PROP(TransportTestProblem, Problem)
{
public:
    typedef Dumux::TestTransportProblem<TTAG(TransportTestProblem)> type;
};

// Set the model properties
SET_PROP(TransportTestProblem, Model)
{
    typedef Dumux::FVSaturation2P<TTAG(TransportTestProblem)> type;
};

// Set the wetting phase
SET_PROP(TransportTestProblem, WettingPhase)
{
private:
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(Scalar)) Scalar;
public:
    typedef Dumux::LiquidPhase<Scalar, Dumux::Unit<Scalar> > type;
};

// Set the non-wetting phase
SET_PROP(TransportTestProblem, NonwettingPhase)
{
private:
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(Scalar)) Scalar;
public:
    typedef Dumux::LiquidPhase<Scalar, Dumux::Unit<Scalar> > type;
};

// Set the spatial parameters
SET_PROP(TransportTestProblem, SpatialParameters)
{
private:
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(Grid)) Grid;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(Scalar)) Scalar;

public:
    typedef Dumux::TestTransportSpatialParams<TypeTag> type;
};

// Disable gravity
SET_BOOL_PROP(TransportTestProblem, EnableGravity, false);

SET_SCALAR_PROP(TransportTestProblem, CFLFactor, 1.0);
}

/*!
 * \ingroup IMPETtests
 *
 * \brief test problem for the explicit transport model
 *
 * A unit "fluid" is injected from the left side into a rectangular 2D
 * domain also this testing fluid. Upper and lower boundary are closed (Neumann = 0),
 * and there is free outflow on the right side.
 *
 * This test solely applies the 2p transport on a given velocity field, without a
 * pressure field being solved.
 *
 * To run the simulation execute the following line in shell:
 * <tt>./test_transport grids/test_transport.dgf 16000</tt>,
 * where the argument defines the simulation endtime.
 */
template<class TypeTag = TTAG(TransportTestProblem)>
class TestTransportProblem: public TransportProblem2P<TypeTag, TestTransportProblem<TypeTag> >
{
    typedef TestTransportProblem<TypeTag> ThisType;
    typedef TransportProblem2P<TypeTag, ThisType> ParentType;
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
    TestTransportProblem(const GridView &gridView, const GlobalPosition lowerLeft = 0, const GlobalPosition upperRight = 0) :
        ParentType(gridView), lowerLeft_(lowerLeft), upperRight_(upperRight)
    {
        GlobalPosition vel(0);
        vel[0] = 1e-5;
        this->variables().velocity() = vel;
        this->variables().initializePotentials(vel);
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
        return "test_transport";
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


    Scalar referencePressure(const GlobalPosition& globalPos, const Element& element) const
    {
        return 1e5; // -> 10°C
    }

    std::vector<Scalar> source(const GlobalPosition& globalPos, const Element& element)
    {
        return std::vector<Scalar>(2, 0.0);
    }

    BoundaryConditions::Flags bctypeSat(const GlobalPosition& globalPos, const Intersection& intersection) const
    {
        if (globalPos[0] < eps_)
            return Dumux::BoundaryConditions::dirichlet;
        else if (globalPos[0] > upperRight_[0] - eps_)
            return Dumux::BoundaryConditions::outflow;
        else
            return Dumux::BoundaryConditions::neumann;
    }

    Scalar dirichletSat(const GlobalPosition& globalPos, const Intersection& intersection) const
    {
        if (globalPos[0] < eps_)
            return 1.0;
        // all other boundaries
        return 0.0;
    }

    std::vector<Scalar> neumann(const GlobalPosition& globalPos, const Intersection& intersection) const
    {
        return std::vector<Scalar>(2, 0.0);
    }

    Scalar initSat(const GlobalPosition& globalPos, const Element& element) const
    {
        return 0.0;
    }

private:
    GlobalPosition lowerLeft_;
    GlobalPosition upperRight_;

    static const Scalar eps_ = 1e-6;
};
} //end namespace

#endif
