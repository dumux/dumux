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
NEW_TYPE_TAG(TransportTestProblem, INHERITS_FROM(DecoupledTwoP, Transport, TestTransportSpatialParams));


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

SET_INT_PROP(TransportTestProblem, VelocityFormulation,
        DecoupledTwoPCommonIndices::velocityTotal);

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
class TestTransportProblem: public TransportProblem2P<TypeTag>
{
    typedef TransportProblem2P<TypeTag> ParentType;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(GridView)) GridView;

    typedef typename GET_PROP_TYPE(TypeTag, PTAG(TwoPIndices)) Indices;

    typedef typename GET_PROP_TYPE(TypeTag, PTAG(FluidSystem)) FluidSystem;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(FluidState)) FluidState;

    typedef typename GET_PROP_TYPE(TypeTag, PTAG(TimeManager)) TimeManager;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(BoundaryTypes)) BoundaryTypes;
    typedef typename GET_PROP(TypeTag, PTAG(SolutionTypes))::PrimaryVariables PrimaryVariables;

    enum
    {
        dim = GridView::dimension, dimWorld = GridView::dimensionworld
    };

    enum
    {
        wPhaseIdx = Indices::wPhaseIdx, nPhaseIdx = Indices::nPhaseIdx,
        satEqIdx = Indices::satEqIdx
    };

    typedef typename GET_PROP_TYPE(TypeTag, PTAG(Scalar)) Scalar;

    typedef typename GridView::Traits::template Codim<0>::Entity Element;
    typedef typename GridView::Intersection Intersection;
    typedef Dune::FieldVector<Scalar, dimWorld> GlobalPosition;

public:
    TestTransportProblem(TimeManager &timeManager, const GridView &gridView) :
        ParentType(timeManager, gridView)
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

    void sourceAtPos(PrimaryVariables &values,const GlobalPosition& globalPos) const
    {
        values = 0;
    }

    /*!
    * \brief Returns the type of boundary condition.
    *
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
                bcTypes.setAllOutflow();
            }
            // all other boundaries
            else
            {
                bcTypes.setAllNeumann();
            }
    }

    //! set dirichlet condition  (saturation [-])
    void dirichletAtPos(PrimaryVariables &values, const GlobalPosition& globalPos) const
    {
        values = 0;
        if (globalPos[0] < eps_)
        {
            values = 1.0;
        }
    }

    //! set neumann condition for phases (flux, [kg/(m^2 s)])
    void neumannAtPos(PrimaryVariables &values, const GlobalPosition& globalPos) const
    {
        values = 0;
    }

    //! return initial solution
    void initialAtPos(PrimaryVariables &values,
            const GlobalPosition &globalPos) const
    {
        values = 0;
    }

private:
    static constexpr Scalar eps_ = 1e-6;
};
} //end namespace

#endif
