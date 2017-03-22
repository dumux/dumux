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
 * \brief test problem for the sequential 2p model
 */
#ifndef DUMUX_TEST_IMPES_ADAPTIVE_PROBLEM_HH
#define DUMUX_TEST_IMPES_ADAPTIVE_PROBLEM_HH


#include <dumux/material/fluidsystems/liquidphase.hh>
#include <dumux/material/components/simpleh2o.hh>
#include <dumux/material/components/lnapl.hh>

#include <dumux/porousmediumflow/2p/sequential/impes/problem.hh>
#include <dumux/porousmediumflow/2p/sequential/diffusion/cellcentered/pressurepropertiesadaptive.hh>
#include <dumux/porousmediumflow/2p/sequential/transport/cellcentered/properties.hh>

#include "test_impesadaptivespatialparams.hh"

#include<dumux/porousmediumflow/2p/sequential/transport/cellcentered/evalcflfluxcoats.hh>

namespace Dumux
{

template<class TypeTag>
class TestIMPESAdaptiveProblem;

//////////
// Specify the properties
//////////
namespace Properties
{
NEW_TYPE_TAG(TestIMPESAdaptiveProblem, INHERITS_FROM(FVPressureTwoPAdaptive, FVTransportTwoP, IMPESTwoPAdaptive, TestIMPESAdaptiveSpatialParams));
NEW_TYPE_TAG(TestIMPESAdaptiveRestartProblem, INHERITS_FROM(TestIMPESAdaptiveProblem));

// Set the grid type
#if HAVE_DUNE_ALUGRID
SET_TYPE_PROP(TestIMPESAdaptiveProblem, Grid, Dune::ALUGrid<2, 2, Dune::cube, Dune::nonconforming>);
#endif

#if HAVE_DUNE_ALUGRID
// reset the GridCreator to the standard one for testing the restart functionality
SET_TYPE_PROP(TestIMPESAdaptiveRestartProblem, GridCreator, GridCreator<TypeTag>);
#endif

// Set the problem property
SET_TYPE_PROP(TestIMPESAdaptiveProblem, Problem, TestIMPESAdaptiveProblem<TypeTag>);

// Set the wetting phase
SET_PROP(TestIMPESAdaptiveProblem, WettingPhase)
{
private:
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
public:
    typedef FluidSystems::LiquidPhase<Scalar, SimpleH2O<Scalar> > type;
};

// Set the non-wetting phase
SET_PROP(TestIMPESAdaptiveProblem, NonwettingPhase)
{
private:
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
public:
    typedef FluidSystems::LiquidPhase<Scalar, SimpleH2O<Scalar> > type;
};

// Enable gravity
SET_BOOL_PROP(TestIMPESAdaptiveProblem, ProblemEnableGravity, false);

//SET_BOOL_PROP(TestIMPESAdaptiveProblem, EnableCompressibility, true);

//SET_TYPE_PROP(TestIMPESAdaptiveProblem, EvalCflFluxFunction, EvalCflFluxCoats<TypeTag>);

SET_SCALAR_PROP(TestIMPESAdaptiveProblem, ImpetCFLFactor, 0.95);
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
template<class TypeTag>
class TestIMPESAdaptiveProblem: public IMPESProblem2P<TypeTag>
{
    typedef IMPESProblem2P<TypeTag> ParentType;
    typedef typename GET_PROP_TYPE(TypeTag, Grid) Grid;
    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;

    typedef typename GET_PROP_TYPE(TypeTag, Indices) Indices;

    typedef typename GET_PROP_TYPE(TypeTag, WettingPhase) WettingPhase;

    typedef typename GET_PROP_TYPE(TypeTag, TimeManager) TimeManager;

    enum
    {
        dim = GridView::dimension, dimWorld = GridView::dimensionworld
    };

    enum
    {
        nPhaseIdx = Indices::nPhaseIdx,
        pwIdx = Indices::pwIdx,
        swIdx = Indices::swIdx,
        eqIdxPress = Indices::pressureEqIdx,
        eqIdxSat = Indices::satEqIdx
    };

    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;

    typedef typename GridView::Traits::template Codim<0>::Entity Element;
    typedef Dune::FieldVector<Scalar, dimWorld> GlobalPosition;

    typedef typename GET_PROP_TYPE(TypeTag, BoundaryTypes) BoundaryTypes;
    typedef typename GET_PROP(TypeTag, SolutionTypes) SolutionTypes;
    typedef typename SolutionTypes::PrimaryVariables PrimaryVariables;
    typedef typename GET_PROP_TYPE(TypeTag, GridCreator) GridCreator;

public:
    TestIMPESAdaptiveProblem(TimeManager &timeManager, const GridView &gridView) :
            ParentType(timeManager, gridView), eps_(1e-6)
    {
        name_ = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, std::string, Problem, Name);

        // Refine the grid provided that no restart occurs. Otherwise, an
        // already refined grid will be read.
        typedef typename GET_PROP(TypeTag, ParameterTree) ParameterTree;
        if (!(ParameterTree::tree().hasKey("Restart")
              || ParameterTree::tree().hasKey("TimeManager.Restart")))
        {
            GridCreator::grid().globalRefine(GET_PARAM_FROM_GROUP(TypeTag, int, GridAdapt, MaxLevel));
        }
        this->setGrid(GridCreator::grid());

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
    const std::string& name() const
    {
        return name_;
    }

    bool shouldWriteRestartFile() const
    {
        return true;
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

    void source(PrimaryVariables &values, const Element& element) const
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
        else if (globalPos[0] > this->bBoxMax()[0] - eps_)
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
            if (GET_PARAM_FROM_GROUP(TypeTag, bool, Problem, EnableGravity))
            {
                Scalar pRef = referencePressureAtPos(globalPos);
                Scalar temp = temperatureAtPos(globalPos);

                values[pwIdx] = (2e5 + (this->bBoxMax()[dim-1] - globalPos[dim-1]) * WettingPhase::density(temp, pRef) * this->gravity().two_norm());
            }
            else
            {
                values[pwIdx] = 2e5;
            }
            values[swIdx] = 0.8;
        }
        else
        {
            values[pwIdx] = 2e5;
            values[swIdx] = 0.2;
        }
    }

//! set neumann condition for phases (flux, [kg/(m^2 s)])
    void neumannAtPos(PrimaryVariables &values, const GlobalPosition& globalPos) const
    {
        values = 0;
        if (globalPos[0] > this->bBoxMax()[0] - eps_)
        {
            values[nPhaseIdx] = 3e-4;
        }
    }
//! return initial solution -> only saturation values have to be given!
    void initial(PrimaryVariables &values, const Element& element) const
    {
        values[pwIdx] = 0;
        values[swIdx] = 0.2;
    }

private:
    const Scalar eps_;
    std::string name_;
};
} //end namespace

#endif // DUMUX_TEST_IMPES_ADAPTIVE_PROBLEM_HH
