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
 *   along with this program. If not, see <http://www.gnu.org/licenses/>.    *
 *****************************************************************************/
/*!
 * \file
 *
 * \brief Soil contamination problem where DNAPL infiltrates a fully
 *        water saturated medium.
 */

#ifndef DUMUX_2PMINC_TEST_PROBLEM_HH
#define DUMUX_2PMINC_TEST_PROBLEM_HH

#if HAVE_UG
#include <dune/grid/uggrid.hh>
#endif

#include <dune/grid/yaspgrid.hh>

#include <dumux/porousmediumflow/2pminc/implicit/model.hh>
#include <dumux/porousmediumflow/2pminc/implicit/volumevariables.hh>
#include <dumux/porousmediumflow/implicit/problem.hh>
#include <dumux/material/components/simpleh2o.hh>
#include <dumux/material/components/dnapl.hh>
#include <dumux/io/cubegridcreator.hh>

#include "2pminctestspatialparams.hh"

namespace Dumux
{

template <typename TypeTag>
class TwoPMincTestProblem;

//////////
// Specify the properties for the lens problem
//////////
namespace Properties
{
NEW_TYPE_TAG(TwoPMincTestProblem, INHERITS_FROM(BoxTwoPMinc, TwoPMincSpatialParams));
NEW_TYPE_TAG(TwoPMincTestBoxProblem, INHERITS_FROM(BoxTwoPMinc, TwoPMincTestProblem));

// Set the grid type
#if HAVE_UG
SET_TYPE_PROP(TwoPMincTestProblem, Grid, Dune::UGGrid<2>);
#else
SET_TYPE_PROP(TwoPMincTestProblem, Grid, Dune::YaspGrid<2>);
#endif

// Set the problem property
SET_TYPE_PROP(TwoPMincTestProblem, Problem, Dumux::TwoPMincTestProblem<TypeTag>);

#define BUFFERSIZE 700
#define PROBLEM_OUTPUT_NAME "newMincProblem"
#define PROBLEM_OUTPUT_NAME_I PROBLEM_OUTPUT_NAME NUM_CONT_CHAR

// Set the wetting phase
SET_PROP(TwoPMincTestProblem, WettingPhase)
{
private:
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
public:
    typedef Dumux::LiquidPhase<Scalar, Dumux::SimpleH2O<Scalar> > type;
};

// Set the non-wetting phase
SET_PROP(TwoPMincTestProblem, NonwettingPhase)
{
private:
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
public:
    typedef Dumux::LiquidPhase<Scalar, Dumux::DNAPL<Scalar> > type;
};

// Set the grid creator
SET_TYPE_PROP(TwoPMincTestProblem, GridCreator, Dumux::CubeGridCreator<TypeTag>);

// Set number of Continua
SET_INT_PROP(TwoPMincTestProblem, NumContinua, 4);//here the number of continua can be set, default is 2
SET_INT_PROP(TwoPMincTestProblem, NumEq, GET_PROP_VALUE(TypeTag, NumContinua) * 2);

// Enable partial reassembly of the jacobian matrix?
SET_BOOL_PROP(TwoPMincTestProblem, ImplicitEnablePartialReassemble, true);

// Enable reuse of Jacobian matrices?
SET_BOOL_PROP(TwoPMincTestProblem, ImplicitEnableJacobianRecycling, true);

// Write the solutions of individual Newton iterations?
SET_BOOL_PROP(TwoPMincTestProblem, NewtonWriteConvergence, false);

// Use forward differences instead of central differences
SET_INT_PROP(TwoPMincTestProblem, ImplicitNumericDifferenceMethod, +1);

// Linear solver settings
SET_TYPE_PROP(TwoPMincTestProblem, LinearSolver, Dumux::ILU0BiCGSTABBackend<TypeTag> );
SET_INT_PROP(TwoPMincTestProblem, LinearSolverVerbosity, 0);
SET_INT_PROP(TwoPMincTestProblem, LinearSolverPreconditionerIterations, 1);
SET_SCALAR_PROP(TwoPMincTestProblem, LinearSolverPreconditionerRelaxation, 1.0);

// Enable gravity
SET_BOOL_PROP(TwoPMincTestProblem, ProblemEnableGravity, true);

NEW_PROP_TAG(BaseProblem);
SET_TYPE_PROP(TwoPMincTestBoxProblem, BaseProblem, ImplicitPorousMediaProblem<TypeTag>);
SET_TYPE_PROP(TwoPMincTestProblem, VolumeVariables, TwoPMincVolumeVariables<TypeTag>);

}

/*!
 * \ingroup TwoPBoxModel
 * \ingroup ImplicitTestProblems
 * \brief Soil contamination problem where DNAPL infiltrates a fully
 *        water saturated medium.
 *
 *  TODO reformulate the documentation according to the actual tested MINC problem
 *
 * The domain is sized 6m times 4m and features a rectangular lens
 * with low permeability which spans from (1 m , 2 m) to (4 m, 3 m)
 * and is surrounded by a medium with higher permeability. Note that
 * this problem is discretized using only two dimensions, so from the
 * point of view of the two-phase model, the depth of the domain
 * implicitly is 1 m everywhere.
 *
 * On the top and the bottom of the domain neumann boundary conditions
 * are used, while dirichlet conditions apply on the left and right
 * boundaries.
 *
 * DNAPL is injected at the top boundary from 3m to 4m at a rate of
 * 0.04 kg/(s m^2), the remaining neumann boundaries are no-flow
 * boundaries.
 *
 * The dirichlet boundaries on the left boundary is the hydrostatic
 * pressure scaled by a factor of 1.125, while on the right side it is
 * just the hydrostatic pressure. The DNAPL saturation on both sides
 * is zero.
 *
 * This problem uses the \ref TwoPBoxModel.
 *
 * This problem should typically be simulated until \f$t_{\text{end}}
 * \approx 20\,000\;s\f$ is reached. A good choice for the initial time step
 * size is \f$t_{\text{inital}} = 250\;s\f$.
 *
 * To run the simulation execute the following line in shell:
 * <tt>./test_box2p -parameterFile test_box2p.input</tt> or
 * <tt>./test_cc2p -parameterFile test_cc2p.input</tt>
 */
template <typename TypeTag >
class TwoPMincTestProblem : public GET_PROP_TYPE(TypeTag, BaseProblem)
{
    typedef typename GET_PROP_TYPE(TypeTag, BaseProblem) ParentType;
    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;

    typedef typename GET_PROP_TYPE(TypeTag, Indices) Indices;

    typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;
    typedef typename GET_PROP_TYPE(TypeTag, WettingPhase) WettingPhase;
    typedef typename GET_PROP_TYPE(TypeTag, NonwettingPhase) NonwettingPhase;

    enum {
        numEq = GET_PROP_VALUE(TypeTag, NumEq),
        numContinua = GET_PROP_VALUE(TypeTag, NumContinua),

        // primary variable indices
        pwIdx = Indices::pwIdx,
        snIdx = Indices::snIdx,

        // equation indices
        contiNEqIdx = Indices::contiNEqIdx,

        // phase indices
        wPhaseIdx = Indices::wPhaseIdx,
        nPhaseIdx = Indices::nPhaseIdx,

        //rock type index
        fractureIdx = 0,
        matrixIdx = 1,

        // world dimension
        dimWorld = GridView::dimensionworld
    };


    typedef typename GET_PROP_TYPE(TypeTag, PrimaryVariables) PrimaryVariables;
    typedef typename GET_PROP_TYPE(TypeTag, BoundaryTypes) BoundaryTypes;
    typedef typename GET_PROP_TYPE(TypeTag, TimeManager) TimeManager;
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef Dune::FieldVector<Scalar, dimWorld> GlobalPosition;
    typedef Dune::FieldVector<int, dimWorld> GridResolution;

public:
    /*!
     * \brief The constructor
     *
     * \param timeManager The time manager
     * \param gridView The grid view
     */
    TwoPMincTestProblem(TimeManager &timeManager,
                const GridView &gridView)
    : ParentType(timeManager, gridView),
      res_(0.0)//
    {
        eps_ = 3e-6;
        temperature_ = 273.15 + 20; // -> 20°C

        name_ = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag,
                std::string,
                Problem,
                Name);
        res_[0] = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag,
                int,
                Grid,
                NumberOfCellsX);
        res_[1] = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag,
                int,
                Grid,
                NumberOfCellsY);
        assert(res_[0] > 0);
        assert(res_[1] > 0);

        this->model().calculateMincGeometricParameters(res_, this->bBoxMin(), this->bBoxMax());
        this->getMincProblemParameters(gridView, res_);
    }

    /*!
     * \name Problem parameters
     */
    // \{

    void getMincProblemParameters(const GridView &gridView, const GridResolution &res)
    {
        distNestedContinua=this->model().getDistNestedContinua();
        volFraction=this->model().getVolFraction();
        interfaceArea=this->model().getInterfaceArea();
    }

    /*!
     * \brief The problem name.
     *
     * This is used as a prefix for files generated by the simulation.
     */
    const char *name() const
    {
        return name_.c_str();
    }

    /*!
     * \brief Called directly after the time integration.
     */
    void postTimeStep()
    {
        // Calculate storage terms
        PrimaryVariables storage;
        this->model().globalStorage(storage);

        // Write mass balance information for rank 0
        if (this->gridView().comm().rank() == 0) {
            std::cout<<"Storage: " << storage << std::endl;
        }
    }

    /*!
     * \brief Returns the temperature within the domain.
     *
     * This problem assumes a uniform temperature of 20 degrees Celsius.
     */
    Scalar temperature() const
    { return temperature_; }


    void sourceAtPos(PrimaryVariables &values,
                const GlobalPosition &globalPos) const
    {
        values = 0.0;
    }

    // \}

    /*!
     * \name Boundary conditions
     */
    // \{

    /*!
     * \brief Specifies which kind of boundary condition should be
     *        used for which equation on a given boundary control volume.
     *
     * \param values The boundary types for the conservation equations
     * \param globalPos The position of the center of the finite volume
     */
    void boundaryTypesAtPos(BoundaryTypes &values,
            const GlobalPosition &globalPos) const
    {
        if (onLeftBoundary_(globalPos) || onRightBoundary_(globalPos)) {
            values.setAllDirichlet();
        }
        else {
            values.setAllNeumann();
        }
    }

    /*!
     * \brief Evaluate the boundary conditions for a dirichlet
     *        control volume.
     *
     * \param values The dirichlet values for the primary variables
     * \param globalPos The center of the finite volume which ought to be set.
     *
     * For this method, the \a values parameter stores primary variables.
     */
    void dirichletAtPos(PrimaryVariables &values,
                        const GlobalPosition &globalPos) const
    {
        typename GET_PROP_TYPE(TypeTag, FluidState) fluidState;//TODO what is this?
        fluidState.setTemperature(temperature_);
        fluidState.setPressure(FluidSystem::wPhaseIdx, /*pressure=*/1e5);
        fluidState.setPressure(FluidSystem::nPhaseIdx, /*pressure=*/1e5);

        Scalar densityW = FluidSystem::density(fluidState, FluidSystem::wPhaseIdx);

        Scalar height = this->bBoxMax()[1] - this->bBoxMin()[1];
        Scalar depth = this->bBoxMax()[1] - globalPos[1];
        Scalar alpha = 1 + 1.5/height;
        Scalar width = this->bBoxMax()[0] - this->bBoxMin()[0];
        Scalar factor = (width*alpha + (1.0 - alpha)*globalPos[0])/width;

        // hydrostatic pressure scaled by alpha for fracture continua
        values[pwIdx] = 1e5 - factor*densityW*this->gravity()[1]*depth;
        values[snIdx] = 0.0;
    }

    /*!
     * \brief Evaluate the boundary conditions for a neumann
     *        boundary segment.
     *
     * \param values The neumann values for the conservation equations [kg / (m^2 *s )]
     * \param globalPos The position of the integration point of the boundary segment.
     *
     * For this method, the \a values parameter stores the mass flux
     * in normal direction of each phase. Negative values mean influx.
     */
    void neumannAtPos(PrimaryVariables &values,
                      const GlobalPosition &globalPos) const
    {
        values = 0.0;
        if (onInlet_(globalPos)) {
            values[contiNEqIdx] = -0.04; // kg / (m * s)
        }
    }
    // \}

    /*!
     * \name Volume terms
     * \{
     */


    /*!
     * \brief Evaluate the initial value for a control volume.
     *
     * \param values The initial values for the primary variables
     * \param globalPos The center of the finite volume which ought to be set.
     *
     * For this method, the \a values parameter stores primary
     * variables.
     */
    void initialAtPos(PrimaryVariables &values,
                      const GlobalPosition &globalPos) const
    {
        typename GET_PROP_TYPE(TypeTag, FluidState) fluidState;
        fluidState.setTemperature(temperature_);
        fluidState.setPressure(FluidSystem::wPhaseIdx, /*pressure=*/1e5);
        fluidState.setPressure(FluidSystem::nPhaseIdx, /*pressure=*/1e5);

        Scalar densityW = FluidSystem::density(fluidState, FluidSystem::wPhaseIdx);

        Scalar depth = this->bBoxMax()[1] - globalPos[1];

        // hydrostatic pressure
        values[pwIdx] = 1e5 - densityW*this->gravity()[1]*depth;
        values[snIdx] = 0.0;
        // Matrix Continua
        for (int nC=1; nC<numContinua; nC++)
        {
            values[2*matrixIdx * nC  + pwIdx] = 1e5 - densityW*this->gravity()[1]*depth;
            values[2*matrixIdx * nC  + snIdx] = 0.01;
        }
    }
    //! \}

private:

    bool onLeftBoundary_(const GlobalPosition &globalPos) const
    {
        return globalPos[0] < this->bBoxMin()[0] + eps_;
    }

    bool onRightBoundary_(const GlobalPosition &globalPos) const
    {
        return globalPos[0] > this->bBoxMax()[0] - eps_;
    }

    bool onLowerBoundary_(const GlobalPosition &globalPos) const
    {
        return globalPos[1] < this->bBoxMin()[1] + eps_;
    }

    bool onUpperBoundary_(const GlobalPosition &globalPos) const
    {
        return globalPos[1] > this->bBoxMax()[1] - eps_;
    }

    bool onInlet_(const GlobalPosition &globalPos) const
    {
        Scalar width = this->bBoxMax()[0] - this->bBoxMin()[0];
        Scalar lambda = (this->bBoxMax()[0] - globalPos[0])/width;
        return onUpperBoundary_(globalPos) && 0.5 < lambda && lambda < 2.0/3.0;
    }

    Scalar temperature_;
    Scalar eps_;
    std::string name_;

    GridResolution res_;
public:
    Dune::FieldVector <Scalar, numContinua> volFraction;
    Dune::FieldVector <Scalar, numContinua> interfaceArea;
    Dune::FieldVector <Scalar, numContinua> distNestedContinua; //distance between two nested continua d[nC+1]=pos[nC+1]-pos[nC]
};
} //end namespace

#endif
