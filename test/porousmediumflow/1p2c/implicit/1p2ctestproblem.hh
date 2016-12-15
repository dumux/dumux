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
/**
 * \file
 * \brief Definition of a problem, for the 1p2c problem:
 * Component transport of nitrogen dissolved in the water phase.
 */
#ifndef DUMUX_1P2C_TEST_PROBLEM_HH
#define DUMUX_1P2C_TEST_PROBLEM_HH

#include <dumux/implicit/box/properties.hh>
#include <dumux/implicit/cellcentered/tpfa/properties.hh>
#include <dumux/implicit/cellcentered/mpfa/properties.hh>
#include <dumux/porousmediumflow/1p2c/implicit/model.hh>
#include <dumux/porousmediumflow/implicit/problem.hh>

#include <dumux/material/fluidsystems/h2on2.hh>
#include "1p2ctestspatialparams.hh"

#define NONISOTHERMAL 0

namespace Dumux
{

template <class TypeTag>
class OnePTwoCTestProblem;

namespace Properties
{
#if NONISOTHERMAL
NEW_TYPE_TAG(OnePTwoCOutflowProblem, INHERITS_FROM(OnePTwoCNI));
NEW_TYPE_TAG(OnePTwoCOutflowBoxProblem, INHERITS_FROM(BoxModel, OnePTwoCOutflowProblem));
NEW_TYPE_TAG(OnePTwoCOutflowCCProblem, INHERITS_FROM(CCTpfaModel, OnePTwoCOutflowProblem));
#else
NEW_TYPE_TAG(OnePTwoCTestProblem, INHERITS_FROM(OnePTwoC));
NEW_TYPE_TAG(OnePTwoCTestBoxProblem, INHERITS_FROM(BoxModel, OnePTwoCTestProblem));
NEW_TYPE_TAG(OnePTwoCTestCCProblem, INHERITS_FROM(CCTpfaModel, OnePTwoCTestProblem));
NEW_TYPE_TAG(OnePTwoCTestCCMpfaProblem, INHERITS_FROM(CCMpfaModel, OnePTwoCTestProblem));
#endif

// Set the grid type
#if HAVE_UG
SET_TYPE_PROP(OnePTwoCTestProblem, Grid, Dune::UGGrid<2>);
#else
SET_TYPE_PROP(OnePTwoCTestProblem, Grid, Dune::YaspGrid<2>);
#endif

// Set the problem property
SET_TYPE_PROP(OnePTwoCTestProblem, Problem, OnePTwoCTestProblem<TypeTag>);

// Set fluid configuration
SET_TYPE_PROP(OnePTwoCTestProblem,
              FluidSystem,
              FluidSystems::H2ON2<typename GET_PROP_TYPE(TypeTag, Scalar), false>);

// Set the spatial parameters
SET_TYPE_PROP(OnePTwoCTestProblem, SpatialParams, OnePTwoCTestSpatialParams<TypeTag>);

// Define whether mole(true) or mass (false) fractions are used
SET_BOOL_PROP(OnePTwoCTestProblem, UseMoles, true);

// Enable velocity output
SET_BOOL_PROP(OnePTwoCTestProblem, VtkAddVelocity, true);

// Disable gravity
SET_BOOL_PROP(OnePTwoCTestProblem, ProblemEnableGravity, false);
}


/*!
 * \ingroup OnePTwoCModel
 * \ingroup ImplicitTestProblems
 *
 * \brief Definition of a problem, for the 1p2c problem:
 * Nitrogen is dissolved in the water phase and
 * is transported with the water flow from the left side to the right.
 *
 * The model domain is 1 m times 1 m with a discretization length of 0.05 m
 * and homogeneous soil properties (\f$ \mathrm{K=10e-10, \Phi=0.4, \tau=0.28}\f$).
 * Initially the domain is filled with pure water.
 *
 * At the left side, a Dirichlet condition defines a nitrogen mole fraction
 * of 0.3 mol/mol.
 * The water phase flows from the left side to the right due to the applied pressure
 * gradient of 1e5 Pa/m. The nitrogen is transported with the water flow
 * and leaves the domain at the right boundary
 * where again Dirichlet boundary conditions are applied. Here, the nitrogen mole
 * fraction is set to 0.0 mol/mol.
 *
 * The model is able to use either mole or mass fractions. The property useMoles can be set to either true or false in the
 * problem file. Make sure that the according units are used in the problem setup. The default setting for useMoles is true.
 *
 * This problem uses the \ref OnePTwoCModel model.
 *
 * To run the simulation execute the following line in shell:
 * <tt>./test_box1p2c -parameterFile ./test_box1p2c.input</tt> or
 * <tt>./test_cc1p2c -parameterFile ./test_cc1p2c.input</tt>
 */
template <class TypeTag>
class OnePTwoCTestProblem : public ImplicitPorousMediaProblem<TypeTag>
{
    using ParentType = ImplicitPorousMediaProblem<TypeTag>;

    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using Indices = typename GET_PROP_TYPE(TypeTag, Indices);
    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
    using FluidSystem = typename GET_PROP_TYPE(TypeTag, FluidSystem);
    using TimeManager = typename GET_PROP_TYPE(TypeTag, TimeManager);
    using BoundaryTypes = typename GET_PROP_TYPE(TypeTag, BoundaryTypes);
    using PrimaryVariables = typename GET_PROP_TYPE(TypeTag, PrimaryVariables);

    // copy some indices for convenience
    enum
    {
        // indices of the primary variables
        pressureIdx = Indices::pressureIdx,
        massOrMoleFracIdx = Indices::massOrMoleFracIdx,

        // indices for the non-isothermal case
#if NONISOTHERMAL
        temperatureIdx = Indices::temperatureIdx,
        energyEqIdx = Indices::energyEqIdx,
#endif

        // indices of the equations
        conti0EqIdx = Indices::conti0EqIdx,
        transportEqIdx = Indices::transportEqIdx
    };

    //! property that defines whether mole or mass fractions are used
    static const bool useMoles = GET_PROP_VALUE(TypeTag, UseMoles);

    static const int dimWorld = GridView::dimensionworld;
    using GlobalPosition = Dune::FieldVector<Scalar, dimWorld>;

public:
    OnePTwoCTestProblem(TimeManager &timeManager, const GridView &gridView)
        : ParentType(timeManager, gridView)
    {
        //initialize fluid system
        FluidSystem::init();

        name_ = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, std::string, Problem, Name);

        // stating in the console whether mole or mass fractions are used
        if(useMoles)
            std::cout<<"problem uses mole fractions"<<std::endl;
        else
            std::cout<<"problem uses mass fractions"<<std::endl;
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

#if !NONISOTHERMAL
    /*!
     * \brief Returns the temperature within the domain [K].
     *
     * This problem assumes a temperature of 20 degrees Celsius.
     */
    Scalar temperature() const
    { return 273.15 + 20; } // in [K]

#endif

    // \}

    /*!
     * \name Boundary conditions
     */
    // \{

    /*!
     * \brief Specifies which kind of boundary condition should be
     *        used for which equation on a given boundary segment.
     *
     * \param globalPos The position for which the bc type should be evaluated
     */
    BoundaryTypes boundaryTypesAtPos(const GlobalPosition &globalPos) const
    {
        BoundaryTypes values;

        if(globalPos[0] < eps_ || globalPos[0] > this->bBoxMax()[0] - eps_)
            values.setAllDirichlet();
        else
            values.setAllNeumann();

        return values;
    }

    /*!
     * \brief Evaluate the boundary conditions for a dirichlet
     *        boundary segment.
     *
     * \param globalPos The position for which the bc type should be evaluated
     */
    PrimaryVariables dirichletAtPos(const GlobalPosition &globalPos) const
    {
        PrimaryVariables values = initial_(globalPos);

        // condition for the N2 molefraction at left boundary
        if (globalPos[0] < eps_)
        {
            values[massOrMoleFracIdx] = 2.0e-5;
#if NONISOTHERMAL
            values[temperatureIdx] = 300.;
#endif
        }

        return values;
    }

    /*!
     * \brief Evaluate the boundary conditions for a Neumann
     *        boundary segment.
     *
     * For this method, the \a priVars parameter stores the mass flux
     * in normal direction of each component. Negative values mean
     * influx.
     *
     * The units must be according to either using mole or mass fractions. (mole/(m^2*s) or kg/(m^2*s))
     */
    PrimaryVariables neumannAtPos(const GlobalPosition& globalPos) const
    { return PrimaryVariables(0.0); }

    // \}

    /*!
     * \name Volume terms
     */
    // \{

    /*!
     * \brief Evaluate the source term for all phases within a given
     *        sub-control-volume.
     *
     * For this method, the \a priVars parameter stores the rate mass
     * of a component is generated or annihilate per volume
     * unit. Positive values mean that mass is created, negative ones
     * mean that it vanishes.
     *
     * The units must be according to either using mole or mass fractions. (mole/(m^3*s) or kg/(m^3*s))
     */
    PrimaryVariables sourceAtPos(const GlobalPosition &globalPos) const
    { return PrimaryVariables(0.0); }

    /*!
     * \brief Evaluate the initial value for a control volume.
     *
     * \param globalPos The position for which the initial condition should be evaluated
     *
     * For this method, the \a values parameter stores primary
     * variables.
     */
    PrimaryVariables initialAtPos(const GlobalPosition &globalPos) const
    { return initial_(globalPos); }

    // \}

private:
    // the internal method for the initial condition
    PrimaryVariables initial_(const GlobalPosition &globalPos) const
    {
        PrimaryVariables priVars;
        priVars[pressureIdx] = 2e5 - 1e5*globalPos[0]; // initial condition for the pressure
        priVars[massOrMoleFracIdx] = 0.0;  // initial condition for the N2 molefraction
#if NONISOTHERMAL
        priVars[temperatureIdx] = 290.;
#endif
        return priVars;
    }

    static constexpr Scalar eps_ = 1e-6;
    std::string name_;
};

} //end namespace Dumux

#endif
