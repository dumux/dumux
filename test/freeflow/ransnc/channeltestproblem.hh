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
 * \ingroup RANSNCTests
 * \brief Channel flow test for the multi-component staggered grid Reynolds-averaged Navier-Stokes model
 */
#ifndef DUMUX_RANS_NC_TEST_PROBLEM_HH
#define DUMUX_RANS_NC_TEST_PROBLEM_HH

#include <dumux/material/fluidsystems/h2oair.hh>

#include <dumux/discretization/staggered/freeflow/properties.hh>

#include <dumux/freeflow/rans/zeroeq/problem.hh>
#include <dumux/freeflow/ransnc/model.hh>

namespace Dumux
{
template <class TypeTag>
class ChannelNCTestProblem;

namespace Properties
{

#if !NONISOTHERMAL
NEW_TYPE_TAG(ChannelNCTestTypeTag, INHERITS_FROM(StaggeredFreeFlowModel, ZeroEqNC));
#else
NEW_TYPE_TAG(ChannelNCTestTypeTag, INHERITS_FROM(StaggeredFreeFlowModel, ZeroEqNCNI));
#endif

NEW_PROP_TAG(FluidSystem);

// Select the fluid system
SET_TYPE_PROP(ChannelNCTestTypeTag, FluidSystem,
              FluidSystems::H2OAir<typename GET_PROP_TYPE(TypeTag, Scalar)>);

SET_INT_PROP(ChannelNCTestTypeTag, PhaseIdx,
             GET_PROP_TYPE(TypeTag, FluidSystem)::wPhaseIdx);

SET_INT_PROP(ChannelNCTestTypeTag, ReplaceCompEqIdx, GET_PROP_VALUE(TypeTag, PhaseIdx));

// Set the grid type
SET_TYPE_PROP(ChannelNCTestTypeTag, Grid,
              Dune::YaspGrid<2, Dune::TensorProductCoordinates<typename GET_PROP_TYPE(TypeTag, Scalar), 2> >);

// Set the problem property
SET_TYPE_PROP(ChannelNCTestTypeTag, Problem, Dumux::ChannelNCTestProblem<TypeTag> );

SET_BOOL_PROP(ChannelNCTestTypeTag, EnableFVGridGeometryCache, true);

SET_BOOL_PROP(ChannelNCTestTypeTag, EnableGridFluxVariablesCache, true);
SET_BOOL_PROP(ChannelNCTestTypeTag, EnableGridVolumeVariablesCache, true);

// Enable gravity
SET_BOOL_PROP(ChannelNCTestTypeTag, UseMoles, true);
}

/*!
 * \ingroup RANSNCTests
 * \brief  Test problem for the one-phase model.
 *
 * Dry air is entering the channel, in 2-D a flat plate, from the left side.
 * In the middle of the inlet, water vapor is injected, which spreads by turbulent diffusion.
 * For the nonisothermal model the bottom has a constant temperature
 * which is \f$ \unit[30]{K} \f$ higher than the initial and inlet temperature.
 */
template <class TypeTag>
class ChannelNCTestProblem : public ZeroEqProblem<TypeTag>
{
    using ParentType = ZeroEqProblem<TypeTag>;
    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using FluidSystem = typename GET_PROP_TYPE(TypeTag, FluidSystem);

    // copy some indices for convenience
    using Indices = typename GET_PROP_TYPE(TypeTag, Indices);
    enum { dimWorld = GridView::dimensionworld };
    enum {
        totalMassBalanceIdx = Indices::totalMassBalanceIdx,
        transportEqIdx = 1-GET_PROP_VALUE(TypeTag, PhaseIdx),
        momentumBalanceIdx = Indices::momentumBalanceIdx,
        pressureIdx = Indices::pressureIdx,
        velocityXIdx = Indices::velocityXIdx,
        velocityYIdx = Indices::velocityYIdx,
#if NONISOTHERMAL
        temperatureIdx = Indices::temperatureIdx,
        energyBalanceIdx = Indices::energyBalanceIdx,
#endif
        transportCompIdx = 1-GET_PROP_VALUE(TypeTag, PhaseIdx)
    };

    using BoundaryTypes = typename GET_PROP_TYPE(TypeTag, BoundaryTypes);

    using FVGridGeometry = typename GET_PROP_TYPE(TypeTag, FVGridGeometry);

    using GlobalPosition = Dune::FieldVector<Scalar, dimWorld>;

    using PrimaryVariables = typename GET_PROP_TYPE(TypeTag, PrimaryVariables);
    using NumEqVector = typename GET_PROP_TYPE(TypeTag, NumEqVector);

    using TimeLoopPtr = std::shared_ptr<CheckPointTimeLoop<Scalar>>;
    using GridVariables = typename GET_PROP_TYPE(TypeTag, GridVariables);
    using SolutionVector = typename GET_PROP_TYPE(TypeTag, SolutionVector);

public:
    ChannelNCTestProblem(std::shared_ptr<const FVGridGeometry> fvGridGeometry)
    : ParentType(fvGridGeometry), eps_(1e-6)
    {
        inletVelocity_ = getParam<Scalar>("Problem.InletVelocity");
        FluidSystem::init();
    }

   /*!
     * \name Problem parameters
     */
    // \{

    bool shouldWriteRestartFile() const
    {
        return false;
    }

   /*!
     * \brief Return the temperature within the domain in [K].
     *
     * This problem assumes a temperature of 10 degrees Celsius.
     */
    Scalar temperature() const
    { return 273.15 + 10; } // 10C

   /*!
     * \brief Return the sources within the domain.
     *
     * \param globalPos The global position
     */
    NumEqVector sourceAtPos(const GlobalPosition &globalPos) const
    {
        return NumEqVector(0.0);
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
     * \param globalPos The position of the center of the finite volume
     */
    BoundaryTypes boundaryTypesAtPos(const GlobalPosition &globalPos) const
    {
        BoundaryTypes values;

        values.setDirichlet(momentumBalanceIdx);
        values.setOutflow(totalMassBalanceIdx);
        values.setOutflow(transportEqIdx);
#if NONISOTHERMAL
        values.setOutflow(energyBalanceIdx);
#endif

        if(isInlet(globalPos))
        {
            values.setDirichlet(momentumBalanceIdx);
            values.setOutflow(totalMassBalanceIdx);
            values.setDirichlet(transportEqIdx);
#if NONISOTHERMAL
            values.setDirichlet(energyBalanceIdx);
#endif
        }
        else if(isOutlet(globalPos))
        {
            values.setOutflow(momentumBalanceIdx);
            values.setDirichlet(totalMassBalanceIdx);
            values.setOutflow(transportEqIdx);
#if NONISOTHERMAL
            values.setOutflow(energyBalanceIdx);
#endif
        }
        else if (isOnWall(globalPos))
        {
            values.setDirichlet(velocityXIdx);
            values.setOutflow(velocityYIdx);
            values.setOutflow(totalMassBalanceIdx);
            values.setOutflow(transportEqIdx);
#if NONISOTHERMAL
            values.setOutflow(energyBalanceIdx);
#endif
        }
        else if (isSymmetry(globalPos))
        {
            values.isSymmetry();
        }

        return values;
    }

   /*!
     * \brief Evaluate the boundary conditions for a Dirichlet control volume.
     *
     * \param globalPos The center of the finite volume which ought to be set.
     */
    PrimaryVariables dirichletAtPos(const GlobalPosition &globalPos) const
    {
        PrimaryVariables values = initialAtPos(globalPos);

        if (isInlet(globalPos)
            && globalPos[1] > 0.4 * this->fvGridGeometry().bBoxMax()[1]
            && globalPos[1] < 0.6 * this->fvGridGeometry().bBoxMax()[1])
        {
            values[transportCompIdx] = 1e-3;
        }
#if NONISOTHERMAL
        if (isOnWall(globalPos))
        {
            values[temperatureIdx] += 30.;
        }
#endif

        return values;
    }

    // \}

   /*!
     * \name Volume terms
     */
    // \{

   /*!
     * \brief Evaluate the initial value for a control volume.
     *
     * \param globalPos The global position
     */
    PrimaryVariables initialAtPos(const GlobalPosition &globalPos) const
    {
        PrimaryVariables values(0.0);
        values[pressureIdx] = 1.0e+5;
        values[transportCompIdx] = 0.0;
#if NONISOTHERMAL
        values[temperatureIdx] = temperature();
#endif

        // block velocity profile
        values[velocityXIdx] = 0.0;
        if (!isOnWall(globalPos))
            values[velocityXIdx] =  inletVelocity_;
        values[velocityYIdx] = 0.0;

        return values;
    }

    // \}

    void setTimeLoop(TimeLoopPtr timeLoop)
    {
        timeLoop_ = timeLoop;
    }

    Scalar time() const
    {
        return timeLoop_->time();
    }

    bool isOnWall(const GlobalPosition& globalPos) const
    {
        return globalPos[1] < eps_;
    }

    bool isSymmetry(const GlobalPosition& globalPos) const
    {
        return globalPos[1] > this->fvGridGeometry().bBoxMax()[0] - eps_;
    }

private:

    bool isInlet(const GlobalPosition& globalPos) const
    {
        return globalPos[0] < eps_;
    }

    bool isOutlet(const GlobalPosition& globalPos) const
    {
        return globalPos[0] > this->fvGridGeometry().bBoxMax()[0] - eps_;
    }

    const Scalar eps_;
    Scalar inletVelocity_;
    TimeLoopPtr timeLoop_;
};
} //end namespace

#endif
