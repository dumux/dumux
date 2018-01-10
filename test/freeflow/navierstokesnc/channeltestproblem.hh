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
 * \ingroup NavierStokesNCTests
 * \brief Channel flow test for the multi-component staggered grid (Navier-)Stokes model
 */
#ifndef DUMUX_CHANNEL_NC_TEST_PROBLEM_HH
#define DUMUX_CHANNEL_NC_TEST_PROBLEM_HH

#include <dumux/material/fluidsystems/liquidphase.hh>
#include <dumux/material/components/simpleh2o.hh>
#include <dumux/material/components/constant.hh>
#include <dumux/material/fluidsystems/h2oair.hh>

#include <dumux/freeflow/navierstokes/problem.hh>
#include <dumux/discretization/staggered/freeflow/properties.hh>

#include <dumux/freeflow/navierstokesnc/model.hh>

namespace Dumux
{
template <class TypeTag>
class ChannelNCTestProblem;

namespace Properties
{

#if !NONISOTHERMAL
NEW_TYPE_TAG(ChannelNCTestProblem, INHERITS_FROM(StaggeredFreeFlowModel, NavierStokesNC));
#else
NEW_TYPE_TAG(ChannelNCTestProblem, INHERITS_FROM(StaggeredFreeFlowModel, NavierStokesNCNI));
#endif

NEW_PROP_TAG(FluidSystem);

// Select the fluid system
SET_TYPE_PROP(ChannelNCTestProblem, FluidSystem,
              FluidSystems::H2OAir<typename GET_PROP_TYPE(TypeTag, Scalar)/*, SimpleH2O<typename GET_PROP_TYPE(TypeTag, Scalar)>, true*/>);

SET_PROP(ChannelNCTestProblem, PhaseIdx)
{
private:
    using FluidSystem = typename GET_PROP_TYPE(TypeTag, FluidSystem);
public:
    static constexpr int value = FluidSystem::wPhaseIdx;
};

SET_INT_PROP(ChannelNCTestProblem, ReplaceCompEqIdx, 0);

// Set the grid type
SET_TYPE_PROP(ChannelNCTestProblem, Grid, Dune::YaspGrid<2>);

// Set the problem property
SET_TYPE_PROP(ChannelNCTestProblem, Problem, Dumux::ChannelNCTestProblem<TypeTag> );

SET_BOOL_PROP(ChannelNCTestProblem, EnableFVGridGeometryCache, true);

SET_BOOL_PROP(ChannelNCTestProblem, EnableGridFluxVariablesCache, true);
SET_BOOL_PROP(ChannelNCTestProblem, EnableGridVolumeVariablesCache, true);

// Enable gravity
SET_BOOL_PROP(ChannelNCTestProblem, UseMoles, true);

// #if ENABLE_NAVIERSTOKES
SET_BOOL_PROP(ChannelNCTestProblem, EnableInertiaTerms, true);
// #else
// SET_BOOL_PROP(ChannelNCTestProblem, EnableInertiaTerms, false);
// #endif
}

/*!
 * \ingroup NavierStokesNCTests
 * \brief  Test problem for the one-phase model.
 * \todo doc me!
 */
template <class TypeTag>
class ChannelNCTestProblem : public NavierStokesProblem<TypeTag>
{
    using ParentType = NavierStokesProblem<TypeTag>;
    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using FluidSystem = typename GET_PROP_TYPE(TypeTag, FluidSystem);

    // copy some indices for convenience
    using Indices = typename GET_PROP_TYPE(TypeTag, Indices);
    enum { dimWorld = GridView::dimensionworld };
    enum {
        massBalanceIdx = Indices::massBalanceIdx,
        transportEqIdx = 1,
        momentumBalanceIdx = Indices::momentumBalanceIdx,
        pressureIdx = Indices::pressureIdx,
        velocityXIdx = Indices::velocityXIdx,
        velocityYIdx = Indices::velocityYIdx,
#if NONISOTHERMAL
        temperatureIdx = Indices::temperatureIdx,
        energyBalanceIdx = Indices::energyBalanceIdx,
#endif
        transportCompIdx = 1/*FluidSystem::wCompIdx*/
    };

    using BoundaryTypes = typename GET_PROP_TYPE(TypeTag, BoundaryTypes);

    using FVGridGeometry = typename GET_PROP_TYPE(TypeTag, FVGridGeometry);

    using GlobalPosition = Dune::FieldVector<Scalar, dimWorld>;

    using PrimaryVariables = typename GET_PROP_TYPE(TypeTag, PrimaryVariables);
    using SourceValues = typename GET_PROP_TYPE(TypeTag, NumEqVector);

    using TimeLoopPtr = std::shared_ptr<CheckPointTimeLoop<Scalar>>;
    using GridVariables = typename GET_PROP_TYPE(TypeTag, GridVariables);
    using SolutionVector = typename GET_PROP_TYPE(TypeTag, SolutionVector);

public:
    ChannelNCTestProblem(std::shared_ptr<const FVGridGeometry> fvGridGeometry)
    : ParentType(fvGridGeometry), eps_(1e-6)
    {
        inletVelocity_ = getParam<Scalar>("Problem.InletVelocity");
        FluidSystem::init();
        deltaP_.resize(this->fvGridGeometry().numCellCenterDofs());
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
    SourceValues sourceAtPos(const GlobalPosition &globalPos) const
    {
        return SourceValues(0.0);
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

        if(isInlet(globalPos))
        {
            values.setDirichlet(momentumBalanceIdx);
            values.setOutflow(massBalanceIdx);
            values.setDirichlet(transportEqIdx);
#if NONISOTHERMAL
            values.setDirichlet(energyBalanceIdx);
#endif
        }
        else if(isOutlet(globalPos))
        {
            values.setOutflow(momentumBalanceIdx);
            values.setDirichlet(massBalanceIdx);
            values.setOutflow(transportEqIdx);
#if NONISOTHERMAL
            values.setOutflow(energyBalanceIdx);
#endif
        }

        else
        {
            // set Dirichlet values for the velocity everywhere
            values.setDirichlet(momentumBalanceIdx);
            values.setOutflow(massBalanceIdx);
            values.setOutflow(transportEqIdx);
#if NONISOTHERMAL
            values.setOutflow(energyBalanceIdx);
#endif
        }

        return values;
    }

   /*!
     * \brief Evaluate the boundary conditions for a dirichlet
     *        control volume.
     *
     * \param globalPos The center of the finite volume which ought to be set.
     */
    PrimaryVariables dirichletAtPos(const GlobalPosition &globalPos) const
    {
        PrimaryVariables values = initialAtPos(globalPos);

        // give the system some time so that the pressure can equilibrate, then start the injection of the tracer
        if(isInlet(globalPos))
        {
            if(time() >= 10.0 || inletVelocity_  < eps_)
            {
                values[transportCompIdx] = 1e-3;
#if NONISOTHERMAL
            values[temperatureIdx] = 293.15;
#endif
            }
        }

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
        PrimaryVariables values;
        values[pressureIdx] = 1.1e+5;
        values[transportCompIdx] = 0.0;
#if NONISOTHERMAL
        values[temperatureIdx] = 283.15;
#endif

        // parabolic velocity profile
        values[velocityXIdx] =  inletVelocity_*(globalPos[1] - this->fvGridGeometry().bBoxMin()[1])*(this->fvGridGeometry().bBoxMax()[1] - globalPos[1])
                               / (0.25*(this->fvGridGeometry().bBoxMax()[1] - this->fvGridGeometry().bBoxMin()[1])*(this->fvGridGeometry().bBoxMax()[1] - this->fvGridGeometry().bBoxMin()[1]));

        values[velocityYIdx] = 0.0;

        return values;
    }

   /*!
     * \brief Adds additional VTK output data to the VTKWriter. Function is called by the output module on every write.
     *
     * \param gridVariables The grid variables
     * \param sol The solution vector
     */
    void calculateDeltaP(const GridVariables& gridVariables, const SolutionVector& sol)
    {
        for (const auto& element : elements(this->fvGridGeometry().gridView()))
        {
            auto fvGeometry = localView(this->fvGridGeometry());
            fvGeometry.bindElement(element);
            for (auto&& scv : scvs(fvGeometry))
            {
                auto ccDofIdx = scv.dofIndex();

                auto elemVolVars = localView(gridVariables.curGridVolVars());
                elemVolVars.bind(element, fvGeometry, sol);

                deltaP_[ccDofIdx] = elemVolVars[scv].pressure() - 1.1e5;
            }
        }
    }

    auto& getDeltaP() const
    { return deltaP_; }

    // \}

    void setTimeLoop(TimeLoopPtr timeLoop)
    {
        timeLoop_ = timeLoop;
    }

    Scalar time() const
    {
        return timeLoop_->time();
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

    bool isWall(const GlobalPosition& globalPos) const
    {
        return globalPos[0] > eps_ || globalPos[0] < this->fvGridGeometry().bBoxMax()[0] - eps_;
    }

    const Scalar eps_;
    Scalar inletVelocity_;
    TimeLoopPtr timeLoop_;
    std::vector<Scalar> deltaP_;
};
} //end namespace

#endif
