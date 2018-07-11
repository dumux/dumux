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

#include <dune/grid/yaspgrid.hh>

#include <dumux/material/components/simpleh2o.hh>
#include <dumux/material/fluidsystems/h2oair.hh>

#include <dumux/freeflow/navierstokes/problem.hh>
#include <dumux/discretization/staggered/freeflow/properties.hh>

#include <dumux/freeflow/compositional/navierstokesncmodel.hh>

namespace Dumux
{
template <class TypeTag>
class ChannelNCTestProblem;

namespace Properties
{

#if !NONISOTHERMAL
NEW_TYPE_TAG(ChannelNCTestTypeTag, INHERITS_FROM(StaggeredFreeFlowModel, NavierStokesNC));
#else
NEW_TYPE_TAG(ChannelNCTestTypeTag, INHERITS_FROM(StaggeredFreeFlowModel, NavierStokesNCNI));
#endif

NEW_PROP_TAG(FluidSystem);

// Select the fluid system
SET_TYPE_PROP(ChannelNCTestTypeTag, FluidSystem,
              FluidSystems::H2OAir<typename GET_PROP_TYPE(TypeTag, Scalar)>);

SET_INT_PROP(ChannelNCTestTypeTag, PhaseIdx,
             GET_PROP_TYPE(TypeTag, FluidSystem)::liquidPhaseIdx);

SET_INT_PROP(ChannelNCTestTypeTag, ReplaceCompEqIdx, GET_PROP_VALUE(TypeTag, PhaseIdx));

// Set the grid type
SET_TYPE_PROP(ChannelNCTestTypeTag, Grid, Dune::YaspGrid<2>);

// Set the problem property
SET_TYPE_PROP(ChannelNCTestTypeTag, Problem, Dumux::ChannelNCTestProblem<TypeTag> );

SET_BOOL_PROP(ChannelNCTestTypeTag, EnableFVGridGeometryCache, true);

SET_BOOL_PROP(ChannelNCTestTypeTag, EnableGridFluxVariablesCache, true);
SET_BOOL_PROP(ChannelNCTestTypeTag, EnableGridVolumeVariablesCache, true);

// Use mole fraction formulation
#if USE_MASS
SET_BOOL_PROP(ChannelNCTestTypeTag, UseMoles, false);
#else
SET_BOOL_PROP(ChannelNCTestTypeTag, UseMoles, true);
#endif

SET_BOOL_PROP(ChannelNCTestTypeTag, EnableInertiaTerms, true);
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

    using BoundaryTypes = typename GET_PROP_TYPE(TypeTag, BoundaryTypes);
    using FluidSystem = typename GET_PROP_TYPE(TypeTag, FluidSystem);
    using FVGridGeometry = typename GET_PROP_TYPE(TypeTag, FVGridGeometry);
    using Indices = typename GET_PROP_TYPE(TypeTag, ModelTraits)::Indices;
    using NumEqVector = typename GET_PROP_TYPE(TypeTag, NumEqVector);
    using PrimaryVariables = typename GET_PROP_TYPE(TypeTag, PrimaryVariables);
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);

    static constexpr auto dimWorld = GET_PROP_TYPE(TypeTag, GridView)::dimensionworld;
    using GlobalPosition = Dune::FieldVector<Scalar, dimWorld>;

    using TimeLoopPtr = std::shared_ptr<CheckPointTimeLoop<Scalar>>;

    static constexpr auto compIdx = 1;
    static constexpr auto transportCompIdx = Indices::conti0EqIdx + compIdx;
    static constexpr auto transportEqIdx = Indices::conti0EqIdx + compIdx;

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

        if(isInlet(globalPos))
        {
            values.setDirichlet(Indices::velocityXIdx);
            values.setDirichlet(Indices::velocityYIdx);
            values.setDirichlet(transportCompIdx);
#if NONISOTHERMAL
            values.setDirichlet(Indices::temperatureIdx);
#endif
        }
        else if(isOutlet(globalPos))
        {
            values.setDirichlet(Indices::pressureIdx);
            values.setOutflow(transportEqIdx);
#if NONISOTHERMAL
            values.setOutflow(Indices::energyBalanceIdx);
#endif
        }
        else
        {
            // set Dirichlet values for the velocity everywhere
            values.setDirichlet(Indices::velocityXIdx);
            values.setDirichlet(Indices::velocityYIdx);
            values.setNeumann(Indices::conti0EqIdx);
            values.setNeumann(transportEqIdx);
#if NONISOTHERMAL
            values.setNeumann(Indices::energyBalanceIdx);
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
                Scalar moleFracTransportedComp = 1e-3;
#if USE_MASS
                Scalar averageMolarMassPhase = moleFracTransportedComp * FluidSystem::molarMass(compIdx)
                                               + (1. - moleFracTransportedComp)  * FluidSystem::molarMass(1-compIdx);
                values[transportCompIdx] = moleFracTransportedComp * FluidSystem::molarMass(compIdx)
                                           / averageMolarMassPhase;
#else
                values[transportCompIdx] = moleFracTransportedComp;
#endif
#if NONISOTHERMAL
            values[Indices::temperatureIdx] = 293.15;
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
        values[Indices::pressureIdx] = 1.1e+5;
        values[transportCompIdx] = 0.0;
#if NONISOTHERMAL
        values[Indices::temperatureIdx] = 283.15;
#endif

        // parabolic velocity profile
        values[Indices::velocityXIdx] =  inletVelocity_*(globalPos[1] - this->fvGridGeometry().bBoxMin()[1])*(this->fvGridGeometry().bBoxMax()[1] - globalPos[1])
                                      / (0.25*(this->fvGridGeometry().bBoxMax()[1] - this->fvGridGeometry().bBoxMin()[1])*(this->fvGridGeometry().bBoxMax()[1] - this->fvGridGeometry().bBoxMin()[1]));

        values[Indices::velocityYIdx] = 0.0;

        return values;
    }

   /*!
     * \brief Adds additional VTK output data to the VTKWriter. Function is called by the output module on every write.
     *
     * \param gridVariables The grid variables
     * \param sol The solution vector
     */
    template<class GridVariables, class SolutionVector>
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
