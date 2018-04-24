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
 * \ingroup RANSTests
 * \brief Pipe flow test for the staggered grid RANS model
 *
 * This test simulates is based on pipe flow experiments by
 * John Laufers experiments in 1954 \cite Laufer1954a.
 */
#ifndef DUMUX_PIPE_LAUFER_PROBLEM_HH
#define DUMUX_PIPE_LAUFER_PROBLEM_HH

#include <dumux/material/fluidsystems/1pgas.hh>
#include <dumux/material/components/air.hh>

#include <dumux/freeflow/turbulenceproperties.hh>
#if LOWREKEPSILON
#include <dumux/freeflow/rans/twoeq/lowrekepsilon/model.hh>
#include <dumux/freeflow/rans/twoeq/lowrekepsilon/problem.hh>
#else
#include <dumux/freeflow/rans/zeroeq/model.hh>
#include <dumux/freeflow/rans/zeroeq/problem.hh>
#endif
#include <dumux/discretization/staggered/freeflow/properties.hh>

namespace Dumux
{
template <class TypeTag>
class PipeLauferProblem;

namespace Properties
{
#if NONISOTHERMAL
NEW_TYPE_TAG(PipeLauferProblem, INHERITS_FROM(StaggeredFreeFlowModel, ZeroEqNI));
#else
#if LOWREKEPSILON
NEW_TYPE_TAG(PipeLauferProblem, INHERITS_FROM(StaggeredFreeFlowModel, LowReKEpsilon));
#else
NEW_TYPE_TAG(PipeLauferProblem, INHERITS_FROM(StaggeredFreeFlowModel, ZeroEq));
#endif
#endif

// the fluid system
SET_PROP(PipeLauferProblem, FluidSystem)
{
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using type = FluidSystems::OnePGas<Scalar, Components::Air<Scalar> >;
};

// Set the grid type
SET_TYPE_PROP(PipeLauferProblem, Grid,
              Dune::YaspGrid<2, Dune::TensorProductCoordinates<typename GET_PROP_TYPE(TypeTag, Scalar), 2> >);

// Set the problem property
SET_TYPE_PROP(PipeLauferProblem, Problem, Dumux::PipeLauferProblem<TypeTag> );

SET_BOOL_PROP(PipeLauferProblem, EnableFVGridGeometryCache, true);

SET_BOOL_PROP(PipeLauferProblem, EnableGridFluxVariablesCache, true);
SET_BOOL_PROP(PipeLauferProblem, EnableGridVolumeVariablesCache, true);
}

/*!
 * \ingroup NavierStokesTests
 * \brief  Test problem for the one-phase (Navier-) Stokes problem in a channel.
 *
 * This test simulates is based on pipe flow experiments by
 * John Laufers experiments in 1954 \cite Laufer1954a.
 */
template <class TypeTag>
#if LOWREKEPSILON
class PipeLauferProblem : public LowReKEpsilonProblem<TypeTag>
{
    using ParentType = LowReKEpsilonProblem<TypeTag>;
#else
class PipeLauferProblem : public ZeroEqProblem<TypeTag>
{
    using ParentType = ZeroEqProblem<TypeTag>;
#endif

    using BoundaryTypes = typename GET_PROP_TYPE(TypeTag, BoundaryTypes);
    using FluidSystem = typename GET_PROP_TYPE(TypeTag, FluidSystem);
    using FluidState = typename GET_PROP_TYPE(TypeTag, FluidState);
    using FVGridGeometry = typename GET_PROP_TYPE(TypeTag, FVGridGeometry);
    using Indices = typename GET_PROP_TYPE(TypeTag, ModelTraits)::Indices;
    using NumEqVector = typename GET_PROP_TYPE(TypeTag, NumEqVector);
    using PrimaryVariables = typename GET_PROP_TYPE(TypeTag, PrimaryVariables);
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);

    static constexpr auto dimWorld = GET_PROP_TYPE(TypeTag, GridView)::dimensionworld;
    using GlobalPosition = Dune::FieldVector<Scalar, dimWorld>;

    using TimeLoopPtr = std::shared_ptr<CheckPointTimeLoop<Scalar>>;

    static const unsigned int phaseIdx = GET_PROP_VALUE(TypeTag, PhaseIdx);

public:
    PipeLauferProblem(std::shared_ptr<const FVGridGeometry> fvGridGeometry)
    : ParentType(fvGridGeometry), eps_(1e-6)
    {
        inletVelocity_ = getParam<Scalar>("Problem.InletVelocity");
        inletTemperature_ = getParam<Scalar>("Problem.InletTemperature", 283.15);
        wallTemperature_ = getParam<Scalar>("Problem.WallTemperature", 323.15);
        sandGrainRoughness_ = getParam<Scalar>("Problem.SandGrainRoughness", 0.0);
        startWithZeroVelocity_ = getParam<bool>("RANS.StartWithZeroVelocity", false);

        FluidSystem::init();
        Dumux::TurbulenceProperties<Scalar, dimWorld, true> turbulenceProperties;
        FluidState fluidState;
        fluidState.setPressure(phaseIdx, 1e5);
        fluidState.setTemperature(inletTemperature_);
        Scalar density = FluidSystem::density(fluidState, phaseIdx);
        Scalar kinematicViscosity = FluidSystem::viscosity(fluidState, phaseIdx) / density;
        Scalar diameter = this->fvGridGeometry().bBoxMax()[1] - this->fvGridGeometry().bBoxMin()[1];
        viscosityTilde_ = turbulenceProperties.viscosityTilde(inletVelocity_, diameter, kinematicViscosity);
        turbulentKineticEnergy_ = turbulenceProperties.turbulentKineticEnergy(inletVelocity_, diameter, kinematicViscosity);
        dissipation_ = turbulenceProperties.dissipation(inletVelocity_, diameter, kinematicViscosity);
        dissipationRate_ = turbulenceProperties.dissipationRate(inletVelocity_, diameter, kinematicViscosity);
        std::cout << std::endl;
    }

   /*!
     * \name Problem parameters
     */
    // \{

    bool isOnWall(const GlobalPosition &globalPos) const
    {
        return globalPos[1] < this->fvGridGeometry().bBoxMin()[1] + eps_
               || globalPos[1] > this->fvGridGeometry().bBoxMax()[1] - eps_;
    }

    Scalar sandGrainRoughnessAtPos(const GlobalPosition &globalPos) const
    {
        return sandGrainRoughness_;
    }

    bool shouldWriteRestartFile() const
    {
        return false;
    }

   /*!
     * \brief Return the temperature [K] within the domain for the isothermal model.
     */
    Scalar temperature() const
    { return inletTemperature_; }

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

        // set Dirichlet values for the velocity and outflow for total mass everywhere
        values.setOutflow(Indices::conti0EqIdx);
        values.setDirichlet(Indices::momentumXBalanceIdx);
        values.setDirichlet(Indices::momentumYBalanceIdx);
        if (isOutlet(globalPos))
        {
            values.setDirichlet(Indices::conti0EqIdx);
            values.setOutflow(Indices::momentumXBalanceIdx);
            values.setOutflow(Indices::momentumYBalanceIdx);
        }

#if NONISOTHERMAL
        values.setDirichlet(Indices::energyBalanceIdx);
        if (isOutlet(globalPos))
        {
            values.setOutflow(Indices::energyBalanceIdx);
        }
#endif

#if LOWREKEPSILON
        values.setDirichlet(Indices::turbulentKineticEnergyIdx);
        values.setDirichlet(Indices::dissipationIdx);
        if (isOutlet(globalPos))
        {
            values.setOutflow(Indices::turbulentKineticEnergyEqIdx);
            values.setOutflow(Indices::dissipationEqIdx);
        }
#endif

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
        PrimaryVariables values(initialAtPos(globalPos));
#if NONISOTHERMAL
        if (time() > 10.0)
        {
            values[Indices::temperatureIdx] = inletTemperature_;
            if (isOnWall(globalPos))
            {
                values[Indices::temperatureIdx] = wallTemperature_;
            }
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
        values[Indices::pressureIdx] = 1.0e+5;
        values[Indices::velocityXIdx] = inletVelocity_;
        if (isOnWall(globalPos))
        {
            values[Indices::velocityXIdx] = 0.0;
        }

#if NONISOTHERMAL
        values[Indices::temperatureIdx] = inletTemperature_;
        if (isOnWall(globalPos))
        {
            values[Indices::temperatureIdx] = wallTemperature_;
        }
#endif

#if LOWREKEPSILON
        if (time() < eps_ && startWithZeroVelocity_)
        {
            values[Indices::velocityXIdx] = 0.0;
        }
        values[Indices::turbulentKineticEnergyEqIdx] = turbulentKineticEnergy_;
        values[Indices::dissipationEqIdx] = dissipation_;
        if (isOnWall(globalPos))
        {
            values[Indices::turbulentKineticEnergyEqIdx] = 0.0;
            values[Indices::dissipationEqIdx] = 0.0;
        }
#endif

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

private:
    bool isInlet(const GlobalPosition& globalPos) const
    {
        return globalPos[0] < this->fvGridGeometry().bBoxMin()[0] + eps_;
    }

    bool isOutlet(const GlobalPosition& globalPos) const
    {
        return globalPos[0] > this->fvGridGeometry().bBoxMax()[0] - eps_;
    }

    Scalar eps_;
    Scalar inletVelocity_;
    Scalar inletTemperature_;
    Scalar wallTemperature_;
    Scalar sandGrainRoughness_;
    bool startWithZeroVelocity_;
    Scalar viscosityTilde_;
    Scalar turbulentKineticEnergy_;
    Scalar dissipation_;
    Scalar dissipationRate_;
    TimeLoopPtr timeLoop_;
};
} //end namespace

#endif
