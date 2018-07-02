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
 * John Laufer's experiments in 1954 \cite Laufer1954a.
 */
#ifndef DUMUX_PIPE_LAUFER_PROBLEM_HH
#define DUMUX_PIPE_LAUFER_PROBLEM_HH

#include <dune/grid/yaspgrid.hh>

#include <dumux/discretization/staggered/freeflow/properties.hh>
#include <dumux/freeflow/turbulenceproperties.hh>
#include <dumux/material/fluidsystems/1pgas.hh>
#include <dumux/material/components/air.hh>

#if LOWREKEPSILON
#include <dumux/freeflow/rans/twoeq/lowrekepsilon/model.hh>
#include <dumux/freeflow/rans/twoeq/lowrekepsilon/problem.hh>
#elif KEPSILON
#include <dumux/freeflow/rans/twoeq/kepsilon/model.hh>
#include <dumux/freeflow/rans/twoeq/kepsilon/problem.hh>
#elif KOMEGA
#include <dumux/freeflow/rans/twoeq/komega/model.hh>
#include <dumux/freeflow/rans/twoeq/komega/problem.hh>
#elif ONEEQ
#include <dumux/freeflow/rans/oneeq/model.hh>
#include <dumux/freeflow/rans/oneeq/problem.hh>
#else
#include <dumux/freeflow/rans/zeroeq/model.hh>
#include <dumux/freeflow/rans/zeroeq/problem.hh>
#endif

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
#elif KEPSILON
NEW_TYPE_TAG(PipeLauferProblem, INHERITS_FROM(StaggeredFreeFlowModel, KEpsilon));
#elif KOMEGA
NEW_TYPE_TAG(PipeLauferProblem, INHERITS_FROM(StaggeredFreeFlowModel, KOmega));
#elif ONEEQ
NEW_TYPE_TAG(PipeLauferProblem, INHERITS_FROM(StaggeredFreeFlowModel, OneEq));
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
#elif KEPSILON
class PipeLauferProblem : public KEpsilonProblem<TypeTag>
{
    using ParentType = KEpsilonProblem<TypeTag>;
#elif KOMEGA
class PipeLauferProblem : public KOmegaProblem<TypeTag>
{
    using ParentType = KOmegaProblem<TypeTag>;
#elif ONEEQ
class PipeLauferProblem : public OneEqProblem<TypeTag>
{
    using ParentType = OneEqProblem<TypeTag>;
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

    using Element = typename FVGridGeometry::GridView::template Codim<0>::Entity;
    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;
    using FVElementGeometry = typename GET_PROP_TYPE(TypeTag, FVGridGeometry)::LocalView;
    using SubControlVolume = typename FVElementGeometry::SubControlVolume;
    using SubControlVolumeFace = typename FVElementGeometry::SubControlVolumeFace;

    using TimeLoopPtr = std::shared_ptr<CheckPointTimeLoop<Scalar>>;

    static const unsigned int phaseIdx = GET_PROP_VALUE(TypeTag, PhaseIdx);
    static constexpr auto dimWorld = FVGridGeometry::GridView::dimensionworld;

public:
    PipeLauferProblem(std::shared_ptr<const FVGridGeometry> fvGridGeometry)
    : ParentType(fvGridGeometry), eps_(1e-6)
    {
        inletVelocity_ = getParam<Scalar>("Problem.InletVelocity");
        inletTemperature_ = getParam<Scalar>("Problem.InletTemperature", 283.15);
        wallTemperature_ = getParam<Scalar>("Problem.WallTemperature", 323.15);
        sandGrainRoughness_ = getParam<Scalar>("Problem.SandGrainRoughness", 0.0);
        startWithZeroVelocity_ = getParam<bool>("Problem.StartWithZeroVelocity", false);

        FluidSystem::init();
        Dumux::TurbulenceProperties<Scalar, dimWorld, true> turbulenceProperties;
        FluidState fluidState;
        fluidState.setPressure(phaseIdx, 1e5);
        fluidState.setTemperature(inletTemperature_);
        Scalar density = FluidSystem::density(fluidState, phaseIdx);
        Scalar kinematicViscosity = FluidSystem::viscosity(fluidState, phaseIdx) / density;
        Scalar diameter = this->fvGridGeometry().bBoxMax()[1] - this->fvGridGeometry().bBoxMin()[1];
        viscosityTilde_ = getParam<Scalar>("Problem.InletViscosityTilde",
                                           turbulenceProperties.viscosityTilde(inletVelocity_, diameter, kinematicViscosity));
        turbulentKineticEnergy_ = getParam<Scalar>("Problem.InletTurbulentKineticEnergy",
                                                   turbulenceProperties.turbulentKineticEnergy(inletVelocity_, diameter, kinematicViscosity));
#if KOMEGA
        dissipation_ = getParam<Scalar>("Problem.InletDissipationRate",
                                        turbulenceProperties.dissipationRate(inletVelocity_, diameter, kinematicViscosity));
#else
        dissipation_ = getParam<Scalar>("Problem.InletDissipation",
                                        turbulenceProperties.dissipation(inletVelocity_, diameter, kinematicViscosity));
#endif
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

        if(isOutlet(globalPos))
        {
            values.setDirichlet(Indices::pressureIdx);

#if NONISOTHERMAL
            values.setOutflow(Indices::energyBalanceIdx);
#endif

#if LOWREKEPSILON || KEPSILON || KOMEGA
            values.setOutflow(Indices::turbulentKineticEnergyEqIdx);
            values.setOutflow(Indices::dissipationEqIdx);
#endif
        }
        else // walls and inflow
        {
            values.setDirichlet(Indices::velocityXIdx);
            values.setDirichlet(Indices::velocityYIdx);

#if NONISOTHERMAL
            values.setDirichlet(Indices::temperatureIdx);
#endif

#if LOWREKEPSILON || KEPSILON || KOMEGA
            values.setDirichlet(Indices::turbulentKineticEnergyIdx);
            values.setDirichlet(Indices::dissipationIdx);
#endif
#if KOMEGA
            // set a fixed dissipation (omega) in one cell
            if (isOnWall(globalPos))
                values.setDirichletCell(Indices::dissipationIdx);
#endif
#if ONEEQ
            values.setDirichlet(Indices::viscosityTildeIdx);
#endif
        }
        return values;
    }

     /*!
      * \brief Evaluate the boundary conditions for a dirichlet values at the boundary.
      *
      * \param element The finite element
      * \param scvf the sub control volume face
      * \note used for cell-centered discretization schemes
      */
    PrimaryVariables dirichlet(const Element &element, const SubControlVolumeFace &scvf) const
    {
        const auto globalPos = scvf.ipGlobal();
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

     /*!
      * \brief Evaluate the boundary conditions for fixed values at cell centers
      *
      * \param element The finite element
      * \param scv the sub control volume
      * \note used for cell-centered discretization schemes
      */
    PrimaryVariables dirichlet(const Element &element, const SubControlVolume &scv) const
    {
        const auto globalPos = scv.center();
        PrimaryVariables values(initialAtPos(globalPos));
#if KOMEGA
        using std::pow;
        unsigned int elementID = this->fvGridGeometry().elementMapper().index(element);
        const auto wallDistance = ParentType::wallDistance_[elementID];
        values[Indices::dissipationEqIdx] = 6.0 * ParentType::kinematicViscosity_[elementID]
                                            / (ParentType::betaOmega() * pow(wallDistance, 2));
#endif
        return values;
    }

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
        if (isOnWall(globalPos)
            || (startWithZeroVelocity_ && time() < eps_))
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

#if LOWREKEPSILON || KEPSILON || KOMEGA
        values[Indices::turbulentKineticEnergyIdx] = turbulentKineticEnergy_;
        values[Indices::dissipationIdx] = dissipation_;
        if (isOnWall(globalPos))
        {
            values[Indices::turbulentKineticEnergyIdx] = 0.0;
            values[Indices::dissipationIdx] = 0.0;
        }
#endif

#if ONEEQ
        values[Indices::viscosityTildeIdx] = viscosityTilde_;
        if (isOnWall(globalPos))
        {
            values[Indices::viscosityTildeIdx] = 0.0;
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
    TimeLoopPtr timeLoop_;
};
} //end namespace

#endif
