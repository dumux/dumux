// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   See the file COPYING for full copying permissions.                      *
 *                                                                           *
 *   This program is free software: you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation, either version 3 of the License, or       *
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
// Create new type tags
namespace TTag {
#if NONISOTHERMAL
struct PipeLauferProblem { using InheritsFrom = std::tuple<ZeroEqNI, StaggeredFreeFlowModel>; };
#else
#if LOWREKEPSILON
struct PipeLauferProblem { using InheritsFrom = std::tuple<LowReKEpsilon, StaggeredFreeFlowModel>; };
#elif KEPSILON
struct PipeLauferProblem { using InheritsFrom = std::tuple<KEpsilon, StaggeredFreeFlowModel>; };
#elif KOMEGA
struct PipeLauferProblem { using InheritsFrom = std::tuple<KOmega, StaggeredFreeFlowModel>; };
#elif ONEEQ
struct PipeLauferProblem { using InheritsFrom = std::tuple<OneEq, StaggeredFreeFlowModel>; };
#else
struct PipeLauferProblem { using InheritsFrom = std::tuple<ZeroEq, StaggeredFreeFlowModel>; };
#endif
#endif
} // end namespace TTag

// the fluid system
template<class TypeTag>
struct FluidSystem<TypeTag, TTag::PipeLauferProblem>
{
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using type = FluidSystems::OnePGas<Scalar, Components::Air<Scalar> >;
};

// Set the grid type
template<class TypeTag>
struct Grid<TypeTag, TTag::PipeLauferProblem>
{ using type = Dune::YaspGrid<2, Dune::TensorProductCoordinates<GetPropType<TypeTag, Properties::Scalar>, 2> >; };

// Set the problem property
template<class TypeTag>
struct Problem<TypeTag, TTag::PipeLauferProblem> { using type = Dumux::PipeLauferProblem<TypeTag> ; };

template<class TypeTag>
struct EnableFVGridGeometryCache<TypeTag, TTag::PipeLauferProblem> { static constexpr bool value = true; };

template<class TypeTag>
struct EnableGridFluxVariablesCache<TypeTag, TTag::PipeLauferProblem> { static constexpr bool value = true; };
template<class TypeTag>
struct EnableGridVolumeVariablesCache<TypeTag, TTag::PipeLauferProblem> { static constexpr bool value = true; };
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

    using BoundaryTypes = GetPropType<TypeTag, Properties::BoundaryTypes>;
    using FluidSystem = GetPropType<TypeTag, Properties::FluidSystem>;
    using FluidState = GetPropType<TypeTag, Properties::FluidState>;
    using FVGridGeometry = GetPropType<TypeTag, Properties::FVGridGeometry>;
    using Indices = typename GetPropType<TypeTag, Properties::ModelTraits>::Indices;
    using NumEqVector = GetPropType<TypeTag, Properties::NumEqVector>;
    using PrimaryVariables = GetPropType<TypeTag, Properties::PrimaryVariables>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;

    using Element = typename FVGridGeometry::GridView::template Codim<0>::Entity;
    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;
    using FVElementGeometry = typename GetPropType<TypeTag, Properties::FVGridGeometry>::LocalView;
    using SubControlVolume = typename FVElementGeometry::SubControlVolume;
    using SubControlVolumeFace = typename FVElementGeometry::SubControlVolumeFace;

    using TimeLoopPtr = std::shared_ptr<CheckPointTimeLoop<Scalar>>;

    static constexpr auto dimWorld = FVGridGeometry::GridView::dimensionworld;

public:
    PipeLauferProblem(std::shared_ptr<const FVGridGeometry> fvGridGeometry)
    : ParentType(fvGridGeometry), eps_(1e-6)
    {
        inletVelocity_ = getParam<Scalar>("Problem.InletVelocity");
        inletTemperature_ = getParam<Scalar>("Problem.InletTemperature", 283.15);
        wallTemperature_ = getParam<Scalar>("Problem.WallTemperature", 323.15);
        sandGrainRoughness_ = getParam<Scalar>("Problem.SandGrainRoughness", 0.0);
#if ONEEQ
        initializationTime_ = getParam<Scalar>("TimeLoop.Initialization", 1.0);
#else
        initializationTime_ = getParam<Scalar>("TimeLoop.Initialization", -1.0);
#endif

        FluidSystem::init();
        Dumux::TurbulenceProperties<Scalar, dimWorld, true> turbulenceProperties;
        FluidState fluidState;
        fluidState.setPressure(0, 1e5);
        fluidState.setTemperature(inletTemperature_);
        Scalar density = FluidSystem::density(fluidState, 0);
        Scalar kinematicViscosity = FluidSystem::viscosity(fluidState, 0) / density;
        Scalar diameter = this->fvGridGeometry().bBoxMax()[1] - this->fvGridGeometry().bBoxMin()[1];
        // ideally the viscosityTilde parameter as inflow for the Spalart-Allmaras model should be zero
        viscosityTilde_ = getParam<Scalar>("Problem.InletViscosityTilde",
                                           1e-3 * turbulenceProperties.viscosityTilde(inletVelocity_, diameter, kinematicViscosity));
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

    bool isOnWallAtPos(const GlobalPosition &globalPos) const
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

        if(isOutlet_(globalPos))
        {
            values.setDirichlet(Indices::pressureIdx);

#if NONISOTHERMAL
            values.setOutflow(Indices::energyEqIdx);
#endif

#if LOWREKEPSILON || KEPSILON || KOMEGA
            values.setOutflow(Indices::turbulentKineticEnergyEqIdx);
            values.setOutflow(Indices::dissipationEqIdx);
#endif

#if ONEEQ
            values.setOutflow(Indices::viscosityTildeIdx);
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

#if ONEEQ
            values.setDirichlet(Indices::viscosityTildeIdx);
#endif
        }
        return values;
    }

#if KOMEGA
    /*!
     * \brief Returns whether a fixed Dirichlet value shall be used at a given cell.
     *
     * \param element The finite element
     * \param fvGeometry The finite-volume geometry
     * \param scv The sub control volume
     */
    template<class Element, class FVElementGeometry, class SubControlVolume>
    bool isDirichletCell(const Element& element,
                         const FVElementGeometry& fvGeometry,
                         const SubControlVolume& scv,
                         int pvIdx) const
    {
        // set a fixed dissipation (omega) for all cells at the wall
        for (const auto& scvf : scvfs(fvGeometry))
            if (isOnWallAtPos(scvf.center()) && pvIdx == Indices::dissipationIdx)
                return true;

        return false;
    }
#endif

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
            if (isOnWallAtPos(globalPos))
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
        unsigned int elementIdx = this->fvGridGeometry().elementMapper().index(element);
        const auto wallDistance = ParentType::wallDistance_[elementIdx];
        values[Indices::dissipationEqIdx] = 6.0 * ParentType::kinematicViscosity_[elementIdx]
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
        values[Indices::velocityXIdx] = time() > initializationTime_
                                        ? inletVelocity_
                                        : time() / initializationTime_ * inletVelocity_;
        if (isOnWallAtPos(globalPos))
        {
            values[Indices::velocityXIdx] = 0.0;
        }

#if NONISOTHERMAL
        values[Indices::temperatureIdx] = inletTemperature_;
        if (isOnWallAtPos(globalPos))
        {
            values[Indices::temperatureIdx] = wallTemperature_;
        }
#endif

#if LOWREKEPSILON || KEPSILON || KOMEGA
        values[Indices::turbulentKineticEnergyIdx] = turbulentKineticEnergy_;
        values[Indices::dissipationIdx] = dissipation_;
        if (isOnWallAtPos(globalPos))
        {
            values[Indices::turbulentKineticEnergyIdx] = 0.0;
            values[Indices::dissipationIdx] = 0.0;
        }
#endif

#if ONEEQ
        values[Indices::viscosityTildeIdx] = viscosityTilde_;
        if (isOnWallAtPos(globalPos))
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
    bool isInlet_(const GlobalPosition& globalPos) const
    {
        return globalPos[0] < this->fvGridGeometry().bBoxMin()[0] + eps_;
    }

    bool isOutlet_(const GlobalPosition& globalPos) const
    {
        return globalPos[0] > this->fvGridGeometry().bBoxMax()[0] - eps_;
    }

    Scalar eps_;
    Scalar inletVelocity_;
    Scalar inletTemperature_;
    Scalar wallTemperature_;
    Scalar sandGrainRoughness_;
    Scalar initializationTime_;
    Scalar viscosityTilde_;
    Scalar turbulentKineticEnergy_;
    Scalar dissipation_;
    TimeLoopPtr timeLoop_;
};
} //end namespace

#endif
