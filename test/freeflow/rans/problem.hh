// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup RANSTests
 * \brief Pipe flow test for the staggered grid RANS model
 *
 * This test simulates pipe flow experiments performed by John Laufer in 1954
 * \cite Laufer1954a.
 */
#ifndef DUMUX_PIPE_LAUFER_PROBLEM_HH
#define DUMUX_PIPE_LAUFER_PROBLEM_HH

#include <dumux/common/properties.hh>
#include <dumux/common/parameters.hh>
#include <dumux/common/timeloop.hh>
#include <dumux/common/numeqvector.hh>

#include <dumux/freeflow/rans/boundarytypes.hh>
#include <dumux/freeflow/rans/problem.hh>
#include <dumux/freeflow/turbulencemodel.hh>
#include <dumux/freeflow/turbulenceproperties.hh>

namespace Dumux {

/*!
 * \ingroup NavierStokesTests
 * \brief  Test problem for the one-phase (Navier-) Stokes problem in a channel.
 *
 * This test simulates is based on pipe flow experiments by
 * John Laufers experiments in 1954 \cite Laufer1954a.
 */
template <class TypeTag>
class PipeLauferProblem : public RANSProblem<TypeTag>
{
    using ParentType = RANSProblem<TypeTag>;

    using FluidSystem = GetPropType<TypeTag, Properties::FluidSystem>;
    using FluidState = GetPropType<TypeTag, Properties::FluidState>;
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using PrimaryVariables = GetPropType<TypeTag, Properties::PrimaryVariables>;
    using NumEqVector = Dumux::NumEqVector<PrimaryVariables>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;

    static constexpr auto dimWorld = GridGeometry::GridView::dimensionworld;
    using ModelTraits = GetPropType<TypeTag, Properties::ModelTraits>;
    using Indices = typename GetPropType<TypeTag, Properties::ModelTraits>::Indices;
    using BoundaryTypes = Dumux::RANSBoundaryTypes<ModelTraits, ModelTraits::numEq()>;
    using Element = typename GridGeometry::GridView::template Codim<0>::Entity;
    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;
    using FVElementGeometry = typename GetPropType<TypeTag, Properties::GridGeometry>::LocalView;
    using SubControlVolume = typename FVElementGeometry::SubControlVolume;
    using SubControlVolumeFace = typename FVElementGeometry::SubControlVolumeFace;

    using TimeLoopPtr = std::shared_ptr<CheckPointTimeLoop<Scalar>>;


public:
    PipeLauferProblem(std::shared_ptr<const GridGeometry> gridGeometry)
    : ParentType(gridGeometry), eps_(1e-6)
    {
        inletVelocity_ = getParam<Scalar>("Problem.InletVelocity");
        inletTemperature_ = getParam<Scalar>("Problem.InletTemperature", 283.15);
        wallTemperature_ = getParam<Scalar>("Problem.WallTemperature", 303.15);

        FluidSystem::init();
        Dumux::TurbulenceProperties<Scalar, dimWorld, true> turbulenceProperties;
        FluidState fluidState;
        fluidState.setPressure(0, 1e5);
        fluidState.setTemperature(this->spatialParams().temperatureAtPos({}));
        Scalar density = FluidSystem::density(fluidState, 0);
        Scalar kinematicViscosity = FluidSystem::viscosity(fluidState, 0) / density;
        Scalar diameter = this->gridGeometry().bBoxMax()[1] - this->gridGeometry().bBoxMin()[1];

        // ideally the viscosityTilde parameter as inflow for the Spalart-Allmaras model should be zero
        viscosityTilde_ = 1e-3 * turbulenceProperties.viscosityTilde(inletVelocity_, diameter, kinematicViscosity);
        turbulentKineticEnergy_ = turbulenceProperties.turbulentKineticEnergy(inletVelocity_, diameter, kinematicViscosity);

        if (ModelTraits::turbulenceModel() == TurbulenceModel::komega || ModelTraits::turbulenceModel() == TurbulenceModel::sst)
            dissipation_ = turbulenceProperties.dissipationRate(inletVelocity_, diameter, kinematicViscosity);
        else
            dissipation_ = turbulenceProperties.dissipation(inletVelocity_, diameter, kinematicViscosity);

        if (ModelTraits::turbulenceModel() == TurbulenceModel::oneeq)
            initializationTime_ = getParam<Scalar>("TimeLoop.Initialization", 1.0);
        else
            initializationTime_ = getParam<Scalar>("TimeLoop.Initialization", -1.0);

        turbulenceModelName_ = turbulenceModelToString(ModelTraits::turbulenceModel());
        std::cout << "Using the "<< turbulenceModelName_ << " Turbulence Model. \n";
        std::cout << std::endl;
    }

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

        // common boundary types for all turbulence models
        if(isOutlet_(globalPos))
            values.setDirichlet(Indices::pressureIdx);
        else
        {
            values.setDirichlet(Indices::velocityXIdx);
            values.setDirichlet(Indices::velocityYIdx);
        }

        if (isLowerWall_(globalPos) || isUpperWall_(globalPos)) //  All walls
            values.setWall();

#if NONISOTHERMAL
        if(isOutlet_(globalPos))
            values.setOutflow(Indices::energyEqIdx);
        else
            values.setDirichlet(Indices::temperatureIdx);
#endif
        // turbulence model-specific boundary types
        setBcTypes_(values, globalPos);

        return values;
    }

    /*!
     * \brief Returns whether a fixed Dirichlet value shall be used at a given cell.
     *
     * \param element The finite element
     * \param fvGeometry The finite-volume geometry
     * \param scv The sub control volume
     * \param pvIdx The primary variable index in the solution vector
     */
    template<class Element, class FVElementGeometry, class SubControlVolume>
    bool isDirichletCell(const Element& element,
                         const FVElementGeometry& fvGeometry,
                         const SubControlVolume& scv,
                         int pvIdx) const
    { return isDirichletCell_(element, fvGeometry, scv, pvIdx); }

     /*!
      * \brief Evaluate the boundary conditions for a dirichlet values at the boundary.
      *
      * \param element The finite element
      * \param scvf the sub control volume face
      * \note used for cell-centered discretization schemes
      */
    PrimaryVariables dirichlet(const Element& element, const SubControlVolumeFace& scvf) const
    {
        const auto globalPos = scvf.ipGlobal();
        PrimaryVariables values(initialAtPos(globalPos));

#if NONISOTHERMAL
        values[Indices::temperatureIdx] =  (isLowerWall_(globalPos) || isUpperWall_(globalPos)) ? wallTemperature_ : inletTemperature_;
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
    PrimaryVariables dirichlet([[maybe_unused]] const Element& element, const SubControlVolume& scv) const
    {
        if constexpr (ModelTraits::turbulenceModel() == TurbulenceModel::kepsilon
                   || ModelTraits::turbulenceModel() == TurbulenceModel::komega
                   || ModelTraits::turbulenceModel() == TurbulenceModel::sst)
            return dirichletTurbulentTwoEq_(element, scv);
        else
        {
            const auto globalPos = scv.center();
            PrimaryVariables values(initialAtPos(globalPos));
            return values;
        }
    }

   /*!
     * \brief Evaluate the initial value for a control volume.
     *
     * \param globalPos The global position
     */
    PrimaryVariables initialAtPos(const GlobalPosition& globalPos) const
    {
        // common initial conditions for all turbulence models
        PrimaryVariables values(0.0);
        values[Indices::pressureIdx] = 1.0e+5;
        values[Indices::velocityXIdx] = time() > initializationTime_
                                        ? inletVelocity_
                                        : time() / initializationTime_ * inletVelocity_;
        if (isLowerWall_(globalPos) || isUpperWall_(globalPos))
            values[Indices::velocityXIdx] = 0.0;
#if NONISOTHERMAL
        values[Indices::temperatureIdx] = (isLowerWall_(globalPos) || isUpperWall_(globalPos)) ? wallTemperature_ : inletTemperature_;
#endif
        // turbulence model-specific initial conditions
        setInitialAtPos_(values, globalPos);
        return values;
    }

    // \}

    void setTimeLoop(TimeLoopPtr timeLoop)
    { timeLoop_ = timeLoop; }

    Scalar time() const
    { return timeLoop_->time(); }

private:
    bool isInlet_(const GlobalPosition& globalPos) const
    { return globalPos[0] < this->gridGeometry().bBoxMin()[0] + eps_; }

    bool isOutlet_(const GlobalPosition& globalPos) const
    { return globalPos[0] > this->gridGeometry().bBoxMax()[0] - eps_; }

    bool isLowerWall_(const GlobalPosition& globalPos) const
    { return globalPos[1] < this->gridGeometry().bBoxMin()[1] + eps_; }

    bool isUpperWall_(const GlobalPosition& globalPos) const
    { return globalPos[1] > this->gridGeometry().bBoxMax()[1] - eps_; }

    //! Initial conditions for the komega, kepsilon and lowrekepsilon turbulence models
    void setInitialAtPos_([[maybe_unused]] PrimaryVariables& values,
                          [[maybe_unused]] const GlobalPosition &globalPos) const
    {
        if constexpr (numTurbulenceEq(ModelTraits::turbulenceModel()) == 0) // zero equation models
            return;
        else if constexpr (numTurbulenceEq(ModelTraits::turbulenceModel()) == 1)  // one equation models
        {
            values[Indices::viscosityTildeIdx] = viscosityTilde_;
            if (isLowerWall_(globalPos) || isUpperWall_(globalPos))
                values[Indices::viscosityTildeIdx] = 0.0;
        }
        else // two equation models
        {
            static_assert(numTurbulenceEq(ModelTraits::turbulenceModel()) == 2, "Only reached by 2eq models");
            values[Indices::turbulentKineticEnergyIdx] = turbulentKineticEnergy_;
            values[Indices::dissipationIdx] = dissipation_;
            if (isLowerWall_(globalPos) || isUpperWall_(globalPos))
            {
                values[Indices::turbulentKineticEnergyIdx] = 0.0;
                values[Indices::dissipationIdx] = 0.0;
            }
        }
    }

    //! Boundary condition types for the one-eq turbulence model
    void setBcTypes_([[maybe_unused]] BoundaryTypes& values,
                     [[maybe_unused]] const GlobalPosition& pos) const
    {
        if constexpr (numTurbulenceEq(ModelTraits::turbulenceModel()) == 0) // zero equation models
            return;
        else if constexpr (numTurbulenceEq(ModelTraits::turbulenceModel()) == 1)  // one equation models
        {
            if(isOutlet_(pos))
                values.setOutflow(Indices::viscosityTildeIdx);
            else // walls and inflow
                values.setDirichlet(Indices::viscosityTildeIdx);
        }
        else // two equation models
        {
            static_assert(numTurbulenceEq(ModelTraits::turbulenceModel()) == 2, "Only reached by 2eq models");
            if(isOutlet_(pos))
            {
                values.setOutflow(Indices::turbulentKineticEnergyEqIdx);
                values.setOutflow(Indices::dissipationEqIdx);
            }
            else
            {
                // walls and inflow
                values.setDirichlet(Indices::turbulentKineticEnergyIdx);
                values.setDirichlet(Indices::dissipationIdx);
            }
        }
    }

    template<class Element, class FVElementGeometry, class SubControlVolume>
    bool isDirichletCell_([[maybe_unused]] const Element& element,
                          const FVElementGeometry& fvGeometry,
                          [[maybe_unused]] const SubControlVolume& scv,
                          const int& pvIdx) const
    {
        if constexpr (ModelTraits::turbulenceModel() == TurbulenceModel::kepsilon)
        {
            const auto eIdx = fvGeometry.gridGeometry().elementMapper().index(element);
            // For the kepsilon model we set fixed values within the matching point and at the wall
            if (this->inNearWallRegion(eIdx))
                return pvIdx == Indices::turbulentKineticEnergyEqIdx || pvIdx == Indices::dissipationEqIdx;
            if (this->isMatchingPoint(eIdx))
                return pvIdx == Indices::dissipationEqIdx;
            return false;
        }
        else if constexpr (ModelTraits::turbulenceModel() == TurbulenceModel::komega ||
                           ModelTraits::turbulenceModel() == TurbulenceModel::sst)
        {
            // For the komega model we set a fixed dissipation (omega) for all cells at the wall
            for (const auto& scvf : scvfs(fvGeometry))
                if (this->boundaryTypes(element, scvf).hasWall() && pvIdx == Indices::dissipationIdx)
                    return true;
            return false;
        }
        else
            return ParentType::isDirichletCell(element, fvGeometry, scv, pvIdx);
    }

    //! Specialization for the kepsilon and komega
    template<class Element, class SubControlVolume>
    PrimaryVariables dirichletTurbulentTwoEq_(const Element& element,
                                              const SubControlVolume& scv) const
    {
        const auto globalPos = scv.center();
        PrimaryVariables values(initialAtPos(globalPos));
        unsigned int  elementIdx = this->gridGeometry().elementMapper().index(element);

        if constexpr (ModelTraits::turbulenceModel() == TurbulenceModel::kepsilon)
        {
            // For the kepsilon model we set a fixed value for the turbulent kinetic energy and the dissipation
            values[Indices::turbulentKineticEnergyEqIdx] = this->turbulentKineticEnergyWallFunction(elementIdx);
            values[Indices::dissipationEqIdx] = this->dissipationWallFunction(elementIdx);
            return values;
        }
        else
        {
            static_assert(ModelTraits::turbulenceModel() == TurbulenceModel::komega
                       || ModelTraits::turbulenceModel() == TurbulenceModel::sst, "Only valid for Komega");
            // For the komega model we set a fixed value for the dissipation
            const auto wallDistance = ParentType::wallDistance(elementIdx);
            using std::pow;
            values[Indices::dissipationEqIdx] = 6.0 * ParentType::kinematicViscosity(elementIdx)
                                                    / (ParentType::betaOmega() * wallDistance * wallDistance);
            return values;
        }
    }

    Scalar eps_;
    Scalar inletVelocity_;
    Scalar inletTemperature_;
    Scalar wallTemperature_;
    Scalar initializationTime_;
    Scalar viscosityTilde_;
    Scalar turbulentKineticEnergy_;
    Scalar dissipation_;
    std::string turbulenceModelName_;
    TimeLoopPtr timeLoop_;
};
} // end namespace Dumux

#endif
