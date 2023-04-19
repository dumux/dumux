// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup RANSNCTests
 * \brief Flat plate test for the multi-component staggered grid Reynolds-averaged Navier-Stokes model.
 */
#ifndef DUMUX_RANS_NC_TEST_PROBLEM_HH
#define DUMUX_RANS_NC_TEST_PROBLEM_HH

#include <dumux/common/properties.hh>
#include <dumux/common/parameters.hh>
#include <dumux/common/timeloop.hh>
#include <dumux/common/numeqvector.hh>

#include <dumux/freeflow/rans/boundarytypes.hh>
#include <dumux/freeflow/turbulencemodel.hh>
#include <dumux/freeflow/turbulenceproperties.hh>
#include <dumux/freeflow/rans/zeroeq/problem.hh>
#include <dumux/freeflow/rans/oneeq/problem.hh>
#include <dumux/freeflow/rans/twoeq/komega/problem.hh>
#include <dumux/freeflow/rans/twoeq/sst/problem.hh>
#include <dumux/freeflow/rans/twoeq/lowrekepsilon/problem.hh>
#include <dumux/freeflow/rans/twoeq/kepsilon/problem.hh>

namespace Dumux {

/*!
 * \ingroup RANSNCTests
 * \brief  Test problem for the one-phase model.
 *
 * Dry air is entering from the left side and flows above a 1-D a flat plate.
 * In the middle of the inlet, water vapor is injected, which spreads by turbulent diffusion.
 * For the non-isothermal model the bottom has a constant temperature
 * which is higher than the initial and inlet temperature.
 */
template <class TypeTag>
class FlatPlateNCTestProblem : public RANSProblem<TypeTag>
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

    static constexpr auto transportEqIdx = Indices::conti0EqIdx + 1;
    static constexpr auto transportCompIdx = Indices::conti0EqIdx + 1;

public:
    FlatPlateNCTestProblem(std::shared_ptr<const GridGeometry> gridGeometry)
    : ParentType(gridGeometry), eps_(1e-6)
    {
        inletVelocity_ = getParam<Scalar>("Problem.InletVelocity", 0.1);
        inletTemperature_ = getParam<Scalar>("Problem.InletTemperature", 283.15);
        wallTemperature_ = getParam<Scalar>("Problem.WallTemperature", 313.15);
        inletMoleFraction_ = getParam<Scalar>("Problem.InletMoleFraction", 1e-3);

        FluidSystem::init();
        Dumux::TurbulenceProperties<Scalar, dimWorld, true> turbulenceProperties;
        FluidState fluidState;
        const auto phaseIdx = 0;
        fluidState.setPressure(phaseIdx, 1e5);
        fluidState.setTemperature(this->spatialParams().temperatureAtPos({}));
        fluidState.setMassFraction(phaseIdx, phaseIdx, 1.0);
        Scalar density = FluidSystem::density(fluidState, phaseIdx);
        Scalar kinematicViscosity = FluidSystem::viscosity(fluidState, phaseIdx) / density;
        Scalar diameter = this->gridGeometry().bBoxMax()[1] - this->gridGeometry().bBoxMin()[1];
        viscosityTilde_ = 1e-3 * turbulenceProperties.viscosityTilde(inletVelocity_, diameter, kinematicViscosity);
        turbulentKineticEnergy_ = turbulenceProperties.turbulentKineticEnergy(inletVelocity_, diameter, kinematicViscosity);
        if (ModelTraits::turbulenceModel() == TurbulenceModel::komega || ModelTraits::turbulenceModel() == TurbulenceModel::sst)
            dissipation_ = turbulenceProperties.dissipationRate(inletVelocity_, diameter, kinematicViscosity);
        else
            dissipation_ = turbulenceProperties.dissipation(inletVelocity_, diameter, kinematicViscosity);
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

        // turbulence model-specific boundary types
        setBcTypes_(values, globalPos);

        if(isInlet_(globalPos))
        {
            values.setDirichlet(Indices::velocityXIdx);
            values.setDirichlet(Indices::velocityYIdx);
            values.setDirichlet(transportCompIdx);
#if NONISOTHERMAL
            values.setDirichlet(Indices::temperatureIdx);
#endif
        }
        else if(isOutlet_(globalPos))
        {
            values.setDirichlet(Indices::pressureIdx);
            values.setOutflow(transportEqIdx);
#if NONISOTHERMAL
            values.setOutflow(Indices::energyEqIdx);
#endif
        }
        else if(isLowerWall_(globalPos))
        {
            values.setWall();
            values.setNeumann(transportEqIdx);
#if NONISOTHERMAL
            values.setDirichlet(Indices::temperatureIdx);
#endif
        }
        else
            values.setAllSymmetry();

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

        if (isInlet_(globalPos)
            && globalPos[1] > 0.4 * this->gridGeometry().bBoxMax()[1]
            && globalPos[1] < 0.6 * this->gridGeometry().bBoxMax()[1])
        {
            values[transportCompIdx] = (time() > 10.0) ? inletMoleFraction_ : 0.0;
        }

#if NONISOTHERMAL
        values[Indices::temperatureIdx] = (isLowerWall_(globalPos) && time() > 10.0) ? wallTemperature_ : inletTemperature_;
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
     * \brief Evaluates the initial value for a control volume.
     *
     * \param globalPos The global position
     */
    PrimaryVariables initialAtPos(const GlobalPosition &globalPos) const
    {
        PrimaryVariables values(0.0);
        values[Indices::pressureIdx] = 1.0e+5;
        values[transportCompIdx] = 0.0;
#if NONISOTHERMAL
        values[Indices::temperatureIdx] = inletTemperature_;
#endif
        // block velocity profile
        values[Indices::velocityXIdx] = 0.0;
        if (!isLowerWall_(globalPos))
            values[Indices::velocityXIdx] =  inletVelocity_;
        values[Indices::velocityYIdx] = 0.0;

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
    { return globalPos[0] < eps_; }

    bool isOutlet_(const GlobalPosition& globalPos) const
    { return globalPos[0] > this->gridGeometry().bBoxMax()[0] - eps_; }

    bool isLowerWall_(const GlobalPosition& globalPos) const
    { return globalPos[1] < eps_; }

    //! Initial conditions for the komega, kepsilon and lowrekepsilon turbulence models
    void setInitialAtPos_([[maybe_unused]] PrimaryVariables& values,
                          [[maybe_unused]] const GlobalPosition &globalPos) const
    {
        if constexpr (numTurbulenceEq(ModelTraits::turbulenceModel()) == 0) // zero equation models
            return;
        else if constexpr (numTurbulenceEq(ModelTraits::turbulenceModel()) == 1)  // one equation models
        {
            values[Indices::viscosityTildeIdx] = viscosityTilde_;
            if (isLowerWall_(globalPos))
                values[Indices::viscosityTildeIdx] = 0.0;
        }
        else // two equation models
        {
            static_assert(numTurbulenceEq(ModelTraits::turbulenceModel()) == 2, "Only reached by 2eq models");
            values[Indices::turbulentKineticEnergyIdx] = turbulentKineticEnergy_;
            values[Indices::dissipationIdx] = dissipation_;
            if (isLowerWall_(globalPos))
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
                           ModelTraits::turbulenceModel() == TurbulenceModel::sst )
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
            static_assert(ModelTraits::turbulenceModel() == TurbulenceModel::komega ||
                          ModelTraits::turbulenceModel() == TurbulenceModel::sst, "Only valid for SST and KOmega Models");
            // For the komega model we set a fixed value for the dissipation
            const auto wallDistance = ParentType::wallDistance(elementIdx);
            using std::pow;
            values[Indices::dissipationEqIdx] = 6.0 * ParentType::kinematicViscosity(elementIdx)
                                                    / (ParentType::betaOmega() * wallDistance * wallDistance);
            return values;
        }
    }

    const Scalar eps_;
    Scalar inletVelocity_;
    Scalar inletMoleFraction_;
    Scalar inletTemperature_;
    Scalar wallTemperature_;
    Scalar viscosityTilde_;
    Scalar turbulentKineticEnergy_;
    Scalar dissipation_;
    TimeLoopPtr timeLoop_;
    std::string turbulenceModelName_;
};
} // end namespace Dumux

#endif
