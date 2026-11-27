// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup NavierStokesTests
 * \brief Turbulent channel flow test using the one-equation (Spalart-Allmaras) RANS model.
 */

#ifndef DUMUX_TEST_RANS_ONEEQ_CHANNEL_PROBLEM_HH
#define DUMUX_TEST_RANS_ONEEQ_CHANNEL_PROBLEM_HH

#include <cmath>
#include <algorithm>

#include <dumux/common/parameters.hh>
#include <dumux/common/properties.hh>

#include <dumux/freeflow/navierstokes/boundarytypes.hh>
#include <dumux/freeflow/navierstokes/momentum/fluxhelper.hh>
#include <dumux/freeflow/navierstokes/scalarfluxhelper.hh>
#include <dumux/freeflow/rans/oneeq/advectiveflux.hh>
#include <dumux/freeflow/turbulenceproperties.hh>

namespace Dumux {

/*!
 * \ingroup NavierStokesTests
 * \brief Test problem for turbulent channel flow using the one-equation (Spalart-Allmaras)
 *        RANS model.
 *
 * Flow from left to right in a two-dimensional channel, with a uniform (top-hat) inflow
 * velocity profile at the inlet (left), a fixed pressure at the outlet (right), and no-slip
 * walls at the top and bottom. The working viscosity ν̃ is fixed at the inlet, vanishes
 * at the walls, and leaves the domain at the outlet.
 */
template<class TypeTag, class BaseProblem>
class RANSOneEqChannelTestProblem : public BaseProblem
{
    using ParentType = BaseProblem;

    using BoundaryTypes = typename ParentType::BoundaryTypes;
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using FVElementGeometry = typename GridGeometry::LocalView;
    using SubControlVolumeFace = typename GridGeometry::SubControlVolumeFace;
    using InitialValues = typename ParentType::InitialValues;
    using DirichletValues = typename ParentType::DirichletValues;
    using BoundaryFluxes = typename ParentType::BoundaryFluxes;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using ModelTraits = GetPropType<TypeTag, Properties::ModelTraits>;

    using Element = typename FVElementGeometry::Element;
    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;

    using CouplingManager = GetPropType<TypeTag, Properties::CouplingManager>;

    static constexpr auto dimWorld = GridGeometry::GridView::dimensionworld;

public:
    using Indices = typename GetPropType<TypeTag, Properties::ModelTraits>::Indices;

    RANSOneEqChannelTestProblem(std::shared_ptr<const GridGeometry> gridGeometry, std::shared_ptr<CouplingManager> couplingManager)
    : ParentType(gridGeometry, couplingManager)
    {
        inletVelocity_ = getParam<Scalar>("Problem.InletVelocity");
        outletPressure_ = getParam<Scalar>("Problem.OutletPressure", 1.0e5);
#if NONISOTHERMAL
        inletTemperature_ = getParam<Scalar>("Problem.InletTemperature", 283.15);
        wallTemperature_ = getParam<Scalar>("Problem.WallTemperature", 303.15);
#endif

        // The inflow value of ν̃ is only needed by the mass problem, which owns the ν̃
        // primary variable (the momentum problem has no fluid state). It is an empirical
        // estimate (Dumux::TurbulenceProperties) from the inlet velocity, the channel height
        // and the kinematic viscosity of the fluid at a reference state.
        if constexpr (!ParentType::isMomentumProblem())
        {
            using FluidSystem = GetPropType<TypeTag, Properties::FluidSystem>;
            using FluidState = GetPropType<TypeTag, Properties::FluidState>;
            FluidState fluidState;
            fluidState.setPressure(0, 1e5);
#if NONISOTHERMAL
            fluidState.setTemperature(inletTemperature_);
#else
            fluidState.setTemperature(283.15);
#endif
            const Scalar density = FluidSystem::density(fluidState, 0);
            const Scalar kinematicViscosity = FluidSystem::viscosity(fluidState, 0)/density;

            const Scalar halfHeight = 0.5*(this->gridGeometry().bBoxMax()[1] - this->gridGeometry().bBoxMin()[1]);

            // ideally ν̃ at the inlet of the Spalart-Allmaras model would be zero, so a small
            // fraction of the empirical estimate is used
            Dumux::TurbulenceProperties<Scalar, dimWorld, true> turbulenceProperties;
            viscosityTilde_ = 1e-3*turbulenceProperties.viscosityTilde(inletVelocity_, 2.0*halfHeight, kinematicViscosity);
        }
    }

    /*!
     * \name Boundary conditions
     */
    // \{

    BoundaryTypes boundaryTypesAtPos(const GlobalPosition& globalPos) const
    {
        BoundaryTypes values;

        if constexpr (ParentType::isMomentumProblem())
        {
            values.setDirichlet(Indices::velocityXIdx);
            values.setDirichlet(Indices::velocityYIdx);

            if (isOutlet_(globalPos))
                values.setAllNeumann();
        }
        else
        {
            values.setAllNeumann();

            // Pressure and ν̃ are Dirichlet at the inlet. At the walls the pressure is
            // Neumann, and the cell-centered scheme does not allow different boundary
            // condition types for different equations on the same face. Therefore ν̃ = 0 at
            // the walls is imposed weakly through the diffusive flux in neumann().
            if (isInlet_(globalPos))
            {
                values.setDirichlet(Indices::pressureIdx);
                values.setDirichlet(Indices::viscosityTildeIdx);
#if NONISOTHERMAL
                values.setDirichlet(Indices::energyEqIdx);
#endif
            }
        }

        return values;
    }

    DirichletValues dirichlet(const Element& element, const SubControlVolumeFace& scvf) const
    {
        const auto& globalPos = scvf.ipGlobal();
        DirichletValues values = initialAtPos(globalPos);

        if constexpr (ParentType::isMomentumProblem())
        {
            // The inflow profile is used for every velocity Dirichlet boundary (inlet and walls)
            values[Indices::velocityXIdx] = inletVelocityProfile_(globalPos[1]);
            values[Indices::velocityYIdx] = 0.0;
        }
        else
        {
            values[Indices::pressureIdx] = this->couplingManager().cellPressure(element, scvf);
            // only queried at the inlet, ν̃ is imposed weakly at the walls (see neumann())
            values[Indices::viscosityTildeIdx] = viscosityTilde_;
#if NONISOTHERMAL
            values[Indices::temperatureIdx] = inletTemperature_;
#endif
        }

        return values;
    }

    template<class ElementVolumeVariables, class ElementFluxVariablesCache>
    BoundaryFluxes neumann(const Element& element,
                           const FVElementGeometry& fvGeometry,
                           const ElementVolumeVariables& elemVolVars,
                           const ElementFluxVariablesCache& elemFluxVarsCache,
                           const SubControlVolumeFace& scvf) const
    {
        BoundaryFluxes values(0.0);

        if constexpr (ParentType::isMomentumProblem())
        {
            using FluxHelper = NavierStokesMomentumBoundaryFlux<typename GridGeometry::DiscretizationMethod>;
            values = FluxHelper::fixedPressureMomentumFlux(*this, fvGeometry, scvf, elemVolVars, elemFluxVarsCache,
                                                            outletPressure_, true /*zeroNormalVelocityGradient*/);
        }
        else
        {
            const auto& globalPos = scvf.ipGlobal();
            if (isOutlet_(globalPos))
            {
                using FluxHelper = NavierStokesScalarBoundaryFluxHelper<AdvectiveFlux<ModelTraits>>;
                values = FluxHelper::scalarOutflowFlux(*this, element, fvGeometry, scvf, elemVolVars);
            }
            else if (isOnWallAtPos(globalPos))
            {
                // Weak imposition of ν̃ = 0 at the wall (see boundaryTypesAtPos()): the
                // two-point diffusive flux towards the wall value, directed out of the cell
                const auto& insideScv = fvGeometry.scv(scvf.insideScvIdx());
                const auto& insideVars = elemVolVars[insideScv];
                const auto distance = (insideScv.dofPosition() - scvf.ipGlobal()).two_norm();
                const auto diffCoeff = (insideVars.viscosity() + insideVars.density()*insideVars.viscosityTilde())/insideVars.sigma();
                values[Indices::viscosityTildeEqIdx] = diffCoeff*insideVars.viscosityTilde()/distance;
#if NONISOTHERMAL
                values[Indices::energyEqIdx] = insideVars.effectiveThermalConductivity()
                    *(insideVars.temperature() - wallTemperature_)/distance;
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

    InitialValues initialAtPos(const GlobalPosition& globalPos) const
    {
        InitialValues values(0.0);

        if constexpr (ParentType::isMomentumProblem())
        {
            values[Indices::velocityXIdx] = inletVelocityProfile_(globalPos[1]);
            values[Indices::velocityYIdx] = 0.0;
        }
        else
        {
            values[Indices::pressureIdx] = outletPressure_;
            values[Indices::viscosityTildeIdx] = isOnWallAtPos(globalPos) ? 0.0 : viscosityTilde_;
#if NONISOTHERMAL
            values[Indices::temperatureIdx] = isOnWallAtPos(globalPos) ? wallTemperature_ : inletTemperature_;
#endif
        }

        return values;
    }

    // \}

    /*!
     * \brief Returns whether the given position lies on a solid (no-slip) wall of the channel.
     * \note Required by Dumux::RANSMomentumProblem for the wall-distance computation.
     */
    bool isOnWallAtPos(const GlobalPosition& globalPos) const
    { return globalPos[1] < this->gridGeometry().bBoxMin()[1] + eps_ || globalPos[1] > this->gridGeometry().bBoxMax()[1] - eps_; }

    Scalar referencePressure(const Element& element,
                             const FVElementGeometry& fvGeometry,
                             const SubControlVolumeFace& scvf) const
    { return outletPressure_; }

private:
    bool isInlet_(const GlobalPosition& globalPos) const
    { return globalPos[0] < this->gridGeometry().bBoxMin()[0] + eps_; }

    bool isOutlet_(const GlobalPosition& globalPos) const
    { return globalPos[0] > this->gridGeometry().bBoxMax()[0] - eps_; }

    //! Uniform (top-hat) inflow profile: zero at the walls, the inlet velocity elsewhere
    Scalar inletVelocityProfile_(const Scalar y) const
    {
        const Scalar yMin = this->gridGeometry().bBoxMin()[1];
        const Scalar yMax = this->gridGeometry().bBoxMax()[1];
        const Scalar distanceToWall = std::min(y - yMin, yMax - y);
        return distanceToWall > eps_ ? inletVelocity_ : 0.0;
    }

    static constexpr Scalar eps_ = 1e-6;
    Scalar inletVelocity_;
    Scalar outletPressure_;
    Scalar viscosityTilde_;
#if NONISOTHERMAL
    Scalar inletTemperature_;
    Scalar wallTemperature_;
#endif
};

} // end namespace Dumux

#endif
