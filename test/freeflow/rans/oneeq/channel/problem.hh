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
 *        RANS model, mirroring test/freeflow/rans/zeroeq/channel's geometry/BCs/fluid, with
 *        the same 1/7 power-law inflow profile and the same graded, wall-refined mesh, but
 *        also setting boundary/initial conditions for the working viscosity ν̃ (Dirichlet at
 *        walls/inlet, outflow at the outlet) - matching
 *        releases/3.10:test/freeflow/rans/problem.hh's PipeLauferProblem (one-eq branch).
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

        // The inflow viscosityTilde_ value is only used on the mass domain (which owns the ν̃
        // primary variable) - the momentum domain has no FluidState property, so this is
        // skipped there. Computed the same way releases/3.10:test/freeflow/rans/problem.hh did:
        // an empirical estimate (Dumux::TurbulenceProperties) from the inlet velocity, a
        // reference length, and the fluid's kinematic viscosity at a reference state.
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

            // ideally the viscosityTilde parameter as inflow for the Spalart-Allmaras model
            // should be zero - matching releases/3.10:test/freeflow/rans/problem.hh exactly.
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

            // Only the inlet gets a true (per-equation-consistent) Dirichlet condition, since
            // there both pressure *and* ν̃ are Dirichlet - no mixing. At walls, pressure must
            // stay Neumann (as everywhere outside the inlet), so ν̃ is *also* kept Neumann there
            // instead of Dirichlet: dumux/assembly/cclocalresidual.hh's evalFlux() rejects a
            // boundary scvf with different BC types across equations ("mixed boundary
            // conditions ... convert Dirichlet BCs to Robin BCs"). The wall condition ν̃=0 is
            // therefore enforced weakly, as a manually-computed Dirichlet-equivalent diffusive
            // flux inside neumann() below - numerically identical to a true Dirichlet BC, just
            // expressed through the Neumann flux interface CCTpfa requires here.
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
            // 1/7 power-law profile, see test/freeflow/rans/zeroeq/channel/problem.hh for the
            // rationale (a physically consistent approximate turbulent mean-velocity shape).
            values[Indices::velocityXIdx] = powerLawProfile_(globalPos[1]);
            values[Indices::velocityYIdx] = 0.0;
        }
        else
        {
            values[Indices::pressureIdx] = this->couplingManager().cellPressure(element, scvf);
            // Only ever queried at the inlet - see boundaryTypesAtPos() (walls are Neumann,
            // enforced weakly instead, see neumann()).
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
                // Weak (Robin-form) equivalent of a true Dirichlet ν̃=0 wall condition - see
                // boundaryTypesAtPos() above for why this can't be a real per-equation
                // Dirichlet BC here. Numerically the same two-point diffusive flux CCTpfa's
                // own Dirichlet-boundary assembly would produce for this equation, computed
                // manually: coeff*(insideValue - 0)/distance*area, directed out of the cell.
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
            values[Indices::velocityXIdx] = powerLawProfile_(globalPos[1]);
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

    Scalar powerLawProfile_(const Scalar y) const
    {
        const Scalar yMin = this->gridGeometry().bBoxMin()[1];
        const Scalar yMax = this->gridGeometry().bBoxMax()[1];
        const Scalar halfHeight = 0.5*(yMax - yMin);
        const Scalar distanceToWall = std::min(y - yMin, yMax - y);

        using std::pow;
        using std::max;
        return inletVelocity_ * pow(max(distanceToWall, Scalar(0.0))/halfHeight, 1.0/7.0);
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
