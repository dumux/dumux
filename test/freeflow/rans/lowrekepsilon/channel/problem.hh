// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup NavierStokesTests
 * \brief Turbulent channel flow test using the two-equation, low-Reynolds-number k-epsilon
 *        RANS model.
 */

#ifndef DUMUX_TEST_RANS_LOWREKEPSILON_CHANNEL_PROBLEM_HH
#define DUMUX_TEST_RANS_LOWREKEPSILON_CHANNEL_PROBLEM_HH

#include <cmath>
#include <algorithm>

#include <dumux/common/parameters.hh>
#include <dumux/common/properties.hh>

#include <dumux/freeflow/navierstokes/boundarytypes.hh>
#include <dumux/freeflow/navierstokes/momentum/fluxhelper.hh>
#include <dumux/freeflow/navierstokes/scalarfluxhelper.hh>
#include <dumux/freeflow/rans/lowrekepsilon/advectiveflux.hh>
#include <dumux/freeflow/turbulenceproperties.hh>

namespace Dumux {

/*!
 * \ingroup NavierStokesTests
 * \brief Test problem for turbulent channel flow using the two-equation, low-Reynolds-number
 *        k-epsilon RANS model, mirroring test/freeflow/rans/komega/channel's geometry/BCs/fluid,
 *        with k and epsilon-tilde playing the role k and omega play there - matching
 *        releases/3.10:test/freeflow/rans/problem.hh, which has no lowrekepsilon-specific
 *        branch at all (confirmed this session): both k and epsilon-tilde get a plain Dirichlet
 *        wall value of zero there, with no internal-cell constraint whatsoever.
 *
 * \note Walls are kept Neumann (not Dirichlet) for pressure, k, *and* epsilon-tilde, exactly as
 *       k is treated in the k-omega test: dumux/assembly/cclocalassembler.hh rejects mixed
 *       Dirichlet/Neumann boundary conditions across equations at the same CCTpfa face, so the
 *       wall values k=0, epsilon-tilde=0 are enforced weakly (a one-sided diffusive/Robin-style
 *       flux) instead - see neumann() below. Unlike k-omega/k-epsilon, no internal Dirichlet
 *       constraint is needed anywhere: the low-Reynolds damping functions (f_mu, D, E) let both
 *       transported quantities be solved by their ordinary PDE all the way to the wall.
 */
template<class TypeTag, class BaseProblem>
class RANSLowReKEpsilonChannelTestProblem : public BaseProblem
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

    RANSLowReKEpsilonChannelTestProblem(std::shared_ptr<const GridGeometry> gridGeometry, std::shared_ptr<CouplingManager> couplingManager)
    : ParentType(gridGeometry, couplingManager)
    {
        inletVelocity_ = getParam<Scalar>("Problem.InletVelocity");
        outletPressure_ = getParam<Scalar>("Problem.OutletPressure", 1.0e5);
#if NONISOTHERMAL
        inletTemperature_ = getParam<Scalar>("Problem.InletTemperature", 283.15);
        wallTemperature_ = getParam<Scalar>("Problem.WallTemperature", 303.15);
#endif

        // k/epsilon inflow estimates are only meaningful (and only compile: FluidState is only
        // defined) on the mass domain, exactly as for the k-omega test's turbulentKineticEnergy_.
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

            Dumux::TurbulenceProperties<Scalar, dimWorld, true> turbulenceProperties;
            turbulentKineticEnergy_ = turbulenceProperties.turbulentKineticEnergy(inletVelocity_, 2.0*halfHeight, kinematicViscosity);
            // plain dissipation() (not dissipationRate()) - this is a k-epsilon-family model,
            // not a k-omega-family one, matching the k-epsilon (wall-function) test's estimator.
            dissipation_ = turbulenceProperties.dissipation(inletVelocity_, 2.0*halfHeight, kinematicViscosity);
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

            if (isInlet_(globalPos))
            {
                values.setDirichlet(Indices::pressureIdx);
                values.setDirichlet(Indices::turbulentKineticEnergyIdx);
                values.setDirichlet(Indices::dissipationIdx);
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
            values[Indices::velocityXIdx] = inletVelocityProfile_(globalPos[1]);
            values[Indices::velocityYIdx] = 0.0;
        }
        else
        {
            values[Indices::pressureIdx] = this->couplingManager().cellPressure(element, scvf);
            // Only ever queried at the inlet - see boundaryTypesAtPos() (walls are Neumann).
            values[Indices::turbulentKineticEnergyIdx] = turbulentKineticEnergy_;
            values[Indices::dissipationIdx] = dissipation_;
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
                // Weak (Robin-form) equivalent of true Dirichlet k=0, epsilon-tilde=0 wall
                // conditions - the CCTpfa mixed-boundary-condition restriction (see class docs
                // above) rules out an actual Dirichlet face BC here, exactly as for k-omega's k.
                // Unlike k-omega/k-epsilon, there is no internal-cell constraint for either
                // quantity in this model - both wall values are enforced only through this flux.
                const auto& insideScv = fvGeometry.scv(scvf.insideScvIdx());
                const auto& insideVars = elemVolVars[insideScv];
                const auto distance = (insideScv.dofPosition() - scvf.ipGlobal()).two_norm();
                const auto diffCoeffK = insideVars.dynamicEddyViscosity()/insideVars.sigmaK() + insideVars.viscosity();
                const auto diffCoeffEpsilon = insideVars.dynamicEddyViscosity()/insideVars.sigmaEpsilon() + insideVars.viscosity();
                values[Indices::turbulentKineticEnergyEqIdx] = diffCoeffK*insideVars.turbulentKineticEnergy()/distance;
                values[Indices::dissipationEqIdx] = diffCoeffEpsilon*insideVars.dissipationTilde()/distance;
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
            values[Indices::turbulentKineticEnergyIdx] = isOnWallAtPos(globalPos) ? 0.0 : turbulentKineticEnergy_;
            values[Indices::dissipationIdx] = isOnWallAtPos(globalPos) ? 0.0 : dissipation_;
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

    //! Uniform (top-hat) inflow profile: zero at the walls, inletVelocity_ everywhere else -
    //! see test/freeflow/rans/zeroeq/channel/problem.hh for the rationale.
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
    Scalar turbulentKineticEnergy_;
    Scalar dissipation_;
#if NONISOTHERMAL
    Scalar inletTemperature_;
    Scalar wallTemperature_;
#endif
};

} // end namespace Dumux

#endif
