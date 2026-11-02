// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup NavierStokesTests
 * \brief Turbulent channel flow test using the two-equation, wall-function k-epsilon RANS model.
 */

#ifndef DUMUX_TEST_RANS_KEPSILON_CHANNEL_PROBLEM_HH
#define DUMUX_TEST_RANS_KEPSILON_CHANNEL_PROBLEM_HH

#include <cmath>
#include <algorithm>

#include <dumux/common/parameters.hh>
#include <dumux/common/properties.hh>

#include <dumux/freeflow/navierstokes/boundarytypes.hh>
#include <dumux/freeflow/navierstokes/momentum/fluxhelper.hh>
#include <dumux/freeflow/navierstokes/scalarfluxhelper.hh>
#include <dumux/freeflow/rans/kepsilon/advectiveflux.hh>
#include <dumux/freeflow/turbulenceproperties.hh>

namespace Dumux {

/*!
 * \ingroup NavierStokesTests
 * \brief Test problem for turbulent channel flow using the two-equation, wall-function k-epsilon
 *        RANS model.
 *
 * Flow from left to right in a two-dimensional channel, with a uniform (top-hat) inflow
 * velocity profile at the inlet (left), a fixed pressure at the outlet (right), and walls
 * at the top and bottom. The turbulent kinetic energy k and the dissipation rate epsilon
 * are fixed at the inlet and leave the domain at the outlet.
 *
 * \note At the walls the pressure, k and epsilon are Neumann. The near-wall values of k and
 *       epsilon are imposed through the internal Dirichlet constraints of
 *       Dumux::KEpsilonMassProblem, see that class for the near-wall region and matching-point
 *       treatment. The momentum is no-slip Dirichlet everywhere except at the outlet and at
 *       the wall face of the matching-point cells, where
 *       Dumux::KEpsilonMomentumProblem::boundaryTypes() switches to Neumann so that the
 *       log-law wall shear flux can be applied (see neumann()).
 */
template<class TypeTag, class BaseProblem>
class RANSKEpsilonChannelTestProblem : public BaseProblem
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

    RANSKEpsilonChannelTestProblem(std::shared_ptr<const GridGeometry> gridGeometry, std::shared_ptr<CouplingManager> couplingManager)
    : ParentType(gridGeometry, couplingManager)
    {
        inletVelocity_ = getParam<Scalar>("Problem.InletVelocity");
        outletPressure_ = getParam<Scalar>("Problem.OutletPressure", 1.0e5);
#if NONISOTHERMAL
        inletTemperature_ = getParam<Scalar>("Problem.InletTemperature", 283.15);
        wallTemperature_ = getParam<Scalar>("Problem.WallTemperature", 303.15);
#endif

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

            // empirical estimates of the inflow values from the inlet velocity, the channel
            // height and the kinematic viscosity
            Dumux::TurbulenceProperties<Scalar, dimWorld, true> turbulenceProperties;
            turbulentKineticEnergy_ = turbulenceProperties.turbulentKineticEnergy(inletVelocity_, 2.0*halfHeight, kinematicViscosity);
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
            // only queried at the inlet, the walls are Neumann
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
            if (isOnWallAtPos(scvf.ipGlobal()))
            {
                // Only reached at the wall face of the matching-point cells, where
                // KEpsilonMomentumProblem::boundaryTypes() sets a Neumann condition
                values = BoundaryFluxes(this->wallFunctionMomentumFlux(element, fvGeometry, elemVolVars, elemFluxVarsCache, scvf));
            }
            else
            {
                using FluxHelper = NavierStokesMomentumBoundaryFlux<typename GridGeometry::DiscretizationMethod>;
                values = FluxHelper::fixedPressureMomentumFlux(*this, fvGeometry, scvf, elemVolVars, elemFluxVarsCache,
                                                                outletPressure_, true /*zeroNormalVelocityGradient*/);
            }
        }
        else
        {
            const auto& globalPos = scvf.ipGlobal();
            if (isOutlet_(globalPos))
            {
                using FluxHelper = NavierStokesScalarBoundaryFluxHelper<AdvectiveFlux<ModelTraits>>;
                values = FluxHelper::scalarOutflowFlux(*this, element, fvGeometry, scvf, elemVolVars);
            }
            // At the walls, k and epsilon have zero flux. Their near-wall values are
            // imposed by the internal Dirichlet constraints of Dumux::KEpsilonMassProblem.
#if NONISOTHERMAL
            else if (isOnWallAtPos(globalPos))
            {
                const auto& insideScv = fvGeometry.scv(scvf.insideScvIdx());
                const auto& insideVars = elemVolVars[insideScv];
                if (this->isMatchingPoint(insideScv.elementIndex()))
                {
                    // Jayatilleke log-law energy wall function at the matching point,
                    // which replaces the conductive wall flux
                    values[Indices::energyEqIdx] = this->wallFunctionEnergyFlux(insideScv.elementIndex(), insideVars, wallTemperature_);
                }
                else
                {
                    // The wall temperature is imposed weakly through a conductive heat flux
                    const auto distance = (insideScv.dofPosition() - scvf.ipGlobal()).two_norm();
                    values[Indices::energyEqIdx] = insideVars.effectiveThermalConductivity()
                        *(insideVars.temperature() - wallTemperature_)/distance;
                }
            }
#endif
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
            values[Indices::dissipationIdx] = dissipation_;
#if NONISOTHERMAL
            values[Indices::temperatureIdx] = isOnWallAtPos(globalPos) ? wallTemperature_ : inletTemperature_;
#endif
        }

        return values;
    }

    // \}

    /*!
     * \brief Returns whether the given position lies on a solid (no-slip) wall of the channel.
     * \note Required by Dumux::RANSMomentumProblem for the wall distance computation, by
     *       Dumux::KEpsilonMassProblem for the near-wall region and matching-point
     *       classification, and by the boundary treatment of this problem.
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
    Scalar turbulentKineticEnergy_;
    Scalar dissipation_;
#if NONISOTHERMAL
    Scalar inletTemperature_;
    Scalar wallTemperature_;
#endif
};

} // end namespace Dumux

#endif
