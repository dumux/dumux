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
 *        RANS model, mirroring test/freeflow/rans/komega/channel's geometry/BCs/fluid, with
 *        additional boundary/initial conditions for k and epsilon, and a wall-function-aware
 *        momentum boundary treatment - matching
 *        releases/3.10:test/freeflow/rans/problem.hh's PipeLauferProblem (kepsilon branch).
 *
 * \note Walls are Neumann (not Dirichlet) for pressure, k, *and* epsilon: the actual near-wall
 *       values are enforced entirely through Dumux::KEpsilonMassProblem's internal (whole-cell)
 *       Dirichlet constraints, not through boundary conditions - see that class and
 *       whatisimplemented.md for the near-wall-region/matching-point mechanism. Momentum is
 *       Dirichlet no-slip everywhere except the outlet *and* the matching-point cell's own wall
 *       face, where Dumux::KEpsilonMomentumProblem::boundaryTypes() switches to Neumann so the
 *       log-law wall-shear flux can be applied instead (see neumann() below).
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

        if constexpr (!ParentType::isMomentumProblem())
        {
            using FluidSystem = GetPropType<TypeTag, Properties::FluidSystem>;
            using FluidState = GetPropType<TypeTag, Properties::FluidState>;
            FluidState fluidState;
            fluidState.setPressure(0, 1e5);
            fluidState.setTemperature(283.15);
            const Scalar density = FluidSystem::density(fluidState, 0);
            const Scalar kinematicViscosity = FluidSystem::viscosity(fluidState, 0)/density;

            const Scalar halfHeight = 0.5*(this->gridGeometry().bBoxMax()[1] - this->gridGeometry().bBoxMin()[1]);

            // Note: plain dissipation(), not dissipationRate() - that estimator is
            // komega/SST-specific (confirmed from the old test problem's own branching).
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
            values[Indices::velocityXIdx] = powerLawProfile_(globalPos[1]);
            values[Indices::velocityYIdx] = 0.0;
        }
        else
        {
            values[Indices::pressureIdx] = this->couplingManager().cellPressure(element, scvf);
            // Only ever queried at the inlet - see boundaryTypesAtPos() (walls are Neumann).
            values[Indices::turbulentKineticEnergyIdx] = turbulentKineticEnergy_;
            values[Indices::dissipationIdx] = dissipation_;
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
                // Only reached at the matching-point cell's own wall face - see
                // KEpsilonMomentumProblem::boundaryTypes(), which is the only place this
                // becomes a Neumann condition in the first place.
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
            // At walls, both k and epsilon are Neumann with zero flux here - their actual
            // near-wall values are enforced entirely via the internal Dirichlet constraints in
            // Dumux::KEpsilonMassProblem, not through this boundary treatment.
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
            values[Indices::turbulentKineticEnergyIdx] = isOnWallAtPos(globalPos) ? 0.0 : turbulentKineticEnergy_;
            values[Indices::dissipationIdx] = dissipation_;
        }

        return values;
    }

    // \}

    /*!
     * \brief Returns whether the given position lies on a solid (no-slip) wall of the channel.
     * \note Required by Dumux::RANSMomentumProblem for the wall-distance computation, by
     *       Dumux::KEpsilonMassProblem for the near-wall-region/matching-point classification,
     *       and directly here (both problem instances) to dispatch boundary treatment.
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
    Scalar turbulentKineticEnergy_;
    Scalar dissipation_;
};

} // end namespace Dumux

#endif
