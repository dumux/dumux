// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup NavierStokesTests
 * \brief Turbulent channel flow test using the zero-equation (mixing-length) RANS model.
 */

#ifndef DUMUX_TEST_RANS_ZEROEQ_CHANNEL_PROBLEM_HH
#define DUMUX_TEST_RANS_ZEROEQ_CHANNEL_PROBLEM_HH

#include <cmath>
#include <algorithm>

#include <dumux/common/parameters.hh>
#include <dumux/common/properties.hh>

#include <dumux/freeflow/navierstokes/boundarytypes.hh>
#include <dumux/freeflow/navierstokes/momentum/fluxhelper.hh>
#include <dumux/freeflow/navierstokes/scalarfluxhelper.hh>
#include <dumux/freeflow/navierstokes/mass/1p/advectiveflux.hh>

namespace Dumux {

/*!
 * \ingroup NavierStokesTests
 * \brief Test problem for turbulent channel flow using the zero-equation RANS model.
 *
 * Flow from left to right in a two-dimensional channel, with a uniform (top-hat) inflow
 * velocity profile at the inlet (left), a fixed pressure at the outlet (right), and no-slip
 * walls at the top and bottom - matching releases/3.10:test/freeflow/rans/problem.hh's
 * PipeLauferProblem inflow exactly, so that Problem.InletVelocity means the same thing (the
 * bulk/mean inflow velocity) here as it did there. The Reynolds number is chosen high enough
 * for the flow to be turbulent (see params.input), so that the zero-equation eddy-viscosity
 * closure implemented in Dumux::ZeroEqProblem (dumux/freeflow/rans/zeroeq/problem.hh) is
 * exercised.
 */
template<class TypeTag, class BaseProblem>
class RANSZeroEqChannelTestProblem : public BaseProblem
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

public:
    using Indices = typename GetPropType<TypeTag, Properties::ModelTraits>::Indices;

    RANSZeroEqChannelTestProblem(std::shared_ptr<const GridGeometry> gridGeometry, std::shared_ptr<CouplingManager> couplingManager)
    : ParentType(gridGeometry, couplingManager)
    {
        inletVelocity_ = getParam<Scalar>("Problem.InletVelocity");
        outletPressure_ = getParam<Scalar>("Problem.OutletPressure", 1.0e5);
#if NONISOTHERMAL
        inletTemperature_ = getParam<Scalar>("Problem.InletTemperature", 283.15);
        wallTemperature_ = getParam<Scalar>("Problem.WallTemperature", 303.15);
#endif
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
            // The uniform inflow profile is used for every velocity-Dirichlet boundary (inlet
            // and walls): it is zero at y=yMin/yMax (walls) and inletVelocity_ everywhere else
            // (inlet), matching releases/3.10's PipeLauferProblem.
            values[Indices::velocityXIdx] = inletVelocityProfile_(globalPos[1]);
            values[Indices::velocityYIdx] = 0.0;
        }
        else
        {
            values[Indices::pressureIdx] = this->couplingManager().cellPressure(element, scvf);
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
            using FluxHelper = NavierStokesScalarBoundaryFluxHelper<AdvectiveFlux<ModelTraits>>;
            if (isOutlet_(globalPos))
                values = FluxHelper::scalarOutflowFlux(*this, element, fvGeometry, scvf, elemVolVars);
#if NONISOTHERMAL
            else if (isOnWallAtPos(globalPos))
            {
                // Weak (Robin-form) equivalent of a true Dirichlet wall-temperature condition -
                // same trick as the other RANS-NI tests use.
                const auto& insideScv = fvGeometry.scv(scvf.insideScvIdx());
                const auto& insideVars = elemVolVars[insideScv];
                const auto distance = (insideScv.dofPosition() - scvf.ipGlobal()).two_norm();
                values[Indices::energyEqIdx] = insideVars.effectiveThermalConductivity()
                    *(insideVars.temperature() - wallTemperature_)/distance;
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
        InitialValues values;

        if constexpr (ParentType::isMomentumProblem())
        {
            values[Indices::velocityXIdx] = inletVelocityProfile_(globalPos[1]);
            values[Indices::velocityYIdx] = 0.0;
        }
        else
        {
            values[Indices::pressureIdx] = outletPressure_;
#if NONISOTHERMAL
            values[Indices::temperatureIdx] = isOnWallAtPos(globalPos) ? wallTemperature_ : inletTemperature_;
#endif
        }

        return values;
    }

    // \}

    /*!
     * \brief Returns whether the given position lies on a solid (no-slip) wall of the channel,
     *        i.e. the top or bottom boundary but not the inlet/outlet.
     * \note Required by Dumux::ZeroEqProblem for the wall-distance computation.
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

    //! Uniform (top-hat) inflow profile: zero at the walls (y=yMin/yMax), inletVelocity_
    //! everywhere else, matching releases/3.10:test/freeflow/rans/problem.hh's
    //! PipeLauferProblem::initialAtPos() exactly, so that Problem.InletVelocity is the bulk
    //! (flux-averaged) inflow velocity, not a centerline/peak value.
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
#if NONISOTHERMAL
    Scalar inletTemperature_;
    Scalar wallTemperature_;
#endif
};

} // end namespace Dumux

#endif
