// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup TwoPTests
 * \brief Buckley-Leverett test problem for incompressible immiscible two-phase flow.
 */
#ifndef DUMUX_TEST_TWOP_BUCKLEYLEVERETT_PROBLEM_HH
#define DUMUX_TEST_TWOP_BUCKLEYLEVERETT_PROBLEM_HH

#include <dumux/common/boundarytypes.hh>
#include <dumux/common/numeqvector.hh>
#include <dumux/common/parameters.hh>
#include <dumux/common/properties.hh>

#include <dumux/porousmediumflow/problem.hh>

namespace Dumux {

template<class TypeTag>
class BuckleyLeverettProblem : public PorousMediumFlowProblem<TypeTag>
{
    using ParentType = PorousMediumFlowProblem<TypeTag>;
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using GridView = typename GridGeometry::GridView;
    using Element = typename GridView::template Codim<0>::Entity;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using FluidSystem = GetPropType<TypeTag, Properties::FluidSystem>;
    using FluidState = GetPropType<TypeTag, Properties::FluidState>;
    using PrimaryVariables = GetPropType<TypeTag, Properties::PrimaryVariables>;
    using NumEqVector = Dumux::NumEqVector<PrimaryVariables>;
    using BoundaryTypes = Dumux::BoundaryTypes<GetPropType<TypeTag, Properties::ModelTraits>::numEq()>;
    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;
    using Indices = typename GetPropType<TypeTag, Properties::ModelTraits>::Indices;

    static constexpr int wettingPhaseIdx = FluidSystem::phase0Idx;
    static constexpr int nonwettingPhaseIdx = FluidSystem::phase1Idx;

public:
    BuckleyLeverettProblem(std::shared_ptr<const GridGeometry> gridGeometry)
    : ParentType(gridGeometry)
    , injectionPressure_(getParam<Scalar>("Problem.InjectionPressure"))
    , referencePressure_(getParam<Scalar>("Problem.ReferencePressure"))
    , totalVelocity_(getParam<Scalar>("Problem.TotalVelocity"))
    {}

    BoundaryTypes boundaryTypesAtPos(const GlobalPosition& globalPos) const
    {
        BoundaryTypes values;

        if (onLeftBoundary_(globalPos))
            values.setAllDirichlet();
        else
            values.setAllNeumann();

        return values;
    }

    PrimaryVariables dirichletAtPos(const GlobalPosition& globalPos) const
    {
        PrimaryVariables values;
        values[Indices::pressureIdx] = injectionPressure_;
        values[Indices::saturationIdx] = residualNonwettingSaturation_(globalPos);
        return values;
    }

    NumEqVector neumannAtPos(const GlobalPosition& globalPos) const
    {
        NumEqVector values(0.0);

        if (onRightBoundary_(globalPos))
            values[Indices::conti0EqIdx + FluidSystem::comp1Idx]
                = totalVelocity_*density_(globalPos, nonwettingPhaseIdx);

        return values;
    }

    PrimaryVariables initialAtPos(const GlobalPosition& globalPos) const
    {
        PrimaryVariables values;
        values[Indices::pressureIdx] = injectionPressure_;
        values[Indices::saturationIdx] = 1.0 - residualWettingSaturation_(globalPos);
        return values;
    }

    Scalar referencePressure() const
    { return referencePressure_; }

    Scalar totalVelocity() const
    { return totalVelocity_; }

private:
    Scalar density_(const GlobalPosition& globalPos, int phaseIdx) const
    {
        FluidState fluidState;
        fluidState.setTemperature(this->spatialParams().temperatureAtPos(globalPos));
        fluidState.setPressure(wettingPhaseIdx, referencePressure_);
        fluidState.setPressure(nonwettingPhaseIdx, referencePressure_);
        return FluidSystem::density(fluidState, phaseIdx);
    }

    Scalar residualWettingSaturation_(const GlobalPosition& globalPos) const
    { return this->spatialParams().fluidMatrixInteractionAtPos(globalPos).pcSwCurve().effToAbsParams().swr(); }

    Scalar residualNonwettingSaturation_(const GlobalPosition& globalPos) const
    { return this->spatialParams().fluidMatrixInteractionAtPos(globalPos).pcSwCurve().effToAbsParams().snr(); }

    bool onLeftBoundary_(const GlobalPosition& globalPos) const
    { return globalPos[0] < this->gridGeometry().bBoxMin()[0] + eps_; }

    bool onRightBoundary_(const GlobalPosition& globalPos) const
    { return globalPos[0] > this->gridGeometry().bBoxMax()[0] - eps_; }

    Scalar injectionPressure_;
    Scalar referencePressure_;
    Scalar totalVelocity_;

    static constexpr Scalar eps_ = 1e-6;
};

} // end namespace Dumux

#endif
