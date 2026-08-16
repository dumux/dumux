// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup TwoPTests
 * \brief Analytical saturation profile for the Buckley-Leverett test.
 */
#ifndef DUMUX_TEST_TWOP_BUCKLEYLEVERETT_ANALYTICSOLUTION_HH
#define DUMUX_TEST_TWOP_BUCKLEYLEVERETT_ANALYTICSOLUTION_HH

#include <cmath>
#include <functional>

#include <dune/common/fvector.hh>
#include <dune/istl/bvector.hh>

#include <dumux/common/parameters.hh>
#include <dumux/nonlinear/findscalarroot.hh>

namespace Dumux {

template<class TypeTag>
class BuckleyLeverettAnalyticSolution
{
    using Problem = GetPropType<TypeTag, Properties::Problem>;
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using GridView = typename GridGeometry::GridView;
    using FluidSystem = GetPropType<TypeTag, Properties::FluidSystem>;
    using FluidState = GetPropType<TypeTag, Properties::FluidState>;
    using Element = typename GridView::template Codim<0>::Entity;
    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;

public:
    using ScalarVector = Dune::BlockVector<Scalar>;

    BuckleyLeverettAnalyticSolution(std::shared_ptr<const Problem> problem)
    : problem_(problem)
    , values_(problem->gridGeometry().gridView().size(0))
    , totalVelocity_(getParam<Scalar>("Problem.TotalVelocity"))
    , refPressure_(getParam<Scalar>("Problem.ReferencePressure"))
    {
        FluidState fluidState;
        fluidState.setTemperature(problem_->spatialParams().temperatureAtPos(GlobalPosition{}));
        fluidState.setPressure(FluidSystem::phase0Idx, refPressure_);
        fluidState.setPressure(FluidSystem::phase1Idx, refPressure_);
        viscosityW_ = FluidSystem::viscosity(fluidState, FluidSystem::phase0Idx);
        viscosityN_ = FluidSystem::viscosity(fluidState, FluidSystem::phase1Idx);

        constPorosity_ = problem_->spatialParams().porosityAtPos(GlobalPosition{});

        const auto fluidMatrixInteraction = problem_->spatialParams().fluidMatrixInteractionAtPos(GlobalPosition{});
        swr_ = fluidMatrixInteraction.pcSwCurve().effToAbsParams().swr();
        snr_ = fluidMatrixInteraction.pcSwCurve().effToAbsParams().snr();

        swLeft_ = 1.0 - snr_;
        swRight_ = swr_;
        swShock_ = findSwShock_();
        characteristicSpeed_ = computeShockSpeed_();

        update(0.0);
    }

    Scalar computeSaturation(Scalar x,
                             Scalar time)
    {
        if(time<=0.0)
            return swr_;
        else
        {
            // self-similar speed coordinate xi
            // analytical solution of Buckley-Leverett depends on the ratio x/time, not x and time separately
            Scalar xi = x/time;
            // dfractionalW_dsw_(1.0-snr_) is zero since dfractionalW_dsw_(1.0-snr_) is zero
            Scalar xi_lower = totalVelocity_/constPorosity_ * dfractionalW_dsw_(1.0-snr_);
            if(xi <= xi_lower)
                return 1.0-snr_;
            else if(xi_lower<xi && xi<characteristicSpeed_)
                return computeRarefactionSw_(xi);
            else
                return swr_;
        }
    }

    void update(Scalar time)
    {
        const auto& gridView = problem_->gridGeometry().gridView();
        for (const auto& element : elements(gridView))
        {
            const auto globalPos = element.geometry().center();
            const auto x = globalPos[0];
            const auto eIdx = problem_->gridGeometry().elementMapper().index(element);

            values_[eIdx] = computeSaturation(x, time);
        }
    }

    const ScalarVector& values() const
    { return values_; }

private:

    Scalar fractionalWFunc_(const Scalar sw)
    {
        const auto fluidMatrixInteraction = problem_->spatialParams().fluidMatrixInteractionAtPos(GlobalPosition{});
        const auto mobilityW = fluidMatrixInteraction.krw(sw)/viscosityW_;
        const auto mobilityN = fluidMatrixInteraction.krn(sw)/viscosityN_;
        Scalar fractionalW = mobilityW/(mobilityW + mobilityN);
        return fractionalW;
    }

    Scalar dfractionalW_dsw_(const Scalar sw)
    {
        const auto fluidMatrixInteraction = problem_->spatialParams().fluidMatrixInteractionAtPos(GlobalPosition{});
        const auto mobilityW = fluidMatrixInteraction.krw(sw)/viscosityW_;
        const auto mobilityN = fluidMatrixInteraction.krn(sw)/viscosityN_;

        const Scalar dkrw_dsw = fluidMatrixInteraction.dkrw_dsw(sw);
        const Scalar dkrn_dsw = fluidMatrixInteraction.dkrn_dsw(sw);

        const Scalar dmobilityW_dsw = dkrw_dsw/viscosityW_;
        const Scalar dmobilityN_dsw = dkrn_dsw/viscosityN_;

        // quotient rule
        const Scalar numerator = dmobilityW_dsw*(mobilityW+mobilityN) - mobilityW*(dmobilityW_dsw + dmobilityN_dsw);
        const Scalar denominator = (mobilityW+mobilityN)*(mobilityW+mobilityN);

        return numerator/denominator;
    }

    Scalar findSwShock_()
    {
        // Welge construction
        std::function<Scalar(Scalar)> tangentToSwr = [this](const Scalar sw)
        {
            Scalar diffQuot = (fractionalWFunc_(sw) - fractionalWFunc_(swr_)) / (sw - swr_);
            return diffQuot;
        };

        std::function<Scalar(Scalar)> residualFunction = [this, &tangentToSwr](const Scalar sw)
        {
            return dfractionalW_dsw_(sw) - tangentToSwr(sw);
        };

        auto swShock = findScalarRootBrent(swRight_+eps_, swLeft_-eps_, residualFunction);

        return swShock;
    }

    Scalar computeShockSpeed_()
    {
        Scalar characteristicSpeed = totalVelocity_/constPorosity_ * dfractionalW_dsw_(swShock_);
        // the shock-front speed selected by Welge’s tangent construction (\cite Welge1952) is the same speed
        // given by the Rankine–Hugoniot jump condition:
        // totalVelocity_/constPorosity_ * (fw(swShock_) - fw(swRight_))/(swShock_ - swRight_)
        return characteristicSpeed;
    }

    Scalar computeRarefactionSw_(Scalar xi)
    {
        auto residual = [this, &xi](Scalar sw)
        {
            return totalVelocity_/constPorosity_ * dfractionalW_dsw_(sw) - xi;
        };

        Scalar sw = findScalarRootBrent(swShock_, swLeft_, residual);

        return sw;
    }

    Scalar eps_ = 1e-12;
    std::shared_ptr<const Problem> problem_;
    ScalarVector values_;
    Scalar totalVelocity_;
    Scalar refPressure_;
    Scalar viscosityW_, viscosityN_;
    Scalar swr_, snr_;
    Scalar constPorosity_;

    Scalar swLeft_, swShock_,  swRight_;
    Scalar characteristicSpeed_;
};

} // end namespace Dumux

#endif
