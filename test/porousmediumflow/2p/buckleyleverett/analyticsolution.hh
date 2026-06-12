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

#include <dune/common/fvector.hh>
#include <dune/istl/bvector.hh>

#include <dumux/common/parameters.hh>

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
    {
        initialize_();
        update(0.0);
    }

    void update(Scalar time)
    {
        const auto& gridView = problem_->gridGeometry().gridView();

        // position of characteristics x = v_t/phi * df_w/dS_w * t.
        Dune::FieldVector<Scalar, pointNum_> frontPosition(0.0);
        for (int i = 0; i < pointNum_; ++i)
            frontPosition[i] = totalVelocity_*time/problem_->spatialParams().constantPorosity()*dfwdsw_[i];

        // initial guess for shock index, gets updated in the following while loop
        int shockIdx = pointNum_/3;
        int shockIdxOld = 0;
        int shockIdxOldOld = 0;
        int frontEndIdx = 0;
        const int frontMaxIdx = firstTurningPoint_(frontPosition);

        bool converged = false;
        while (!converged)
        {
            // avoid multiple solutions by balancing mass ahead and behind the shock front
            Scalar areaAhead = 0.0;
            for (int i = 0; i <= shockIdx - 1; ++i)
                areaAhead += (satVec_[i] - swr_ + satVec_[i + 1] - swr_)*0.5*(frontPosition[i + 1] - frontPosition[i]);

            Scalar areaBehind = 0.0;
            for (int i = shockIdx; i <= frontMaxIdx - 1; ++i)
                areaBehind += (satVec_[frontMaxIdx] - satVec_[i] + satVec_[frontMaxIdx] - satVec_[i + 1])*0.5*(frontPosition[i + 1] - frontPosition[i]);

            auto x = frontPosition[frontMaxIdx];
            frontEndIdx = frontMaxIdx;
            while (frontEndIdx + 1 < pointNum_ && x > frontPosition[shockIdx])
            {
                ++frontEndIdx;
                x = frontPosition[frontEndIdx];
            }

            for (int i = frontMaxIdx; i <= frontEndIdx - 1; ++i)
                areaBehind += (satVec_[i] - satVec_[frontMaxIdx] + satVec_[i + 1] - satVec_[frontMaxIdx])*0.5*(frontPosition[i] - frontPosition[i + 1]);

            shockIdxOldOld = shockIdxOld;
            shockIdxOld = shockIdx;

            if (std::abs(areaAhead) > std::abs(areaBehind))
                --shockIdx;
            else if (std::abs(areaAhead) < std::abs(areaBehind))
                ++shockIdx;

            // converged when index oscillates between two neighboring values
            converged = shockIdx == shockIdxOldOld;
        }

        for (const auto& element : elements(gridView))
        {
            const auto globalPos = element.geometry().center();
            const auto eIdx = problem_->gridGeometry().elementMapper().index(element);

            if (globalPos[0] > frontPosition[frontEndIdx])
            {
                values_[eIdx] = swr_;
                continue;
            }

            int nextIdx = 0;
            for (int i = intervalNum_; i >= 0; --i)
            {
                if (globalPos[0] < frontPosition[i])
                {
                    nextIdx = i;
                    break;
                }
            }

            values_[eIdx] = satVec_[nextIdx];
        }
    }

    const ScalarVector& values() const
    { return values_; }

private:
    static constexpr int intervalNum_ = 1000;
    static constexpr int pointNum_ = intervalNum_ + 1;

    void initialize_()
    {
        GlobalPosition pos(0.0);
        const auto fluidMatrixInteraction = problem_->spatialParams().fluidMatrixInteractionAtPos(pos);
        swr_ = fluidMatrixInteraction.pcSwCurve().effToAbsParams().swr();
        snr_ = fluidMatrixInteraction.pcSwCurve().effToAbsParams().snr();

        satVec_ = swr_;
        for (int i = 1; i < pointNum_; ++i)
            satVec_[i] = satVec_[i - 1] + (1.0 - snr_ - swr_)/intervalNum_;

        FluidState fluidState;
        Scalar referencePressure = getParam<Scalar>("Problem.ReferencePressure");
        fluidState.setTemperature(problem_->spatialParams().temperatureAtPos(pos));
        fluidState.setPressure(FluidSystem::phase0Idx, referencePressure);
        fluidState.setPressure(FluidSystem::phase1Idx, referencePressure);

        const auto viscosityW = FluidSystem::viscosity(fluidState, FluidSystem::phase0Idx);
        const auto viscosityN = FluidSystem::viscosity(fluidState, FluidSystem::phase1Idx);

        std::cout << "ViscosityN: " << viscosityN << std::endl;

        fractionalW_ = 0.0;
        for (int i = 0; i < pointNum_; ++i)
        {
            const auto mobilityW = fluidMatrixInteraction.krw(satVec_[i])/viscosityW;
            const auto mobilityN = fluidMatrixInteraction.krn(satVec_[i])/viscosityN;
            fractionalW_[i] = mobilityW/(mobilityW + mobilityN);
        }

        dfwdsw_ = 0.0;
        for (int i = 1; i < intervalNum_; ++i)
            dfwdsw_[i] = (fractionalW_[i + 1] - fractionalW_[i - 1])/(satVec_[i + 1] - satVec_[i - 1]);
    }

    int firstTurningPoint_(const Dune::FieldVector<Scalar, pointNum_>& values) const
    {
        for (int i = 0; i < intervalNum_; ++i)
            if (values[i] > values[i + 1])
                return i;

        return intervalNum_;
    }

    std::shared_ptr<const Problem> problem_;
    ScalarVector values_;
    Scalar totalVelocity_;
    Scalar swr_;
    Scalar snr_;

    Dune::FieldVector<Scalar, pointNum_> satVec_;
    Dune::FieldVector<Scalar, pointNum_> fractionalW_;
    Dune::FieldVector<Scalar, pointNum_> dfwdsw_;
};

} // end namespace Dumux

#endif
