// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup TwoPNCTests
 * \brief   Material law for a surfactant model.
 */
#ifndef DUMUX_TEST_2P3C_SURFACTANT_MATERIALLAW_HH
#define DUMUX_TEST_2P3C_SURFACTANT_MATERIALLAW_HH

#include <algorithm>
#include <dune/common/math.hh>

namespace Dumux::FluidMatrix {

/*!
 * \brief Material law for a surfactant model.This is based on the surfactant model described
 *        in K. Jørgensen's master thesis (NTNU, 2013) http://hdl.handle.net/11250/240038
 */
template <class Scalar>
class SurfactantPcKrSw
{
    // surfactant-related values
    struct SurfactantParams { double srw, sro, ww, wo; };

public:
    SurfactantPcKrSw(const Scalar Ncv, const Scalar sMinKr, const Scalar sMaxKr, const Scalar srwSurf, const Scalar sroSurf)
    : sMinKr_(sMinKr), sMaxKr_(sMaxKr), srwSurf_(srwSurf), sroSurf_(sroSurf)
    {
        using std::clamp; using std::log;
        params_.ww = 1.0 - clamp(log(Ncv / 1e-3) / log(8e-7 / 1e-3), 0.0, 1.0);
        params_.wo = 1.0 - clamp(log(Ncv / 1e-4) / log(2e-6 / 1e-4), 0.0, 1.0);

        const Scalar srwNoSurf = sMinKr;
        const Scalar sroNoSurf = 1.0 - sMaxKr;
        params_.srw = (1 - params_.ww) * srwNoSurf + params_.ww * srwSurf;
        params_.sro = (1 - params_.wo) * sroNoSurf + params_.wo * sroSurf;
    }

    Scalar pc(Scalar sw) const
    {
        return 0.0;
    }

    Scalar krw(Scalar sw) const
    {
        using std::clamp;

        const auto swStar = (sw - params_.srw) / (1.0 - params_.sro - params_.srw);
        const auto swMapped = sMinKr_ + (sMaxKr_ - sMinKr_) * swStar;
        const auto krwSurf = Dune::power(clamp(swStar, 0.0, 1.0), 2);

        return params_.ww * krwSurf + (1.0 - params_.ww) * krwNoSurfactant_(clamp(swMapped, sMinKr_, sMaxKr_));
    }

    Scalar krn(Scalar sw) const
    {
        using std::clamp;

        const auto swStar = (sw - params_.srw) / (1.0 - params_.sro - params_.srw);
        const auto swMapped = sMinKr_ + (sMaxKr_ - sMinKr_) * swStar;
        const auto kroSurf = Dune::power(clamp(1.0 - swStar, 0.0, 1.0), 2);

        return params_.wo * kroSurf + (1.0 - params_.wo) * krnNoSurfactant_(clamp(swMapped, sMinKr_, sMaxKr_));
    }

private:
    Scalar krnNoSurfactant_(Scalar sw) const
    {
        using std::clamp;

        const auto clampedS = std::clamp((sw - sMinKr_)/(sMaxKr_ - sMinKr_), 0.0, 1.0);
        return 0.14 * Dune::power(1.0 - clampedS, 2);
    }

    Scalar krwNoSurfactant_ (Scalar sw) const
    {
        using std::clamp;

        return 0.02 * Dune::power(clamp((sw - sMinKr_) / (sMaxKr_ - sMinKr_), 0.0, 1.0), 2);
    }

    SurfactantParams params_;
    Scalar sMinKr_, sMaxKr_;
    Scalar srwSurf_, sroSurf_;
};

} // end namespace Dumux::FluidMatrix

#endif
