// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup SurfaceRunoff
 * \brief Green-Ampt infiltration with a soil store of finite depth.
 */
#ifndef DUMUX_FREEFLOW_SHALLOWWATER_SURFACERUNOFF_GREENAMPT_HH
#define DUMUX_FREEFLOW_SHALLOWWATER_SURFACERUNOFF_GREENAMPT_HH

#include <algorithm>
#include <cmath>
#include <limits>

#include <dune/common/exceptions.hh>

namespace Dumux::SurfaceRunoff {

/*!
 * \ingroup SurfaceRunoff
 * \brief Green-Ampt infiltration into a soil column of finite depth.
 *
 * A sharp wetting front descends through soil of uniform initial moisture, driven by
 * gravity and by the suction `psi` ahead of it. The capacity
 * `f = Ks*(1 + psi*dTheta/F)` starts unbounded and decays towards `Ks` as the wetted
 * depth `F` grows.
 *
 * Classic Green-Ampt assumes a semi-infinite column, so `F` grows without bound and the
 * only way to produce runoff is for the rainfall intensity to exceed `Ks`. Daily rainfall
 * fields average intensity over 24 h and rarely do that, so the column is given a depth:
 * once `F` reaches `soilDepth*dTheta` the store is full and everything runs off. That is
 * saturation excess rather than infiltration excess, and on a wet catchment under a
 * multi-day storm it is the mechanism that actually generates the flood.
 */
template<class Scalar>
class GreenAmptSoil
{
public:
    GreenAmptSoil(const Scalar conductivity,
                  const Scalar suction,
                  const Scalar moistureDeficit,
                  const Scalar soilDepth)
    : conductivity_(conductivity)
    , suctionTimesDeficit_(suction*moistureDeficit)
    , capacity_(soilDepth*moistureDeficit)
    {
        if (conductivity < 0.0 || suction < 0.0 || moistureDeficit < 0.0 || soilDepth < 0.0)
            DUNE_THROW(Dune::InvalidStateException, "Green-Ampt parameters must be non-negative");
    }

    //! Wetted depth at which the store is full and no more water can enter
    Scalar capacity() const
    { return capacity_; }

    /*!
     * \brief Water entering the soil over one time step, as a depth.
     *
     * \param wettedDepth cumulative infiltration F at the start of the step
     * \param supplyRate water available at the surface, rainfall plus any ponding
     * \param dt the time step
     *
     * The return value never exceeds `supplyRate*dt` nor the remaining store, so it is
     * safe to subtract from the surface without driving the water depth negative.
     */
    Scalar infiltration(const Scalar wettedDepth, const Scalar supplyRate, const Scalar dt) const
    {
        using std::min;
        const auto remaining = capacity_ - wettedDepth;
        if (remaining <= 0.0 || supplyRate <= 0.0 || dt <= 0.0)
            return 0.0;

        return min(remaining, min(supplyRate*dt, unlimitedInfiltration_(wettedDepth, supplyRate, dt)));
    }

    /*!
     * \brief Time from dry until the surface ponds under constant rainfall.
     *
     * Before ponding the soil takes everything that falls, so the process is limited by
     * the rainfall and not by the capacity: `F = i*t`, and ponding begins once `F` reaches
     * `Fp`. Rainfall at or below `Ks` never ponds, the capacity never dropping that far.
     */
    Scalar timeToPonding(const Scalar rainfallRate) const
    {
        if (rainfallRate <= conductivity_)
            return std::numeric_limits<Scalar>::infinity();
        return wettedDepthAtPonding(rainfallRate)/rainfallRate;
    }

    //! Wetted depth at which the capacity has fallen to the rainfall rate
    Scalar wettedDepthAtPonding(const Scalar rainfallRate) const
    {
        if (rainfallRate <= conductivity_)
            return std::numeric_limits<Scalar>::infinity();
        return suctionTimesDeficit_*conductivity_/(rainfallRate - conductivity_);
    }

    //! Infiltration capacity at a given wetted depth
    Scalar capacityRate(const Scalar wettedDepth) const
    {
        if (wettedDepth <= 0.0)
            return std::numeric_limits<Scalar>::infinity();
        return conductivity_*(1.0 + suctionTimesDeficit_/wettedDepth);
    }

private:
    //! Green-Ampt without the finite store, split at the moment ponding begins
    Scalar unlimitedInfiltration_(const Scalar wettedDepth, const Scalar supplyRate, const Scalar dt) const
    {
        // a sealed surface ponds at once and passes nothing
        if (conductivity_ <= 0.0)
            return 0.0;

        // with no suction the front is gravity-driven and the capacity is Ks throughout
        if (suctionTimesDeficit_ <= 0.0)
            return std::min(supplyRate, conductivity_)*dt;

        if (supplyRate <= conductivity_)
            return supplyRate*dt;

        const auto pondingDepth = wettedDepthAtPonding(supplyRate);
        if (wettedDepth >= pondingDepth)
            return ponded_(wettedDepth, dt) - wettedDepth;

        // the surface only ponds part way through the step, if at all
        const auto timeToPond = (pondingDepth - wettedDepth)/supplyRate;
        if (timeToPond >= dt)
            return supplyRate*dt;

        return ponded_(pondingDepth, dt - timeToPond) - wettedDepth;
    }

    /*!
     * \brief Wetted depth after ponded infiltration for a time dt, starting from F0.
     *
     * Solves `F - F0 - psi*dTheta*ln((psi*dTheta + F)/(psi*dTheta + F0)) = Ks*dt`, which is
     * the cumulative form of Green-Ampt shifted to start at F0 rather than at zero.
     */
    Scalar ponded_(const Scalar startDepth, const Scalar dt) const
    {
        using std::abs, std::log;
        const auto& s = suctionTimesDeficit_;
        const auto target = conductivity_*dt;

        // the residual is increasing and convex in F, so Newton converges from either side
        auto depth = startDepth + target;
        for (int i = 0; i < maxIterations_; ++i)
        {
            const auto residual = depth - startDepth - s*log((s + depth)/(s + startDepth)) - target;
            const auto derivative = depth/(s + depth);
            const auto increment = residual/derivative;
            depth -= increment;
            if (abs(increment) < tolerance_*(s + depth))
                return depth;
        }

        DUNE_THROW(Dune::MathError, "Green-Ampt cumulative infiltration did not converge");
    }

    static constexpr int maxIterations_ = 100;
    static constexpr Scalar tolerance_ = 1e-12;

    Scalar conductivity_, suctionTimesDeficit_, capacity_;
};

} // end namespace Dumux::SurfaceRunoff

#endif
