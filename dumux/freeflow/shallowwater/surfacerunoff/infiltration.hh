// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup SurfaceRunoff
 * \brief Rainfall losses on the surface: canopy interception and Green-Ampt infiltration.
 */
#ifndef DUMUX_FREEFLOW_SHALLOWWATER_SURFACERUNOFF_INFILTRATION_HH
#define DUMUX_FREEFLOW_SHALLOWWATER_SURFACERUNOFF_INFILTRATION_HH

#include <algorithm>
#include <cstddef>
#include <utility>
#include <vector>

#include <dune/common/exceptions.hh>

#include "greenampt.hh"

namespace Dumux::SurfaceRunoff {

/*!
 * \ingroup SurfaceRunoff
 * \brief The rainfall a surface loses before it can run off, one store per degree of freedom.
 *
 * Canopy interception and infiltration act in series: the canopy store fills first, and the
 * throughfall together with part of the ponded water is offered to the Green-Ampt soil column.
 *
 * The loss rates are computed before a time step and held constant over it. The stores advance
 * in `commit()` with the time step size that was actually taken, which may be smaller than the
 * one the rates were computed for, so that the stores gain exactly what the surface lost.
 */
template<class Scalar>
class SurfaceLosses
{
    using Soil = GreenAmptSoil<Scalar>;

public:
    /*!
     * \brief One uniform soil and canopy everywhere.
     * \param numDofs number of surface degrees of freedom
     * \param soil the Green-Ampt soil column
     * \param interceptionDepth canopy storage capacity, as a depth
     */
    SurfaceLosses(const std::size_t numDofs, const Soil& soil, const Scalar interceptionDepth = 0.0)
    : soil_(numDofs, soil)
    , interceptionCapacity_(numDofs, interceptionDepth)
    { resizeStores_(numDofs); }

    /*!
     * \brief A soil and a canopy store per degree of freedom.
     */
    SurfaceLosses(std::vector<Soil>&& soil, std::vector<Scalar>&& interceptionDepth)
    : soil_(std::move(soil))
    , interceptionCapacity_(std::move(interceptionDepth))
    {
        if (soil_.size() != interceptionCapacity_.size())
            DUNE_THROW(Dune::InvalidStateException,
                       "Got " << soil_.size() << " soil columns but "
                       << interceptionCapacity_.size() << " interception depths");
        resizeStores_(soil_.size());
    }

    //! Water already held in the soil column at the start, as a depth
    void setInitialWettedDepth(const Scalar depth)
    { std::fill(wettedDepth_.begin(), wettedDepth_.end(), depth); }

    /*!
     * \brief Water already held in each soil column at the start, as a depth, capped at the
     *        capacity of the column
     */
    void setInitialWettedDepth(const std::vector<Scalar>& depths)
    {
        if (depths.size() != wettedDepth_.size())
            DUNE_THROW(Dune::InvalidStateException,
                       "Got " << depths.size() << " initial wetted depths but "
                       << wettedDepth_.size() << " soil columns");
        for (std::size_t i = 0; i < depths.size(); ++i)
        {
            if (depths[i] < 0.0)
                DUNE_THROW(Dune::InvalidStateException,
                           "Initial wetted depth at column " << i << " is negative");
            wettedDepth_[i] = std::min(depths[i], soil_[i].capacity());
        }
    }

    /*!
     * \brief Fix the loss rates at one degree of freedom for the step about to be taken.
     *
     * \param dofIdx the degree of freedom
     * \param rainfallRate rain falling on the canopy
     * \param pondedDepth water standing on the surface
     * \param dt the time step
     * \param pondDrawdown fraction of the pond offered to the soil within one step
     *
     * Only a fraction of the pond is offered, since the surface flux over the step is
     * computed from the same depth.
     */
    void update(const std::size_t dofIdx,
                const Scalar rainfallRate,
                const Scalar pondedDepth,
                const Scalar dt,
                const Scalar pondDrawdown)
    {
        if (dt <= 0.0)
            return;

        using std::max, std::min;
        const auto canopyRoom = max(Scalar(0.0), interceptionCapacity_[dofIdx] - interceptionStored_[dofIdx]);
        interceptionRate_[dofIdx] = min(max(Scalar(0.0), rainfallRate), canopyRoom/dt);

        const auto throughfall = max(Scalar(0.0), rainfallRate) - interceptionRate_[dofIdx];
        const auto supply = throughfall + pondDrawdown*max(Scalar(0.0), pondedDepth)/dt;
        infiltrationRate_[dofIdx] = soil_[dofIdx].infiltration(wettedDepth_[dofIdx], supply, dt)/dt;
    }

    //! Total rate the surface loses, canopy plus soil
    Scalar lossRate(const std::size_t dofIdx) const
    { return interceptionRate_[dofIdx] + infiltrationRate_[dofIdx]; }

    Scalar infiltrationRate(const std::size_t dofIdx) const
    { return infiltrationRate_[dofIdx]; }

    Scalar interceptionRate(const std::size_t dofIdx) const
    { return interceptionRate_[dofIdx]; }

    //! Move both stores by what the step actually took
    void commit(const Scalar dt)
    {
        for (std::size_t i = 0; i < wettedDepth_.size(); ++i)
        {
            wettedDepth_[i] += infiltrationRate_[i]*dt;
            interceptionStored_[i] += interceptionRate_[i]*dt;
        }
    }

    //! Water held in the soil column, as a depth
    const std::vector<Scalar>& wettedDepth() const
    { return wettedDepth_; }

    //! Water held on the canopy, as a depth
    const std::vector<Scalar>& interceptionStored() const
    { return interceptionStored_; }

    std::size_t size() const
    { return soil_.size(); }

private:
    void resizeStores_(const std::size_t numDofs)
    {
        wettedDepth_.assign(numDofs, 0.0);
        interceptionStored_.assign(numDofs, 0.0);
        infiltrationRate_.assign(numDofs, 0.0);
        interceptionRate_.assign(numDofs, 0.0);
    }

    std::vector<Soil> soil_;
    std::vector<Scalar> interceptionCapacity_;
    std::vector<Scalar> wettedDepth_, interceptionStored_;
    std::vector<Scalar> infiltrationRate_, interceptionRate_;
};

} // end namespace Dumux::SurfaceRunoff

#endif
