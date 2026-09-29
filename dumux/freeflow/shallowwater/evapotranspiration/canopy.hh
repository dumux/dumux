// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup Evapotranspiration
 * \brief Splitting a reference evapotranspiration between a canopy and the ground below it.
 */
#ifndef DUMUX_FREEFLOW_SHALLOWWATER_EVAPOTRANSPIRATION_CANOPY_HH
#define DUMUX_FREEFLOW_SHALLOWWATER_EVAPOTRANSPIRATION_CANOPY_HH

#include <algorithm>
#include <cmath>

#include <dune/common/exceptions.hh>

namespace Dumux::Evapotranspiration {

/*!
 * \ingroup Evapotranspiration
 * \brief A canopy that partitions a reference evapotranspiration between transpiration and
 *        the ground by Beer's law on the leaf area index
 *
 * For a deciduous canopy, the leaf area index ramps linearly over `transitionDays` after
 * leaf-on and before leaf-off.
 */
template<class Scalar>
class DeciduousCanopy
{
public:
    /*!
     * \param maxLeafAreaIndex leaf area index of the closed canopy [-]
     * \param extinctionCoefficient Beer's law extinction for the canopy [-]
     * \param leafOnDay day of year at which leaf-out begins
     * \param leafOffDay day of year at which leaf fall is complete
     * \param transitionDays length of the leaf-out and leaf-fall ramps [d]
     * \param woodyAreaIndex area index of stems and branches [-]
     */
    DeciduousCanopy(Scalar maxLeafAreaIndex,
                    Scalar extinctionCoefficient,
                    Scalar leafOnDay,
                    Scalar leafOffDay,
                    Scalar transitionDays,
                    Scalar woodyAreaIndex = 0.0)
    : maxLeafAreaIndex_(maxLeafAreaIndex)
    , extinctionCoefficient_(extinctionCoefficient)
    , leafOn_(leafOnDay)
    , leafOff_(leafOffDay)
    , transition_(transitionDays)
    , woodyAreaIndex_(woodyAreaIndex)
    , deciduous_(true)
    {
        if (!(leafOff_ > leafOn_))
            DUNE_THROW(Dune::RangeError, "a deciduous canopy needs leaf-off (" << leafOff_
                       << ") after leaf-on (" << leafOn_ << "); a season that wraps the new "
                       "year is not supported");
        if (!(leafOff_ <= 365.25))
            DUNE_THROW(Dune::RangeError, "leaf-off on day " << leafOff_ << " is past the end of "
                       "the year, so leaf fall would be cut short by the wrap");
        if (!(transition_ > 0.0))
            DUNE_THROW(Dune::RangeError, "the leaf-out ramp must be longer than zero days");
    }

    //! An evergreen canopy, or one whose phenology is not known.
    DeciduousCanopy(Scalar maxLeafAreaIndex, Scalar extinctionCoefficient,
                    Scalar woodyAreaIndex = 0.0)
    : maxLeafAreaIndex_(maxLeafAreaIndex)
    , extinctionCoefficient_(extinctionCoefficient)
    , leafOn_(0.0)
    , leafOff_(0.0)
    , transition_(1.0)
    , woodyAreaIndex_(woodyAreaIndex)
    , deciduous_(false)
    {}

    //! Leaf area index on a given day of year.
    Scalar leafAreaIndex(Scalar dayOfYear) const
    {
        if (!deciduous_)
            return maxLeafAreaIndex_;

        using std::clamp;
        using std::min;
        const auto out = clamp((dayOfYear - leafOn_)/transition_, Scalar(0.0), Scalar(1.0));
        const auto fall = clamp((leafOff_ - dayOfYear)/transition_, Scalar(0.0), Scalar(1.0));
        return maxLeafAreaIndex_*min(out, fall);
    }

    //! The share of the reference rate intercepted by the leaves and available for
    //! transpiration, \f$1 - e^{-k\,\mathrm{LAI}}\f$
    Scalar transpiringFraction(Scalar dayOfYear) const
    {
        using std::exp;
        return 1.0 - exp(-extinctionCoefficient_*leafAreaIndex(dayOfYear));
    }

    /*!
     * \brief The share that reaches the ground, \f$e^{-k(\mathrm{LAI}+\mathrm{WAI})}\f$.
     *
     * Stems and branches keep intercepting after leaf fall. The woody share is neither
     * transpired nor passed to the ground, so for a positive woody area index this and
     * `transpiringFraction` sum to less than one.
     */
    Scalar groundFraction(Scalar dayOfYear) const
    {
        using std::exp;
        return exp(-extinctionCoefficient_*(leafAreaIndex(dayOfYear) + woodyAreaIndex_));
    }

private:
    Scalar maxLeafAreaIndex_;
    Scalar extinctionCoefficient_;
    Scalar leafOn_;
    Scalar leafOff_;
    Scalar transition_;
    Scalar woodyAreaIndex_;
    bool deciduous_;
};

//! Day of year counted from zero at `time` (in seconds) of a run that began on day
//! `startDayOfYear` counted from one
template<class Scalar>
Scalar dayOfYear(Scalar time, Scalar startDayOfYear)
{
    using std::fmod;
    return fmod(startDayOfYear - 1.0 + time/86400.0, 365.25);
}

} // end namespace Dumux::Evapotranspiration

#endif
