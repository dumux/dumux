// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup ShallowWaterTests
 * \brief Bed and roughness of Wooding's V-catchment with the channel resolved in 2D.
 */
#ifndef DUMUX_TEST_FREEFLOW_SHALLOWWATER_LONGWAVE_WOODING_CATCHMENT_HH
#define DUMUX_TEST_FREEFLOW_SHALLOWWATER_LONGWAVE_WOODING_CATCHMENT_HH

#include <cmath>

#include <dumux/common/parameters.hh>

namespace Dumux::Wooding {

/*!
 * \ingroup ShallowWaterTests
 * \brief Geometry and roughness shared by the two models solved on this catchment.
 *
 * The channel occupies a strip of Problem.ChannelWidth centred on Problem.ChannelX and
 * falls along y at Problem.ChannelSlope towards the outlet at y = 0.
 *
 * The planes fall along y at Problem.PlaneAlongSlope, which defaults to the channel slope
 * but is the analytic solution's assumption only when it is zero: a plane that descends
 * towards the outlet is a second, far wider and smoother conveyance path down the valley,
 * and then the channel routes only a fraction of the discharge. Levelling the planes costs
 * the constant incision depth, since a channel that keeps falling under a level plane has
 * to cut deeper the further downstream it goes, so Problem.RiverDepth becomes the incision
 * at the channel head rather than everywhere.
 */
template<class Scalar>
class Catchment
{
public:
    Catchment()
    : channelX_(getParam<Scalar>("Problem.ChannelX"))
    , halfWidth_(0.5*getParam<Scalar>("Problem.ChannelWidth"))
    , channelLength_(getParam<Scalar>("Problem.ChannelLength"))
    , crossSlope_(getParam<Scalar>("Problem.CrossSlope"))
    , channelSlope_(getParam<Scalar>("Problem.ChannelSlope"))
    , planeAlongSlope_(getParam<Scalar>("Problem.PlaneAlongSlope", channelSlope_))
    , riverDepth_(getParam<Scalar>("Problem.RiverDepth"))
    , planeManningN_(getParam<Scalar>("Problem.PlaneManningN"))
    , channelManningN_(getParam<Scalar>("Problem.ChannelManningN"))
    {}

    //! distance from the channel edge, zero inside the channel
    template<class GlobalPosition>
    Scalar distanceToChannel(const GlobalPosition& globalPos) const
    {
        using std::abs, std::max;
        return max(Scalar(0.0), abs(globalPos[0] - channelX_) - halfWidth_);
    }

    template<class GlobalPosition>
    bool inChannel(const GlobalPosition& globalPos) const
    { return distanceToChannel(globalPos) <= 0.0; }

    template<class GlobalPosition>
    Scalar bedElevation(const GlobalPosition& globalPos) const
    {
        const auto distance = distanceToChannel(globalPos);
        const auto plane = planeAlongSlope_*globalPos[1];
        if (distance > 0.0)
            return plane + crossSlope_*distance;
        return plane - riverDepth_
               - (channelSlope_ - planeAlongSlope_)*(channelLength_ - globalPos[1]);
    }

    template<class GlobalPosition>
    Scalar manningN(const GlobalPosition& globalPos) const
    { return inChannel(globalPos) ? channelManningN_ : planeManningN_; }

private:
    Scalar channelX_, halfWidth_, channelLength_;
    Scalar crossSlope_, channelSlope_, planeAlongSlope_, riverDepth_;
    Scalar planeManningN_, channelManningN_;
};

} // end namespace Dumux::Wooding

#endif
