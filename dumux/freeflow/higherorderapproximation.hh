// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   See the file COPYING for full copying permissions.                      *
 *                                                                           *
 *   This program is free software: you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation, either version 2 of the License, or       *
 *   (at your option) any later version.                                     *
 *                                                                           *
 *   This program is distributed in the hope that it will be useful,         *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of          *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the            *
 *   GNU General Public License for more details.                            *
 *                                                                           *
 *   You should have received a copy of the GNU General Public License       *
 *   along with this program.  If not, see <http://www.gnu.org/licenses/>.   *
 *****************************************************************************/
/** \file
  * \brief This file contains different higher order methods for approximating the velocity.
  */

#ifndef DUMUX_HIGHER_ORDER_VELOCITY_APPROXIMATION_HH
#define DUMUX_HIGHER_ORDER_VELOCITY_APPROXIMATION_HH

#include <iostream>

#include<dune/common/fvector.hh>

namespace Dumux {

/**
  * \brief This file contains different higher order methods for approximating the velocity.
  */
template<class Scalar>
class HigherOrderApproximation
{
public:
    /**
      * \brief Upwind Method
      */
    Scalar upwind(const Scalar downstreamVelocity,
                            const Scalar upstreamVelocity
                            const Scalar density) const
    {
        return upstreamVelocity * density;
    }

    /**
      * \brief Central Differencing Method
      */
    Scalar centralDifference(const Scalar downstreamVelocity,
                            const Scalar upstreamVelocity
                            const Scalar density) const
    {
        return (0.5 * (upstreamVelocity + downstreamVelocity)) * density;
    }


    /**
      * \brief Linear Upwind Method
      */
    Scalar linearUpwind(const Scalar upstreamVelocity,
                             const Scalar previousVelocity,
                             const Scalar upstreamToDownstreamDistance,
                             const Scalar previousToUpstreamDistance,
                             const Scalar density) const
    {
        Scalar zeroOrder = upstreamVelocity;
        Scalar firstOrder = -1.0 * ((upstreamVelocity - previousVelocity) / previousToUpstreamDistance) * ( upstreamToDownstreamDistance / -2.0);
        return (zeroOrder + firstOrder) * density;
    }

    /**
      * \brief QUICK upwinding Scheme: Quadratic Upstream Interpolation for Convective Kinematics
      */
    Scalar upwindQUICK(const Scalar downstreamVelocity,
                       const Scalar upstreamVelocity,
                       const Scalar previousVelocity,
                       const Scalar upstreamToDownstreamDistance,
                       const Scalar previousToUpstreamDistance,
                       const Scalar density) const
    {
        Scalar normalDistance = (previousToUpstreamDistance + upstreamToDownstreamDistance) / 2.0;
        Scalar zeroOrder = upstreamVelocity;
        Scalar firstOrder = ((downstreamVelocity - upstreamVelocity) / 2.0);
        Scalar secondOrder = -(((downstreamVelocity - upstreamVelocity) / upstreamToDownstreamDistance) - ((upstreamVelocity - previousVelocity) / previousToUpstreamDistance))
                  * std::pow(upstreamToDownstreamDistance , 2 ) / (8.0 * normalDistance);
        return (zeroOrder + firstOrder + secondOrder) * density;
    }

};
} // end namespace Dumux

#endif // DUMUX_TURBULENCE_PROPERTIES_HH
