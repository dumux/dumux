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

#include <cmath>
#include <functional>
#include <iostream>
#include <string>

#include <dumux/common/exceptions.hh>

namespace Dumux {

//! \brief Available Tvd approaches
enum class TvdApproach
{
    uniform = 1,
    li      = 2,
    hou     = 3
};

//! \biraf Available differencing schemes
enum class DifferencingScheme
{
    vanleer   = 1,
    vanalbada = 2,
    minmod    = 3,
    superbee  = 4,
    umist     = 5,
    mclimiter = 6,
    wahyd     = 7
};

/**
  * \brief This file contains different higher order methods for approximating the velocity.
  */
template<class Scalar>
class HigherOrderApproximation
{
public:
    HigherOrderApproximation(const std::string& paramGroup = "")
    {
        tvdApproach_ = static_cast<TvdApproach>(getParamFromGroup<int>(paramGroup, "Discretization.TvdApproach", 1));
        differencingScheme_ = static_cast<DifferencingScheme>(getParamFromGroup<int>(paramGroup, "Discretization.DifferencingScheme", 1));

        switch (differencingScheme_)
        {
            case DifferencingScheme::vanleer :
            {
                limiter_ = this->vanleer;
                break;
            }
            case DifferencingScheme::vanalbada :
            {
                if (tvdApproach_ == TvdApproach::hou)
                    DUNE_THROW(ParameterException, "\nDifferencing scheme " << static_cast<int>(DifferencingScheme::vanalbada) <<
                        " (Van Albada) is not implemented for the Tvd approach " << static_cast<int>(TvdApproach::hou) <<" (Hou).");
                limiter_ = this->vanalbada;
                break;
            }
            case DifferencingScheme::minmod :
            {
                limiter_ = this->minmod;
                break;
            }
            case DifferencingScheme::superbee :
            {
                limiter_ = this->superbee;
                break;
            }
            case DifferencingScheme::umist :
            {
                limiter_ = this->umist;
                break;
            }
            case DifferencingScheme::mclimiter :
            {
                limiter_ = this->mclimiter;
                break;
            }
            case DifferencingScheme::wahyd :
            {
                limiter_ = this->wahyd;
                break;
            }
            default:
            {
                DUNE_THROW(ParameterException, "\nDifferencing scheme " << static_cast<int>(differencingScheme_) <<
                    " is not implemented.\n");
                break;
            }
        }
    }

    /**
      * \brief Upwind Method
      */
    Scalar upwind(const Scalar downstreamVelocity,
                  const Scalar upstreamVelocity,
                  const Scalar density) const
    {
        return upstreamVelocity * density;
    }

    /**
      * \brief Central Differencing Method
      */
    Scalar centralDifference(const Scalar downstreamVelocity,
                             const Scalar upstreamVelocity,
                             const Scalar density) const
    {
        return (0.5 * (upstreamVelocity + downstreamVelocity)) * density;
    }

    /**
      * \brief Linear Upwind Method
      */
    Scalar linearUpwind(const Scalar upstreamVelocity,
                        const Scalar upUpstreamVelocity,
                        const Scalar upstreamToDownstreamDistance,
                        const Scalar upUpstreamToUpstreamDistance,
                        const Scalar density) const
    {
        Scalar zeroOrder = upstreamVelocity;
        Scalar firstOrder = -1.0 * ((upstreamVelocity - upUpstreamVelocity) / upUpstreamToUpstreamDistance) * ( upstreamToDownstreamDistance / -2.0);
        return (zeroOrder + firstOrder) * density;
    }

    /**
      * \brief QUICK upwinding Scheme: Quadratic Upstream Interpolation for Convective Kinematics
      */
    Scalar upwindQUICK(const Scalar downstreamVelocity,
                       const Scalar upstreamVelocity,
                       const Scalar upUpstreamVelocity,
                       const Scalar upstreamToDownstreamDistance,
                       const Scalar upUpstreamToUpstreamDistance,
                       const Scalar density) const
    {
        Scalar normalDistance = (upUpstreamToUpstreamDistance + upstreamToDownstreamDistance) / 2.0;
        Scalar zeroOrder = upstreamVelocity;
        Scalar firstOrder = ((downstreamVelocity - upstreamVelocity) / 2.0);
        Scalar secondOrder = -(((downstreamVelocity - upstreamVelocity) / upstreamToDownstreamDistance) - ((upstreamVelocity - upUpstreamVelocity) / upUpstreamToUpstreamDistance))
                           * upstreamToDownstreamDistance * upstreamToDownstreamDistance / (8.0 * normalDistance);
        return (zeroOrder + firstOrder + secondOrder) * density;
    }

    /**
      * \brief Tvd Scheme: Total Variation Diminuishing
      */
    Scalar tvd(const Scalar downstreamVelocity,
               const Scalar upstreamVelocity,
               const Scalar upUpstreamVelocity,
               const Scalar density) const
    {
        const Scalar ratio = (upstreamVelocity - upUpstreamVelocity) / (downstreamVelocity - upstreamVelocity);

        // If the velocity field is uniform (like at the first newton step) we get a NaN
        if(ratio > 0.0 && std::isfinite(ratio))
        {
            const Scalar secondOrderTerm = 0.5 * limiter_(ratio, 2.0) * (downstreamVelocity - upstreamVelocity);
            return density * (upstreamVelocity + secondOrderTerm);
        }
        else
            return density * upstreamVelocity;
    }

    /**
      * \brief Tvd Scheme: Total Variation Diminuishing
      *
      * This functions manages the non uniformities of the grid according to [Li, Liao 2007].
      * It tries to reconstruct the value for the velocity at the upstream-upstream point
      * if the grid was uniform.
      */
    Scalar tvd(const Scalar downstreamVelocity,
               const Scalar upstreamVelocity,
               const Scalar upUpstreamVelocity,
               const Scalar upstreamToDownstreamDistance,
               const Scalar upUpstreamToUpstreamDistance,
               const bool selfIsUpstream,
               const Scalar density) const
    {
        // I need the information of selfIsUpstream to get the correct sign because upUpstreamToUpstreamDistance is always positive
        const Scalar upUpstreamGradient = (upstreamVelocity - upUpstreamVelocity) / upUpstreamToUpstreamDistance * selfIsUpstream;

        // Distance between the upUpstream node and the position where it should be if the grid were uniform.
        const Scalar correctionDistance = upUpstreamToUpstreamDistance - upstreamToDownstreamDistance;
        const Scalar reconstrutedUpUpstreamVelocity = upUpstreamVelocity + upUpstreamGradient * correctionDistance;
        const Scalar ratio = (upstreamVelocity - reconstrutedUpUpstreamVelocity) / (downstreamVelocity - upstreamVelocity);

        // If the velocity field is uniform (like at the first newton step) we get a NaN
        if(ratio > 0 && std::isfinite(ratio))
        {
            const Scalar secondOrderTerm = 0.5 * limiter_(ratio, 2.0) * (downstreamVelocity - upstreamVelocity);
            return density * (upstreamVelocity + secondOrderTerm);
        }
        else
            return density * upstreamVelocity;
    }

    /**
     * \brief Tvd Scheme: Total Variation Diminuishing
     *
     * This functions manages the non uniformities of the grid according to [Hou, Simons, Hinkelmann 2007].
     * It should behave better then the Li's version in very stretched grids.
     */
    Scalar tvd(const Scalar downstreamVelocity,
               const Scalar upstreamVelocity,
               const Scalar upUpstreamVelocity,
               const Scalar upstreamToDownstreamDistance,
               const Scalar upUpstreamToUpstreamDistance,
               const Scalar downstreamStaggeredCellSize,
               const Scalar density) const
    {
        const Scalar ratio = (upstreamVelocity - upUpstreamVelocity) / (downstreamVelocity - upstreamVelocity)
                           * upstreamToDownstreamDistance / upUpstreamToUpstreamDistance;

        // If the velocity field is uniform (like at the first newton step) we get a NaN
        if(ratio > 0.0 && std::isfinite(ratio))
        {
            const Scalar upstreamStaggeredCellSize = 0.5 * (upstreamToDownstreamDistance + upUpstreamToUpstreamDistance);
            const Scalar R = (upstreamStaggeredCellSize + downstreamStaggeredCellSize) / upstreamStaggeredCellSize;
            const Scalar secondOrderTerm = limiter_(ratio, R) / R * (downstreamVelocity - upstreamVelocity);
            return density * (upstreamVelocity + secondOrderTerm);
        }
        else
            return density * upstreamVelocity;
    }

    /**
      * \brief Van Leer flux limiter function [Van Leer 1974]
      *
      * With R != 2 is the modified Van Leer flux limiter function [Hou, Simons, Hinkelmann 2007]
      */
    static Scalar vanleer(const Scalar r, const Scalar R)
    {
        return R * r / (R - 1.0 + r);
    }

    /**
      * \brief Van Albada flux limiter function [Van Albada et al. 1982]
      */
    static Scalar vanalbada(const Scalar r, const Scalar R)
    {
        return r * (r + 1.0) / (1.0 + r * r);
    }

    /**
      * \brief MinMod flux limiter function [Roe 1985]
      */
    static Scalar minmod(const Scalar r, const Scalar R)
    {
        return std::min(r, 1.0);
    }

    /**
      * \brief SUPERBEE flux limiter function [Roe 1985]
      *
      * With R != 2 is the modified SUPERBEE flux limiter function [Hou, Simons, Hinkelmann 2007]
      */
    static Scalar superbee(const Scalar r, const Scalar R)
    {
        return std::max(std::min(R * r, 1.0), std::min(r, R));
    }

    /**
      * \brief UMIST flux limiter function [Lien and Leschziner 1993]
      */
    static Scalar umist(const Scalar r, const Scalar R)
    {
        return std::min({R * r, (r * (5.0 - R) + R - 1.0) / 4.0, (r * (R - 1.0) + 5.0 - R) / 4.0, R});
    }

    /*
     * \brief Monotonized-Central limiter [Van Leer 1977]
     */
    static Scalar mclimiter(const Scalar r, const Scalar R)
    {
        return std::min({R * r, (r + 1.0) / 2.0, R});
    }

    /**
      * \brief WAHYD Scheme [Hou, Simons, Hinkelmann 2007];
      */
    static Scalar wahyd(const Scalar r, const Scalar R)
    {
        return r > 1 ? std::min((r + R * r * r) / (R + r * r), R)
                     : vanleer(r, R);
    }

    //! Returns the Tvd approach
    const TvdApproach& tvdApproach() const
    {
        return tvdApproach_;
    }

    //! Returns the differencing scheme
    const DifferencingScheme& differencingScheme() const
    {
        return differencingScheme_;
    }

private:
    TvdApproach tvdApproach_;
    DifferencingScheme differencingScheme_;

    std::function<Scalar(const Scalar, const Scalar)> limiter_;
};

} // end namespace Dumux

#endif // DUMUX_HIGHER_ORDER_VELOCITY_APPROXIMATION_HH
