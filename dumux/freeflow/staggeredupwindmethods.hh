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
/*!
 * \file
 * \ingroup FreeflowModels
 * \brief This file contains different higher order methods for approximating the velocity.
 */

#ifndef DUMUX_UPWINDING_METHODS_HH
#define DUMUX_UPWINDING_METHODS_HH

#include <cmath>
#include <functional>
#include <iostream>
#include <string>

#include <dumux/common/exceptions.hh>
#include <dumux/common/parameters.hh>

namespace Dumux {

/*!
 * \ingroup FreeflowModels
 * \brief Available Tvd approaches
 */
enum class TvdApproach
{
    none, uniform, li, hou
};

/*!
 * \ingroup FreeflowModels
 * \brief Available differencing schemes
 */
enum class DifferencingScheme
{
    none, vanleer, vanalbada, minmod, superbee, umist, mclimiter, wahyd
};

/*!
 * \ingroup FreeflowModels
 * \brief This file contains different higher order methods for approximating the velocity.
 */
template<class Scalar, int upwindSchemeOrder>
class StaggeredUpwindMethods
{
public:
    StaggeredUpwindMethods(const std::string& paramGroup = "")
    {
        upwindWeight_ = getParamFromGroup<Scalar>(paramGroup, "Flux.UpwindWeight");

        if (upwindSchemeOrder > 1)
        {
            // Read the runtime parameters
            tvdApproach_ = tvdApproachFromString(getParamFromGroup<std::string>(paramGroup, "Flux.TvdApproach", "Uniform"));
            differencingScheme_ = differencingSchemeFromString(getParamFromGroup<std::string>(paramGroup, "Flux.DifferencingScheme", "Minmod"));

            // Assign the limiter_ depending on the differencing scheme
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
                        DUNE_THROW(ParameterException, "\nDifferencing scheme (Van Albada) is not implemented for the Tvd approach (Hou).");
                    else
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
                    DUNE_THROW(ParameterException, "\nDifferencing scheme " << static_cast<std::string>(differencingSchemeToString(differencingScheme_)) <<
                        " is not implemented.\n");  // should never be reached
                    break;
                }
            }

            if (!hasParamInGroup(paramGroup, "Flux.TvdApproach"))
            {
                std::cout << "No TvdApproach specified. Defaulting to the Uniform method." << "\n";
                std::cout << "Other available TVD approaches for uniform (and nonuniform) grids are as follows: \n"
                          << "  " << tvdApproachToString(TvdApproach::uniform) << ": assumes a Uniform cell size distribution\n"
                          << "  " << tvdApproachToString(TvdApproach::li) << ": Li's approach for nonuniform cell sizes\n"
                          << "  " << tvdApproachToString(TvdApproach::hou) << ": Hou's approach for nonuniform cell sizes \n";
                std::cout << "Each approach can be specified as written above in the Flux group under the title TvdApproach in your input file. \n";
            }
            if (!hasParamInGroup(paramGroup, "Flux.DifferencingScheme"))
            {

                std::cout << "No DifferencingScheme specified. Defaulting to the Minmod scheme." << "\n";
                std::cout << "Other available Differencing Schemes are as follows: \n"
                          << "  " << differencingSchemeToString(DifferencingScheme::vanleer) << ": The Vanleer flux limiter\n"
                          << "  " << differencingSchemeToString(DifferencingScheme::vanalbada) << ": The Vanalbada flux limiter\n"
                          << "  " << differencingSchemeToString(DifferencingScheme::minmod) << ": The Minmod flux limiter\n"
                          << "  " << differencingSchemeToString(DifferencingScheme::superbee) << ": The Superbee flux limiter\n"
                          << "  " << differencingSchemeToString(DifferencingScheme::umist) << ": The Umist flux limiter\n"
                          << "  " << differencingSchemeToString(DifferencingScheme::mclimiter) << ": The Mclimiter flux limiter\n"
                          << "  " << differencingSchemeToString(DifferencingScheme::wahyd) << ": The Wahyd flux limiter";
                std::cout << "Each scheme can be specified as written above in the Flux group under the variable DifferencingScheme in your input file. \n";
            }

            std::cout << "Using the tvdApproach \"" << tvdApproachToString(tvdApproach_)
                      << "\" and the differencing Scheme \" " << differencingSchemeToString(differencingScheme_) << "\" \n";
        }
        else
        {
            // If the runtime parameters are not specified we will use upwind
            tvdApproach_ = TvdApproach::none;
            differencingScheme_ = DifferencingScheme::none;
        }
    }

    /**
     * \brief Convenience function to convert user input given as std::string
     *        to the corresponding enum class used for choosing the TVD Approach
     */
    TvdApproach tvdApproachFromString(const std::string& tvd)
    {
        if (tvd == "Uniform") return TvdApproach::uniform;
        if (tvd == "Li") return TvdApproach::li;
        if (tvd == "Hou") return TvdApproach::hou;
        DUNE_THROW(ParameterException, "\nThis tvd approach : \"" << tvd << "\" is not implemented.\n"
                                       << "The available TVD approaches for uniform (and nonuniform) grids are as follows: \n"
                                       << tvdApproachToString(TvdApproach::uniform) << ": assumes a uniform cell size distribution\n"
                                       << tvdApproachToString(TvdApproach::li) << ": li's approach for nonuniform cell sizes\n"
                                       << tvdApproachToString(TvdApproach::hou) << ": hou's approach for nonuniform cell sizes");
    }

    /**
     * \brief return the name of the TVD approach
     */
    std::string tvdApproachToString(TvdApproach tvd)
    {
        switch (tvd)
        {
            case TvdApproach::uniform: return "Uniform";
            case TvdApproach::li: return "Li";
            case TvdApproach::hou: return "Hou";
            default: return "Invalid"; // should never be reached
        }
    }

    /**
     * \brief Convenience function to convert user input given as std::string
     *        to the corresponding enum class used for choosing the Discretization Method
     */
    DifferencingScheme differencingSchemeFromString(const std::string& differencingScheme)
    {
        if (differencingScheme == "Vanleer") return DifferencingScheme::vanleer;
        if (differencingScheme == "Vanalbada") return DifferencingScheme::vanalbada;
        if (differencingScheme == "Minmod") return DifferencingScheme::minmod;
        if (differencingScheme == "Superbee") return DifferencingScheme::superbee;
        if (differencingScheme == "Umist") return DifferencingScheme::umist;
        if (differencingScheme == "Mclimiter") return DifferencingScheme::mclimiter;
        if (differencingScheme == "Wahyd") return DifferencingScheme::wahyd;
        DUNE_THROW(ParameterException, "\nThis differencing scheme: \"" << differencingScheme << "\" is not implemented.\n"
                                       << "The available differencing schemes are as follows: \n"
                                       << differencingSchemeToString(DifferencingScheme::vanleer) << ": The Vanleer flux limiter\n"
                                       << differencingSchemeToString(DifferencingScheme::vanalbada) << ": The VanAlbada flux limiter\n"
                                       << differencingSchemeToString(DifferencingScheme::minmod) << ": The Min-Mod flux limiter\n"
                                       << differencingSchemeToString(DifferencingScheme::superbee) << ": The SuperBEE flux limiter\n"
                                       << differencingSchemeToString(DifferencingScheme::umist) << ": The UMist flux limiter\n"
                                       << differencingSchemeToString(DifferencingScheme::mclimiter) << ": The McLimiter flux limiter\n"
                                       << differencingSchemeToString(DifferencingScheme::wahyd) << ": The Wahyd flux limiter");
    }

    /**
     * \brief return the name of the Discretization Method
     */
    std::string differencingSchemeToString(DifferencingScheme differencingScheme)
    {
        switch (differencingScheme)
        {
            case DifferencingScheme::vanleer: return "Vanleer";
            case DifferencingScheme::vanalbada: return "Vanalbada";
            case DifferencingScheme::minmod: return "Minmod";
            case DifferencingScheme::superbee: return "Superbee";
            case DifferencingScheme::umist: return "Umist";
            case DifferencingScheme::mclimiter: return "Mclimiter";
            case DifferencingScheme::wahyd: return "Wahyd";
            default: return "Invalid"; // should never be reached
        }
    }

    /**
      * \brief Upwind Method
      */
    Scalar upwind(const Scalar downstreamMomentum,
                  const Scalar upstreamMomentum) const
    {
        return (upwindWeight_ * upstreamMomentum + (1.0 - upwindWeight_) * downstreamMomentum);
    }

    /**
     * \brief Tvd Scheme: Total Variation Diminishing
     *
     */
    Scalar tvd(const std::array<Scalar,3>& momenta,
               const std::array<Scalar,3>& distances,
               const bool selfIsUpstream,
               const TvdApproach tvdApproach) const
    {
        Scalar momentum = 0.0;
        switch(tvdApproach)
        {
            case TvdApproach::uniform :
            {
                momentum += tvdUniform(momenta, distances, selfIsUpstream);
                break;
            }
            case TvdApproach::li :
            {
                momentum += tvdLi(momenta, distances, selfIsUpstream);
                break;
            }
            case TvdApproach::hou :
            {
                momentum += tvdHou(momenta, distances, selfIsUpstream);
                break;
            }
            default:
            {
                DUNE_THROW(ParameterException, "\nThis Tvd Approach is not implemented.\n");
                break;
            }
        }
        return momentum;
    }

    /**
     * \brief Tvd Scheme: Total Variation Diminishing
     *
     * This function assumes the cell size distribution to be uniform.
     */
    Scalar tvdUniform(const std::array<Scalar,3>& momenta,
                      const std::array<Scalar,3>& distances,
                      const bool selfIsUpstream) const
    {
        using std::isfinite;
        const Scalar downstreamMomentum = momenta[0];
        const Scalar upstreamMomentum = momenta[1];
        const Scalar upUpstreamMomentum = momenta[2];
        const Scalar ratio = (upstreamMomentum - upUpstreamMomentum) / (downstreamMomentum - upstreamMomentum);

        // If the velocity field is uniform (like at the first newton step) we get a NaN
        if(ratio > 0.0 && isfinite(ratio))
        {
            const Scalar secondOrderTerm = 0.5 * limiter_(ratio, 2.0) * (downstreamMomentum - upstreamMomentum);
            return (upstreamMomentum + secondOrderTerm);
        }
        else
            return upstreamMomentum;
    }

    /**
      * \brief Tvd Scheme: Total Variation Diminishing
      *
      * This function manages the non uniformities of the grid according to [Li, Liao 2007].
      * It tries to reconstruct the value for the velocity at the upstream-upstream point
      * if the grid was uniform.
      */
    Scalar tvdLi(const std::array<Scalar,3>& momenta,
                 const std::array<Scalar,3>& distances,
                 const bool selfIsUpstream) const
    {
        using std::isfinite;
        const Scalar downstreamMomentum = momenta[0];
        const Scalar upstreamMomentum = momenta[1];
        const Scalar upUpstreamMomentum = momenta[2];
        const Scalar upstreamToDownstreamDistance = distances[0];
        const Scalar upUpstreamToUpstreamDistance = distances[1];

        // selfIsUpstream is required to get the correct sign because upUpstreamToUpstreamDistance is always positive
        const Scalar upUpstreamGradient = (upstreamMomentum - upUpstreamMomentum) / upUpstreamToUpstreamDistance * selfIsUpstream;

        // Distance between the upUpstream node and the position where it should be if the grid were uniform.
        const Scalar correctionDistance = upUpstreamToUpstreamDistance - upstreamToDownstreamDistance;
        const Scalar reconstrutedUpUpstreamVelocity = upUpstreamMomentum + upUpstreamGradient * correctionDistance;
        const Scalar ratio = (upstreamMomentum - reconstrutedUpUpstreamVelocity) / (downstreamMomentum - upstreamMomentum);

        // If the velocity field is uniform (like at the first newton step) we get a NaN
        if(ratio > 0.0 && isfinite(ratio))
        {
            const Scalar secondOrderTerm = 0.5 * limiter_(ratio, 2.0) * (downstreamMomentum - upstreamMomentum);
            return (upstreamMomentum + secondOrderTerm);
        }
        else
            return upstreamMomentum;
    }

    /**
     * \brief Tvd Scheme: Total Variation Diminishing
     *
     * This function manages the non uniformities of the grid according to [Hou, Simons, Hinkelmann 2007].
     * It should behave better then the Li's version in very stretched grids.
     */
    Scalar tvdHou(const std::array<Scalar,3>& momenta,
                  const std::array<Scalar,3>& distances,
                  const bool selfIsUpstream) const
    {
        using std::isfinite;
        const Scalar downstreamMomentum = momenta[0];
        const Scalar upstreamMomentum = momenta[1];
        const Scalar upUpstreamMomentum = momenta[2];
        const Scalar upstreamToDownstreamDistance = distances[0];
        const Scalar upUpstreamToUpstreamDistance = distances[1];
        const Scalar downstreamStaggeredCellSize = distances[2];
        const Scalar ratio = (upstreamMomentum - upUpstreamMomentum) / (downstreamMomentum - upstreamMomentum)
                           * upstreamToDownstreamDistance / upUpstreamToUpstreamDistance;

        // If the velocity field is uniform (like at the first newton step) we get a NaN
        if(ratio > 0.0 && isfinite(ratio))
        {
            const Scalar upstreamStaggeredCellSize = 0.5 * (upstreamToDownstreamDistance + upUpstreamToUpstreamDistance);
            const Scalar R = (upstreamStaggeredCellSize + downstreamStaggeredCellSize) / upstreamStaggeredCellSize;
            const Scalar secondOrderTerm = limiter_(ratio, R) / R * (downstreamMomentum - upstreamMomentum);
            return (upstreamMomentum + secondOrderTerm);
        }
        else
            return upstreamMomentum;
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
        using std::min;
        return min(r, 1.0);
    }

    /**
      * \brief SUPERBEE flux limiter function [Roe 1985]
      *
      * With R != 2 is the modified SUPERBEE flux limiter function [Hou, Simons, Hinkelmann 2007]
      */
    static Scalar superbee(const Scalar r, const Scalar R)
    {
        using std::min;
        using std::max;
        return max(min(R * r, 1.0), min(r, R));
    }

    /**
      * \brief UMIST flux limiter function [Lien and Leschziner 1993]
      */
    static Scalar umist(const Scalar r, const Scalar R)
    {
        using std::min;
        return min({R * r, (r * (5.0 - R) + R - 1.0) / 4.0, (r * (R - 1.0) + 5.0 - R) / 4.0, R});
    }

    /*
     * \brief Monotonized-Central limiter [Van Leer 1977]
     */
    static Scalar mclimiter(const Scalar r, const Scalar R)
    {
        using std::min;
        return min({R * r, (r + 1.0) / 2.0, R});
    }

    /**
      * \brief WAHYD Scheme [Hou, Simons, Hinkelmann 2007];
      */
    static Scalar wahyd(const Scalar r, const Scalar R)
    {
        using std::min;
        return r > 1 ? min((r + R * r * r) / (R + r * r), R)
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
    Scalar upwindWeight_;

    std::function<Scalar(const Scalar, const Scalar)> limiter_;
};

} // end namespace Dumux

#endif
