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
 * \ingroup Fluidmatrixinteractions
 * \brief This is a policy for 2p material laws how to convert absolute to relative
 *        saturations and vice versa.
 *
 */
#ifndef DUMUX_MATERIAL_FLUIDMATRIX_TWOP_EFF_TO_ABS_DEFAULT_POLICY_HH
#define DUMUX_MATERIAL_FLUIDMATRIX_TWOP_EFF_TO_ABS_DEFAULT_POLICY_HH

#include <dune/common/float_cmp.hh>
#include <dumux/common/parameters.hh>

namespace Dumux::FluidMatrix {

/*!
 * \ingroup Fluidmatrixinteractions
 *
 * \brief This is a policy for 2p material laws how to convert absolute to relative
 *        saturations and vice versa.
 *
 *        Material laws (like VanGenuchten or BrooksCorey) are defined for effective saturations.
 *        The numeric calculations however are performed with absolute saturations. The policy class converts
 *        the saturations. This allows for changing the calculation of the effective
 *        saturations easily, as this is subject of discussion / may be problem specific.
 */
class TwoPEffToAbsDefaultPolicy
{
public:
    /*!
     * \brief The parameter type
     * \tparam Scalar The scalar type
     * \note The efftoabs policy need two parameters: \f$\mathrm{S_{w,r}}, \mathrm{S_{n,r}}\f$.
     *       For the respective formulas check out the description of the free function.
     */
    template<class Scalar>
    struct Params
    {
        Params(const Scalar swr = 0.0, const Scalar snr = 0.0)
        : swr_(swr), snr_(snr)
        {}

        /*!
         * \brief Return the residual wetting saturation.
         */
        Scalar swr() const
        { return swr_; }

        /*!
         * \brief Set the residual wetting saturation.
         */
        void setSwr(Scalar v)
        { swr_ = v; }

        /*!
         * \brief Return the residual non-wetting saturation.
         */
        Scalar snr() const
        { return snr_; }

        /*!
         * \brief Set the residual non-wetting saturation.
         */
        void setSnr(Scalar v)
        { snr_ = v; }

        bool operator== (const Params& p) const
        {
            return Dune::FloatCmp::eq(swr(), p.swr(), 1e-6)
                   && Dune::FloatCmp::eq(snr(), p.snr(), 1e-6);
        }
    private:
        Scalar swr_;
        Scalar snr_;
    };

    /*!
     * \brief Construct from a subgroup from the global parameter tree
     * \note This will give you nice error messages if a mandatory parameter is missing
     */
    template<class Scalar>
    static Params<Scalar> makeParams(const std::string& paramGroup)
    {
        Params<Scalar> params;
        params.setSwr(getParamFromGroup<Scalar>(paramGroup, "Swr", 0.0));
        params.setSnr(getParamFromGroup<Scalar>(paramGroup, "Snr", 0.0));
        return params;
    }

    /*!
     * \brief Convert an absolute wetting saturation to an effective one.
     *
     * \param sw Absolute saturation of the wetting phase \f$\mathrm{[{S}_w]}\f$.
     * \param params A container object that is populated with the appropriate coefficients for the respective law.
     *                  Therefore, in the (problem specific) spatialParameters  first, the material law is chosen,
     *                  and then the params container is constructed accordingly. Afterwards the values are set there, too.
     * \return Effective saturation of the wetting phase.
     */
    template<class Scalar>
    static Scalar swToSwe(const Scalar sw, const Params<Scalar>& params)
    {
        return (sw - params.swr())/(1.0 - params.swr() - params.snr());
    }

    /*!
     * \brief Convert an effective wetting saturation to an absolute one.
     *
     * \param swe Effective saturation of the non-wetting phase \f$\mathrm{[\overline{S}_n]}\f$.
     * \param params A container object that is populated with the appropriate coefficients for the respective law.
     *                  Therefore, in the (problem specific) spatialParameters  first, the material law is chosen,
     *                  and then the params container is constructed accordingly. Afterwards the values are set there, too.
     * \return Absolute saturation of the non-wetting phase.
     */
    template<class Scalar>
    static Scalar sweToSw(const Scalar swe, const Params<Scalar>& params)
    {
        return swe*(1.0 - params.swr() - params.snr()) + params.swr();
    }

    /*!
     * \brief Derivative of the effective saturation w.r.t. the absolute saturation.
     *
     * \param params A container object that is populated with the appropriate coefficients for the respective law.
     *                  Therefore, in the (problem specific) spatialParameters  first, the material law is chosen,
     *                  and then the params container is constructed accordingly. Afterwards the values are set there, too.
     * \return Derivative of the effective saturation w.r.t. the absolute saturation.
     */
    template<class Scalar>
    static Scalar dswe_dsw(const Params<Scalar>& params)
    {
        return 1.0/(1.0 - params.swr() - params.snr());
    }

    /*!
     * \brief Derivative of the absolute saturation w.r.t. the effective saturation.
     *
     * \param params A container object that is populated with the appropriate coefficients for the respective law.
     *                  Therefore, in the (problem specific) spatialParameters  first, the material law is chosen,
     *                  and then the params container is constructed accordingly. Afterwards the values are set there, too.
     * \return Derivative of the absolute saturation w.r.t. the effective saturation.
     */
    template<class Scalar>
    static Scalar dsw_dswe(const Params<Scalar>& params)
    {
        return 1.0 - params.swr() - params.snr();
    }
};

} // end namespace Dumux::FluidMatrix

#endif
