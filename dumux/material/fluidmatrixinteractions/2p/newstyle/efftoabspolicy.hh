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
 *
 * \brief This is a policy for 2p material laws how to convert absolute to relative
 *        saturations and vice versa.
 *
 */
#ifndef DUMUX_TWOP_EFF_TO_ABS_POLICY_HH
#define DUMUX_TWOP_EFF_TO_ABS_POLICY_HH

namespace Dumux
{

/*!
 * \ingroup fluidmatrixinteractionsparams
 *
 * \brief A default implementation of the parameters for the adapter
 *        class to convert material laws from effective to absolute
 *        saturations.
 */
template <class Scalar>
class TwoPEffToAbsPolicyParams
{
public:

    /*!
     * \brief Constructor sets default
     */
    TwoPEffToAbsPolicyParams()
    { swr_ = snr_ = 0; }

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

private:
    Scalar swr_;
    Scalar snr_;
};

/*!
 * \ingroup fluidmatrixinteractionslaws
 *
 * \brief This is a policy for 2p material laws how to convert absolute to relative
 *        saturations and vice versa.
 *
 *        Material laws (like VanGenuchten or BrooksCorey) are defined for effective saturations.
 *        The numeric calculations however are performed with absolute saturations. The policy class converts
 *        the saturations. This allows for changing the calculation of the effective
 *        saturations easily, as this is subject of discussion / may be problem specific.
 *
 *        The EffToAbsPolicy is a template parameter for the actual material laws (C++ policy concept)
 */
template <class Scalar>
class TwoPEffToAbsPolicy
{
public:
    using Params = TwoPEffToAbsPolicyParams<Scalar>;
    /*!
     * \brief Convert an absolute wetting saturation to an effective one.
     *
     * \param sw Absolute saturation of the wetting phase \f$\mathrm{[{S}_w]}\f$.
     * \param params A container object that is populated with the appropriate coefficients for the respective law.
     *                  Therefore, in the (problem specific) spatialParameters  first, the material law is chosen,
     *                  and then the params container is constructed accordingly. Afterwards the values are set there, too.
     * \return Effective saturation of the wetting phase.
     */
    template<class Params>
    static Scalar swToSwe(const Params &params, const Scalar sw)
    {
        return (sw - params.swr())/(1. - params.swr() - params.snr());
    }

    /*!
     * \brief Convert an absolute non-wetting saturation to an effective one.
     *
     * \param sn Absolute saturation of the non-wetting phase \f$\mathrm{[{S}_n]}\f$.
     * \param params A container object that is populated with the appropriate coefficients for the respective law.
     *                  Therefore, in the (problem specific) spatialParameters  first, the material law is chosen,
     *                  and then the params container is constructed accordingly. Afterwards the values are set there, too.
     * \return Effective saturation of the non-wetting phase.
     */
    template<class Params>
    static Scalar snToSne(const Params &params, const Scalar sn)
    {
        return (sn - params.snr())/(1. - params.swr() - params.snr());
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
    template<class Params>
    static Scalar sweToSw(const Params &params, const Scalar swe)
    {
        return swe*(1. - params.swr() - params.snr()) + params.swr();
    }

    /*!
     * \brief Derivative of the effective saturation w.r.t. the absolute saturation.
     *
     * \param params A container object that is populated with the appropriate coefficients for the respective law.
     *                  Therefore, in the (problem specific) spatialParameters  first, the material law is chosen,
     *                  and then the params container is constructed accordingly. Afterwards the values are set there, too.
     * \return Derivative of the effective saturation w.r.t. the absolute saturation.
     */
    template<class Params>
    static Scalar dswe_dsw(const Params &params)
    {
        return 1.0/(1. - params.swr() - params.snr());
    }

    /*!
     * \brief Derivative of the absolute saturation w.r.t. the effective saturation.
     *
     * \param params A container object that is populated with the appropriate coefficients for the respective law.
     *                  Therefore, in the (problem specific) spatialParameters  first, the material law is chosen,
     *                  and then the params container is constructed accordingly. Afterwards the values are set there, too.
     * \return Derivative of the absolute saturation w.r.t. the effective saturation.
     */
    template<class Params>
    static Scalar dsw_dswe(const Params &params)
    {
        return 1. - params.swr() - params.snr();
    }
};

} // end namespace Dumux

#endif
