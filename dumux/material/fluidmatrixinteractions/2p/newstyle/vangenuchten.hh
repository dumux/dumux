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
 * \brief   Implementation of the capillary pressure and
 *          relative permeability <-> saturation relations according to van Genuchten.
 */
#ifndef DUMUX_VAN_GENUCHTEN_HH
#define DUMUX_VAN_GENUCHTEN_HH

#include <algorithm>
#include <cmath>
#include <cassert>

#include "efftoabspolicy.hh"
#include "vangenuchtenregularization.hh"
#include "vangenuchtenparams.hh"

namespace Dumux
{
/*!
 * \ingroup fluidmatrixinteractionslaws
 *
 * \brief Implementation of the van Genuchten capillary pressure <->
 *        saturation relation.
 * \tparam Scalar the scalar type
 * \tparam RegularizationPolicy the policy on what to do for low and high saturations
 *         should implement the methods:
 * \tparam EffToAbsPolicy the policy how to convert absolute to effective saturations
 *         and vice versa. Should implement the methods: swToSwe, SweToSw, dswe_dsw, dsw_dswe.
 *
 * \see VanGenuchtenParams
 */
template <class Scalar,
          class RegularizationPolicyType = VanGenuchtenRegularization<Scalar>,
          class EffToAbsPolicyType = TwoPEffToAbsPolicy<Scalar>>
class VanGenuchten
{
    using ThisType = VanGenuchten<Scalar, RegularizationPolicyType, EffToAbsPolicyType>;
public:
    //! export our parameter class and the policies
    using RegularizationPolicy = RegularizationPolicyType;
    using EffToAbsPolicy = EffToAbsPolicyType;
    using Params = VanGenuchtenParams<Scalar, ThisType>;

    /*!
     * \brief The capillary pressure-saturation curve according to van Genuchten.
     * \note The regularization depends on the regularization policy
     * \param sw Absolute saturation of the wetting phase \f$\mathrm{S_w}\f$
     */
    static Scalar pc(const Params &params, const Scalar sw)
    {
        //! regularization for low saturations
        if (RegularizationPolicy::underLowSwThreshold_pc(params, sw))
            return RegularizationPolicy::pcForLowSw(params, sw);

        //! regularization for high saturations
        else if (RegularizationPolicy::overHighSwThreshold_pc(params, sw))
            return RegularizationPolicy::pcForHighSw(params, sw);

        //! no regularization if saturation is in good range
        else
            return pcRaw(params, EffToAbsPolicy::swToSwe(params, sw));
    }

    /*!
     * \brief The regularized saturation-capillary pressure curve according to van Genuchten.
     * \note The regularization depends on the regularization policy
     * \param pc The capillary pressure \f$\mathrm{p_C}\f$ in \f$\mathrm{[Pa]}\f$
     */
    static Scalar sw(const Params &params, const Scalar pc)
    {
        //! regularization for low pcs
        if (RegularizationPolicy::underLowPcThreshold_sw(params, pc))
            return RegularizationPolicy::swForLowPc(params, pc);

        //! regularization for high pcs
        else if (RegularizationPolicy::overHighPcThreshold_sw(params, pc))
            return RegularizationPolicy::swForHighPc(params, pc);

        //! no regularization if pc is in good range
        else
            return EffToAbsPolicy::sweToSw(params, sweRaw(params, pc));
    }

    /*!
     * \brief The regularized partial derivative of the capillary
     *        pressure w.r.t. the effective saturation according to van Genuchten.
     * \note The regularization depends on the regularization policy
     * \param sw Absolute saturation of the wetting phase \f$\mathrm{S_w}\f$
     */
    // static Scalar dpc_dsw(const Params &params, Scalar sw)
    // {
    //     auto swe = EffToAbsPolicy::swToSwe(params, sw);

    //     //! regularization for low saturations
    //     if (swe < RegularizationPolicy::lowSweThreshold())
    //         return RegularizationPolicy::dpc_dswLowSwe(params, swe);

    //     //! regularization for high saturations
    //     else if (swe > RegularizationPolicy::highSweThreshold())
    //         return RegularizationPolicy::dpc_dswHighSwe(params, swe);

    //     //! no regularization if saturation is in good range
    //     else
    //         return dpc_dsweRaw(params, swe)*EffToAbsPolicy::dswe_dsw(params);
    // }

    /*!
     * \brief The partial derivative of the effective
     *        saturation to the capillary pressure according to van Genuchten.
     *
     * \param pc Capillary pressure \f$\mathrm{p_C}\f$ in \f$\mathrm{[Pa]}\f$
     * \param params A container object that is populated with the appropriate coefficients for the respective law.
     *                  Therefore, in the (problem specific) spatialParameters  first, the material law is chosen, and then the params container
     *                  is constructed accordingly. Afterwards the values are set there, too.
     *
     * \note Instead of undefined behaviour if pc is not in the valid range, we return a valid number,
     *       by clamping the input.
     */
    // static Scalar dsw_dpc(const Params &params, Scalar pc)
    // {
    //     using std::pow;
    //     using std::max;

    //     pc = max(pc, 0.0); // the equation below is undefined for negative pcs

    //     const Scalar powAlphaPc = pow(params.vgAlpha()*pc, params.vgn());
    //     const Scalar dswe_dpc = -pow(powAlphaPc + 1, -params.vgm()-1)*params.vgm()*powAlphaPc/pc*params.vgn();

    //     return dswe_dpc*EffToAbsPolicy::dsw_dswe(params);
    // }

    /*!
     * \brief The regularized relative permeability for the wetting phase of
     *        the medium implied by van Genuchten's parameterization.
     * \note The regularization depends on the regularization policy
     * \param sw Absolute saturation of the wetting phase \f$\mathrm{S_w}\f$
     */
    static Scalar krw(const Params &params, const Scalar sw)
    {
        //! regularization for high saturations
        if (RegularizationPolicy::overHighSwThreshold_krw(params, sw))
            return RegularizationPolicy::krwForHighSw(params, sw);

        //! no regularization if saturation is in good range
        else
            return krwRaw(params, EffToAbsPolicy::swToSwe(params, sw));
    }

    /*!
     * \brief The regularized relative permeability for the non-wetting phase of
     *        the medium implied by van Genuchten's parameterization.
     * \note The regularization depends on the regularization policy
     * \param sw Absolute saturation of the wetting phase \f$\mathrm{S_w}\f$
     */
    static Scalar krn(const Params &params, const Scalar sw)
    {
        //! regularization for low saturations
        if (RegularizationPolicy::underLowSwThreshold_krn(params, sw))
            return RegularizationPolicy::krnForLowSw(params, sw);

        //! no regularization if saturation is in good range
        else
            return krnRaw(params, EffToAbsPolicy::swToSwe(params, sw));
    }

    /*!
     * \brief The capillary pressure-saturation curve according to van Genuchten.
     *
     * Van Genuchten's empirical capillary pressure <-> saturation
     * function is given by
     * \f$\mathrm{
     p_C = (\overline{S}_w^{-1/m} - 1)^{1/n}/\alpha
     }\f$
     * \param swe Effective saturation of the wetting phase \f$\mathrm{\overline{S}_w}\f$
     * \param params A container object that is populated with the appropriate coefficients for the respective law.
     *                  Therefore, in the (problem specific) spatialParameters  first, the material law is chosen, and then the params container
     *                  is constructed accordingly. Afterwards the values are set there, too.
     * \note Instead of undefined behaviour if swe is not in the valid range, we return a valid number,
     *       by clamping the input. Note that for pc(swe = 0.0) = inf, have a look at RegularizedVanGenuchten if this is a problem.
     */
    static Scalar pcRaw(const Params &params, Scalar swe)
    {
        // the equation below is only defined for 0.0 <= sw <= 1.0
        using std::pow; using std::min; using std::max;
        swe = min(max(swe, 0.0), 1.0);

        return pow(pow(swe, -1.0/params.vgm()) - 1, 1.0/params.vgn())/params.vgAlpha();
    }

    /*!
     * \brief The saturation-capillary pressure curve according to van Genuchten.
     *
     * This is the inverse of the capillary pressure-saturation curve:
     * \f$\mathrm{
     \overline{S}_w = {p_C}^{-1} = ((\alpha p_C)^n + 1)^{-m}
     }\f$
     *
     * \param pc Capillary pressure \f$\mathrm{p_C}\f$ in \f$\mathrm{[Pa]}\f$
     * \param params A container object that is populated with the appropriate coefficients for the respective law.
     *                  Therefore, in the (problem specific) spatialParameters  first, the material law is chosen, and then the params container
     *                  is constructed accordingly. Afterwards the values are set there, too.
     * \return          The effective saturation of the wetting phase \f$\mathrm{\overline{S}_w}\f$
     * \note Instead of undefined behaviour if pc is not in the valid range, we return a valid number,
     *       i.e. sw(pc < 0.0) = 0.0, by clamping the input to the physical bounds.
     */
    static Scalar sweRaw(const Params &params, Scalar pc)
    {
        using std::pow;
        using std::max;

        // the equation below is undefined for negative pcs
        pc = max(pc, 0.0);

        return pow(pow(params.vgAlpha()*pc, params.vgn()) + 1, -params.vgm());
    }

    /*!
     * \brief The partial derivative of the capillary
     *        pressure w.r.t. the effective saturation according to van Genuchten.
     *
     * This is equivalent to
     * \f$\mathrm{
     \frac{\partial p_C}{\partial \overline{S}_w} =
     -\frac{1}{\alpha} (\overline{S}_w^{-1/m} - 1)^{1/n - }
     \overline{S}_w^{-1/m} / \overline{S}_w / m
     }\f$
     *
     * \param swe Effective saturation of the wetting phase \f$\mathrm{\overline{S}_w}\f$
     * \param params A container object that is populated with the appropriate coefficients for the respective law.
     *                  Therefore, in the (problem specific) spatialParameters  first, the material law is chosen, and then the params container
     *                  is constructed accordingly. Afterwards the values are set there, too.
     *
     * \note Instead of undefined behaviour if swe is not in the valid range, we return a valid number,
     *       by clamping the input.
     */
    static Scalar dpc_dsweRaw(const Params &params, Scalar swe)
    {
        using std::pow;
        using std::min;
        using std::max;

        // the equation below is only defined for 0.0 <= sw <= 1.0
        swe = min(max(swe, 0.0), 1.0);

        const Scalar powSwe = pow(swe, -1/params.vgm());
        const Scalar dpc_dswe = - 1.0/params.vgAlpha() * pow(powSwe - 1, 1.0/params.vgn() - 1)/params.vgn()
                                      * powSwe/swe/params.vgm();

        return dpc_dswe;
    }

    /*!
     * \brief The relative permeability for the wetting phase of
     *        the medium implied by van Genuchten's
     *        parameterization.
     *
     * \param swe The mobile saturation of the wetting phase.
     * \param params A container object that is populated with the appropriate coefficients for the respective law.
     *                  Therefore, in the (problem specific) spatialParameters  first, the material law is chosen, and then the params container
     *                  is constructed accordingly. Afterwards the values are set there, too.
     *
     * \note Instead of undefined behaviour if pc is not in the valid range, we return a valid number,
     *       by clamping the input.
     */
    static Scalar krwRaw(const Params &params, Scalar swe)
    {
        using std::pow;
        using std::sqrt;
        using std::min;
        using std::max;

        swe = min(max(swe, 0.0), 1.0); // the equation below is only defined for 0.0 <= sw <= 1.0

        const Scalar r = 1.0 - pow(1.0 - pow(swe, 1.0/params.vgm()), params.vgm());
        return sqrt(swe)*r*r;
    }

    /*!
     * \brief The relative permeability for the non-wetting phase
     *        of the medium implied by van Genuchten's
     *        parameterization.
     *
     * \param swe The mobile saturation of the wetting phase.
     * \param params A container object that is populated with the appropriate coefficients for the respective law.
     *                  Therefore, in the (problem specific) spatialParameters  first, the material law is chosen, and then the params container
     *                  is constructed accordingly. Afterwards the values are set there, too.
     *
     * \note Instead of undefined behaviour if pc is not in the valid range, we return a valid number,
     *       by clamping the input.
     */
    static Scalar krnRaw(const Params &params, Scalar swe)
    {
        using std::pow;
        using std::cbrt;
        using std::min;
        using std::max;

        swe = min(max(swe, 0.0), 1.0); // the equation below is only defined for 0.0 <= sw <= 1.0

        return cbrt(1 - swe) * pow(1 - pow(swe, 1.0/params.vgm()), 2*params.vgm());
    }

    /*!
     * \brief The derivative of the relative permeability for the
     *        wetting phase in regard to the wetting saturation of the
     *        medium implied by the van Genuchten parameterization.
     *
     * \param swe The mobile saturation of the wetting phase.
     * \param params A container object that is populated with the appropriate coefficients for the respective law.
     *                  Therefore, in the (problem specific) spatialParameters  first, the material law is chosen, and then the params container
     *                  is constructed accordingly. Afterwards the values are set there, too.
     *
     * \note Instead of undefined behaviour if pc is not in the valid range, we return a valid number,
     *       by clamping the input.
     */
    static Scalar dkrw_dsweRaw(const Params &params, Scalar swe)
    {
        using std::pow;
        using std::sqrt;
        using std::min;
        using std::max;

        swe = min(max(swe, 0.0), 1.0); // the equation below is only defined for 0.0 <= sw <= 1.0

        const Scalar x = 1.0 - pow(swe, 1.0/params.vgm());
        const Scalar xToM = pow(x, params.vgm());
        const Scalar dkrw_dswe = (1.0 - xToM)/sqrt(swe) * ( (1.0 - xToM)/2 + 2*xToM*(1.0-x)/x );

        return dkrw_dswe;
    }

    /*!
     * \brief The derivative of the relative permeability for the
     *        non-wetting phase in regard to the wetting saturation of
     *        the medium as implied by the van Genuchten
     *        parameterization.
     *
     * \param swe The mobile saturation of the wetting phase.
     * \param params A container object that is populated with the appropriate coefficients for the respective law.
     *                  Therefore, in the (problem specific) spatialParameters  first, the material law is chosen, and then the params container
     *                  is constructed accordingly. Afterwards the values are set there, too.
     *
     * \note Instead of undefined behaviour if pc is not in the valid range, we return a valid number,
     *       by clamping the input.
     */
    static Scalar dkrn_dsweRaw(const Params &params, Scalar swe)
    {
        using std::pow;
        using std::min;
        using std::max;

        swe = min(max(swe, 0.0), 1.0); // the equation below is only defined for 0.0 <= sw <= 1.0

        const Scalar x = pow(swe, 1.0/params.vgm());
        const Scalar dkrn_dswe = -pow(1.0 - x, 2*params.vgm()) * pow(1.0 - swe, -2.0/3) * (1.0/3 + 2*x/swe);

        return dkrn_dswe;
    }

    /*!
     * \brief The capillary pressure at Swe = 1.0 also called end point capillary pressure
     *
     * \param params A container object that is populated with the appropriate coefficients for the respective law.
     *                  Therefore, in the (problem specific) spatialParameters  first, the material law is chosen, and then the params container
     *                  is constructed accordingly. Afterwards the values are set there, too.
     */
    static Scalar endPointPc(const Params &params)
    { return 0.0; }
};

} // end namespace Dumux

#endif
