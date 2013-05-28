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
#ifndef VAN_GENUCHTEN_HH
#define VAN_GENUCHTEN_HH

#include "vangenuchtenparams.hh"

#include <algorithm>
#include <cmath>
#include <cassert>

namespace Dumux
{
/*!
 * \ingroup fluidmatrixinteractionslaws
 *
 * \brief Implementation of the van Genuchten capillary pressure <->
 *        saturation relation. This class bundles the "raw" curves
 *        as static members and doesn't concern itself converting
 *        absolute to effective saturations and vice versa.
 *
 * For general info: EffToAbsLaw
 *
 * \see VanGenuchtenParams
 */
template <class ScalarT, class ParamsT = VanGenuchtenParams<ScalarT> >
class VanGenuchten
{
public:
    typedef ParamsT     Params;
    typedef typename    Params::Scalar Scalar;

    /*!
     * \brief The capillary pressure-saturation curve according to van Genuchten.
     *
     * Van Genuchten's empirical capillary pressure <-> saturation
     * function is given by
     * \f[
     p_C = (\overline{S}_w^{-1/m} - 1)^{1/n}/\alpha
     \f]
     * \param Swe       Effective saturation of the wetting phase \f$\overline{S}_w\f$
     * \param params    A container object that is populated with the appropriate coefficients for the respective law.
     *                  Therefore, in the (problem specific) spatialParameters  first, the material law is chosen, and then the params container
     *                  is constructed accordingly. Afterwards the values are set there, too.
     */
    static Scalar pc(const Params &params, Scalar Swe)
    {
        assert(0 <= Swe && Swe <= 1);
        return pow(pow(Swe, -1.0/params.vgm()) - 1, 1.0/params.vgn())/params.vgAlpha();
    }

    DUNE_DEPRECATED_MSG("use pc() (uncapitalized 'c') instead")
    static Scalar pC(const Params &params, Scalar Swe)
    {
        return pc(params, Swe);
    }

    /*!
     * \brief The saturation-capillary pressure curve according to van Genuchten.
     *
     * This is the inverse of the capillary pressure-saturation curve:
     * \f[
     \overline{S}_w = {p_C}^{-1} = ((\alpha p_C)^n + 1)^{-m}
     \f]
     *
     * \param pC        Capillary pressure \f$p_C\f$
     * \param params    A container object that is populated with the appropriate coefficients for the respective law.
     *                  Therefore, in the (problem specific) spatialParameters  first, the material law is chosen, and then the params container
     *                  is constructed accordingly. Afterwards the values are set there, too.
     * \return          The effective saturation of the wetting phase \f$\overline{S}_w\f$
     */
    static Scalar sw(const Params &params, Scalar pC)
    {
        assert(pC >= 0);

        return pow(pow(params.vgAlpha()*pC, params.vgn()) + 1, -params.vgm());
    }

    DUNE_DEPRECATED_MSG("use sw() (uncapitalized 's') instead")
    static Scalar Sw(const Params &params, Scalar pC)
    {
        return sw(params, pC);
    }

    /*!
     * \brief The partial derivative of the capillary
     *        pressure w.r.t. the effective saturation according to van Genuchten.
     *
     * This is equivalent to
     * \f[
     \frac{\partial p_C}{\partial \overline{S}_w} =
     -\frac{1}{\alpha} (\overline{S}_w^{-1/m} - 1)^{1/n - }
     \overline{S}_w^{-1/m} / \overline{S}_w / m
     \f]
     *
     * \param Swe       Effective saturation of the wetting phase \f$\overline{S}_w\f$
     * \param params    A container object that is populated with the appropriate coefficients for the respective law.
     *                  Therefore, in the (problem specific) spatialParameters  first, the material law is chosen, and then the params container
     *                  is constructed accordingly. Afterwards the values are set there, too.
    */
    static Scalar dpc_dsw(const Params &params, Scalar Swe)
    {
        assert(0 <= Swe && Swe <= 1);

        Scalar powSwe = pow(Swe, -1/params.vgm());
        return - 1/params.vgAlpha() * pow(powSwe - 1, 1/params.vgn() - 1)/params.vgn()
            * powSwe/Swe/params.vgm();
    }

    DUNE_DEPRECATED_MSG("use dpc_dsw() (uncapitalized 'c', 's') instead")
    static Scalar dpC_dSw(const Params &params, Scalar Swe)
    {
        return dpc_dsw(params, Swe);
    }

    /*!
     * \brief The partial derivative of the effective
     *        saturation to the capillary pressure according to van Genuchten.
     *
     * \param pC        Capillary pressure \f$p_C\f$
     * \param params    A container object that is populated with the appropriate coefficients for the respective law.
     *                  Therefore, in the (problem specific) spatialParameters  first, the material law is chosen, and then the params container
     *                  is constructed accordingly. Afterwards the values are set there, too.
     */
    static Scalar dsw_dpc(const Params &params, Scalar pC)
    {
        assert(pC >= 0);

        Scalar powAlphaPc = pow(params.vgAlpha()*pC, params.vgn());
        return -pow(powAlphaPc + 1, -params.vgm()-1)*
            params.vgm()*powAlphaPc/pC*params.vgn();
    }

    DUNE_DEPRECATED_MSG("use dsw_dpc() (uncapitalized 's', 'c') instead")
    static Scalar dSw_dpC(const Params &params, Scalar pC)
    {
        return dsw_dpc(params, pC);
    }

    /*!
     * \brief The relative permeability for the wetting phase of
     *        the medium implied by van Genuchten's
     *        parameterization.
     *
     * \param Swe        The mobile saturation of the wetting phase.
     * \param params    A container object that is populated with the appropriate coefficients for the respective law.
     *                  Therefore, in the (problem specific) spatialParameters  first, the material law is chosen, and then the params container
     *                  is constructed accordingly. Afterwards the values are set there, too.     */
    static Scalar krw(const Params &params, Scalar Swe)
    {
        assert(0 <= Swe && Swe <= 1);

        Scalar r = 1. - pow(1 - pow(Swe, 1/params.vgm()), params.vgm());
        return sqrt(Swe)*r*r;
    };

    /*!
     * \brief The derivative of the relative permeability for the
     *        wetting phase in regard to the wetting saturation of the
     *        medium implied by the van Genuchten parameterization.
     *
     * \param Swe       The mobile saturation of the wetting phase.
     * \param params    A container object that is populated with the appropriate coefficients for the respective law.
     *                  Therefore, in the (problem specific) spatialParameters  first, the material law is chosen, and then the params container
     *                  is constructed accordingly. Afterwards the values are set there, too.
     */
    static Scalar dkrw_dsw(const Params &params, Scalar Swe)
    {
        assert(0 <= Swe && Swe <= 1);

        const Scalar x = 1 - std::pow(Swe, 1.0/params.vgm());
        const Scalar xToM = std::pow(x, params.vgm());
        return (1 - xToM)/std::sqrt(Swe) * ( (1 - xToM)/2 + 2*xToM*(1-x)/x );
    };

    DUNE_DEPRECATED_MSG("use dkrw_dsw() (uncapitalized 's') instead")
    static Scalar dkrw_dSw(const Params &params, Scalar Swe)
    {
        return dkrw_dsw(params, Swe);
    }


    /*!
     * \brief The relative permeability for the non-wetting phase
     *        of the medium implied by van Genuchten's
     *        parameterization.
     *
     * \param Swe        The mobile saturation of the wetting phase.
     * \param params    A container object that is populated with the appropriate coefficients for the respective law.
     *                  Therefore, in the (problem specific) spatialParameters  first, the material law is chosen, and then the params container
     *                  is constructed accordingly. Afterwards the values are set there, too.
     */
    static Scalar krn(const Params &params, Scalar Swe)
    {
        assert(0 <= Swe && Swe <= 1);

        return
            pow(1 - Swe, 1.0/3) *
            pow(1 - pow(Swe, 1/params.vgm()), 2*params.vgm());
    };

    /*!
     * \brief The derivative of the relative permeability for the
     *        non-wetting phase in regard to the wetting saturation of
     *        the medium as implied by the van Genuchten
     *        parameterization.
     *
     * \param Swe        The mobile saturation of the wetting phase.
     * \param params    A container object that is populated with the appropriate coefficients for the respective law.
     *                  Therefore, in the (problem specific) spatialParameters  first, the material law is chosen, and then the params container
     *                  is constructed accordingly. Afterwards the values are set there, too.
     */
    static Scalar dkrn_dsw(const Params &params, Scalar Swe)
    {
        assert(0 <= Swe && Swe <= 1);

        const Scalar x = std::pow(Swe, 1.0/params.vgm());
        return
            -std::pow(1 - x, 2*params.vgm())
            *std::pow(1 - Swe, -2/3)
            *(1.0/3 + 2*x/Swe);
    }

    DUNE_DEPRECATED_MSG("use dkrn_dsw() (uncapitalized 's') instead")
    static Scalar dkrn_dSw(const Params &params, Scalar Swe)
    {
        return dkrn_dsw(params, Swe);
    }

};
}

#endif
