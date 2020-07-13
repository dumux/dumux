// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   See the file COPYING for full copying permissions.                      *
 *                                                                           *
 *   This program is free software: you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation, either version 3 of the License, or       *
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
 * \brief Implementation of the capillary pressure and
 * relative permeability <-> saturation relations according to Joekar-Niasar et al., 2010.
 *
 */
#ifndef DUMUX_PNM_LOCAL_RULES_HH
#define DUMUX_PNM_LOCAL_RULES_HH

#include "porenetworklocalrulesparams.hh"

#include <algorithm>

namespace Dumux
{

/*!
 * \ingroup fluidmatrixinteractionslaws
 *
 * \brief Implementation of the pore network capillary pressure <->
 *        saturation relation. This class bundles the "raw" curves
 *        as static members and doesn't concern itself converting
 *        absolute to effective saturations and vice versa.
 *        References: Joekar-Niasar et al., 2010
 *
 *\see PNMLocalRulesParams
 */
// template <class ScalarT, class ParamsT = PNMLocalRulesParams<ScalarT>>
template <class ScalarT, class ParamsT >
class PNMLocalRules
{

public:
    using Params = ParamsT;
    using Scalar = ScalarT;

    /*!
     * \brief The capillary pressure-saturation curve according to Joekar-Niasar et al., 2010.
     *
     * Empirical  capillary pressure <-> saturation
     * function is given by
     *
     *  \f$\mathrm{
     *  p_C = \frac{2*\sigma}{R(1-e^{-6.83S})}
     *  }\f$
     *
     * \param sw Saturation of the wetting phase \f$\mathrm{[\overline{S}_w]}\f$
     * \param params A container object that is populated with the appropriate coefficients for the respective law.
     *                  Therefore, in the (problem specific) spatialParameters  first, the material law is chosen,
     *                  and then the params container is constructed accordingly. Afterwards the values are set there, too.
     * \param poreRadius The pore body radius
     * \return Capillary pressure calculated by Joekar-Niasar et al., 2010 constitutive relation.
     */
    static Scalar pc(const Params& params, const Scalar sw)
    {
        assert(0 <= sw && sw <= 1);
        const Scalar poreRadius = params.poreRadius();
        const Scalar sigma = params.surfaceTension() ;
        return 2*sigma / (poreRadius*(1 - std::exp(-6.83*sw))) ;
    }

    /*!
    * \brief The curvature in a plane perpendicular to the direction of the throat
    *
    * \param params A container object that is populated with the appropriate coefficients for the respective law.
    *                  Therefore, in the (problem specific) spatialParameters  first, the material law is chosen,
    *                  and then the params container is constructed accordingly. Afterwards the values are set there, too.
    * \param pc The capillary pressure
    */
    static Scalar curvatureRadius(const Params &params, const Scalar pc)
    {
        const Scalar sigma = params.sigma();
        const Scalar theta = params.theta();
        const Scalar cosTheta = std::cos(theta);
        const Scalar sinTheta = std::sin(theta);
        return sigma / pc * (cosTheta - sinTheta);
    }

     /*!
     * \brief The minimum wetting-phase saturation of a pore body
     *
     * \param params A container object that is populated with the appropriate coefficients for the respective law.
     *                  Therefore, in the (problem specific) spatialParameters  first, the material law is chosen,
     *                  and then the params container is constructed accordingly. Afterwards the values are set there, too.
     * \param poreRadius The pore body radius
     * \param pcMin The smaller of either the global maximum pressure difference or the smallest entry capillary pressure of all throats not yet invaded
     */
    static Scalar swMin(const Params &params, const Scalar poreRadius, const Scalar pcMin)
    {
        const Scalar sigma = params.sigma() ;
        Scalar result = - 1/6.83* std::log(1 - 1/poreRadius * 2*sigma / pcMin);
        if (result != result)
         return 0.0;
        else
            return result;
    }

     /*!
     * \brief The wetting-phase saturation of a pore body
     *
     * \param params A container object that is populated with the appropriate coefficients for the respective law.
     *                  Therefore, in the (problem specific) spatialParameters  first, the material law is chosen,
     *                  and then the params container is constructed accordingly. Afterwards the values are set there, too.
     * \param poreRadius The pore body radius
     * \param pcMin The capillary pressure
     */
    static Scalar sw(const Params& params, const Scalar pc)
    {
        const Scalar poreRadius = params.poreRadius();
        const Scalar sigma = params.surfaceTension();
        using std::log;
        return  - 1/6.83* log(1 - 1/poreRadius * 2*sigma / pc);
    }

    /*!
    * \brief The partial derivative of the capillary
    *        pressure w.r.t. the wetting phase saturation.
    *
    *
    * \param sw Saturation of the wetting phase \f$\mathrm{[\overline{S}_w]}\f$
    * \param params A container object that is populated with the appropriate coefficients for the respective law.
    */
    static Scalar dpc_dsw(const Params &params, const Scalar sw)
    {
        assert(0 <= sw && sw <= 1);
        using std::exp;
        const Scalar sigma = params.surfaceTension();
        const Scalar poreRadius = params.poreRadius();
        const Scalar e = exp(6.83*sw);
        return -(13.66*sigma*e) / (poreRadius*(e-1.0)*(e-1.0));
    }

     /*!
     * \brief The partial derivative of the capillary
     *        pressure w.r.t. the effective saturation.
     *
     *
     * \param swe Effective saturation of the wetting phase \f$\mathrm{[\overline{S}_w]}\f$
     * \param params A container object that is populated with the appropriate coefficients for the respective law.
     *                  Therefore, in the (problem specific) spatialParameters  first, the material law is chosen, and then the params container
     *                  is constructed accordingly. Afterwards the values are set there, too.
     * \param poreRadius The pore body radius
     */
    static Scalar dpc_dsw(const Params &params, const Scalar poreRadius, const Scalar sw)
    {
        assert(0 <= sw && sw <= 1);
        using std::exp;
        const Scalar sigma = params.surfaceTension();
        const Scalar e = exp(6.83*sw);
        return -(13.66*sigma*e) / (poreRadius*(e-1.0)*(e-1.0));
    }


    template<typename... Args>
    static Scalar krw(Args&&... args)
    {
        DUNE_THROW(Dune::NotImplemented, "krw for PNM is not implemented.");
    }

    template<typename... Args>
    static Scalar krn(Args&&... args)
    {
        DUNE_THROW(Dune::NotImplemented, "krw for PNM is not implemented.");
    }

};
}

#endif // DUMUX_PNM_LOCAL_RULES_HH
