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
 * \ingroup Fluidmatrixinteractions
 * \brief Pore-local pc-Sw curves for cubic pore bodies.
 */
#ifndef DUMUX_PNM_2P_LOCAL_RULES_FOR_CUBE_HH
#define DUMUX_PNM_2P_LOCAL_RULES_FOR_CUBE_HH

#include <cmath>
#include <dumux/porenetworkflow/common/poreproperties.hh>
#include "../baselocalrules.hh"

namespace Dumux
{

/*!
 * \ingroup Fluidmatrixinteractions
 * \brief Implementation of the simplified pore-local capillary pressure-saturation curve
 *        according to Joekar-Niasar et al., 2010.
 */
template<class ScalarT>
struct TwoPLocalRulesCubeJoekarNiasar : public TwoPLocalRulesBase
{
    using Scalar = ScalarT;
    using Params = TwoPLocalRulesBase::Params<Scalar>;

    static constexpr bool supportsMultipleGeometries()
    { return false; }

    /*!
     * \brief The capillary pressure-saturation curve according to Joekar-Niasar et al., 2010.
     *
     *  \f$\mathrm{
     *  p_C = \frac{2*\sigma}{R(1-e^{-6.83S})}
     *  }\f$
     *
     * \param sw Saturation of the wetting phase \f$\mathrm{S_{w,i}}\f$ at pore \f$i\f$
     * \param params The parameters container
     */
    static Scalar pc(const Params& params, const Scalar sw)
    {
        assert(0 <= sw && sw <= 1);
        assert(params.shape == Pore::Shape::cube);
        const Scalar poreRadius = params.poreRadius;
        const Scalar sigma = params.surfaceTension ;
        // TODO incorporate contact angle!!!
        return 2*sigma / (poreRadius*(1 - std::exp(-6.83*sw))) ;
    }

    /*!
     * \brief The wetting-phase saturation of a pore body
     *
     * \param pc The capillary pressure \f$\mathrm{p_{c,i}}\f$ at pore \f$i\f$
     * \param params The parameters container
     */
    static Scalar sw(const Params& params, const Scalar pc)
    {
        assert(params.shape == Pore::Shape::cube);
        const Scalar poreRadius = params.poreRadius;
        const Scalar sigma = params.surfaceTension;
        using std::log;
        return  - 1/6.83* log(1 - 1/poreRadius * 2*sigma / pc);
    }

    /*!
     * \brief The partial derivative of the capillary
     *        pressure w.r.t. the wetting phase saturation.
     *
     *
     * \param sw Saturation of the wetting phase \f$\mathrm{S_{w,i}}\f$ at pore \f$i\f$
     * \param params A container object that is populated with the appropriate coefficients for the respective law.
     */
    static Scalar dpc_dsw(const Params& params, const Scalar sw)
    {
        assert(0 <= sw && sw <= 1);
        assert(params.shape == Pore::Shape::cube);
        using std::exp;
        const Scalar sigma = params.surfaceTension;
        const Scalar poreRadius = params.poreRadius;
        const Scalar e = exp(6.83*sw);
        return -(13.66*sigma*e) / (poreRadius*(e-1.0)*(e-1.0));
    }

    /*!
     * \brief DOCU
     *
     *
     * \param pc The capillary pressure \f$\mathrm{p_{c,i}}\f$ at pore \f$i\f$
     * \param params A container object that is populated with the appropriate coefficients for the respective law.
     */
    static Scalar dsw_dpc(const Params& params, const Scalar pc)
    {
        return 0; // TODO
    }
};

}

#endif
