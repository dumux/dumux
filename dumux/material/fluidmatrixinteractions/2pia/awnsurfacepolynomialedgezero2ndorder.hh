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
 *   along with this program. If not, see <http://www.gnu.org/licenses/>.    *
 *****************************************************************************/
/*!
 * \file
 * \ingroup Fluidmatrixinteractions
 * \brief Specification of a function relating volume specific interfacial area to capillary pressure and saturation.
 * This parametrization is a second order polynomial which is zero for saturations of zero and one.
 */
#ifndef AWN_SURFACE_POLYNOMIAL_EDGE_ZERO_2ND_ORDER_HH
#define AWN_SURFACE_POLYNOMIAL_EDGE_ZERO_2ND_ORDER_HH

#include "awnsurfacepolynomialedgezero2ndorderparams.hh"

#include <dune/common/exceptions.hh>

#include <algorithm>
#include <math.h>
#include <assert.h>

#warning "This header is deprecated. Removal after 3.3. Use new material laws."

namespace Dumux {

/*!
 * \ingroup Fluidmatrixinteractions
 * \brief Implementation of the polynomial of second order relating
 *        specific interfacial  area to wetting phase saturation and capillary pressure.
 */
template <class ParamsT>
class [[deprecated("Use new material laws and FluidMatrix::InterfacialAreaolynomialEdgeZero2ndOrder instead!")]] AwnSurfacePolynomialEdgeZero2ndOrder
{
public:
    using Params = ParamsT;
    using Scalar = typename Params::Scalar;

    /*!
     * \brief The awn surface
     * the suggested (as estimated from pore network models) awn surface:
     * \f$[
     a_{wn} = a_{1} (S_{wr}-S_{w})(1.-S_{w}) + a_2 (S_{wr}-S_{w})(1.-S_{w}) p_{c} + a_{3} (S_{wr}-S_{w})(1.-S_{w}) p_{c}^2
     \f$]
     * \param  params parameter container for the coefficients of the surface
     * \param  Sw Effective saturation of the wetting phase
     * \param  pc Capillary pressure in \f$\mathrm{[Pa]}\f$
     */
    static Scalar interfacialArea(const Params &params, const Scalar Sw, const Scalar pc)
    {
        const Scalar a1     = params.a1();
        const Scalar a2     = params.a2();
        const Scalar a3     = params.a3();
        const Scalar Swr    = params.Swr();
        const Scalar factor = (Swr-Sw)*(1.-Sw) ;

        const Scalar aAlphaBeta = a1*factor + a2*factor*pc + a3*factor*pc*pc;
        return aAlphaBeta;
    }


    /*!
     * \brief the derivative of specific interfacial area function w.r.t. capillary pressure
     *
     * \param  params parameter container for the coefficients of the surface
     * \param  Sw Effective saturation of the wetting phase
     * \param  pc Capillary pressure in \f$\mathrm{[Pa]}\f$
     */
    static Scalar dawn_dpc (const Params &params, const Scalar Sw, const Scalar pc)
    {
        const Scalar Swr    = params.Swr();
        const Scalar a1     = params.a1();
        const Scalar a2     = params.a2();
        const Scalar a3     = params.a3();
        const Scalar value =  a2*(Swr-Sw)*(1-Sw) + a3*(Swr-Sw)*(1-Sw)*pc;
        return value;
    }

    /*!
     * \brief the derivative of specific interfacial area function w.r.t. saturation
     *
     * \param  params parameter container for the coefficients of the surface
     * \param  Sw Effective saturation of the wetting phase
     * \param  pc Capillary pressure in \f$\mathrm{[Pa]}\f$
     */
    static Scalar dawn_dsw (const Params &params, const Scalar Sw, const Scalar pc)
    {
        const Scalar Swr                = params.Swr();
        const Scalar a1     = params.a1();
        const Scalar a2     = params.a2();
        const Scalar a3     = params.a3();
        const Scalar derivativeFactor   = (Sw-1.)+(Sw-Swr);
        const Scalar value = a1 * derivativeFactor + a2 * derivativeFactor * pc + a3 * derivativeFactor * pc*pc ;
        return value;
    }

};
} // namespace Dumux

#endif
