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
 * This function is exponential.
 */
#ifndef AWN_SURFACE_EXP_FCT_HH
#define AWN_SURFACE_EXP_FCT_HH

#include "awnsurfaceexpfctparams.hh"

#include <dune/common/exceptions.hh>

#include <algorithm>
#include <cmath>
#include <assert.h>

#warning "This header is deprecated. Removal after 3.3. Use new material laws."

namespace Dumux {

/*!
 * \ingroup Fluidmatrixinteractions
 * \brief Implementation of the exponential function relating
 *        specific interfacial area to wetting phase saturation and c
 *        apillary pressure as suggested by Nuske(2009) (Diploma thesis) \cite nuske2009 .
 */
template <class ParamsT>
class [[deprecated("Use new material laws and FluidMatrix::InterfacialAreaExponential instead!")]] AwnSurfaceExpFct
{
public:
    using Params = ParamsT;
    using Scalar = typename Params::Scalar;

    /*!
     * \brief The interfacial area surface
     *
     * the suggested (as estimated from pore network models) awn surface:
     * \f$\mathrm{
        a_{wn} = a_1 * (S_{wr}-S_w) .* (1-S_w) + a_2 * (S_{wr}-S_w) * (1-S_w) * \exp( a_3 * p_c) ;
     }\f$
     * \param  params parameter container for the coefficients of the surface
     * \param  Sw Effective saturation of the wetting phase
     * \param  pc Capillary pressure in \f$\mathrm{[Pa]}\f$
     */
    static Scalar interfacialArea(const Params & params, const Scalar Sw, const Scalar pc)
    {
        const Scalar a1 = params.a1();
        const Scalar a2 = params.a2();
        const Scalar a3 = params.a3();
        const Scalar Swr = params.Swr();
        using std::exp;
        const Scalar aAlphaBeta = a1 * (Swr-Sw) * (1-Sw) + a2 * (Swr-Sw) * (1-Sw) * exp( a3 * pc) ;
        return aAlphaBeta;
    }

    /*! \brief the derivative of specific interfacial area function w.r.t. capillary pressure
     *
     * \param  params parameter container for the coefficients of the surface
     * \param  Sw Effective saturation of the wetting phase
     * \param  pc Capillary pressure in \f$\mathrm{[Pa]}\f$
     */
    static Scalar dawn_dpc (const Params & params, const Scalar Sw, const Scalar pc)
    {
        const Scalar a2 = params.a2();
        const Scalar a3 = params.a3();
        const Scalar Swr = params.Swr();
        using std::exp;
        const Scalar value = a2 * a3 * (Swr-Sw) * (1-Sw) * exp(a3*pc);
        return value;
    }

    /*! \brief the derivative of specific interfacial area function w.r.t. saturation
     *
     * \param  params parameter container for the coefficients of the surface
     * \param  Sw Effective saturation of the wetting phase
     * \param  pc Capillary pressure in \f$\mathrm{[Pa]}\f$
     */
    static Scalar dawn_dsw (const Params & params, const Scalar Sw, const Scalar pc)
    {
        Scalar value;
        Scalar a1 = params.a1();
        Scalar a2 = params.a2();
        Scalar a3 = params.a3();
        Scalar Swr = params.Swr();
        using std::exp;
        value = - a1 *( 1+Swr-2*Sw ) - a2 * exp(a3*pc) * ( 1+Swr-2*Sw );
        return value;
    }

};
} // namespace Dumux

#endif
