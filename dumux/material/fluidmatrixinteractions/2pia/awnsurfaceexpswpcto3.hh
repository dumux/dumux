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
 * This function is of third order in pc.
 */
#ifndef AWN_SURFACE_EXP_SW_PC_TO_3
#define AWN_SURFACE_EXP_SW_PC_TO_3

#include "awnsurfaceexpswpcto3params.hh"

#include <dune/common/exceptions.hh>

#include <algorithm>
#include <cmath>
#include <assert.h>

#warning "This header is deprecated. Removal after 3.3. Use new material laws."

namespace Dumux {

/*!
 * \ingroup Fluidmatrixinteractions
 * \brief Implementation of a exponential function relating
 * specific interfacial  area to wetting phase saturation and capillary pressure.
 */
template <class ScalarT, class ParamsT =AwnSurfaceExpSwPcTo3Params<ScalarT> >
class [[deprecated("Use new material laws and FluidMatrix::InterfacialAreaExponentialCubic instead!")]] AwnSurfaceExpSwPcTo3
{
public:
    using Params = ParamsT;
    using Scalar = typename Params::Scalar;

    /*!
     * \brief The awn surface
     *
     * the suggested (as estimated from pore network models) interfacial area surface:
     * \f$\mathrm{
        a_{wn} = a_1  e^{a_2 * S_w }  + a_3 * p_c^3 ;
     }\f$
     * \param  params parameter container for the coefficients of the surface
     * \param  Sw Effective saturation of the wetting phase
     * \param  pc Capillary pressure in \f$\mathrm{[Pa]}\f$
     */
    static Scalar interfacialArea(const Params & params, const Scalar Sw, const Scalar pc)
    {
        // TODO think about awn surface for relative saturation
        const Scalar a1 = params.a1();
        const Scalar a2 = params.a2();
        const Scalar a3 = params.a3();

        using std::exp;
        const Scalar aAlphaBeta = a1 * exp( a2 * Sw) + a3 * pc * pc * pc ;
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
        DUNE_THROW(Dune::NotImplemented, __FILE__ << "  dawndpc()");
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
        DUNE_THROW(Dune::NotImplemented, __FILE__ << "  dawndSw()");
    }

};
} // end namespace Dumux

#endif
