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
 * This parametrization uses a maximum value of capillary pressure.
 */
#ifndef AWN_SURFACE_PCMAX_FCT_HH
#define AWN_SURFACE_PCMAX_FCT_HH

#include "awnsurfacepcmaxfctparams.hh"
#include <dune/common/exceptions.hh>

#warning "This header is deprecated. Removal after 3.3. Use new material laws."

namespace Dumux {

/*!
 * \ingroup Fluidmatrixinteractions
 * \brief Implementation of an interfacial area surface.
 *
 * The idea here is to make sure that interfacial area be zero for Sw=1 and pc=pcmax.
 * Because imbibition may have bigger awn than drainage demanding that awn be zero for the whole range of Sw=0 does not fit.
 * However, as we are in primary drainage, demanding that awn=0 for Sw=0 is the same as demanding that awn=0 for pc=pcmax,
 * because Sw=0 cannnot be reached on different ways. We only need to define pcmax ...
 *
 */

template <class ScalarT, class ParamsT =AwnSurfacePcMaxFctParams<ScalarT> >
class [[deprecated("Use new material laws and FluidMatrix::InterfacialAreaPcMax instead!")]] AwnSurfacePcMaxFct
{
public:
    using Params = ParamsT;
    using Scalar = typename Params::Scalar;

    /*!
     * \brief The awn surface
     *
     * the suggested (as estimated from pore network models) awn surface:
     * \f$\mathrm{
     a_{wn} = a_{1} (p_{c { \sf max} } - p_{c}) (1-S_{w}) +  a_{2}  (p_{c {\sf max} } -p_{c})^2 (1-S_{w}) +  a_{3} (p_{c {\sf max} }- p_{c}) (1-S_{w})^2
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
        const Scalar pcMax = params.pcMax() ;

        const Scalar  aAlphaBeta = a1 * (pcMax-pc) * (1.-Sw) +  a2*(pcMax-pc)*(pcMax-pc) * (1.-Sw) +  a3 * (pcMax-pc)*(1-Sw)*(1-Sw);
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
