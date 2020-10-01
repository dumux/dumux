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
 * \brief This material law takes a material law (interfacial area surface)
 *        defined for effective saturations and converts it to a material law defined on
 *        absolute saturations.
 */
#ifndef DUMUX_EFF_TO_ABS_LAW_IA_HH
#define DUMUX_EFF_TO_ABS_LAW_IA_HH

#warning "This header is deprecated. Removal after 3.3. Use new material laws."

#include <dumux/material/fluidmatrixinteractions/2p/efftoabslaw.hh>
#include "efftoabslawiaparams.hh"

namespace Dumux {

/*!
 * \ingroup Fluidmatrixinteractions
 * \brief This material law takes a material law (interfacial area surface) defined for effective
 *        saturations and converts it to a material (interfacial area surface) law defined on absolute
 *        saturations.
 *
 *        The idea: "material laws" (like VanGenuchten or BrooksCorey) are defined for effective saturations.
 *        The numeric calculations however are performed with absolute saturations. The EffToAbsLaw class gets
 *        the "material laws" actually used as well as the corresponding parameter container as template arguments.
 *
 *        Subsequently, the desired function (pc, sw... ) of the actually used "material laws" are called but with the
 *        saturations already converted from absolute to effective.
 *
 *        This approach makes sure that in the "material laws" only effective saturations are considered, which makes sense,
 *        as these laws only deal with effective saturations. This also allows for changing the calculation of the effective
 *        saturations easily, as this is subject of discussion may be problem specific.
 *
 *        Additionally, handing over effective saturations to the "material laws" in stead of them calculating effective
 *        saturations prevents accidently "converting twice".
 *
 *        This boils down to:
 *        - the actual material laws (linear, VanGenuchten...) do not need to deal with any kind of conversion
 *        - the definition of the material law in the spatial parameters is not really intuitive, but using it is:
 *          Hand in values, get back values, do not deal with conversion.
 */
template <class EffLawIAT,
          class MaterialAbsParamsT ,
          class InterfacialAreaAbsParamsT = typename EffLawIAT::Params>
class EffToAbsLawIA
{
    using EffLawIA = EffLawIAT;
    using MaterialParams = MaterialAbsParamsT;

public:
    using Params = InterfacialAreaAbsParamsT;
    using Scalar = typename MaterialParams::Scalar;

    /*!
     * \brief The interfacial area relation
     * \param sw Absolute saturation of the wetting phase \f$\mathrm{{S}_w}\f$.
     * \param iaParams parameter container for the interfacial area
     * \param pc Capillary pressure in \f$\mathrm{[Pa]}\f$
     * \param params parameter container for the saturation/materials
     */
    static Scalar interfacialArea(const Params & iaParams,
                                  const MaterialParams & params,
                                  const Scalar sw,
                                  const Scalar pc)
    {
        return EffLawIA::interfacialArea(iaParams,
                                         swToSwe(params, sw),
                                         pc);
    }

protected:
    /*!
     * \brief Convert an absolute wetting saturation to an effective one.
     *
     * \param sw Absolute saturation of the wetting phase \f$\mathrm{{S}_w}\f$.
     * \param params A container object that is populated with the appropriate coefficients for the respective law.
     *                  Therefore, in the (problem specific) spatialParameters first, the material law is chosen, and then the params container
     *                  is constructed accordingly. Afterwards the values are set there, too.
     * \return          Effective saturation of the wetting phase.
     */
    static Scalar swToSwe(const MaterialParams & params, Scalar sw)
    {
        return (sw - params.swr())/(1. - params.swr() - params.snr());
    }
};
} // end namespace Dumux

#endif
