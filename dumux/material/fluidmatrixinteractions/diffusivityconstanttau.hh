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
 * \brief   Relation for the saturation-dependent effective diffusion coefficient
 */
#ifndef DIFFUSIVITY_CONSTANT_TAU_HH
#define DIFFUSIVITY_CONSTANT_TAU_HH

#include <dumux/common/parameters.hh>
#include <dumux/common/basicproperties.hh>

namespace Dumux
{
/*!
 * \ingroup fluidmatrixinteractionslaws
 *
 * \brief Relation for the saturation-dependent effective diffusion coefficient
 *
 * The material law is:
 * \f[
 *  D_\text{eff,pm} = \phi * S_w * \tau * D
 * \f]
 *
 * with a constant tau.
 */
template<class TypeTag, class Scalar>
class DiffusivityConstantTau
{
public:
    /*!
     * \brief Returns the effective diffusion coefficient \f$\mathrm{[m^2/s]}\f$ based
     *        on a constant tortuosity value
     *
     * \param porosity The porosity
     * \param saturation The saturation of the phase
     * \param diffCoeff The diffusion coefficient of the phase in \f$\mathrm{[m^2/s]}\f$
     */
    static Scalar effectiveDiffusivity(const Scalar porosity,
                                       const Scalar saturation,
                                       const Scalar diffCoeff)

    {
        Scalar tau = GET_RUNTIME_PARAM(TypeTag, Scalar, tau);

        return porosity * saturation * tau * diffCoeff;
    }
};
}
#endif
