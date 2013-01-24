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
 * \brief   Relation for the saturation-dependent effective thermal conductivity
 */
#ifndef THERMALCONDUCTIVITY_SOMERTON_HH
#define THERMALCONDUCTIVITY_SOMERTON_HH

#include <algorithm>

namespace Dumux
{
/*!
 * \ingroup fluidmatrixinteractionslaws
 *
 * \brief Relation for the saturation-dependent effective thermal conductivity
 *
 *  The Somerton method computes the thermal conductivity of dry and the wet soil material
 *  and uses a root function of the wetting saturation to compute the
 *  effective thermal conductivity for a two-phase fluidsystem. It is assumed, that the
 *  non-wetting phase does not contribute to the conduction due to a large contrast
 *  in the fluid conductivities.
 */
template<class Scalar>
class ThermalConductivitySomerton
{
public:
    /*!
     * \brief Returns the effective thermal conductivity \f$[W/m^2]\f$ after Somerton (1974).
     *
     * The material law is:
     * \f[
     l_eff = l_solid + (l_wet - l_solid)
     \f]
     *
     * \param Sw The saturation of the wetting phase
     * \param lambdaW the thermal conductivity of the wetting phase
     * \param lambdaN the thermal conductivity of the non-wetting phase
     * \param lambdaSolid the thermal conductivity of the solid phase
     * \param porosity The porosity
     *
     * \return Effective thermal conductivity \f$[W/m^2]\f$ after Somerton (1974)
     */
    static Scalar effectiveThermalConductivity(const Scalar Sw,
                                               const Scalar lambdaW,
                                               const Scalar lamdaN,
                                               const Scalar lambdaSolid,
                                               const Scalar porosity)
    {
        const Scalar satW = std::max<Scalar>(0.0, Sw);
        const Scalar lSat = std::pow(lambdaSolid, (1.0 - porosity)) * std::pow(lambdaW, porosity);
        const Scalar lDry = std::pow(lambdaSolid, (1.0 - porosity));

        return lDry + std::sqrt(satW) * (lDry - lSat);
    }

};
}
#endif
