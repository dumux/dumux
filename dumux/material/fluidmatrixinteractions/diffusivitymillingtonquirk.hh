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
 * \brief   Relation for the saturation-dependent effective diffusion coefficient
 */
#ifndef DUMUX_MATERIAL_DIFFUSIVITY_MILLINGTON_QUIRK_HH
#define DUMUX_MATERIAL_DIFFUSIVITY_MILLINGTON_QUIRK_HH

#include <cmath>

namespace Dumux {

/*!
 * \ingroup Fluidmatrixinteractions
 * \brief Relation for the saturation-dependent effective diffusion coefficient
 *
 * The material law is:
 * \f[
 *  D_\text{eff,pm} = \phi * S_w * \tau * D
 * \f]
 *
 * with
 * \f[
 *  \tau = \frac{1}{\phi^2} * \left(\phi S_w\right)^{7/3}
 * \f]
 *
 * after Millington and Quirk 1961: <i>Permeability of porous solids</i> \cite millington1961
 * and Helmig 1997: <i>Multiphase Flow and Transport Processes in the Subsurface: A Contribution
 * to the Modeling of Hydrosystems</i>, page 129 \cite helmig1997
 */
template<class Scalar>
class DiffusivityMillingtonQuirk
{
public:
    /*!
     * \brief Returns the effective diffusion coefficient \f$\mathrm{[m^2/s]}\f$ after Millington Quirk.
     *
     * \param porosity The porosity
     * \param saturation The saturation of the phase
     * \param diffCoeff The diffusion coefficient of the phase \f$\mathrm{[m^2/s]}\f$
     */
    [[deprecated("Signature deprecated. Use signature with volume variables!")]]
    static Scalar effectiveDiffusivity(const Scalar porosity,
                                       const Scalar saturation,
                                       const Scalar diffCoeff)
    {
        // instead of D_eff,pm = phi * Sw * 1/phi^2 * (phi * Sw)^(7/3) * D
        // we calculate the more efficient
        // D_eff,pm = phi * Sw^3 * cubicroot(phi * Sw) * D

        using std::cbrt;
        return porosity * (saturation * saturation * saturation)
               * cbrt(porosity * saturation) * diffCoeff;
    }

    /*!
     * \brief Returns the effective diffusion coefficient \f$\mathrm{[m^2/s]}\f$ after Millington Quirk.
     *
     * \param volVars The Volume Variables
     * \param phaseIdx the index of the phase
     * \param compIdx the component index
     */
    template<class VolumeVariables>
    static Scalar effectiveDiffusivity(const VolumeVariables& volVars,
                                       const Scalar diffCoeff,
                                       const int phaseIdx)
    {
        // instead of D_eff,pm = phi * Sw * 1/phi^2 * (phi * Sw)^(7/3) * D
        // we calculate the more efficient
        // D_eff,pm = phi * Sw^3 * cubicroot(phi * Sw) * D

        using std::cbrt;
        return volVars.porosity() * (volVars.saturation(phaseIdx) * volVars.saturation(phaseIdx) * volVars.saturation(phaseIdx))
               * cbrt(volVars.porosity() * volVars.saturation(phaseIdx)) * diffCoeff;
    }
};
}
#endif
