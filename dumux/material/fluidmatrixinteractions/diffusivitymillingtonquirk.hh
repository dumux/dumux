// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup Fluidmatrixinteractions
 * \brief   Relation for the saturation-dependent effective diffusion coefficient
 */
#ifndef DUMUX_MATERIAL_DIFFUSIVITY_MILLINGTON_QUIRK_HH
#define DUMUX_MATERIAL_DIFFUSIVITY_MILLINGTON_QUIRK_HH

#include <cmath>
#include <algorithm>

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
     * \param volVars The Volume Variables
     * \param phaseIdx the index of the phase
     * \param compIdxI the component index i
     * \param compIdxJ the component index j
     */
    template<class VolumeVariables>
    static Scalar effectiveDiffusionCoefficient(const VolumeVariables& volVars,
                                                const int phaseIdx,
                                                const int compIdxI,
                                                const int compIdxJ)
    {
        // instead of D_eff,pm = phi * Sw * 1/phi^2 * (phi * Sw)^(7/3) * D
        // we calculate the more efficient
        // D_eff,pm = phi * Sw^3 * cubicroot(phi * Sw) * D

        using std::cbrt;
        using std::max;
        const Scalar diffCoeff = volVars.diffusionCoefficient(phaseIdx, compIdxI, compIdxJ);
        const Scalar porosity = volVars.porosity();
        const Scalar sat = max<Scalar>(volVars.saturation(phaseIdx), 0.0);
        return porosity * (sat*sat*sat) * cbrt(porosity * sat) * diffCoeff;
    }

};
}
#endif
