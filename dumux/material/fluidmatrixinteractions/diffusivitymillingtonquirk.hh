// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
#ifndef DUMUX_MATERIAL_DIFFUSIVITY_MILLINGTON_QUIRK_HH
#define DUMUX_MATERIAL_DIFFUSIVITY_MILLINGTON_QUIRK_HH

#include <cmath>
#include <algorithm>

namespace Dumux {

/*!
* \addtogroup EffectiveDiffusivity
* \copydetails Dumux::DiffusivityMillingtonQuirk
*/

/*!
 * \ingroup EffectiveDiffusivity
 * \brief Relation for the effective diffusion coefficient after Millington and Quirk
 *
 * ### Millington Quirk
 *
 * For `DiffusivityMillingtonQuirk`,
 * the tortuosity coefficient is estimated after Millington and Quirk (1961) \cite millington1961
 * by
 * \f[
 *  \tau = \frac{1}{\phi^2} \left(\phi S_\alpha\right)^{7/3}.
 * \f]
 *
 * See also Helmig (1997) \cite helmig1997 [page 129].
 */
template<class Scalar>
class DiffusivityMillingtonQuirk
{
public:
    /*!
     * \brief Returns the effective diffusion coefficient (\f$\mathrm{m^2/s}\f$)
     *
     * Returns the effective diffusion coefficient (\f$\mathrm{m^2/s}\f$)
     * of component \f$ \kappa \f$ (index `compIdxI`) in phase \f$ \alpha \f$ computed as
     * \f$ D^\kappa_{\text{eff},\alpha} = \phi S_\alpha \tau D^\kappa_\alpha \f$.
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
        // instead of D_eff = phi S tau D = phi S 1/phi^2 (phi S)^(7/3) D
        // we implement more efficiently D_eff = phi S^3 cubicroot(phi S) D
        using std::cbrt;
        using std::max;
        const Scalar diffCoeff = volVars.diffusionCoefficient(phaseIdx, compIdxI, compIdxJ);
        const Scalar porosity = volVars.porosity();
        const Scalar sat = max<Scalar>(volVars.saturation(phaseIdx), 0.0);
        return porosity * (sat*sat*sat) * cbrt(porosity * sat) * diffCoeff;
    }
};

} // end namespace Dumux

#endif
