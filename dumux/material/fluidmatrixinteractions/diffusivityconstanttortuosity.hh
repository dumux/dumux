// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
#ifndef DUMUX_MATERIAL_DIFFUSIVITY_CONSTANT_TORTUOSITY_HH
#define DUMUX_MATERIAL_DIFFUSIVITY_CONSTANT_TORTUOSITY_HH

#include <dumux/common/parameters.hh>

namespace Dumux {

/*!
* \addtogroup EffectiveDiffusivity
* \copydetails Dumux::DiffusivityConstantTortuosity
*/

/*!
 * \ingroup EffectiveDiffusivity
 * \brief Relation for the effective diffusion coefficient with constant tortuosity
 *
 * ### Constant Tortuosity
 *
 * For `DiffusivityConstantTortuosity`, \f$ \tau = \text{const.} \f$,
 * with default value 0.5, empirically obtained by Carman \cite carman1937.
 * The value can be changed at runtime by setting parameter
 * `"SpatialParams.Tortuosity"`. This will change the value of \f$ \tau \f$,
 * and therefore the effective diffusion coefficient wherever
 * the function `effectiveDiffusionCoefficient` is used.
 */
template<class Scalar>
class DiffusivityConstantTortuosity
{
public:
    /*!
     * \brief Returns the effective diffusion coefficient (\f$\mathrm{m^2/s}\f$)
     *
     * Returns the effective diffusion coefficient (\f$\mathrm{m^2/s}\f$)
     * of component \f$ \kappa \f$ (index `compIdxI`) in phase \f$ \alpha \f$ based
     * on a constant tortuosity coefficient: \f$ D^\kappa_{\text{eff},\alpha} = \phi S_\alpha \tau D^\kappa_\alpha \f$.
     *
     * \param volVars The volume variables
     * \param phaseIdx the index of phase \f$ \alpha \f$
     * \param compIdxI the component index i (the component diffusing in phase \f$ \alpha \f$)
     * \param compIdxJ the component index j (the main component of phase \f$ \alpha \f$)
     */
    template<class VolumeVariables>
    static Scalar effectiveDiffusionCoefficient(const VolumeVariables& volVars,
                                                const int phaseIdx,
                                                const int compIdxI,
                                                const int compIdxJ)
    {
        static const Scalar tau = getParam<Scalar>("SpatialParams.Tortuosity", 0.5);
        const Scalar diffCoeff = volVars.diffusionCoefficient(phaseIdx, compIdxI, compIdxJ);
        return volVars.porosity() * volVars.saturation(phaseIdx) * tau * diffCoeff;
    }

};
} // end namespace Dumux

#endif
