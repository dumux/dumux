// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
/*!
 * \file
 * \ingroup TwoPTwoCTests
 * \brief Effective gas diffusion coefficient of the heat-pipe reference solution of Huang et al. (2015).
 */
#ifndef DUMUX_HEATPIPE_DIFFUSIVITY_HH
#define DUMUX_HEATPIPE_DIFFUSIVITY_HH

namespace Dumux {

/*!
 * \ingroup TwoPTwoCTests
 * \brief Constant macroscopic diffusion coefficient of the reference code of Huang et al. (2015).
 *
 * The binary diffusion coefficient of the fluid system is used as the effective
 * coefficient of the porous medium, without scaling by porosity, saturation or tortuosity.
 */
template<class Scalar>
class HeatPipeDiffusivity
{
public:
    template<class VolumeVariables>
    static Scalar effectiveDiffusionCoefficient(const VolumeVariables& volVars,
                                                const int phaseIdx,
                                                const int compIdxI,
                                                const int compIdxJ)
    { return volVars.diffusionCoefficient(phaseIdx, compIdxI, compIdxJ); }
};

} // end namespace Dumux

#endif
