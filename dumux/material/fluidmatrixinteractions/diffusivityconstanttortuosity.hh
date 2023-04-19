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
#ifndef DUMUX_MATERIAL_DIFFUSIVITY_CONSTANT_TORTUOSITY_HH
#define DUMUX_MATERIAL_DIFFUSIVITY_CONSTANT_TORTUOSITY_HH

#include <dumux/common/parameters.hh>

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
 * with a constant tau.
 *
 * The default value is 0.5, empirically obtained in Carman 1937:
 * <i>Fluid flow through granular beds</i> \cite carman1937
 * Additionally, Bear 1972 \cite bear1972 mentions values 0.4 and in the
 * range of 0.56 to 0.8.
 */
template<class Scalar>
class DiffusivityConstantTortuosity
{
public:
    /*!
     * \brief Returns the effective diffusion coefficient \f$\mathrm{[m^2/s]}\f$ based
     *        on a constant tortuosity value
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
        static const Scalar tau = getParam<Scalar>("SpatialParams.Tortuosity", 0.5);
        const Scalar diffCoeff = volVars.diffusionCoefficient(phaseIdx, compIdxI, compIdxJ);
        return volVars.porosity() * volVars.saturation(phaseIdx) * tau * diffCoeff;
    }

};
} // end namespace Dumux

#endif
