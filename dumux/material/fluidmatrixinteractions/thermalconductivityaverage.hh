// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup EffectiveHeatConductivity
 * \brief Effective thermal conductivity based on weighted arithmetic average
 */

#ifndef DUMUX_MATERIAL_THERMALCONDUCTIVITY_AVERAGED_HH
#define DUMUX_MATERIAL_THERMALCONDUCTIVITY_AVERAGED_HH

#include <algorithm>

namespace Dumux {

/*!
 * \addtogroup EffectiveHeatConductivity
 * \copydetails Dumux::ThermalConductivityAverage
 */

/*!
 * \ingroup EffectiveHeatConductivity
 * \brief Effective thermal conductivity based on weighted arithmetic average
 *
 * ### Average (multiple fluid phases, one solid phase)
 *
 * The effective thermal conductivity of `ThermalConductivityAverage`
 * is calculated as a weighted arithmetic average of the thermal
 * conductivities of the solid and the fluid phases. The weights are determined by the volume
 * fraction the phase occupies. Denoting the volume fractions by \f$ n_\alpha \f$, we have
 * \f[
 * \lambda_\text{eff} = \sum_\alpha \lambda_\alpha n_\alpha / \sum_\alpha n_\alpha,
 * \f]
 * summing over both fluid and solid phases. With the porosity \f$ \phi \f$ as
 * the sum of all fluid volume fractions, we can equivalently write
 * \f[
 * \lambda_\text{eff} = \lambda_\text{s} (1-\phi) + \lambda_\text{f} \phi,
 * \f]
 * where \f$ \lambda_\text{s} \f$ is the thermal conductivity of the solid phase,
 * and the effective thermal conductivity of the liquid phases is computed as
 * an arithmetic average weighted with the fluid saturations.
 */
template<class Scalar>
class ThermalConductivityAverage
{
public:
    /*!
     * \brief Effective thermal conductivity in \f$\mathrm{W/(m K)}\f$
     * \param volVars volume variables
     * \return Effective thermal conductivity in \f$\mathrm{W/(m K)}\f$
     */
    template<class VolumeVariables>
    static Scalar effectiveThermalConductivity(const VolumeVariables& volVars)
    {
        constexpr int numFluidPhases = VolumeVariables::numFluidPhases();

        // Get the thermal conductivities and the porosity from the volume variables
        Scalar lambdaFluid = 0.0;
        for (int phaseIdx = 0; phaseIdx < numFluidPhases; ++phaseIdx)
            lambdaFluid += volVars.fluidThermalConductivity(phaseIdx)*volVars.saturation(phaseIdx);

        const Scalar lambdaSolid = volVars.solidThermalConductivity();
        const Scalar porosity = volVars.porosity();

        return lambdaSolid*(1-porosity) + lambdaFluid*porosity;
    }
};

} // end namespace Dumux

#endif
