// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup Fluidmatrixinteractions
 * \brief Reation for a simple effective thermal conductivity
 */
#ifndef DUMUX_MATERIAL_THERMALCONDUCTIVITY_AVERAGE_HH
#define DUMUX_MATERIAL_THERMALCONDUCTIVITY_AVERAGE_HH

#include <algorithm>

namespace Dumux {

/*!
 * \ingroup Fluidmatrixinteractions
 * \brief Relation for a simple effective thermal conductivity
 */
template<class Scalar>
class ThermalConductivityAverage
{
public:
    /*!
     * \brief Relation for a simple effective thermal conductivity \f$\mathrm{[W/(m K)]}\f$
     *
     * \param volVars volume variables
     * \return Effective thermal conductivity \f$\mathrm{[W/(m K)]}\f$
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
