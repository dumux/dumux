// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup EffectiveHeatConductivity
 * \brief Relation for the saturation-dependent effective thermal conductivity
 */

#ifndef DUMUX_MATERIAL_FLUIDMATRIX_THERMALCONDUCTIVITY_SIMPLE_FLUID_LUMPING_HH
#define DUMUX_MATERIAL_FLUIDMATRIX_THERMALCONDUCTIVITY_SIMPLE_FLUID_LUMPING_HH

#include <assert.h>
#include <algorithm>
#warning "This header is deprecated and will be removed after 3.9. Use ThermalConductivityAverage"

namespace Dumux {

/*!
 * \ingroup EffectiveHeatConductivity
 * \brief Relation for the saturation-dependent effective thermal conductivity
 * \deprecated This does the same as `ThermalConductivityAverage` but only works for two fluid phases
 */
template<class Scalar>
class [[deprecated("Use ThermalConductivityAverage. Will be removed after 3.9.")]] ThermalConductivitySimpleFluidLumping
{
public:
    /*!
     * \brief Effective thermal conductivity in \f$\mathrm{W/(m K)}\f$ for two phases
     * \param volVars volume variables
     * \return Effective thermal conductivity in \f$\mathrm{W/(m K)}\f$ for two phases
     */
    template<class VolumeVariables>
    static Scalar effectiveThermalConductivity(const VolumeVariables& volVars)
    {
        using FluidSystem = typename VolumeVariables::FluidSystem;
        const Scalar sw = volVars.saturation(FluidSystem::phase0Idx);
        const Scalar lambdaW = volVars.fluidThermalConductivity(FluidSystem::phase0Idx);
        const Scalar lambdaN = volVars.fluidThermalConductivity(FluidSystem::phase1Idx);
        const Scalar lambdaSolid = volVars.solidThermalConductivity();
        const Scalar porosity = volVars.porosity();

        return effectiveThermalConductivity_(sw, lambdaW, lambdaN, lambdaSolid, porosity);
    }

private:
    /*!
     * \brief Effective thermal conductivity in \f$\mathrm{W/(m K)}\f$ for two phases
     *
     * \param sw The saturation of the wetting phase
     * \param lambdaW The thermal conductivity of the wetting phase in \f$\mathrm{W/(m K)}\f$
     * \param lambdaN The thermal conductivity of the nonwetting phase in \f$\mathrm{W/(m K)}\f$
     * \param lambdaSolid The thermal conductivity of the solid phase in \f$\mathrm{W/(m K)}\f$
     * \param porosity The porosity
     *
     * \return Effective thermal conductivity in \f$\mathrm{W/(m K)}\f$ for two phases
     */
    static Scalar effectiveThermalConductivity_(const Scalar sw,
                                                const Scalar lambdaW,
                                                const Scalar lambdaN,
                                                const Scalar lambdaSolid,
                                                const Scalar porosity)
    {
        // Franz Lindner / Shi & Wang 2011
        using std::max;
        const Scalar satW = max<Scalar>(0.0, sw);
        return porosity * ( (1. - satW) * lambdaN + satW * lambdaW ) + (1.0 - porosity) * lambdaSolid ; ; // arithmetic
    }
};

} // end namespace Dumux

#endif
