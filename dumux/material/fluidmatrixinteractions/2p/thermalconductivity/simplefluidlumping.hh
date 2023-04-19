// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup Fluidmatrixinteractions
 * \brief   Relation for the saturation-dependent effective thermal conductivity
 */
#ifndef DUMUX_MATERIAL_FLUIDMATRIX_THERMALCONDUCTIVITY_SIMPLE_FLUID_LUMPING_HH
#define DUMUX_MATERIAL_FLUIDMATRIX_THERMALCONDUCTIVITY_SIMPLE_FLUID_LUMPING_HH

#include <assert.h>
#include <algorithm>

namespace Dumux {

/*!
 * \ingroup Fluidmatrixinteractions
 * \brief   Relation for the saturation-dependent effective thermal conductivity
 */
template<class Scalar>
class ThermalConductivitySimpleFluidLumping
{
public:
    /*!
     * \brief Effective thermal conductivity \f$\mathrm{[W/(m K)]}\f$
     *
     * \param volVars volume variables
     * \return effective thermal conductivity \f$\mathrm{[W/(m K)]}\f$
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
     * \brief Returns the effective thermal conductivity \f$\mathrm{[W/(m K)]}\f$.
     *
     * \param sw The saturation of the wetting phase
     * \param lambdaW The thermal conductivity of the wetting phase in \f$\mathrm{[W/(m K)]}\f$
     * \param lambdaN The thermal conductivity of the nonwetting phase in \f$\mathrm{[W/(m K)]}\f$
     * \param lambdaSolid The thermal conductivity of the solid phase in \f$\mathrm{[W/(m K)]}\f$
     * \param porosity The porosity
     *
     * \return Effective thermal conductivity of the fluid phases
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
