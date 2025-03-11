// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup EffectiveHeatConductivity
 * \brief Effective thermal conductivity after Somerton
 */

#ifndef DUMUX_MATERIAL_FLUIDMATRIX_THERMALCONDUCTIVITY_SOMERTON_THREE_P_HH
#define DUMUX_MATERIAL_FLUIDMATRIX_THERMALCONDUCTIVITY_SOMERTON_THREE_P_HH

#include <algorithm>
#include <cmath>

namespace Dumux {

/*!
 * \addtogroup EffectiveHeatConductivity
 * \copydetails Dumux::ThermalConductivitySomertonThreeP
 */

/*!
 * \ingroup EffectiveHeatConductivity
 * \brief Effective thermal conductivity after Somerton
 *
 * ### Somerton (three fluid phases)
 *
 * The Somerton method \cite somerton1974 computes the thermal conductivity of dry and the wet soil material.
 * It is extended here to a three phase system of water (w), NAPL (n) and gas (g).
 * It uses a root function of the water and NAPL saturation to compute the
 * effective thermal conductivity for a three-phase fluidsystem. The individual thermal
 * conductivities are calculated as geometric mean of the thermal conductivity of the porous
 * material and of the respective fluid phase.
 *
 * The effective thermal conductivity of `ThermalConductivitySomertonThreeP` is given by
 * \f[
 * \lambda_\text{eff} = \lambda_\text{g,eff} + \sqrt{S_\text{w}} \left(\lambda_\text{w,eff} - \lambda_\text{g,eff}\right) +
 * \sqrt{S_\text{n}} \left(\lambda_\text{n,eff} - \lambda_\text{g,eff}\right)
 * \f]
 *
 * with \f$ S_\text{w} \f$ the water saturation,
 * \f$ S_\text{n} \f$ the NAPL saturation, the effective phase saturations given by
 * \f$ \lambda_{\alpha,\text{eff}} = (\lambda_\text{s})^{\left(1-\phi\right)} (\lambda_\alpha)^\phi, \alpha \in \{\text{w,n,g}\}\f$
 * (geometric mean) and \f$ \lambda_\text{s} \f$ is the thermal conductivity of the solid phase.
 */
template<class Scalar>
class ThermalConductivitySomertonThreeP
{
public:
    /*!
     * \brief Effective thermal conductivity in \f$\mathrm{W/(m K)}\f$ for three phases
     * \param volVars volume variables
     * \return Effective thermal conductivity in \f$\mathrm{W/(m K)}\f$ for three phases
     */
    template<class VolumeVariables>
    static Scalar effectiveThermalConductivity(const VolumeVariables& volVars)
    {
        using FluidSystem = typename VolumeVariables::FluidSystem;

        const Scalar sw = volVars.saturation(FluidSystem::wPhaseIdx);
        const Scalar sn = volVars.saturation(FluidSystem::nPhaseIdx);
        const Scalar lambdaW = volVars.fluidThermalConductivity(FluidSystem::wPhaseIdx);
        const Scalar lambdaN = volVars.fluidThermalConductivity(FluidSystem::nPhaseIdx);
        const Scalar lambdaG = volVars.fluidThermalConductivity(FluidSystem::gPhaseIdx);
        const Scalar lambdaSolid = volVars.solidThermalConductivity();
        const Scalar porosity = volVars.porosity();

        return effectiveThermalConductivity(sw, sn, lambdaW, lambdaN, lambdaG, lambdaSolid, porosity);
    }

    /*!
     * \brief Effective thermal conductivity in \f$\mathrm{W/(m K)}\f$ for three phases
     *
     * \param sw The saturation of the wetting phase
     * \param sn The saturation of the nonwetting phase
     * \param lambdaW The thermal conductivity of the water phase in \f$\mathrm{W/(m K)}\f$
     * \param lambdaN The thermal conductivity of the NAPL phase in \f$\mathrm{W/(m K)}\f$
     * \param lambdaG The thermal conductivity of the gas phase in \f$\mathrm{W/(m K)}\f$
     * \param lambdaSolid The thermal conductivity of the solid phase in \f$\mathrm{W/(m K)}\f$
     * \param porosity The porosity
     *
     * \return Effective thermal conductivity in \f$\mathrm{W/(m K)}\f$ for three phases
     */
    static Scalar effectiveThermalConductivity(const Scalar sw,
                                               const Scalar sn,
                                               const Scalar lambdaW,
                                               const Scalar lambdaN,
                                               const Scalar lambdaG,
                                               const Scalar lambdaSolid,
                                               const Scalar porosity)
    {
        using std::max;
        using std::pow;
        using std::sqrt;
        const Scalar satW = max<Scalar>(0.0, sw);
        const Scalar satN = max<Scalar>(0.0, sn);

        // porosity weighted geometric mean
        const Scalar lSw = pow(lambdaSolid, (1.0 - porosity)) * pow(lambdaW, porosity);
        const Scalar lSn = pow(lambdaSolid, (1.0 - porosity)) * pow(lambdaN, porosity);
        const Scalar lSg = pow(lambdaSolid, (1.0 - porosity)) * pow(lambdaG, porosity);
        const Scalar lambdaEff = lSg + sqrt(satW) * (lSw - lSg) + sqrt(satN) * (lSn -lSg);

        return lambdaEff;

    }
};


} // end namespace Dumux

#endif
