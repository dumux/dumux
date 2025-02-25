// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
#ifndef DUMUX_MATERIAL_FLUIDMATRIX_THERMALCONDUCTIVITY_SOMERTON_TWO_P_HH
#define DUMUX_MATERIAL_FLUIDMATRIX_THERMALCONDUCTIVITY_SOMERTON_TWO_P_HH

#include <algorithm>
#include <cmath>

namespace Dumux {

/*!
 * \addtogroup EffectiveHeatConductivity
 * \copydetails Dumux::ThermalConductivitySomertonTwoP
 */

/*!
 * \ingroup EffectiveHeatConductivity
 * \brief Effective thermal conductivity after Somerton
 *
 * ### Somerton (two fluid phases)
 *
 * The Somerton method \cite somerton1974 computes the thermal conductivity of dry and the wet soil material.
 * It uses a root function of the water saturation to compute the
 * effective thermal conductivity for a two-phase fluidsystem. The individual thermal
 * conductivities are calculated as geometric mean of the thermal conductivity of the porous
 * material and of the respective fluid phase.
 *
 * The effective thermal conductivity of `ThermalConductivitySomertonTwoP` is given by
 * \f[
 * \lambda_\text{eff} = \lambda_\text{g,eff} + \sqrt{S_\text{w}} \left(\lambda_\text{w,eff} - \lambda_\text{g,eff}\right)
 * \f]
 *
 * with \f$ S_\text{w} \f$ the water saturation,
 * \f$ S_\text{n} \f$ the NAPL saturation, the effective phase saturations given by
 * \f$ \lambda_{\alpha,\text{eff}} = (\lambda_\text{s})^{\left(1-\phi\right)} (\lambda_\alpha)^\phi, \alpha \in \lbrace\text{w,n,g}\rbrace \f$
 * (geometric mean) and \f$ \lambda_\text{s} \f$ is the thermal conductivity of the solid phase.
 * The effective conductivity \f$ \lambda_\text{g,eff} \f$ corresponds to dry conditions, whereas the
 * effective conductivity \f$ \lambda_\text{g,eff} \f$ corresponds to wet conditions.
 */
template<class Scalar>
class ThermalConductivitySomertonTwoP
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
        static_assert(FluidSystem::numPhases == 2, "ThermalConductivitySomertonTwoP only works for two-phase fluid systems!");
        static_assert((FluidSystem::isGas(0) && !FluidSystem::isGas(1)) || (!FluidSystem::isGas(0) && FluidSystem::isGas(1)),
                     "ThermalConductivitySomertonTwoP only works if one phase is gaseous and one is liquid!");

        constexpr int liquidPhaseIdx = FluidSystem::isGas(0) ? 1 : 0;
        constexpr int gasPhaseIdx = FluidSystem::isGas(0) ? 0 : 1;

        const Scalar satLiquid = volVars.saturation(liquidPhaseIdx);
        const Scalar lambdaLiquid = volVars.fluidThermalConductivity(liquidPhaseIdx);
        const Scalar lambdaGas = volVars.fluidThermalConductivity(gasPhaseIdx);
        const Scalar lambdaSolid = volVars.solidThermalConductivity();
        const Scalar porosity = volVars.porosity();

        return effectiveThermalConductivity_(satLiquid, lambdaLiquid, lambdaGas, lambdaSolid, porosity);
    }

private:
    /*!
     * \brief Effective thermal conductivity in \f$\mathrm{W/(m K)}\f$ for two phases
     *
     * \param satLiquid The saturation of the liquid phase
     * \param lambdaLiquid The thermal conductivity of the liquid phase in \f$\mathrm{W/(m K)}\f$
     * \param lambdaGas The thermal conductivity of the gas phase in \f$\mathrm{W/(m K)}\f$
     * \param lambdaSolid The thermal conductivity of the solid phase in \f$\mathrm{W/(m K)}\f$
     * \param porosity The porosity
     *
     * \brief Effective thermal conductivity in \f$\mathrm{W/(m K)}\f$ for two phases
     */
    static Scalar effectiveThermalConductivity_(const Scalar satLiquid,
                                                const Scalar lambdaLiquid,
                                                const Scalar lambdaGas,
                                                const Scalar lambdaSolid,
                                                const Scalar porosity)
    {
        using std::max;
        using std::pow;
        using std::sqrt;
        const Scalar satLiquidPhysical = max<Scalar>(0.0, satLiquid);
        // geometric mean, using ls^(1-p)*l^p = ls*(l/ls)^p
        const Scalar lambdaSaturated = lambdaSolid * pow(lambdaLiquid / lambdaSolid, porosity);
        const Scalar lambdaDry = lambdaSolid * pow(lambdaGas / lambdaSolid, porosity);

        return lambdaDry + sqrt(satLiquidPhysical) * (lambdaSaturated - lambdaDry);
    }
};

#ifndef DOXYGEN
template<class Scalar>
using ThermalConductivitySomerton [[deprecated("Use ThermalConductivitySomertonTwoP. Will be removed after 3.9.")]] = ThermalConductivitySomertonTwoP<Scalar>;
#endif

} // end namespace Dumux

#endif
