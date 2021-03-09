// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   See the file COPYING for full copying permissions.                      *
 *                                                                           *
 *   This program is free software: you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation, either version 3 of the License, or       *
 *   (at your option) any later version.                                     *
 *                                                                           *
 *   This program is distributed in the hope that it will be useful,         *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of          *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the            *
 *   GNU General Public License for more details.                            *
 *                                                                           *
 *   You should have received a copy of the GNU General Public License       *
 *   along with this program.  If not, see <http://www.gnu.org/licenses/>.   *
 *****************************************************************************/
/*!
 * \file
 * \ingroup Fluidmatrixinteractions
 * \brief   Relation for the saturation-dependent effective thermal conductivity
 */
#ifndef DUMUX_MATERIAL_FLUIDMATRIX_THERMALCONDUCTIVITY_SOMERTON_HH
#define DUMUX_MATERIAL_FLUIDMATRIX_THERMALCONDUCTIVITY_SOMERTON_HH

#include <algorithm>
#include <cmath>

namespace Dumux {

/*!
 * \ingroup Fluidmatrixinteractions
 * \brief Relation for the saturation-dependent effective thermal conductivity
 *
 *  The Somerton method computes the thermal conductivity of dry and the wet soil material
 *  and uses a root function of the wetting saturation to compute the
 *  effective thermal conductivity for a two-phase fluidsystem. The individual thermal
 *  conductivities are calculated as geometric mean of the thermal conductivity of the porous
 *  material and of the respective fluid phase.
 *
 * The material law is:
 * \f$\mathrm{
 \lambda_\text{eff} = \lambda_{\text{dry}} + \sqrt{(S_w)} \left(\lambda_\text{wet} - \lambda_\text{dry}\right)
 }\f$
 *
 * with
 * \f$\mathrm{
 \lambda_\text{wet} = \lambda_{solid}^{\left(1-\phi\right)}*\lambda_w^\phi
 }\f$
 * and
 *
 * \f$\mathrm{
 \lambda_\text{dry} = \lambda_{solid}^{\left(1-\phi\right)}*\lambda_n^\phi.
 }\f$
 *
 */
template<class Scalar>
class ThermalConductivitySomerton
{
public:
    /*!
     * \brief effective thermal conductivity \f$\mathrm{[W/(m K)]}\f$ after Somerton (1974) \cite somerton1974 <BR>
     *
     * \param volVars volume variables
     * \return effective thermal conductivity \f$\mathrm{[W/(m K)]}\f$ after Somerton (1974) \cite somerton1974 <BR>
     *
     * This gives an interpolation of the effective thermal conductivities of a porous medium
     * filled with the nonwetting phase and a porous medium filled with the wetting phase.
     * These two effective conductivities are computed as geometric mean of the solid and the
     * fluid conductivities and interpolated with the square root of the wetting saturation.
     * See f.e. Ebigbo, A.: Thermal Effects of Carbon Dioxide Sequestration in the Subsurface, Diploma thesis \cite ebigbo2005 .
     */
    template<class VolumeVariables>
    static Scalar effectiveThermalConductivity(const VolumeVariables& volVars)
    {
        using FluidSystem = typename VolumeVariables::FluidSystem;
        static_assert(FluidSystem::numPhases == 2, "ThermalConductivitySomerton only works for two-phase fluid systems!");
        static_assert((FluidSystem::isGas(0) && !FluidSystem::isGas(1)) || (!FluidSystem::isGas(0) && FluidSystem::isGas(1)),
                     "ThermalConductivitySomerton only works if one phase is gaseous and one is liquid!");

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
     * \brief effective thermal conductivity \f$\mathrm{[W/(m K)]}\f$ after Somerton (1974) \cite somerton1974 <BR>
     *
     * \param satLiquid The saturation of the liquid phase
     * \param lambdaLiquid The thermal conductivity of the liquid phase in \f$\mathrm{[W/(m K)]}\f$
     * \param lambdaGas The thermal conductivity of the gas phase in \f$\mathrm{[W/(m K)]}\f$
     * \param lambdaSolid The thermal conductivity of the solid phase in \f$\mathrm{[W/(m K)]}\f$
     * \param porosity The porosity
     * \param rhoSolid The density of solid phase in \f$\mathrm{[kg/m^3]}\f$
     *
     * \return effective thermal conductivity \f$\mathrm{[W/(m K)]}\f$ after Somerton (1974) \cite somerton1974
     */
    static Scalar effectiveThermalConductivity_(const Scalar satLiquid,
                                                const Scalar lambdaLiquid,
                                                const Scalar lambdaGas,
                                                const Scalar lambdaSolid,
                                                const Scalar porosity,
                                                const Scalar rhoSolid = 0.0 /*unused*/)
    {
        using std::max;
        using std::pow;
        using std::sqrt;
        const Scalar satLiquidPhysical = max<Scalar>(0.0, satLiquid);
        // geometric mean, using ls^(1-p)*l^p = ls*(l/ls)^p
        const Scalar lSat = lambdaSolid * pow(lambdaLiquid / lambdaSolid, porosity);
        const Scalar lDry = lambdaSolid * pow(lambdaGas / lambdaSolid, porosity);

        return lDry + sqrt(satLiquidPhysical) * (lSat - lDry);
    }
};

} // end namespace Dumux

#endif
