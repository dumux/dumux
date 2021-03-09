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
#ifndef DUMUX_MATERIAL_THERMALCONDUCTIVITY_SOMERTON_3P_HH
#define DUMUX_MATERIAL_THERMALCONDUCTIVITY_SOMERTON_3P_HH

#include <algorithm>
#include <cmath>

namespace Dumux {

/*!
 * \ingroup Fluidmatrixinteractions
 * \brief Relation for the saturation-dependent effective thermal conductivity
 *
 *  The Somerton method computes the thermal conductivity of dry and the wet soil material.
 *  It is extended here to a three phase system of water (w), NAPL (n) and gas (g).
 *  It uses a root function of the water and NAPL saturation to compute the
 *  effective thermal conductivity for a three-phase fluidsystem. The individual thermal
 *  conductivities are calculated as geometric mean of the thermal conductivity of the porous
 *  material and of the respective fluid phase.
 *
 * The material law is:
 * \f[
 \lambda_\text{eff} = \lambda_\text{g,eff} + \sqrt{(S_w)} \left(\lambda_\text{w,eff} - \lambda_\text{g,eff}\right) +
 \sqrt{(S_n)} \left(\lambda0_\text{n,eff} - \lambda_\text{g,eff}\right)
 \f]
 *
 * with
 * \f[
 \lambda_\text{w,eff} = \lambda_{solid}^{\left(1-\phi\right)}*\lambda_w^\phi
 \f]
 * and
 *
 * \f[
 \lambda0_\text{n,eff} = \lambda_{solid}^{\left(1-\phi\right)}*\lambda_n^\phi.
 \f]
 *
 * * \f[
 \lambda_\text{g,eff} = \lambda_{solid}^{\left(1-\phi\right)}*\lambda_g^\phi.
 \f]
 */
template<class Scalar>
class ThermalConductivitySomerton
{
public:
    /*!
     * \brief effective thermal conductivity \f$\mathrm{[W/(m K)]}\f$ after Somerton (1974) extended for a three phase system
     *
     * \param volVars volume variables
     *
     * \return effective thermal conductivity \f$\mathrm{[W/(m K)]}\f$ after Somerton (1974)
     *
     * This gives an interpolation of the effective thermal conductivities of a porous medium
     * filled with the water phase (w), a NAPL phase (n) and a gas phase (g).
     * These two effective conductivities are computed as geometric mean of the solid and the
     * fluid conductivities and interpolated with the square root of the wetting saturation.
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
     * \brief effective thermal conductivity \f$\mathrm{[W/(m K)]}\f$ after Somerton (1974)
     *
     * \param sw The saturation of the wetting phase
     * \param sn The saturation of the nonwetting phase
     * \param lambdaW The thermal conductivity of the water phase in \f$\mathrm{[W/(m K)]}\f$
     * \param lambdaN The thermal conductivity of the NAPL phase in \f$\mathrm{[W/(m K)]}\f$
     * \param lambdaG The thermal conductivity of the gas phase in \f$\mathrm{[W/(m K)]}\f$
     * \param lambdaSolid The thermal conductivity of the solid phase in \f$\mathrm{[W/(m K)]}\f$
     * \param porosity The porosity
     *
     * \return effective thermal conductivity \f$\mathrm{[W/(m K)]}\f$ after Somerton (1974)
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
