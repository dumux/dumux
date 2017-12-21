// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   See the file COPYING for full copying permissions.                      *
 *                                                                           *
 *   This program is free software: you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation, either version 2 of the License, or       *
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
#ifndef THERMALCONDUCTIVITY_SOMERTON_3P_HH
#define THERMALCONDUCTIVITY_SOMERTON_3P_HH

#include <algorithm>
#include <cmath>

namespace Dumux
{

struct SimpleThreePIndices
{
    static const int wPhaseIdx = 0;
    static const int nPhaseIdx = 1;
    static const int gPhaseIdx = 2;
};

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
template<class Scalar, class Indices = SimpleThreePIndices>
class ThermalConductivitySomerton
{
public:
    /*!
     * \brief effective thermal conductivity \f$\mathrm{[W/(m K)]}\f$ after Somerton (1974) extended for a three phase system
     *
     * \param volVars volume variables
     * \param spatialParams spatial parameters
     * \param element element (to be passed to spatialParams)
     * \param fvGeometry fvGeometry (to be passed to spatialParams)
     * \param scvIdx scvIdx (to be passed to spatialParams)
     *
     * \return effective thermal conductivity \f$\mathrm{[W/(m K)]}\f$ after Somerton (1974)
     *
     * This gives an interpolation of the effective thermal conductivities of a porous medium
     * filled with the water phase (w), a NAPL phase (n) and a gas phase (g).
     * These two effective conductivities are computed as geometric mean of the solid and the
     * fluid conductivities and interpolated with the square root of the wetting saturation.
     */
    template<class VolumeVariables, class SpatialParams, class Element, class FVGeometry, class SubControlVolume>
    static Scalar effectiveThermalConductivity(const VolumeVariables& volVars,
                                               const SpatialParams& spatialParams,
                                               const Element& element,
                                               const FVGeometry& fvGeometry,
                                               const SubControlVolume& scv)
    {
        Scalar sw = volVars.saturation(Indices::wPhaseIdx);
        Scalar sn = volVars.saturation(Indices::nPhaseIdx);
        Scalar lambdaW = volVars.fluidThermalConductivity(Indices::wPhaseIdx);
        Scalar lambdaN = volVars.fluidThermalConductivity(Indices::nPhaseIdx);
        Scalar lambdaG = volVars.fluidThermalConductivity(Indices::gPhaseIdx);
        Scalar lambdaSolid = volVars.solidThermalConductivity();
        Scalar porosity = volVars.porosity();

        return effectiveThermalConductivity(sw, sn, lambdaW, lambdaN, lambdaG, lambdaSolid, porosity);
    }

    /*!
     * \brief effective thermal conductivity \f$\mathrm{[W/(m K)]}\f$ after Somerton (1974)
     *
     * \param sw The saturation of the wetting phase
     * \param sn The saturation of the non-wetting phase
     * \param lambdaW the thermal conductivity of the water phase in \f$\mathrm{[W/(m K)]}\f$
     * \param lambdaN the thermal conductivity of the NAPL phase in \f$\mathrm{[W/(m K)]}\f$
     * \param lambdaG the thermal conductivity of the gas phase in \f$\mathrm{[W/(m K)]}\f$
     * \param lambdaSolid the thermal conductivity of the solid phase in \f$\mathrm{[W/(m K)]}\f$
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

//        const Scalar lSw = 1.8; //pow(lambdaSolid, (1.0 - porosity)) * pow(lambdaW, porosity);
//        const Scalar lSn = 0.65; //pow(lambdaSolid, (1.0 - porosity)) * pow(lambdaN, porosity);
//        const Scalar lSg = 0.35; //pow(lambdaSolid, (1.0 - porosity)) * pow(lambdaG, porosity);
        // porosity weighted geometric mean
        Scalar lSw = pow(lambdaSolid, (1.0 - porosity)) * pow(lambdaW, porosity);
        Scalar lSn = pow(lambdaSolid, (1.0 - porosity)) * pow(lambdaN, porosity);
        Scalar lSg = pow(lambdaSolid, (1.0 - porosity)) * pow(lambdaG, porosity);
        Scalar lambdaEff = lSg + sqrt(satW) * (lSw - lSg) + sqrt(satN) * (lSn -lSg);

        return lambdaEff;

    }
};
}
#endif
