 // -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
#ifndef DUMUX_MATERIAL_FLUIDMATRIX_THERMALCONDUCTIVITY_JOHANSEN_HH
#define DUMUX_MATERIAL_FLUIDMATRIX_THERMALCONDUCTIVITY_JOHANSEN_HH

#include <cmath>
#include <algorithm>

namespace Dumux {

/*!
 * \addtogroup EffectiveHeatConductivity
 * \copydetails Dumux::ThermalConductivityJohansen
 */

/*!
 * \ingroup EffectiveHeatConductivity
 * \brief Relation for the saturation-dependent effective thermal conductivity
 *
 * ### Johansen (two fluid phases)
 *
 * `ThermalConductivityJohansen` \cite johansen1977 computes the thermal conductivity of dry and the
 * wet soil material and interpolates using the Kersten number. The effective wet conductivity
 * is based on a geometric average and the effective dry conductivity is based on a semi-emprical
 * relation and fitted to medium quartz sand.
 *
 * The effective thermal conductivity is given by
 * \f[
 * \lambda_\text{eff} = \lambda_{\text{dry}} + \text{Ke} \left(\lambda_\text{wet} - \lambda_\text{dry}\right), \quad
 * \lambda_\text{wet} = \lambda_\text{s}^{\left(1-\phi\right)} \lambda_\text{w}^\phi, \quad
 * \lambda_\text{dry} = \frac{0.135 \rho_\text{s} \phi + 64.7}{\rho_\text{s} - 0.947 \rho_\text{s} \phi},
 * \f]
 * where \f$ \phi \f$ is the porosity, \f$ \lambda_\alpha \f$ is the thermal conductivity
 * of phase \f$ \alpha \f$, \f$ \rho_\text{s} \f$ denotes the density of the solid phase, and the
 * Kersten number is given by \f$ \text{Ke} = (\kappa S_\text{w})/(1 + (1-\kappa) S_\text{w}) \f$,
 * with the wetting phase saturation \f$ S_w \f$ and a fitting parameter \f$ \kappa = 15.6 \f$
 * for medium quartz sand.
 */
template<class Scalar>
class ThermalConductivityJohansen
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
        static_assert(FluidSystem::numPhases == 2, "ThermalConductivitySomerton only works for two-phase fluid systems!");
        // TODO: there should be an assertion that the indices are correct and 0 is actually the wetting phase!

        const Scalar sw = volVars.saturation(volVars.wettingPhase());
        const Scalar lambdaW = volVars.fluidThermalConductivity(volVars.wettingPhase());
        const Scalar lambdaN = volVars.fluidThermalConductivity(1-volVars.wettingPhase());
        const Scalar lambdaSolid = volVars.solidThermalConductivity();
        const Scalar porosity = volVars.porosity();
        const Scalar rhoSolid = volVars.solidDensity();

        return effectiveThermalConductivity_(sw, lambdaW, lambdaN, lambdaSolid, porosity, rhoSolid);
    }

private:
    /*!
     * \brief Effective thermal conductivity in \f$\mathrm{W/(m K)}\f$ for two phases
     *
     * \param Sw The saturation of the wetting phase
     * \param lambdaW The thermal conductivity of the wetting phase in \f$\mathrm{W/(m K)}\f$
     * \param lambdaN The thermal conductivity of the nonwetting phase in \f$\mathrm{W/(m K)}\f$
     * \param lambdaSolid The thermal conductivity of the solid phase in \f$\mathrm{W/(m K)}\f$
     * \param porosity The porosity
     * \param rhoSolid The density of solid phase in \f$\mathrm{kg/m^3}\f$
     *
     * \return Effective thermal conductivity in \f$\mathrm{W/(m K)}\f$ for two phases
     */
    static Scalar effectiveThermalConductivity_(const Scalar Sw,
                                                const Scalar lambdaW,
                                                const Scalar lambdaN,
                                                const Scalar lambdaSolid,
                                                const Scalar porosity,
                                                const Scalar rhoSolid)
    {
        using std::max;
        const Scalar satW = max<Scalar>(0.0, Sw);

        const Scalar kappa = 15.6; // fitted to medium quartz sand
        const Scalar rhoBulk = rhoSolid*porosity;

        using std::pow;

        const Scalar lambdaSaturated = lambdaSolid * pow(lambdaW / lambdaSolid, porosity);
        const Scalar lambdaDry = (0.135*rhoBulk + 64.7)/(rhoSolid - 0.947*rhoBulk);
        const Scalar Ke = (kappa*satW)/(1+(kappa-1)*satW);// Kersten number, equation 13

        return lambdaDry + Ke * (lambdaSaturated - lambdaDry); // equation 14
    }
};

} // end namespace Dumux

#endif
