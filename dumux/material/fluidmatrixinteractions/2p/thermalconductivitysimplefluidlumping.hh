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
#ifndef DUMUX_MATERIAL_THERMALCONDUCTIVITY_SIMPLE_FLUID_LUMPING_HH
#define DUMUX_MATERIAL_THERMALCONDUCTIVITY_SIMPLE_FLUID_LUMPING_HH

#include <assert.h>
#include <algorithm>

#include <dune/common/deprecated.hh>

namespace Dumux {

/*!
 * \ingroup Fluidmatrixinteractions
 * \brief   Relation for the saturation-dependent effective thermal conductivity
 * \todo This shouldn't depend on TypeTag!!
 */
template<class Scalar, int numEnergyEquationsFluid>
class ThermalConductivitySimpleFluidLumping
{

public:
    /*!
     * \brief effective thermal conductivity \f$\mathrm{[W/(m K)]}\f$
     */
    template<class VolumeVariables, class SpatialParams, class Element, class FVGeometry>
    DUNE_DEPRECATED_MSG("Signature deprecated. Use signature with volume variables only!")
    static Scalar effectiveThermalConductivity(const VolumeVariables& volVars,
                                               const SpatialParams& spatialParams,
                                               const Element& element,
                                               const FVGeometry& fvGeometry,
                                               const typename FVGeometry::SubControlVolume& scv)
    {
        return effectiveThermalConductivity(volVars);
    }

    /*!
     * \brief Effective thermal conductivity \f$\mathrm{[W/(m K)]}\f$
     *
     * \param volVars volume variables
     * \return effective thermal conductivity \f$\mathrm{[W/(m K)]}\f$
     * \todo TODO: Fix this law for changing wettability
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

        return effectiveThermalConductivity(sw, lambdaW, lambdaN, lambdaSolid, porosity);
    }

    /*!
     * \brief Returns the effective thermal conductivity \f$\mathrm{[W/(m K)]}\f$.
     *
     * \param sw The saturation of the wetting phase
     * \param lambdaW The thermal conductivity of the wetting phase in \f$\mathrm{[W/(m K)]}\f$
     * \param lambdaN The thermal conductivity of the non-wetting phase in \f$\mathrm{[W/(m K)]}\f$
     * \param lambdaSolid The thermal conductivity of the solid phase in \f$\mathrm{[W/(m K)]}\f$
     * \param porosity The porosity
     * \param rhoSolid The density of the solid phase in \f$\mathrm{[kg/m^3]}\f$
     *
     * \return Effective thermal conductivity of the fluid phases
     */
    static Scalar effectiveThermalConductivity(const Scalar sw,
                                               const Scalar lambdaW,
                                               const Scalar lambdaN,
                                               const Scalar lambdaSolid,
                                               const Scalar porosity,
                                               const Scalar rhoSolid = 0.0 /*unused*/)
    {
        assert(numEnergyEquationsFluid != 2) ;

        // Franz Lindner / Shi & Wang 2011
        using std::max;
        const Scalar satW = max<Scalar>(0.0, sw);

        const Scalar kfeff = porosity *((1.-satW)*lambdaN + satW*lambdaW) ; // arithmetic

        Scalar keff ;

        if (numEnergyEquationsFluid == 1){ // solid dealed with individually (extra balance equation)
            keff = kfeff ;
        }
        else {
            const Scalar kseff = (1.0-porosity)  * lambdaSolid ;
            keff = kfeff  + kseff;
        }

        return keff ;
    }
};
} // end namespace Dumux
#endif
