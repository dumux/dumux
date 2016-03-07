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
 *
 * \brief   Relation for the saturation-dependent effective thermal conductivity
 */
#ifndef THERMALCONDUCTIVITY_SIMPLE_FLUID_LUMPING_HH
#define THERMALCONDUCTIVITY_SIMPLE_FLUID_LUMPING_HH

#include <algorithm>
#include <dumux/porousmediumflow/mpnc/implicit/properties.hh>

namespace Dumux
{

/*!
 * \ingroup fluidmatrixinteractionslaws
 * \brief   Relation for the saturation-dependent effective thermal conductivity
 */
template<class TypeTag>
class ThermalConductivitySimpleFluidLumping
{
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, Indices) Indices;
    enum { numEnergyEquations = GET_PROP_VALUE(TypeTag, NumEnergyEquations)};


public:
    /*!
     * \brief effective thermal conductivity \f$\mathrm{[W/(m K)]}\f$
     *
     * \param volVars volume variables
     * \param spatialParams spatial parameters
     * \param element element (to be passed to spatialParams)
     * \param fvGeometry fvGeometry (to be passed to spatialParams)
     * \param scvIdx scvIdx (to be passed to spatialParams)
     *
     * \return effective thermal conductivity \f$\mathrm{[W/(m K)]}\f$
     */
    template<class VolumeVariables, class SpatialParams, class Element, class FVGeometry>
    static Scalar effectiveThermalConductivity(const VolumeVariables& volVars,
                                               const SpatialParams& spatialParams,
                                               const Element& element,
                                               const FVGeometry& fvGeometry,
                                               int scvIdx)
    {
        Scalar sw = volVars.saturation(Indices::wPhaseIdx);
        Scalar lambdaW = volVars.fluidThermalConductivity(Indices::wPhaseIdx);
        Scalar lambdaN = volVars.fluidThermalConductivity(Indices::nPhaseIdx);
        Scalar lambdaSolid = volVars.solidThermalConductivity();
        Scalar porosity = volVars.porosity();

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
        assert(numEnergyEquations != 3) ;

        // Franz Lindner / Shi & Wang 2011
        const Scalar satW = std::max<Scalar>(0.0, sw);

        const Scalar kfeff = porosity *((1.-satW)*lambdaN + satW*lambdaW) ; // arithmetic

        Scalar keff ;

        if (numEnergyEquations == 2){ // solid dealed with individually (extra balance equation)
            keff = kfeff ;
        }
        else {
            const Scalar kseff = (1.0-porosity)  * lambdaSolid ;
            keff = kfeff  + kseff;
        }

        return keff ;
    }
};
}
#endif
