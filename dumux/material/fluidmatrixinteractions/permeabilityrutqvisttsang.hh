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
 * \brief   Relation for the effective permeability
 */
#ifndef PERMEABILITYRUTQVISTSTANG_HH
#define PERMEABILITYRUTQVISTSTANG_HH

#include <algorithm>
#include <cmath>

namespace Dumux
{

/*!
 * \ingroup Fluidmatrixinteractions
 * \brief Relation for the effective permeability
 *
 *  After Rutqvist and Tsang (2002) \cite rutqvist2002, the effective permeability can be
 *  calculated from the intrinsic permeability \f$ k_{\text{0}} \f$, the initial porosity
 *  \f$ \phi_{\text{0}} \f$ and the effective porosity \f$ \phi_{\text{eff}} \f$ as used
 *  in Darcis, M.: Coupling Models of Different Complexity for the Simulation of CO2
 *  Storage in Deep Saline Aquifers, PhD thesis \cite darcis2013.
 *  The calculation is shown here for one diagonal entry of the permeability tensor:
 *
 * \f$ k_\text{eff} = k_{\text{0}} \cdot \text{exp} \left[22.2\left(\phi_\text{eff}/ \phi_\text{0} -1 \right) \right] \f$
 *
 */
template<class Scalar>
class PermeabilityRutqvistTsang
{
public:
    /*!
     * \brief effective permeability tensor \f$\mathrm{[m^{2})]}\f$ after Rutqvist and Tsang (2002) \cite rutqvist2002 <BR>
     *
     * \param volVars volume variables
     * \param spatialParams spatial parameters
     * \param element element (to be passed to spatialParams)
     * \param fvGeometry fvGeometry (to be passed to spatialParams)
     * \param scvIdx control volumne
     *
     * \return effective permeability tensor \f$\mathrm{[m^{2})]}\f$ after Rutqvist and Tsang (2002) \cite rutqvist2002 <BR>
     *
     * This calculates the entries of effective permeability tensor from the intrinsic
     * permeability, the initial and the effective porosity as used in Darcis, M.:
     * Coupling Models of Different Complexity for the Simulation of CO2 Storage in
     * Deep Saline Aquifers, PhD thesis \cite darcis2013 .
     */
    template<class VolumeVariables, class SpatialParams, class Element, class FVGeometry>
    static auto effectivePermeability(const VolumeVariables& volVars,
                                        const SpatialParams& spatialParams,
                                        const Element& element,
                                        const FVGeometry& fvGeometry,
                                               int scvIdx)
    {
        Scalar effPorosity = volVars.effPorosity;
        Scalar initialPorosity = spatialParams.porosity(element, fvGeometry, scvIdx);

        // evaluate effective permeabilities for nodes i and j based on the
        // effective porosities, the initial porosities and the initial permeabilities
        Scalar exponent =
                22.2
                        * (effPorosity / initialPorosity - 1);

        auto Keff
            = spatialParams.intrinsicPermeability(element, fvGeometry, scvIdx);
        using std::exp;
        Keff *= exp(exponent);

        return Keff;
    }

};
}
#endif
