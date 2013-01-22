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
#ifndef SOMERTON_HH
#define SOMERTON_HH

#include <algorithm>

namespace Dumux
{
/*!
 * \ingroup fluidmatrixinteractionslaws
 *
 * \brief Relation for the saturation-dependent effective thermal conductivity
 *
 *  The Somerton method computes the thermal conductivity of dry and the wet soil material
 *  and uses a root function of the wetting saturation to compute the
 *  effective thermal conductivity for a two-phase fluidsystem.
 */
template<class TypeTag>
class Somerton
{
public:
    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
    typedef typename GridView::template Codim<0>::Entity Element;
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, ElementVolumeVariables) ElementVolumeVariables;
    typedef typename GET_PROP_TYPE(TypeTag, SpatialParams) SpatialParams;
    typedef typename GET_PROP_TYPE(TypeTag, FVElementGeometry) FVElementGeometry;

    typedef typename GET_PROP_TYPE(TypeTag, Indices) Indices;
    enum { wPhaseIdx = Indices::wPhaseIdx };


    /*!
     * \brief Returns the effective thermal conductivity \f$[W/m^2]\f$ after Somerton (1974).
     *
     * The material law is:
     * \f[
     l_eff = l_solid + (l_wet - l_solid)
     \f]
     *
     * \param element The finite element
     * \param elemVolVars The volume variables on the element
     * \param fvGeometry The finite volume geometry
     * \param spatialParams The spatial parameters
     * \param scvIdx The local index of the sub-control volume where
     *                    the effective thermal conductivity is computed
     *
     * \return Effective thermal conductivity \f$[W/m^2]\f$ after Somerton (1974)
     */
    static Scalar effectiveThermalConductivity(const Element &element,
                                               const ElementVolumeVariables &elemVolVars,
                                               const FVElementGeometry &fvGeometry,
                                               const SpatialParams &spatialParams,
                                               const int scvIdx)
    {
        const Scalar lambdaSolid = spatialParams.thermalConductivitySolid(element, fvGeometry, scvIdx);
        const Scalar porosity = spatialParams.porosity(element, fvGeometry, scvIdx);

        const Scalar Sw = std::max<Scalar>(0.0, elemVolVars[scvIdx].saturation(wPhaseIdx));
        const Scalar lWater = elemVolVars[scvIdx].thermalConductivity(wPhaseIdx);

        const Scalar lSat = std::pow(lambdaSolid, (1.0 - porosity)) * std::pow(lWater, porosity);
        const Scalar lDry = std::pow(lambdaSolid, (1.0 - porosity));

        return lDry + std::sqrt(Sw) * (lDry - lSat);
    }
};
}
#endif
