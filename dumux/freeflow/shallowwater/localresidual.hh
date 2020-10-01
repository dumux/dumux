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
 * \ingroup ShallowWaterModel
 * \copydoc Dumux::ShallowWaterResidual
 */
#ifndef DUMUX_FREEFLOW_SHALLOW_WATER_LOCAL_RESIDUAL_HH
#define DUMUX_FREEFLOW_SHALLOW_WATER_LOCAL_RESIDUAL_HH

#include <dumux/common/properties.hh>

namespace Dumux{

/*!
 * \ingroup ShallowWaterModel
 * \brief Element-wise calculation of the residual for the shallow water equations
 */
template<class TypeTag>
class ShallowWaterResidual
: public GetPropType<TypeTag, Properties::BaseLocalResidual>
{
    using ParentType = GetPropType<TypeTag, Properties::BaseLocalResidual>;
    using Problem = GetPropType<TypeTag, Properties::Problem>;
    using GridView = typename GetPropType<TypeTag, Properties::GridGeometry>::GridView;
    using NumEqVector = GetPropType<TypeTag, Properties::NumEqVector>;
    using VolumeVariables = GetPropType<TypeTag, Properties::VolumeVariables>;
    using ElementVolumeVariables = typename GetPropType<TypeTag, Properties::GridVolumeVariables>::LocalView;
    using ElementFluxVariablesCache = typename GetPropType<TypeTag, Properties::GridFluxVariablesCache>::LocalView;
    using FVElementGeometry = typename GetPropType<TypeTag, Properties::GridGeometry>::LocalView;
    using SubControlVolume = typename FVElementGeometry::SubControlVolume;
    using SubControlVolumeFace = typename FVElementGeometry::SubControlVolumeFace;
    using Element = typename GridView::template Codim<0>::Entity;
    using Indices = typename GetPropType<TypeTag, Properties::ModelTraits>::Indices;
    using FluxVariables = GetPropType<TypeTag, Properties::FluxVariables>;

public:

    using ParentType::ParentType;

    /*!
     * \brief Evaluate the rate of change of all conservation
     *        quantites (e.g. mass, momentum) within a sub-control
     *        volume of a finite volume element.
     * \param problem The problem
     * \param scv The sub control volume
     * \param volVars The current or previous volVars
     * \note This function should not include the source and sink terms.
     * \note The volVars can be different to allow computing
     *       the implicit euler time derivative here
     */
    NumEqVector computeStorage(const Problem& problem,
                               const SubControlVolume& scv,
                               const VolumeVariables& volVars) const
    {
        // partial time derivative of the phase mass
        NumEqVector storage(0.0);
        storage[Indices::massBalanceIdx] = volVars.waterDepth();
        storage[Indices::momentumXBalanceIdx] = volVars.waterDepth() * volVars.velocity(0);
        storage[Indices::momentumYBalanceIdx] = volVars.waterDepth() * volVars.velocity(1);
        return storage;
    }

       /*!
     * \brief Evaluate the mass flux over a face of a sub control volume
     *
     * \param problem The problem
     * \param element The current element.
     * \param fvGeometry The finite-volume geometry
     * \param elemVolVars The volume variables of the current element
     * \param scvf The sub control volume face to compute the flux on
     * \param elemFluxVarsCache the flux variable cache for the element stencil
    */
    NumEqVector computeFlux(const Problem& problem,
                            const Element& element,
                            const FVElementGeometry& fvGeometry,
                            const ElementVolumeVariables& elemVolVars,
                            const SubControlVolumeFace& scvf,
                            const ElementFluxVariablesCache& elemFluxVarsCache) const
    {
        NumEqVector flux(0.0);
        FluxVariables fluxVars;
        flux += fluxVars.advectiveFlux(problem, element, fvGeometry, elemVolVars, scvf);
        flux += fluxVars.diffusiveFlux(problem, element, fvGeometry, elemVolVars, scvf);
        return flux;
    }
};
} // end namespace Dumux

#endif
