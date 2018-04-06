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
 * \ingroup SweModel
 * \copydoc Dumux::SweResidual
 */
#ifndef DUMUX_SWE_LOCAL_RESIDUAL_HH
#define DUMUX_SWE_LOCAL_RESIDUAL_HH

#include <dumux/common/properties.hh>
#include <dumux/discretization/methods.hh>
#include <dumux/assembly/fvlocalresidual.hh>
#include <dune/common/hybridutilities.hh>

namespace Dumux
{


/*!
 * \ingroup SweModel
 * \brief Element-wise calculation of the residual for the shallow water equations
 */
template<class TypeTag>
class SweResidual : public GET_PROP_TYPE(TypeTag, BaseLocalResidual)
{
    using ParentType = typename GET_PROP_TYPE(TypeTag, BaseLocalResidual);
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using Problem = typename GET_PROP_TYPE(TypeTag, Problem);
    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
    using ResidualVector = typename GET_PROP_TYPE(TypeTag, NumEqVector);
    using VolumeVariables = typename GET_PROP_TYPE(TypeTag, VolumeVariables);
    using ElementVolumeVariables = typename GET_PROP_TYPE(TypeTag, ElementVolumeVariables);
    using FluxVariables = typename GET_PROP_TYPE(TypeTag, FluxVariables);
    using ElementFluxVariablesCache = typename GET_PROP_TYPE(TypeTag, ElementFluxVariablesCache);
    using FVElementGeometry = typename GET_PROP_TYPE(TypeTag, FVGridGeometry)::LocalView;
    using SubControlVolume = typename FVElementGeometry::SubControlVolume;
    using SubControlVolumeFace = typename FVElementGeometry::SubControlVolumeFace;
    using Element = typename GridView::template Codim<0>::Entity;
    using Indices = typename GET_PROP_TYPE(TypeTag, Indices);
    using AdvectionType = typename GET_PROP_TYPE(TypeTag, AdvectionType);


    enum {
         // grid and world dimension
        dim = GridView::dimension,
        dimWorld = GridView::dimensionworld,

        // indices for primary variables
        massBalanceIdx = Indices::massBalanceIdx,
        momentumXBalanceIdx = Indices::momentumXBalanceIdx,
        momentumYBalanceIdx = Indices::momentumYBalanceIdx,
        waterdepthIdx = Indices::waterdepthIdx,
        velocityXIdx = Indices::velocityXIdx,
        velocityYIdx = Indices::velocityYIdx
    };


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
    ResidualVector computeStorage(const Problem& problem,
                                  const SubControlVolume& scv,
                                  const VolumeVariables& volVars) const
    {
        // partial time derivative of the phase mass
        ResidualVector storage(0.0);
        storage[massBalanceIdx] = volVars.getH();
        storage[momentumXBalanceIdx] = volVars.getH() * volVars.getU();
        storage[momentumYBalanceIdx] = volVars.getH() * volVars.getV();

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
     * \param elemFluxVarsCache The cache related to flux compuation
     */
    ResidualVector computeFlux(const Problem& problem,
                               const Element& element,
                               const FVElementGeometry& fvGeometry,
                               const ElementVolumeVariables& elemVolVars,
                               const SubControlVolumeFace& scvf,
                               const ElementFluxVariablesCache& elemFluxVarsCache) const
    {
        ResidualVector flux(0.0);

        auto numFlux = AdvectionType::flux(problem, element, fvGeometry, elemVolVars, scvf, elemFluxVarsCache);
        flux = numFlux;// + turbFlux;

        return flux;
    }
};
}

#endif   // DUMUX_SWE_LOCAL_RESIDUAL_HH
