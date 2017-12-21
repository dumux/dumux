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
 * \brief This file contains the class which is required to calculate
 *        heat conduction fluxes with Fourier's law for cell-centered MPFA models.
 */
#ifndef DUMUX_DISCRETIZATION_CC_MPFA_FOURIERS_LAW_HH
#define DUMUX_DISCRETIZATION_CC_MPFA_FOURIERS_LAW_HH

#include <dumux/common/properties.hh>
#include <dumux/discretization/methods.hh>

namespace Dumux
{
// forward declaration
template<class TypeTag, DiscretizationMethods discMethod>
class FouriersLawImplementation;

/*!
 * \ingroup Mpfa
 * \brief Specialization of Fourier's Law for the CCMpfa method.
 */
template <class TypeTag>
class FouriersLawImplementation<TypeTag, DiscretizationMethods::CCMpfa>
{
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using Problem = typename GET_PROP_TYPE(TypeTag, Problem);
    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
    using Element = typename GridView::template Codim<0>::Entity;
    using SubControlVolumeFace = typename GET_PROP_TYPE(TypeTag, SubControlVolumeFace);
    using FVElementGeometry = typename GET_PROP_TYPE(TypeTag, FVElementGeometry);
    using ElementVolumeVariables = typename GET_PROP_TYPE(TypeTag, ElementVolumeVariables);
    using ElementFluxVarsCache = typename GET_PROP_TYPE(TypeTag, ElementFluxVariablesCache);
    using FluxVariablesCache = typename GET_PROP_TYPE(TypeTag, FluxVariablesCache);
    using ThermalConductivityModel = typename GET_PROP_TYPE(TypeTag, ThermalConductivityModel);

    // Always use the dynamic type for vectors (compatibility with the boundary)
    using PrimaryInteractionVolume = typename GET_PROP_TYPE(TypeTag, PrimaryInteractionVolume);
    using CoefficientVector = typename PrimaryInteractionVolume::Traits::DynamicVector;
    using DataHandle = typename PrimaryInteractionVolume::Traits::DataHandle;

    static constexpr int dim = GridView::dimension;
    static constexpr int dimWorld = GridView::dimensionworld;
    static constexpr int energyEqIdx = GET_PROP_TYPE(TypeTag, Indices)::energyEqIdx;

    //! Class that fills the cache corresponding to mpfa Darcy's Law
    class MpfaFouriersLawCacheFiller
    {
    public:
        //! Function to fill an MpfaDarcysLawCache of a given scvf
        //! This interface has to be met by any cache filler class for heat conduction quantities
        template<class FluxVariablesCacheFiller>
        static void fill(FluxVariablesCache& scvfFluxVarsCache,
                         const Problem& problem,
                         const Element& element,
                         const FVElementGeometry& fvGeometry,
                         const ElementVolumeVariables& elemVolVars,
                         const SubControlVolumeFace& scvf,
                         const FluxVariablesCacheFiller& fluxVarsCacheFiller)
        {
          // get interaction volume from the flux vars cache filler & upate the cache
          if (fvGeometry.fvGridGeometry().vertexUsesSecondaryInteractionVolume(scvf.vertexIndex()))
              scvfFluxVarsCache.updateHeatConduction(fluxVarsCacheFiller.secondaryInteractionVolume(),
                                                     fluxVarsCacheFiller.dataHandle(),
                                                     scvf);
          else
              scvfFluxVarsCache.updateHeatConduction(fluxVarsCacheFiller.primaryInteractionVolume(),
                                                     fluxVarsCacheFiller.dataHandle(),
                                                     scvf);
        }
    };

    //! The cache used in conjunction with the mpfa Fourier's Law
    class MpfaFouriersLawCache
    {
        // We always use the dynamic types here to be compatible on the boundary
        using Stencil = typename PrimaryInteractionVolume::Traits::DynamicGlobalIndexContainer;
        using DirichletDataContainer = typename PrimaryInteractionVolume::DirichletDataContainer;
    public:
        // export filler type
        using Filler = MpfaFouriersLawCacheFiller;

        // update cached objects for heat conduction
        template<class InteractionVolume>
        void updateHeatConduction(const InteractionVolume& iv, const DataHandle& dataHandle, const SubControlVolumeFace &scvf)
        {
            const auto& localFaceData = iv.getLocalFaceData(scvf);

            // update the quantities that are equal for all phases
            heatConductionSwitchFluxSign_ = localFaceData.isOutside();
            heatConductionVolVarsStencil_ = &dataHandle.volVarsStencil();
            heatConductionDirichletData_ = &dataHandle.dirichletData();

            // the transmissibilities on surface grids have to be obtained from the outside
            if (dim == dimWorld)
                heatConductionTij_ = &dataHandle.T()[localFaceData.ivLocalScvfIndex()];
            else
                heatConductionTij_ = localFaceData.isOutside() ?
                                     &dataHandle.outsideTij()[localFaceData.ivLocalOutsideScvfIndex()] :
                                     &dataHandle.T()[localFaceData.ivLocalScvfIndex()];
        }

        //! Returns the stencil for heat conduction flux computation on an scvf
        const Stencil& heatConductionVolVarsStencil() const { return *heatConductionVolVarsStencil_; }

        //! Returns the transmissibilities associated with the volume variables
        const CoefficientVector& heatConductionTij() const { return *heatConductionTij_; }

        //! On faces that are "outside" w.r.t. a face in the interaction volume,
        //! we have to take the negative value of the fluxes, i.e. multiply by -1.0
        bool heatConductionSwitchFluxSign() const { return heatConductionSwitchFluxSign_; }

        //! Returns the data on dirichlet boundary conditions affecting
        //! the flux computation on this face
        const DirichletDataContainer& heatConductionDirichletData() const { return *heatConductionDirichletData_; }

    private:
        bool heatConductionSwitchFluxSign_;
        const Stencil* heatConductionVolVarsStencil_;
        const CoefficientVector* heatConductionTij_;
        const DirichletDataContainer* heatConductionDirichletData_;
    };

public:
    // state the discretization method this implementation belongs to
    static const DiscretizationMethods myDiscretizationMethod = DiscretizationMethods::CCMpfa;

    // state the type for the corresponding cache and its filler
    using Cache = MpfaFouriersLawCache;

    static Scalar flux(const Problem& problem,
                       const Element& element,
                       const FVElementGeometry& fvGeometry,
                       const ElementVolumeVariables& elemVolVars,
                       const SubControlVolumeFace& scvf,
                       const ElementFluxVarsCache& elemFluxVarsCache)
    {
        // prepare computations
        Scalar flux(0.0);
        unsigned int i = 0;
        const auto& fluxVarsCache = elemFluxVarsCache[scvf];
        const auto& tij = fluxVarsCache.heatConductionTij();

        // calculate Tij*tj
        for (const auto volVarIdx : fluxVarsCache.heatConductionVolVarsStencil())
            flux += tij[i++]*elemVolVars[volVarIdx].temperature();

        // add contributions from dirichlet BCs
        for (const auto& d : fluxVarsCache.heatConductionDirichletData())
            flux += tij[i++]*elemVolVars[d.volVarIndex()].temperature();

        // return overall resulting flux
        return fluxVarsCache.heatConductionSwitchFluxSign() ? -flux : flux;
    }
};

} // end namespace Dumux

#endif
