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
 * \ingroup CCMpfaFlux
 * \brief Fourier's law for cell-centered finite volume schemes with multi-point flux approximation
 */
#ifndef DUMUX_DISCRETIZATION_CC_MPFA_FOURIERS_LAW_HH
#define DUMUX_DISCRETIZATION_CC_MPFA_FOURIERS_LAW_HH

#include <dune/common/dynvector.hh>
#include <dune/common/dynmatrix.hh>

#include <dumux/common/properties.hh>
#include <dumux/discretization/method.hh>

namespace Dumux {

//! forward declaration of the method-specific implementation
template<class TypeTag, DiscretizationMethod discMethod>
class FouriersLawImplementation;

/*!
 * \ingroup CCMpfaFlux
 * \brief Fourier's law for cell-centered finite volume schemes with two-point flux approximation
 */
template <class TypeTag>
class FouriersLawImplementation<TypeTag, DiscretizationMethod::ccmpfa>
{
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using Problem = GetPropType<TypeTag, Properties::Problem>;
    using GridView = typename GetPropType<TypeTag, Properties::GridGeometry>::GridView;
    using Element = typename GridView::template Codim<0>::Entity;

    static constexpr int dim = GridView::dimension;
    static constexpr int dimWorld = GridView::dimensionworld;

    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using FVElementGeometry = typename GridGeometry::LocalView;
    using SubControlVolumeFace = typename FVElementGeometry::SubControlVolumeFace;
    using ElementVolumeVariables = typename GetPropType<TypeTag, Properties::GridVolumeVariables>::LocalView;
    using ElementFluxVarsCache = typename GetPropType<TypeTag, Properties::GridFluxVariablesCache>::LocalView;
    using FluxVariablesCache = GetPropType<TypeTag, Properties::FluxVariablesCache>;

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
          if (fvGeometry.gridGeometry().vertexUsesSecondaryInteractionVolume(scvf.vertexIndex()))
              scvfFluxVarsCache.updateHeatConduction(fluxVarsCacheFiller.secondaryInteractionVolume(),
                                                     fluxVarsCacheFiller.secondaryIvLocalFaceData(),
                                                     fluxVarsCacheFiller.secondaryIvDataHandle());
          else
              scvfFluxVarsCache.updateHeatConduction(fluxVarsCacheFiller.primaryInteractionVolume(),
                                                     fluxVarsCacheFiller.primaryIvLocalFaceData(),
                                                     fluxVarsCacheFiller.primaryIvDataHandle());
        }
    };

    //! The cache used in conjunction with the mpfa Fourier's Law
    class MpfaFouriersLawCache
    {
        using DualGridNodalIndexSet = GetPropType<TypeTag, Properties::DualGridNodalIndexSet>;
        using Stencil = typename DualGridNodalIndexSet::NodalGridStencilType;

        static constexpr bool considerSecondaryIVs = GridGeometry::MpfaHelper::considerSecondaryIVs();
        using PrimaryDataHandle = typename ElementFluxVarsCache::PrimaryIvDataHandle::HeatConductionHandle;
        using SecondaryDataHandle = typename ElementFluxVarsCache::SecondaryIvDataHandle::HeatConductionHandle;

        //! sets the pointer to the data handle (overload for secondary data handles)
        template< bool doSecondary = considerSecondaryIVs, std::enable_if_t<doSecondary, int> = 0 >
        void setHandlePointer_(const SecondaryDataHandle& dataHandle)
        { secondaryHandlePtr_ = &dataHandle; }

        //! sets the pointer to the data handle (overload for primary data handles)
        void setHandlePointer_(const PrimaryDataHandle& dataHandle)
        { primaryHandlePtr_ = &dataHandle; }

    public:
        // export filler type
        using Filler = MpfaFouriersLawCacheFiller;

        /*!
         * \brief Update cached objects (transmissibilities).
         *        This is used for updates with primary interaction volumes.
         *
         * \param iv The interaction volume this scvf is embedded in
         * \param localFaceData iv-local info on this scvf
         * \param dataHandle Transmissibility matrix & gravity data of this iv
         */
         template<class IV, class LocalFaceData, class DataHandle>
         void updateHeatConduction(const IV& iv,
                                   const LocalFaceData& localFaceData,
                                   const DataHandle& dataHandle)
        {
            switchFluxSign_ = localFaceData.isOutsideFace();
            stencil_ = &iv.stencil();
            setHandlePointer_(dataHandle.heatConductionHandle());
        }

        //! The stencil corresponding to the transmissibilities (primary type)
        const Stencil& heatConductionStencil() const { return *stencil_; }

        //! The corresponding data handles
        const PrimaryDataHandle& heatConductionPrimaryDataHandle() const { return *primaryHandlePtr_; }
        const SecondaryDataHandle& heatConductionSecondaryDataHandle() const { return *secondaryHandlePtr_; }

        //! Returns whether or not this scvf is an "outside" face in the scope of the iv.
        bool heatConductionSwitchFluxSign() const { return switchFluxSign_; }

    private:
        bool switchFluxSign_;

        //! pointers to the corresponding iv-data handles
        const PrimaryDataHandle* primaryHandlePtr_;
        const SecondaryDataHandle* secondaryHandlePtr_;

        //! The stencil, i.e. the grid indices j
        const Stencil* stencil_;
    };

public:
    // state the discretization method this implementation belongs to
    static const DiscretizationMethod discMethod = DiscretizationMethod::ccmpfa;

    // state the type for the corresponding cache and its filler
    using Cache = MpfaFouriersLawCache;

    /*!
     * \brief Returns the heat flux within the porous medium
     *        (in J/s) across the given sub-control volume face.
     * \note This law assumes thermal equilibrium between the fluid
     *       and solid phases, and uses an effective thermal conductivity
     *       for the overall aggregate.
     */
    static Scalar flux(const Problem& problem,
                       const Element& element,
                       const FVElementGeometry& fvGeometry,
                       const ElementVolumeVariables& elemVolVars,
                       const SubControlVolumeFace& scvf,
                       const ElementFluxVarsCache& elemFluxVarsCache)
    {
        const auto& fluxVarsCache = elemFluxVarsCache[scvf];

        // forward to the private function taking the iv data handle
        if (fluxVarsCache.usesSecondaryIv())
            return flux_(problem, fluxVarsCache, fluxVarsCache.heatConductionSecondaryDataHandle());
        else
            return flux_(problem, fluxVarsCache, fluxVarsCache.heatConductionPrimaryDataHandle());
    }

private:
    template< class Problem, class FluxVarsCache, class DataHandle >
    static Scalar flux_(const Problem& problem,
                        const FluxVarsCache& cache,
                        const DataHandle& dataHandle)
    {
        const bool switchSign = cache.heatConductionSwitchFluxSign();

        const auto localFaceIdx = cache.ivLocalFaceIndex();
        const auto idxInOutside = cache.indexInOutsideFaces();
        const auto& Tj = dataHandle.uj();
        const auto& tij = dim == dimWorld ? dataHandle.T()[localFaceIdx]
                                          : (!switchSign ? dataHandle.T()[localFaceIdx]
                                                         : dataHandle.tijOutside()[localFaceIdx][idxInOutside]);
        Scalar scvfFlux = tij*Tj;

        // switch the sign if necessary
        if (switchSign)
            scvfFlux *= -1.0;

        return scvfFlux;
    }
};

} // end namespace Dumux

#endif
