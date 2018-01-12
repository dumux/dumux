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
 * \ingroup CCMpfaDiscretization
 * \brief Fourier's law for cell-centered finite volume schemes with multi-point flux approximation
 */
#ifndef DUMUX_DISCRETIZATION_CC_MPFA_FOURIERS_LAW_HH
#define DUMUX_DISCRETIZATION_CC_MPFA_FOURIERS_LAW_HH

#include <dune/common/dynvector.hh>
#include <dune/common/dynmatrix.hh>

#include <dumux/common/properties.hh>
#include <dumux/discretization/methods.hh>

#include <dumux/discretization/methods.hh>

namespace Dumux
{
//! forward declaration of the method-specific implementation
template<class TypeTag, DiscretizationMethods discMethod>
class FouriersLawImplementation;

/*!
* \ingroup CCMpfaDiscretization
* \brief Fourier's law for cell-centered finite volume schemes with two-point flux approximation
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
                                                     fluxVarsCacheFiller.secondaryIvLocalFaceData(),
                                                     fluxVarsCacheFiller.secondaryIvDataHandle(),
                                                     scvf);
          else
              scvfFluxVarsCache.updateHeatConduction(fluxVarsCacheFiller.primaryInteractionVolume(),
                                                     fluxVarsCacheFiller.primaryIvLocalFaceData(),
                                                     fluxVarsCacheFiller.primaryIvDataHandle(),
                                                     scvf);
        }
    };

    //! The cache used in conjunction with the mpfa Fourier's Law
    class MpfaFouriersLawCache
    {
        using DualGridNodalIndexSet = typename GET_PROP_TYPE(TypeTag, DualGridNodalIndexSet);
        using Stencil = typename DualGridNodalIndexSet::GridStencilType;

        using MpfaHelper = typename GET_PROP_TYPE(TypeTag, MpfaHelper);
        static constexpr bool considerSecondaryIVs = MpfaHelper::considerSecondaryIVs();

        // In the current implementation of the flux variables cache we cannot make a
        // disctinction between dynamic (e.g. mpfa-o unstructured) and static (e.g.mpfa-l)
        // matrix and vector types, as currently the cache class can only be templated
        // by a type tag (and there can only be one). Currently, pointers to both the
        // primary and secondary iv data is stored. Before accessing it has to be checked
        // whether or not the scvf is embedded in a secondary interaction volume.
        using PrimaryInteractionVolume = typename GET_PROP_TYPE(TypeTag, PrimaryInteractionVolume);
        using PrimaryIvLocalFaceData = typename PrimaryInteractionVolume::Traits::LocalFaceData;
        using PrimaryIvDataHandle = typename PrimaryInteractionVolume::Traits::DataHandle;
        using PrimaryIvVector = typename PrimaryInteractionVolume::Traits::Vector;
        using PrimaryIvMatrix = typename PrimaryInteractionVolume::Traits::Matrix;
        using PrimaryIvTij = typename PrimaryIvMatrix::row_type;

        using SecondaryInteractionVolume = typename GET_PROP_TYPE(TypeTag, SecondaryInteractionVolume);
        using SecondaryIvLocalFaceData = typename SecondaryInteractionVolume::Traits::LocalFaceData;
        using SecondaryIvDataHandle = typename SecondaryInteractionVolume::Traits::DataHandle;
        using SecondaryIvVector = typename SecondaryInteractionVolume::Traits::Vector;
        using SecondaryIvMatrix = typename SecondaryInteractionVolume::Traits::Matrix;
        using SecondaryIvTij = typename SecondaryIvMatrix::row_type;

        static constexpr int dim = GridView::dimension;
        static constexpr int dimWorld = GridView::dimensionworld;

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
         * \param scvf The sub-control volume face
         */
        void updateHeatConduction(const PrimaryInteractionVolume& iv,
                                  const PrimaryIvLocalFaceData& localFaceData,
                                  const PrimaryIvDataHandle& dataHandle,
                                  const SubControlVolumeFace &scvf)
        {
            stencil_ = &iv.stencil();
            switchFluxSign_ = localFaceData.isOutside();

            // store pointer to the temperature vector of this iv
            primaryTj_ = &dataHandle.temperatures();

            const auto ivLocalIdx = localFaceData.ivLocalScvfIndex();
            if (dim == dimWorld)
                primaryTij_ = &dataHandle.heatConductionT()[ivLocalIdx];
            else
                primaryTij_ = localFaceData.isOutside() ? &dataHandle.heatConductionTout()[ivLocalIdx][localFaceData.scvfLocalOutsideScvfIndex()]
                                                        : &dataHandle.heatConductionT()[ivLocalIdx];
        }

        /*!
         * \brief Update cached objects (transmissibilities).
         *        This is used for updates with secondary interaction volumes.
         *
         * \param iv The interaction volume this scvf is embedded in
         * \param localFaceData iv-local info on this scvf
         * \param dataHandle Transmissibility matrix & gravity data of this iv
         * \param scvf The sub-control volume face
         */
        template< bool doSecondary = considerSecondaryIVs, std::enable_if_t<doSecondary, int > = 0 >
        void updateHeatConduction(const SecondaryInteractionVolume& iv,
                                  const SecondaryIvLocalFaceData& localFaceData,
                                  const SecondaryIvDataHandle& dataHandle,
                                  const SubControlVolumeFace &scvf)
        {
            stencil_ = &iv.stencil();
            switchFluxSign_ = localFaceData.isOutside();

            // store pointer to the temperature vector of this iv
            secondaryTj_ = &dataHandle.temperatures();

            const auto ivLocalIdx = localFaceData.ivLocalScvfIndex();
            if (dim == dimWorld)
                secondaryTij_ = &dataHandle.heatConductionT()[ivLocalIdx];
            else
                secondaryTij_ = localFaceData.isOutside() ? &dataHandle.heatConductionTout()[ivLocalIdx][localFaceData.scvfLocalOutsideScvfIndex()]
                                                          : &dataHandle.heatConductionT()[ivLocalIdx];
        }

        //! Coefficients for the cell (& Dirichlet) unknowns in flux expressions (primary type)
        const PrimaryIvTij& heatConductionTijPrimaryIv() const { return *primaryTij_; }

        //! Coefficients for the cell (& Dirichlet) unknowns in flux expressions (secondary type)
        const SecondaryIvTij& heatConductionTijSecondaryIv() const { return *secondaryTij_; }

        //! The stencil corresponding to the transmissibilities (primary type)
        const Stencil& heatConductionStencil() const { return *stencil_; }

        //! The cell (& Dirichlet) temperatures within this interaction volume (primary type)
        const PrimaryIvVector& temperaturesPrimaryIv() const { return *primaryTj_; }

        //! The cell (& Dirichlet) temperatures within this interaction volume (secondary type)
        const SecondaryIvVector& temperaturesSecondaryIv() const { return *secondaryTj_; }

        //! In the interaction volume-local system of eq we have one unknown per face.
        //! On scvfs on this face, but in "outside" (neighbor) elements of it, we have
        //! to take the negative value of the fluxes due to the flipped normal vector.
        //! This function returns whether or not this scvf is an "outside" face in the iv.
        bool heatConductionSwitchFluxSign() const { return switchFluxSign_; }

    private:
        bool switchFluxSign_;

        //! The stencil, i.e. the grid indices j
        const Stencil* stencil_;

        //! The transmissibilities such that f = Tij*Tj
        const PrimaryIvTij* primaryTij_;
        const SecondaryIvTij* secondaryTij_;

        //! The interaction-volume wide temperature Tj
        const PrimaryIvVector* primaryTj_;
        const SecondaryIvVector* secondaryTj_;
    };

public:
    // state the discretization method this implementation belongs to
    static const DiscretizationMethods myDiscretizationMethod = DiscretizationMethods::CCMpfa;

    // state the type for the corresponding cache and its filler
    using Cache = MpfaFouriersLawCache;

    //! Compute the conductive flux across an scvf
    static Scalar flux(const Problem& problem,
                       const Element& element,
                       const FVElementGeometry& fvGeometry,
                       const ElementVolumeVariables& elemVolVars,
                       const SubControlVolumeFace& scvf,
                       const ElementFluxVarsCache& elemFluxVarsCache)
    {
        const auto& fluxVarsCache = elemFluxVarsCache[scvf];

        // compute Tij*tj
        Scalar flux;
        if (fluxVarsCache.usesSecondaryIv())
        {
            const auto& tij = fluxVarsCache.heatConductionTijSecondaryIv();
            const auto& Tj = fluxVarsCache.temperaturesSecondaryIv();
            flux = tij*Tj;
        }
        else
        {
            const auto& tij = fluxVarsCache.heatConductionTijPrimaryIv();
            const auto& Tj = fluxVarsCache.temperaturesPrimaryIv();
            flux = tij*Tj;
        }

        if (fluxVarsCache.heatConductionSwitchFluxSign())
            flux *= -1.0;

        return flux;
    }
};

} // end namespace Dumux

#endif
