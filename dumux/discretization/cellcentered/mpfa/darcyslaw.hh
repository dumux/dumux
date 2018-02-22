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
 * \brief Darcy's Law for cell-centered finite volume schemes
 *        with multi-point flux approximation.
 */
#ifndef DUMUX_DISCRETIZATION_CC_MPFA_DARCYS_LAW_HH
#define DUMUX_DISCRETIZATION_CC_MPFA_DARCYS_LAW_HH

#include <dune/common/dynvector.hh>
#include <dune/common/dynmatrix.hh>

#include <dumux/common/properties.hh>
#include <dumux/common/parameters.hh>
#include <dumux/discretization/methods.hh>

namespace Dumux
{
//! forward declaration of the method-specific implementation
template<class TypeTag, DiscretizationMethod discMethod>
class DarcysLawImplementation;

/*!
 * \ingroup CCMpfaDiscretization
 * \brief Darcy's law for cell-centered finite volume schemes with multi-point flux approximation.
 */
template<class TypeTag>
class DarcysLawImplementation<TypeTag, DiscretizationMethod::ccmpfa>
{
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using Problem = typename GET_PROP_TYPE(TypeTag, Problem);
    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
    using Element = typename GridView::template Codim<0>::Entity;

    using FVGridGeometry = typename GET_PROP_TYPE(TypeTag, FVGridGeometry);
    using FVElementGeometry = typename FVGridGeometry::LocalView;
    using SubControlVolumeFace = typename FVElementGeometry::SubControlVolumeFace;
    using ElementVolumeVariables = typename GET_PROP_TYPE(TypeTag, ElementVolumeVariables);
    using ElementFluxVariablesCache = typename GET_PROP_TYPE(TypeTag, ElementFluxVariablesCache);
    using FluxVariablesCache = typename GET_PROP_TYPE(TypeTag, FluxVariablesCache);

    //! Class that fills the cache corresponding to mpfa Darcy's Law
    class MpfaDarcysLawCacheFiller
    {
    public:
        //! Function to fill an MpfaDarcysLawCache of a given scvf
        //! This interface has to be met by any advection-related cache filler class
        template<class FluxVariablesCacheFiller>
        static void fill(FluxVariablesCache& scvfFluxVarsCache,
                         const Problem& problem,
                         const Element& element,
                         const FVElementGeometry& fvGeometry,
                         const ElementVolumeVariables& elemVolVars,
                         const SubControlVolumeFace& scvf,
                         const FluxVariablesCacheFiller& fluxVarsCacheFiller)
        {
            // get interaction volume related data from the filler class & upate the cache
            if (fvGeometry.fvGridGeometry().vertexUsesSecondaryInteractionVolume(scvf.vertexIndex()))
                scvfFluxVarsCache.updateAdvection(fluxVarsCacheFiller.secondaryInteractionVolume(),
                                                  fluxVarsCacheFiller.secondaryIvLocalFaceData(),
                                                  fluxVarsCacheFiller.secondaryIvDataHandle(),
                                                  scvf);
            else
                scvfFluxVarsCache.updateAdvection(fluxVarsCacheFiller.primaryInteractionVolume(),
                                                  fluxVarsCacheFiller.primaryIvLocalFaceData(),
                                                  fluxVarsCacheFiller.primaryIvDataHandle(),
                                                  scvf);
        }
    };

    //! The cache used in conjunction with the mpfa Darcy's Law
    class MpfaDarcysLawCache
    {
        static constexpr int dim = GridView::dimension;
        static constexpr int dimWorld = GridView::dimensionworld;
        static constexpr int numPhases = GET_PROP_VALUE(TypeTag, NumPhases);

        using DualGridNodalIndexSet = typename GET_PROP_TYPE(TypeTag, DualGridNodalIndexSet);
        using Stencil = typename DualGridNodalIndexSet::NodalGridStencilType;

        using MpfaHelper = typename FVGridGeometry::MpfaHelper;
        static constexpr bool considerSecondaryIVs = MpfaHelper::considerSecondaryIVs();

        using PrimaryInteractionVolume = typename GET_PROP_TYPE(TypeTag, PrimaryInteractionVolume);
        using PrimaryIvLocalFaceData = typename PrimaryInteractionVolume::Traits::LocalFaceData;
        using PrimaryIvDataHandle = typename ElementFluxVariablesCache::PrimaryIvDataHandle;
        using PrimaryIvCellVector = typename PrimaryInteractionVolume::Traits::MatVecTraits::CellVector;
        using PrimaryIvTij = typename PrimaryInteractionVolume::Traits::MatVecTraits::TMatrix::row_type;

        using SecondaryInteractionVolume = typename GET_PROP_TYPE(TypeTag, SecondaryInteractionVolume);
        using SecondaryIvLocalFaceData = typename SecondaryInteractionVolume::Traits::LocalFaceData;
        using SecondaryIvDataHandle = typename ElementFluxVariablesCache::SecondaryIvDataHandle;
        using SecondaryIvCellVector = typename SecondaryInteractionVolume::Traits::MatVecTraits::CellVector;
        using SecondaryIvTij = typename SecondaryInteractionVolume::Traits::MatVecTraits::TMatrix::row_type;

    public:
        // export the filler type
        using Filler = MpfaDarcysLawCacheFiller;

        /*!
         * \brief Update cached objects (transmissibilities and gravity).
         *        This is used for updates with primary interaction volumes.
         *
         * \param iv The interaction volume this scvf is embedded in
         * \param localFaceData iv-local info on this scvf
         * \param dataHandle Transmissibility matrix & gravity data of this iv
         * \param scvf The sub-control volume face
         */
        void updateAdvection(const PrimaryInteractionVolume& iv,
                             const PrimaryIvLocalFaceData& localFaceData,
                             const PrimaryIvDataHandle& dataHandle,
                             const SubControlVolumeFace &scvf)
        {
            switchFluxSign_ = localFaceData.isOutside();
            stencil_ = &iv.stencil();

            for (unsigned int pIdx = 0; pIdx < numPhases; ++pIdx)
                primaryPj_[pIdx] = &dataHandle.pressures(pIdx);

            static const bool enableGravity = getParamFromGroup<bool>(GET_PROP_VALUE(TypeTag, ModelParameterGroup), "Problem.EnableGravity");
            const auto ivLocalIdx = localFaceData.ivLocalScvfIndex();

            // standard grids
            if (dim == dimWorld)
            {
                primaryTij_ = &dataHandle.advectionT()[ivLocalIdx];
                if (enableGravity)
                    for (unsigned int phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx)
                        g_[phaseIdx] = dataHandle.gravity(phaseIdx)[ivLocalIdx];
            }
            // surface grids
            else
            {
                if (!localFaceData.isOutside())
                {
                    primaryTij_ = &dataHandle.advectionT()[ivLocalIdx];
                    if (enableGravity)
                        for (unsigned int phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx)
                            g_[phaseIdx] = dataHandle.gravity(phaseIdx)[ivLocalIdx];
                }
                else
                {
                    const auto idxInOutsideFaces = localFaceData.scvfLocalOutsideScvfIndex();
                    primaryTij_ = &dataHandle.advectionTout()[ivLocalIdx][idxInOutsideFaces];
                    if (enableGravity)
                        for (unsigned int phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx)
                            g_[phaseIdx] = dataHandle.gravityOutside(phaseIdx)[ivLocalIdx][idxInOutsideFaces];
                }
            }
        }

        /*!
         * \brief Update cached objects (transmissibilities and gravity).
         *        This is used for updates with secondary interaction volumes.
         *
         * \param iv The interaction volume this scvf is embedded in
         * \param localFaceData iv-local info on this scvf
         * \param dataHandle Transmissibility matrix & gravity data of this iv
         * \param scvf The sub-control volume face
         */
        template< bool doSecondary = considerSecondaryIVs, std::enable_if_t<doSecondary, int > = 0 >
        void updateAdvection(const SecondaryInteractionVolume& iv,
                             const SecondaryIvLocalFaceData& localFaceData,
                             const SecondaryIvDataHandle& dataHandle,
                             const SubControlVolumeFace &scvf)
        {
            switchFluxSign_ = localFaceData.isOutside();
            stencil_ = &iv.stencil();

            for (unsigned int pIdx = 0; pIdx < numPhases; ++pIdx)
                secondaryPj_[pIdx] = &dataHandle.pressures(pIdx);

            static const bool enableGravity = getParamFromGroup<bool>(GET_PROP_VALUE(TypeTag, ModelParameterGroup), "Problem.EnableGravity");
            const auto ivLocalIdx = localFaceData.ivLocalScvfIndex();

            // standard grids
            if (dim == dimWorld)
            {
                secondaryTij_ = &dataHandle.advectionT()[ivLocalIdx];
                if (enableGravity)
                    for (unsigned int phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx)
                        g_[phaseIdx] = dataHandle.gravity(phaseIdx)[ivLocalIdx];
            }
            // surface grids
            else
            {
                if (!localFaceData.isOutside())
                {
                    secondaryTij_ = &dataHandle.advectionT()[ivLocalIdx];
                    if (enableGravity)
                        for (unsigned int phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx)
                            g_[phaseIdx] = dataHandle.gravity(phaseIdx)[ivLocalIdx];
                }
                else
                {
                    const auto idxInOutsideFaces = localFaceData.scvfLocalOutsideScvfIndex();
                    secondaryTij_ = &dataHandle.advectionTout()[ivLocalIdx][idxInOutsideFaces];
                    if (enableGravity)
                        for (unsigned int phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx)
                            g_[phaseIdx] = dataHandle.gravityOutside(phaseIdx)[ivLocalIdx][idxInOutsideFaces];
                }
            }
        }

        //! The stencil corresponding to the transmissibilities (primary type)
        const Stencil& advectionStencil() const { return *stencil_; }

        //! Coefficients for the cell (& Dirichlet) unknowns in flux expressions (primary type)
        const PrimaryIvTij& advectionTijPrimaryIv() const { return *primaryTij_; }

        //! Coefficients for the cell (& Dirichlet) unknowns in flux expressions (secondary type)
        const SecondaryIvTij& advectionTijSecondaryIv() const { return *secondaryTij_; }

        //! The cell (& maybe Dirichlet) pressures within this interaction volume (primary type)
        const PrimaryIvCellVector& pressuresPrimaryIv(unsigned int phaseIdx) const { return *primaryPj_[phaseIdx]; }

        //! The cell (& maybe Dirichlet) pressures within this interaction volume (secondary type)
        const SecondaryIvCellVector& pressuresSecondaryIv(unsigned int phaseIdx) const { return *secondaryPj_[phaseIdx]; }

        //! The gravitational acceleration for a phase on this scvf
        Scalar gravity(unsigned int phaseIdx) const { return g_[phaseIdx]; }

        //! In the interaction volume-local system of eq we have one unknown per face.
        //! On scvfs on this face, but in "outside" (neighbor) elements of it, we have
        //! to take the negative value of the fluxes due to the flipped normal vector.
        //! This function returns whether or not this scvf is an "outside" face in the iv.
        bool advectionSwitchFluxSign() const { return switchFluxSign_; }

    private:
        bool switchFluxSign_;

        //! The stencil, i.e. the grid indices j
        const Stencil* stencil_;

        //! The transmissibilities such that f = Tij*pj
        const PrimaryIvTij* primaryTij_;
        const SecondaryIvTij* secondaryTij_;

        //! The interaction-volume wide phase pressures pj
        std::array<const PrimaryIvCellVector*, numPhases> primaryPj_;
        std::array<const SecondaryIvCellVector*, numPhases> secondaryPj_;

        //! Gravitational flux contribution on this face
        std::array< Scalar, numPhases > g_;
    };

public:
    // state the discretization method this implementation belongs to
    static const DiscretizationMethod discMethod = DiscretizationMethod::ccmpfa;

    // export the type for the corresponding cache
    using Cache = MpfaDarcysLawCache;

    //! Compute the advective flux across an scvf
    static Scalar flux(const Problem& problem,
                       const Element& element,
                       const FVElementGeometry& fvGeometry,
                       const ElementVolumeVariables& elemVolVars,
                       const SubControlVolumeFace& scvf,
                       const unsigned int phaseIdx,
                       const ElementFluxVariablesCache& elemFluxVarsCache)
    {
        const auto& fluxVarsCache = elemFluxVarsCache[scvf];

        // compute t_ij*p_j
        Scalar scvfFlux;
        if (fluxVarsCache.usesSecondaryIv())
        {
            const auto& tij = fluxVarsCache.advectionTijSecondaryIv();
            const auto& pj = fluxVarsCache.pressuresSecondaryIv(phaseIdx);
            scvfFlux = tij*pj;
        }
        else
        {
            const auto& tij = fluxVarsCache.advectionTijPrimaryIv();
            const auto& pj = fluxVarsCache.pressuresPrimaryIv(phaseIdx);
            scvfFlux = tij*pj;
        }

        // maybe add gravitational acceleration
        static const bool enableGravity = getParamFromGroup<bool>(GET_PROP_VALUE(TypeTag, ModelParameterGroup), "Problem.EnableGravity");
        if (enableGravity)
            scvfFlux += fluxVarsCache.gravity(phaseIdx);

        // switch the sign if necessary
        if (fluxVarsCache.advectionSwitchFluxSign())
            scvfFlux *= -1.0;

        return scvfFlux;
    }
};

} // end namespace

#endif
