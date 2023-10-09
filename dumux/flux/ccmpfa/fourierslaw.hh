// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
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
#include <dumux/discretization/cellcentered/mpfa/dualgridindexset.hh>

namespace Dumux {

//! forward declaration of the method-specific implementation
template<class TypeTag, class DiscretizationMethod>
class FouriersLawImplementation;

/*!
 * \ingroup CCMpfaFlux
 * \brief Fourier's law for cell-centered finite volume schemes with multi-point flux approximation
 */
template <class TypeTag>
class FouriersLawImplementation<TypeTag, DiscretizationMethods::CCMpfa>
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
    using GridFluxVariablesCache = GetPropType<TypeTag, Properties::GridFluxVariablesCache>;
    using ElementFluxVarsCache = typename GridFluxVariablesCache::LocalView;
    using FluxVariablesCache = typename GridFluxVariablesCache::FluxVariablesCache;

    //! Class that fills the cache corresponding to mpfa Darcy's Law
    class MpfaFouriersLawCacheFiller
    {
    public:
        //! Function to fill an MpfaDarcysLawCache of a given scvf
        //! This interface has to be met by any cache filler class for heat conduction quantities
        template<class FluxVariablesCacheFiller>
        static void fill(FluxVariablesCache& scvfFluxVarsCache,
                         const Problem&,
                         const Element&,
                         const FVElementGeometry&,
                         const ElementVolumeVariables&,
                         const SubControlVolumeFace&,
                         const FluxVariablesCacheFiller& fluxVarsCacheFiller)
        {
            scvfFluxVarsCache.updateHeatConduction(
                fluxVarsCacheFiller.fluxesPtr(),
                fluxVarsCacheFiller.heatConductionFluxId()
            );
        }
    };

    //! The cache used in conjunction with the mpfa Fourier's Law
    class MpfaFouriersLawCache
    {
        using Fluxes = CCMpfaFluxes<GridGeometry, Scalar>;
        using FluxId = typename Fluxes::FluxId;

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
         void updateHeatConduction(const Fluxes* fluxesPtr, FluxId id)
        {
            fluxes_ = fluxesPtr;
            id_ = std::move(id);
        }

        const Fluxes& heatConductionFluxes() const
        { return *fluxes_; }

        FluxId heatConductionId() const
        { return id_; }

    private:
        const Fluxes* fluxes_;
        FluxId id_;
    };

public:
    using DiscretizationMethod = DiscretizationMethods::CCMpfa;
    // state the discretization method this implementation belongs to
    static constexpr DiscretizationMethod discMethod{};

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
        return fluxVarsCache.heatConductionFluxes().computeFluxFor(
            fluxVarsCache.heatConductionId(),
            scvf
        );
    }
};

} // end namespace Dumux

#endif
