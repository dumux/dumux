// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup CCMpfaFlux
 * \brief Darcy's Law for cell-centered finite volume schemes
 *        with multi-point flux approximation.
 */
#ifndef DUMUX_DISCRETIZATION_CC_MPFA_DARCYS_LAW_HH
#define DUMUX_DISCRETIZATION_CC_MPFA_DARCYS_LAW_HH

#include <dune/common/dynvector.hh>
#include <dune/common/dynmatrix.hh>

#include <dumux/common/properties.hh>
#include <dumux/common/parameters.hh>

#include <dumux/discretization/method.hh>
#include <dumux/discretization/cellcentered/mpfa/dualgridindexset.hh>
#include <dumux/discretization/cellcentered/mpfa/interactionvolume.hh>

namespace Dumux {

//! forward declaration of the method-specific implementation
template<class TypeTag, class DiscretizationMethod>
class DarcysLawImplementation;

/*!
 * \ingroup CCMpfaFlux
 * \brief Darcy's law for cell-centered finite volume schemes
 *        with multi-point flux approximation.
 */
template<class TypeTag>
class DarcysLawImplementation<TypeTag, DiscretizationMethods::CCMpfa>
{
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using Problem = GetPropType<TypeTag, Properties::Problem>;
    using GridView = typename GetPropType<TypeTag, Properties::GridGeometry>::GridView;
    using Element = typename GridView::template Codim<0>::Entity;

    static constexpr int dim = GridView::dimension;
    static constexpr int dimWorld = GridView::dimensionworld;
    static constexpr int numPhases = GetPropType<TypeTag, Properties::ModelTraits>::numFluidPhases();

    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using FVElementGeometry = typename GridGeometry::LocalView;
    using SubControlVolumeFace = typename FVElementGeometry::SubControlVolumeFace;
    using ElementVolumeVariables = typename GetPropType<TypeTag, Properties::GridVolumeVariables>::LocalView;

    //! Class that fills the cache corresponding to mpfa Darcy's Law
    class MpfaDarcysLawCacheFiller
    {
    public:
        //! Function to fill an MpfaDarcysLawCache of a given scvf
        //! This interface has to be met by any advection-related cache filler class
        template<class FluxVariablesCache, class FluxVariablesCacheFiller>
        static void fill(FluxVariablesCache& scvfFluxVarsCache,
                         const Problem&,
                         const Element&,
                         const FVElementGeometry&,
                         const ElementVolumeVariables&,
                         const SubControlVolumeFace&,
                         const FluxVariablesCacheFiller& fluxVarsCacheFiller)
        {
            scvfFluxVarsCache.updateAdvection(
                fluxVarsCacheFiller.fluxesPtr(),
                fluxVarsCacheFiller.advectionFluxIds()
            );
        }
    };

    //! The cache used in conjunction with the mpfa Darcy's Law
    class MpfaDarcysLawCache
    {
        using Fluxes = CCMpfaFluxes<GridGeometry, Scalar>;
        using FluxId = typename Fluxes::FluxId;

    public:
        using Filler = MpfaDarcysLawCacheFiller;

        /*!
         * \brief Update cached objects (i.e. ptrs to flux computation instance).
         * \param fluxes The flux computation instance
         * \param ids Flux ids for advective fluxes
         */
        void updateAdvection(const Fluxes* fluxes, const std::array<FluxId, numPhases>& ids)
        {
            fluxes_ = fluxes;
            phaseFluxIds_ = ids;
        }

        const Fluxes& advectionFluxes() const
        { return *fluxes_; }

        FluxId advectionId(int phaseIdx) const
        { return phaseFluxIds_[phaseIdx]; }

    private:
        const Fluxes* fluxes_{nullptr};
        std::array<FluxId, numPhases> phaseFluxIds_;
    };

public:
    using DiscretizationMethod = DiscretizationMethods::CCMpfa;
    // state the discretization method this implementation belongs to
    static constexpr DiscretizationMethod discMethod{};

    // export the type for the corresponding cache
    using Cache = MpfaDarcysLawCache;

    /*!
     * \brief Returns the advective flux of a fluid phase
     *        across the given sub-control volume face.
     * \note This assembles the term
     *       \f$-|\sigma| \mathbf{n}^T \mathbf{K} \left( \nabla p - \rho \mathbf{g} \right)\f$,
     *       where \f$|\sigma|\f$ is the area of the face and \f$\mathbf{n}\f$ is the outer
     *       normal vector. Thus, the flux is given in N*m, and can be converted
     *       into a volume flux (m^3/s) or mass flux (kg/s) by applying an upwind scheme
     *       for the mobility or the product of density and mobility, respectively.
     */
    template<class ElementFluxVariablesCache>
    static Scalar flux(const Problem& problem,
                       const Element& element,
                       const FVElementGeometry& fvGeometry,
                       const ElementVolumeVariables& elemVolVars,
                       const SubControlVolumeFace& scvf,
                       const unsigned int phaseIdx,
                       const ElementFluxVariablesCache& elemFluxVarsCache)
    {
        const auto& fluxVarsCache = elemFluxVarsCache[scvf];
        return fluxVarsCache.advectionFluxes().computeFluxFor(
            fluxVarsCache.advectionId(phaseIdx),
            scvf
        );
    }
};

} // end namespace

#endif
