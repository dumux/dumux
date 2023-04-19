// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
* \file
* \ingroup CCTpfaFlux
* \brief Fourier's law for cell-centered finite volume schemes with two-point flux approximation
*/
#ifndef DUMUX_DISCRETIZATION_CC_TPFA_FOURIERS_LAW_HH
#define DUMUX_DISCRETIZATION_CC_TPFA_FOURIERS_LAW_HH

#include <dumux/common/parameters.hh>
#include <dumux/common/properties.hh>

#include <dumux/discretization/method.hh>
#include <dumux/discretization/extrusion.hh>
#include <dumux/discretization/cellcentered/tpfa/computetransmissibility.hh>

namespace Dumux {

// forward declaration
template<class TypeTag, class DiscretizationMethod>
class FouriersLawImplementation;

/*!
 * \ingroup CCTpfaFlux
 * \brief Fourier's law for cell-centered finite volume schemes with two-point flux approximation
 */
template <class TypeTag>
class FouriersLawImplementation<TypeTag, DiscretizationMethods::CCTpfa>
{
    using Implementation = FouriersLawImplementation<TypeTag, DiscretizationMethods::CCTpfa>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using Problem = GetPropType<TypeTag, Properties::Problem>;
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using FVElementGeometry = typename GridGeometry::LocalView;
    using SubControlVolumeFace = typename GridGeometry::SubControlVolumeFace;
    using Extrusion = Extrusion_t<GridGeometry>;
    using GridView = typename GridGeometry::GridView;
    using ElementVolumeVariables = typename GetPropType<TypeTag, Properties::GridVolumeVariables>::LocalView;
    using VolumeVariables = typename ElementVolumeVariables::VolumeVariables;
    using Element = typename GridView::template Codim<0>::Entity;
    using GridFluxVariablesCache = GetPropType<TypeTag, Properties::GridFluxVariablesCache>;
    using ElementFluxVarsCache = typename GridFluxVariablesCache::LocalView;
    using FluxVariablesCache = typename GridFluxVariablesCache::FluxVariablesCache;

    static const int dim = GridView::dimension;
    static const int dimWorld = GridView::dimensionworld;

    //! Class that fills the cache corresponding to tpfa Fick's Law
    class TpfaFouriersLawCacheFiller
    {
    public:
        //! Function to fill a TpfaFicksLawCache of a given scvf
        //! This interface has to be met by any diffusion-related cache filler class
        template<class FluxVariablesCacheFiller>
        static void fill(FluxVariablesCache& scvfFluxVarsCache,
                         const Problem& problem,
                         const Element& element,
                         const FVElementGeometry& fvGeometry,
                         const ElementVolumeVariables& elemVolVars,
                         const SubControlVolumeFace& scvf,
                         const FluxVariablesCacheFiller& fluxVarsCacheFiller)
        {
            scvfFluxVarsCache.updateHeatConduction(problem, element, fvGeometry, elemVolVars, scvf);
        }
    };

    //! Class that caches the transmissibility
    class TpfaFouriersLawCache
    {
    public:
        using Filler = TpfaFouriersLawCacheFiller;

        void updateHeatConduction(const Problem& problem,
                                  const Element& element,
                                  const FVElementGeometry& fvGeometry,
                                  const ElementVolumeVariables& elemVolVars,
                                  const SubControlVolumeFace &scvf)
        {
            tij_ = Implementation::calculateTransmissibility(problem, element, fvGeometry, elemVolVars, scvf);
        }

        const Scalar& heatConductionTij() const
        { return tij_; }

    private:
        Scalar tij_;
    };

public:
    using DiscretizationMethod = DiscretizationMethods::CCTpfa;
    //! state the discretization method this implementation belongs to
    static constexpr DiscretizationMethod discMethod{};

    //! export the type for the corresponding cache
    using Cache = TpfaFouriersLawCache;

    /*!
     * \brief Returns the heat flux within the porous medium
     *        (in J/s) across the given sub-control volume face.
     * \note This law assumes thermal equilibrium between the fluid
     *       and solid phases, and uses an effective thermal conductivity
     *       for the overall aggregate.
     *       This overload allows to explicitly specify the inside and outside volume variables
     *       which can be useful to evaluate conductive fluxes at boundaries with given outside values.
     *       This only works if scvf.numOutsideScv() == 1.
     *
     */
    static Scalar flux(const Problem& problem,
                       const Element& element,
                       const FVElementGeometry& fvGeometry,
                       const VolumeVariables& insideVolVars,
                       const VolumeVariables& outsideVolVars,
                       const SubControlVolumeFace& scvf,
                       const ElementFluxVarsCache& elemFluxVarsCache)
    {
        if constexpr (isMixedDimensional_)
            if (scvf.numOutsideScv() != 1)
                DUNE_THROW(Dune::Exception, "This flux overload requires scvf.numOutsideScv() == 1");

        // heat conductivities are always solution dependent (?)
        Scalar tij = elemFluxVarsCache[scvf].heatConductionTij();

        // get the inside/outside temperatures
        const auto tInside = insideVolVars.temperature();
        const auto tOutside = outsideVolVars.temperature();

        return tij*(tInside - tOutside);
    }

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
        // heat conductivities are always solution dependent (?)
        Scalar tij = elemFluxVarsCache[scvf].heatConductionTij();

        // get the inside/outside temperatures
        const auto tInside = elemVolVars[scvf.insideScvIdx()].temperature();
        const auto tOutside = scvf.numOutsideScvs() == 1 ? elemVolVars[scvf.outsideScvIdx()].temperature()
                              : branchingFacetTemperature_(problem, element, fvGeometry, elemVolVars, elemFluxVarsCache, scvf, tInside, tij);

        return tij*(tInside - tOutside);
    }

    //! Compute transmissibilities
    static Scalar calculateTransmissibility(const Problem& problem,
                                            const Element& element,
                                            const FVElementGeometry& fvGeometry,
                                            const ElementVolumeVariables& elemVolVars,
                                            const SubControlVolumeFace& scvf)
    {
        Scalar tij;

        const auto insideScvIdx = scvf.insideScvIdx();
        const auto& insideScv = fvGeometry.scv(insideScvIdx);
        const auto& insideVolVars = elemVolVars[insideScvIdx];

        const auto insideLambda = insideVolVars.effectiveThermalConductivity();
        const Scalar ti = computeTpfaTransmissibility(fvGeometry, scvf, insideScv, insideLambda, insideVolVars.extrusionFactor());

        // for the boundary (dirichlet) or at branching points we only need ti
        if (scvf.boundary() || scvf.numOutsideScvs() > 1)
        {
            tij = Extrusion::area(fvGeometry, scvf)*ti;
        }
        // otherwise we compute a tpfa harmonic mean
        else
        {
            const auto outsideScvIdx = scvf.outsideScvIdx();
            const auto& outsideScv = fvGeometry.scv(outsideScvIdx);
            const auto& outsideVolVars = elemVolVars[outsideScvIdx];

            const auto outsideLambda = outsideVolVars.effectiveThermalConductivity();
            Scalar tj;
            if constexpr (dim == dimWorld)
                // assume the normal vector from outside is anti parallel so we save flipping a vector
                tj = -1.0*computeTpfaTransmissibility(fvGeometry, scvf, outsideScv, outsideLambda, outsideVolVars.extrusionFactor());
            else
                tj = computeTpfaTransmissibility(fvGeometry, fvGeometry.flipScvf(scvf.index()), outsideScv, outsideLambda, outsideVolVars.extrusionFactor());

            // check for division by zero!
            if (ti*tj <= 0.0)
                tij = 0;
            else
                tij = Extrusion::area(fvGeometry, scvf)*(ti * tj)/(ti + tj);
        }

        return tij;
    }

private:

    //! compute the temperature at branching facets for network grids
    static Scalar branchingFacetTemperature_(const Problem& problem,
                                             const Element& element,
                                             const FVElementGeometry& fvGeometry,
                                             const ElementVolumeVariables& elemVolVars,
                                             const ElementFluxVarsCache& elemFluxVarsCache,
                                             const SubControlVolumeFace& scvf,
                                             Scalar insideTemperature,
                                             Scalar insideTi)
    {
        Scalar sumTi(insideTi);
        Scalar sumTempTi(insideTi*insideTemperature);

        for (unsigned int i = 0; i < scvf.numOutsideScvs(); ++i)
        {
            const auto outsideScvIdx = scvf.outsideScvIdx(i);
            const auto& outsideVolVars = elemVolVars[outsideScvIdx];
            const auto& flippedScvf = fvGeometry.flipScvf(scvf.index(), i);
            const auto& outsideFluxVarsCache = elemFluxVarsCache[flippedScvf];

            auto outsideTi = outsideFluxVarsCache.heatConductionTij();
            sumTi += outsideTi;
            sumTempTi += outsideTi*outsideVolVars.temperature();
        }
        return sumTempTi/sumTi;
    }

    static constexpr bool isMixedDimensional_ = static_cast<int>(GridView::dimension) < static_cast<int>(GridView::dimensionworld);
};

} // end namespace Dumux

#endif
