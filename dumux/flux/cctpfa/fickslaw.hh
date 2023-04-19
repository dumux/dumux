// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup CCTpfaFlux
 * \brief Fick's law for cell-centered finite volume schemes with two-point flux approximation
 */
#ifndef DUMUX_DISCRETIZATION_CC_TPFA_FICKS_LAW_HH
#define DUMUX_DISCRETIZATION_CC_TPFA_FICKS_LAW_HH

#include <dune/common/fvector.hh>

#include <dumux/common/parameters.hh>
#include <dumux/common/properties.hh>

#include <dumux/discretization/method.hh>
#include <dumux/discretization/extrusion.hh>
#include <dumux/discretization/cellcentered/tpfa/computetransmissibility.hh>
#include <dumux/flux/fickiandiffusioncoefficients.hh>

#include <dumux/flux/referencesystemformulation.hh>

namespace Dumux {

// forward declaration
template<class TypeTag, class DiscretizationMethod, ReferenceSystemFormulation referenceSystem>
class FicksLawImplementation;

/*!
 * \ingroup CCTpfaFlux
 * \brief Fick's law for cell-centered finite volume schemes with two-point flux approximation
 */
template <class TypeTag, ReferenceSystemFormulation referenceSystem>
class FicksLawImplementation<TypeTag, DiscretizationMethods::CCTpfa, referenceSystem>
{
    using Implementation = FicksLawImplementation<TypeTag, DiscretizationMethods::CCTpfa, referenceSystem>;
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
    using ElementFluxVariablesCache = typename GridFluxVariablesCache::LocalView;
    using FluxVariablesCache = typename GridFluxVariablesCache::FluxVariablesCache;
    using BalanceEqOpts = GetPropType<TypeTag, Properties::BalanceEqOpts>;

    using ModelTraits = GetPropType<TypeTag, Properties::ModelTraits>;
    static const int dim = GridView::dimension;
    static const int dimWorld = GridView::dimensionworld;
    static const int numPhases = ModelTraits::numFluidPhases();
    static const int numComponents = ModelTraits::numFluidComponents();

    using ComponentFluxVector = Dune::FieldVector<Scalar, numComponents>;

    //! Class that fills the cache corresponding to tpfa Fick's Law
    class TpfaFicksLawCacheFiller
    {
    public:
        //! Function to fill a TpfaFicksLawCache of a given scvf
        //! This interface has to be met by any diffusion-related cache filler class
        template<class FluxVariablesCacheFiller>
        static void fill(FluxVariablesCache& scvfFluxVarsCache,
                         unsigned int phaseIdx, unsigned int compIdx,
                         const Problem& problem,
                         const Element& element,
                         const FVElementGeometry& fvGeometry,
                         const ElementVolumeVariables& elemVolVars,
                         const SubControlVolumeFace& scvf,
                         const FluxVariablesCacheFiller& fluxVarsCacheFiller)
        {
            scvfFluxVarsCache.updateDiffusion(problem, element, fvGeometry, elemVolVars, scvf, phaseIdx, compIdx);
        }
    };

    //! Class that caches the transmissibility
    class TpfaFicksLawCache
    {
    public:
        using Filler = TpfaFicksLawCacheFiller;

        void updateDiffusion(const Problem& problem,
                             const Element& element,
                             const FVElementGeometry& fvGeometry,
                             const ElementVolumeVariables& elemVolVars,
                             const SubControlVolumeFace &scvf,
                             const unsigned int phaseIdx,
                             const unsigned int compIdx)
        {
            tij_[phaseIdx][compIdx] = Implementation::calculateTransmissibility(problem, element, fvGeometry, elemVolVars, scvf, phaseIdx, compIdx);
        }

        const Scalar& diffusionTij(unsigned int phaseIdx, unsigned int compIdx) const
        { return tij_[phaseIdx][compIdx]; }

    private:
        std::array< std::array<Scalar, numComponents>, numPhases> tij_;
    };

public:
    using DiscretizationMethod = DiscretizationMethods::CCTpfa;
    //! state the discretization method this implementation belongs to
    static constexpr DiscretizationMethod discMethod{};

    //! Return the reference system
    static constexpr ReferenceSystemFormulation referenceSystemFormulation()
    { return referenceSystem; }

    //! state the type for the corresponding cache and its filler
    using Cache = TpfaFicksLawCache;

    using DiffusionCoefficientsContainer = FickianDiffusionCoefficients<Scalar, numPhases, numComponents>;

    /*!
     * \brief Returns the diffusive fluxes of all components within
     *        a fluid phase across the given sub-control volume face.
     *        The computed fluxes are given in mole/s or kg/s, depending
     *        on the template parameter ReferenceSystemFormulation.
     *
     * \note This overload allows to explicitly specify the inside and outside volume variables
     *       which can be useful to evaluate diffusive fluxes at boundaries with given outside values.
     *       This only works if scvf.numOutsideScv() == 1.
     */
    static ComponentFluxVector flux(const Problem& problem,
                                    const Element& element,
                                    const FVElementGeometry& fvGeometry,
                                    const VolumeVariables& insideVolVars,
                                    const VolumeVariables& outsideVolVars,
                                    const SubControlVolumeFace& scvf,
                                    const int phaseIdx,
                                    const ElementFluxVariablesCache& elemFluxVarsCache)
    {
        if constexpr (isMixedDimensional_)
            if (scvf.numOutsideScv() != 1)
                DUNE_THROW(Dune::Exception, "This flux overload requires scvf.numOutsideScv() == 1");

        // helper lambda to get the outside mole or mass fraction
        const auto getOutsideX = [&](const Scalar xInside, const Scalar tij, const int phaseIdx, const int compIdx)
        {
            return massOrMoleFraction(outsideVolVars, referenceSystem, phaseIdx, compIdx);
        };

        // helper lambda to get the averaged density at the scvf
        const auto getRho = [&](const int phaseIdx, const Scalar rhoInside, const Scalar rhoOutside)
        {
            return 0.5*(rhoInside + rhoOutside);
        };

        return flux_(element, fvGeometry, insideVolVars, outsideVolVars, elemFluxVarsCache, scvf, phaseIdx, getOutsideX, getRho);
    }


    /*!
     * \brief Returns the diffusive fluxes of all components within
     *        a fluid phase across the given sub-control volume face.
     *        The computed fluxes are given in mole/s or kg/s, depending
     *        on the template parameter ReferenceSystemFormulation.
     */
    static ComponentFluxVector flux(const Problem& problem,
                                    const Element& element,
                                    const FVElementGeometry& fvGeometry,
                                    const ElementVolumeVariables& elemVolVars,
                                    const SubControlVolumeFace& scvf,
                                    const int phaseIdx,
                                    const ElementFluxVariablesCache& elemFluxVarsCache)
    {
        // get inside/outside volume variables
        const auto& insideVolVars = elemVolVars[scvf.insideScvIdx()];
        const auto& outsideVolVars = elemVolVars[scvf.outsideScvIdx()];

        // helper lambda to get the outside mole or mass fraction
        const auto getOutsideX = [&](const Scalar xInside, const Scalar tij, const int phaseIdx, const int compIdx)
        {
            const Scalar massOrMoleFractionOutside = massOrMoleFraction(outsideVolVars, referenceSystem, phaseIdx, compIdx);
            if constexpr (isMixedDimensional_)
            {
                return scvf.numOutsideScvs() == 1 ? massOrMoleFractionOutside
                                                : branchingFacetX_(problem, element, fvGeometry, elemVolVars,
                                                                   elemFluxVarsCache, scvf, xInside, tij, phaseIdx, compIdx);
            }
            else
                return massOrMoleFractionOutside;
        };

        // helper lambda to get the averaged density at the scvf
        const auto getRho = [&](const int phaseIdx, const Scalar rhoInside, const Scalar rhoOutside)
        {
            if constexpr (isMixedDimensional_)
            {
                return scvf.numOutsideScvs() == 1 ? 0.5*(rhoInside + rhoOutside)
                                                  : branchingFacetDensity_(elemVolVars, scvf, phaseIdx, rhoInside);
            }
            else
                return 0.5*(rhoInside + rhoOutside);
        };

        return flux_(element, fvGeometry, insideVolVars, outsideVolVars, elemFluxVarsCache, scvf, phaseIdx, getOutsideX, getRho);
    }

    //! compute diffusive transmissibilities
    static Scalar calculateTransmissibility(const Problem& problem,
                                            const Element& element,
                                            const FVElementGeometry& fvGeometry,
                                            const ElementVolumeVariables& elemVolVars,
                                            const SubControlVolumeFace& scvf,
                                            const int phaseIdx, const int compIdx)
    {


        const auto insideScvIdx = scvf.insideScvIdx();
        const auto& insideScv = fvGeometry.scv(insideScvIdx);
        const auto& insideVolVars = elemVolVars[insideScvIdx];
        const auto getDiffCoeff = [&](const auto& vv)
        {
            using FluidSystem = typename ElementVolumeVariables::VolumeVariables::FluidSystem;
            if constexpr (FluidSystem::isTracerFluidSystem())
                return vv.effectiveDiffusionCoefficient(0, 0, compIdx);
            else
                return vv.effectiveDiffusionCoefficient(phaseIdx, FluidSystem::getMainComponent(phaseIdx), compIdx);
        };

        const auto insideDiffCoeff = getDiffCoeff(insideVolVars);

        const Scalar ti = computeTpfaTransmissibility(fvGeometry, scvf, insideScv, insideDiffCoeff, insideVolVars.extrusionFactor());

        // for the boundary (dirichlet) or at branching points we only need ti
        Scalar tij;
        if (scvf.boundary() || scvf.numOutsideScvs() > 1)
            tij = Extrusion::area(fvGeometry, scvf)*ti;

        // otherwise we compute a tpfa harmonic mean
        else
        {
            const auto outsideScvIdx = scvf.outsideScvIdx();
            const auto& outsideScv = fvGeometry.scv(outsideScvIdx);
            const auto& outsideVolVars = elemVolVars[outsideScvIdx];
            const auto outsideDiffCoeff = getDiffCoeff(outsideVolVars);

            Scalar tj;
            if constexpr (dim == dimWorld)
                // assume the normal vector from outside is anti parallel so we save flipping a vector
                tj = -1.0*computeTpfaTransmissibility(fvGeometry, scvf, outsideScv, outsideDiffCoeff, outsideVolVars.extrusionFactor());
            else
                tj = computeTpfaTransmissibility(fvGeometry, fvGeometry.flipScvf(scvf.index()),
                                                 outsideScv,
                                                 outsideDiffCoeff,
                                                 outsideVolVars.extrusionFactor());

            // check if we are dividing by zero!
            if (ti*tj <= 0.0)
                tij = 0;
            else
                tij = Extrusion::area(fvGeometry, scvf)*(ti * tj)/(ti + tj);
        }

        return tij;
    }

private:
    template<class OutsideFractionHelper, class DensityHelper>
    static ComponentFluxVector flux_(const Element& element,
                                     const FVElementGeometry& fvGeometry,
                                     const VolumeVariables& insideVolVars,
                                     const VolumeVariables& outsideVolVars,
                                     const ElementFluxVariablesCache& elemFluxVarsCache,
                                     const SubControlVolumeFace& scvf,
                                     const int phaseIdx,
                                     const OutsideFractionHelper& getOutsideX,
                                     const DensityHelper& getRho)
    {
        ComponentFluxVector componentFlux(0.0);
        for (int compIdx = 0; compIdx < numComponents; compIdx++)
        {
            using FluidSystem = typename VolumeVariables::FluidSystem;
            if constexpr (!FluidSystem::isTracerFluidSystem())
                if (compIdx == FluidSystem::getMainComponent(phaseIdx))
                    continue;

            // diffusion tensors are always solution dependent
            const Scalar tij = elemFluxVarsCache[scvf].diffusionTij(phaseIdx, compIdx);

            // the inside and outside mass/mole fractions fractions
            const Scalar xInside = massOrMoleFraction(insideVolVars, referenceSystem, phaseIdx, compIdx);
            const Scalar xOutside = getOutsideX(xInside, tij, phaseIdx, compIdx);

            const Scalar rhoInside = massOrMolarDensity(insideVolVars, referenceSystem, phaseIdx);
            const Scalar rhoOutside = massOrMolarDensity(outsideVolVars, referenceSystem, phaseIdx);

            const Scalar rho = getRho(phaseIdx, rhoInside, rhoOutside);

            componentFlux[compIdx] = rho*tij*(xInside - xOutside);
            if constexpr (!FluidSystem::isTracerFluidSystem())
                if (BalanceEqOpts::mainComponentIsBalanced(phaseIdx))
                    componentFlux[FluidSystem::getMainComponent(phaseIdx)] -= componentFlux[compIdx];
        }

        return componentFlux;
    }

    //! compute the mole/mass fraction at branching facets for network grids
    static Scalar branchingFacetX_(const Problem& problem,
                                   const Element& element,
                                   const FVElementGeometry& fvGeometry,
                                   const ElementVolumeVariables& elemVolVars,
                                   const ElementFluxVariablesCache& elemFluxVarsCache,
                                   const SubControlVolumeFace& scvf,
                                   const Scalar insideX, const Scalar insideTi,
                                   const int phaseIdx, const int compIdx)
    {
        Scalar sumTi(insideTi);
        Scalar sumXTi(insideTi*insideX);

        for (unsigned int i = 0; i < scvf.numOutsideScvs(); ++i)
        {
            const auto outsideScvIdx = scvf.outsideScvIdx(i);
            const auto& outsideVolVars = elemVolVars[outsideScvIdx];
            const Scalar massOrMoleFractionOutside = massOrMoleFraction(outsideVolVars, referenceSystem, phaseIdx, compIdx);
            const auto& flippedScvf = fvGeometry.flipScvf(scvf.index(), i);

            const Scalar outsideTi = elemFluxVarsCache[flippedScvf].diffusionTij(phaseIdx, compIdx);
            sumTi += outsideTi;
            sumXTi += outsideTi*massOrMoleFractionOutside;
        }

        return sumTi > 0 ? sumXTi/sumTi : 0;
    }

    //! compute the density at branching facets for network grids as arithmetic mean
    static Scalar branchingFacetDensity_(const ElementVolumeVariables& elemVolVars,
                                         const SubControlVolumeFace& scvf,
                                         const int phaseIdx,
                                         const Scalar insideRho)
    {
        Scalar rho(insideRho);
        for (unsigned int i = 0; i < scvf.numOutsideScvs(); ++i)
        {
            const auto outsideScvIdx = scvf.outsideScvIdx(i);
            const auto& outsideVolVars = elemVolVars[outsideScvIdx];
            const Scalar rhoOutside = massOrMolarDensity(outsideVolVars, referenceSystem, phaseIdx);
            rho += rhoOutside;
        }
        return rho/(scvf.numOutsideScvs()+1);
    }

    static constexpr bool isMixedDimensional_ = static_cast<int>(GridView::dimension) < static_cast<int>(GridView::dimensionworld);
};

} // end namespace Dumux

#endif
