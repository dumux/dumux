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
 * \brief The global object of flux var caches
 */
#ifndef DUMUX_DISCRETIZATION_CCMPFA_GLOBAL_FLUXVARSCACHE_HH
#define DUMUX_DISCRETIZATION_CCMPFA_GLOBAL_FLUXVARSCACHE_HH

#include <dumux/implicit/properties.hh>

namespace Dumux
{

/*!
 * \ingroup ImplicitModel
 * \brief Base class for the flux variables cache vector, we store one cache per face
 */
template<class TypeTag, bool EnableGlobalFluxVariablesCache>
class CCMpfaGlobalFluxVariablesCache;

/*!
 * \ingroup ImplicitModel
 * \brief Spezialization when caching globally
 */
template<class TypeTag>
class CCMpfaGlobalFluxVariablesCache<TypeTag, true>
{
    // the local class needs access to the problem
    friend typename GET_PROP_TYPE(TypeTag, ElementFluxVariablesCache);
    using Problem = typename GET_PROP_TYPE(TypeTag, Problem);
    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
    using IndexType = typename GridView::IndexSet::IndexType;
    using FluxVariablesCache = typename GET_PROP_TYPE(TypeTag, FluxVariablesCache);
    using ElementFluxVariablesCache = typename GET_PROP_TYPE(TypeTag, ElementFluxVariablesCache);
    using SubControlVolume = typename GET_PROP_TYPE(TypeTag, SubControlVolume);
    using SubControlVolumeFace = typename GET_PROP_TYPE(TypeTag, SubControlVolumeFace);
    using InteractionVolume = typename GET_PROP_TYPE(TypeTag, InteractionVolume);
    using BoundaryInteractionVolume = typename GET_PROP_TYPE(TypeTag, BoundaryInteractionVolume);
    using FVElementGeometry = typename GET_PROP_TYPE(TypeTag, FVElementGeometry);
    using VolumeVariables = typename GET_PROP_TYPE(TypeTag, VolumeVariables);
    using ElementVolumeVariables = typename GET_PROP_TYPE(TypeTag, ElementVolumeVariables);

    using Element = typename GridView::template Codim<0>::Entity;

    static const bool advection = GET_PROP_VALUE(TypeTag, EnableAdvection);
    static const bool diffusion = GET_PROP_VALUE(TypeTag, EnableMolecularDiffusion);
    static const bool energy = GET_PROP_VALUE(TypeTag, EnableEnergyBalance);
    static const int numEq = GET_PROP_VALUE(TypeTag, NumEq);

public:
    // When global caching is enabled, precompute transmissibilities and stencils for all the scv faces
    void update(Problem& problem)
    {
        problemPtr_ = &problem;
        const auto& globalFvGeometry = problem.model().globalFvGeometry();
        fluxVarsCache_.resize(globalFvGeometry.numScvf());
        for (const auto& element : elements(problem.gridView()))
        {
            // Prepare the geometries within the elements of the stencil
            auto fvGeometry = localView(globalFvGeometry);
            fvGeometry.bind(element);

            auto elemVolVars = localView(problem.model().curGlobalVolVars());
            elemVolVars.bind(element, fvGeometry, problem.model().curSol());

            // prepare all the caches of the scvfs inside the corresponding interaction volume
            for (auto&& scvf : scvfs(fvGeometry))
            {
                if (!fluxVarsCache_[scvf.index()].isUpdated())
                    fillFluxVarCache(element, fvGeometry, elemVolVars, scvf);
            }
        }
    }

    /*!
     * \brief Return a local restriction of this global object
     *        The local object is only functional after calling its bind/bindElement method
     *        This is a free function that will be found by means of ADL
     */
    friend inline ElementFluxVariablesCache localView(const CCMpfaGlobalFluxVariablesCache& global)
    { return ElementFluxVariablesCache(global); }

private:
    // functions to fill the flux var caches in the case of pure advection
    template<class T = TypeTag>
    typename std::enable_if<advection && !diffusion && !energy>::type
    fillFluxVarCache(const Element& element,
                     const FVElementGeometry& fvGeometry,
                     const ElementVolumeVariables& elemVolVars,
                     const SubControlVolumeFace& scvf)
    {
        // lambda function to get the permeability tensor
        const auto* prob = &problem_();
        auto permFunction = [prob](const Element& element, const VolumeVariables& volVars, const SubControlVolume& scv)
                            { return prob->spatialParams().intrinsicPermeability(scv, volVars); };

        // update the flux var caches for this scvf
        if (problem_().model().globalFvGeometry().scvfTouchesBoundary(scvf))
        {
            // we assume phaseIdx = eqIdx here for purely advective problems
            for (unsigned int eqIdx = 0; eqIdx < numEq; ++eqIdx)
            {
                auto boundarySeed = problem_().model().globalFvGeometry().boundaryInteractionVolumeSeed(scvf, eqIdx);
                BoundaryInteractionVolume iv(boundarySeed, problem_(), fvGeometry, elemVolVars);
                iv.solveLocalSystem(permFunction);

                // lambda function defining the upwind factor of the advective flux
                auto advectionUpwindFunction = [eqIdx](const VolumeVariables& volVars) { return volVars.density(eqIdx)*volVars.mobility(eqIdx); };
                iv.assembleNeumannFluxes(advectionUpwindFunction, eqIdx);

                // update flux variables cache
                fluxVarsCache_[scvf.index()].updateBoundaryAdvection(problem_(), element, fvGeometry, elemVolVars, scvf, iv, eqIdx);

                // update flux variable caches of the other scvfs of the interaction volume
                for (const auto& scvfIdx : iv.globalScvfs())
                {
                    if (scvfIdx != scvf.index())
                    {
                        const auto& scvfJ = fvGeometry.scvf(scvfIdx);
                        const auto elementJ = problem_().model().globalFvGeometry().element(scvfJ.insideScvIdx());
                        fluxVarsCache_[scvfIdx].updateBoundaryAdvection(problem_(), elementJ, fvGeometry, elemVolVars, scvfJ, iv, eqIdx);
                        if (eqIdx == numEq - 1)
                            fluxVarsCache_[scvfIdx].setUpdated();
                    }
                }
            }
        }
        else
        {
            auto seed = problem_().model().globalFvGeometry().interactionVolumeSeed(scvf);
            InteractionVolume iv(seed, problem_(), fvGeometry, elemVolVars);
            iv.solveLocalSystem(permFunction);

            // update flux variables cache
            fluxVarsCache_[scvf.index()].updateInnerAdvection(problem_(), element, fvGeometry, elemVolVars, scvf, iv);

            // update flux variable caches of the other scvfs of the interaction volume
            for (const auto& scvfIdx : iv.globalScvfs())
            {
                if (scvfIdx != scvf.index())
                {
                    const auto& scvfJ = fvGeometry.scvf(scvfIdx);
                    const auto elementJ = problem_().model().globalFvGeometry().element(scvfJ.insideScvIdx());
                    fluxVarsCache_[scvfIdx].updateInnerAdvection(problem_(), elementJ, fvGeometry, elemVolVars, scvfJ, iv);
                    fluxVarsCache_[scvfIdx].setUpdated();
                }
            }
        }

        // the flux var cache has been updated
        fluxVarsCache_[scvf.index()].setUpdated();
    }


    // access operators in the case of caching
    const FluxVariablesCache& operator [](IndexType scvfIdx) const
    { return fluxVarsCache_[scvfIdx]; }

    FluxVariablesCache& operator [](IndexType scvfIdx)
    { return fluxVarsCache_[scvfIdx]; }

    const Problem& problem_() const
    { return *problemPtr_; }

    const Problem* problemPtr_;

    std::vector<FluxVariablesCache> fluxVarsCache_;
    std::vector<IndexType> globalScvfIndices_;
};

/*!
 * \ingroup ImplicitModel
 * \brief Spezialization when not using global caching
 */
template<class TypeTag>
class CCMpfaGlobalFluxVariablesCache<TypeTag, false>
{
    // the local class needs access to the problem
    friend typename GET_PROP_TYPE(TypeTag, ElementFluxVariablesCache);
    using Problem = typename GET_PROP_TYPE(TypeTag, Problem);
    using ElementFluxVariablesCache = typename GET_PROP_TYPE(TypeTag, ElementFluxVariablesCache);

public:
    // When global flux variables caching is disabled, we don't need to update the cache
    void update(Problem& problem)
    { problemPtr_ = &problem; }

    /*!
     * \brief Return a local restriction of this global object
     *        The local object is only functional after calling its bind/bindElement method
     *        This is a free function that will be found by means of ADL
     */
    friend inline ElementFluxVariablesCache localView(const CCMpfaGlobalFluxVariablesCache& global)
    { return ElementFluxVariablesCache(global); }

private:

    const Problem& problem_() const
    { return *problemPtr_; }

    const Problem* problemPtr_;
};

} // end namespace

#endif
