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
#include <dumux/discretization/cellcentered/mpfa/fluxvariablescachefiller.hh>

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
    // the local class need access to the problem
    friend typename GET_PROP_TYPE(TypeTag, ElementFluxVariablesCache);
    // the filler needs access to the operators
    friend CCMpfaFluxVariablesCacheFiller<TypeTag>;

    using Problem = typename GET_PROP_TYPE(TypeTag, Problem);
    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
    using IndexType = typename GridView::IndexSet::IndexType;
    using Element = typename GridView::template Codim<0>::Entity;
    using FluxVariablesCache = typename GET_PROP_TYPE(TypeTag, FluxVariablesCache);
    using ElementFluxVariablesCache = typename GET_PROP_TYPE(TypeTag, ElementFluxVariablesCache);
    using SubControlVolume = typename GET_PROP_TYPE(TypeTag, SubControlVolume);
    using SubControlVolumeFace = typename GET_PROP_TYPE(TypeTag, SubControlVolumeFace);
    using InteractionVolume = typename GET_PROP_TYPE(TypeTag, InteractionVolume);
    using BoundaryInteractionVolume = typename GET_PROP_TYPE(TypeTag, BoundaryInteractionVolume);
    using DataHandle = typename BoundaryInteractionVolume::Traits::DataHandle;
    using FVElementGeometry = typename GET_PROP_TYPE(TypeTag, FVElementGeometry);
    using VolumeVariables = typename GET_PROP_TYPE(TypeTag, VolumeVariables);
    using ElementVolumeVariables = typename GET_PROP_TYPE(TypeTag, ElementVolumeVariables);
    using FluxVariablesCacheFiller = CCMpfaFluxVariablesCacheFiller<TypeTag>;

    static constexpr int dim = GridView::dimension;
    static constexpr int dimWorld = GridView::dimensionworld;

public:
    // When global caching is enabled, precompute transmissibilities and stencils for all the scv faces
    void update(Problem& problem)
    {
        problemPtr_ = &problem;

        // clear data
        clear_();

        // reserve memory estimate for caches, interaction volumes and corresponding data
        const auto& globalFvGeometry = problem.model().globalFvGeometry();
        fluxVarsCache_.resize(globalFvGeometry.numScvf());

        const auto numInnerIVs = globalFvGeometry.numInteractionVolumeSeeds();
        const auto numBoundaryIVs = globalFvGeometry.numBoundaryInteractionVolumeSeeds();
        interactionVolumes_.reserve(numInnerIVs);
        boundaryInteractionVolumes_.reserve(numBoundaryIVs);
        ivDataHandles_.reserve(numInnerIVs);
        boundaryIvDataHandles_.reserve(numBoundaryIVs);

        // instantiate helper class to fill the caches
        FluxVariablesCacheFiller filler(problem);
        for (const auto& element : elements(problem.gridView()))
        {
            //! TODO The contexts have to be more general
            //! and the coupling manager should never have a local state
            prepareContext_(problem, element);

            // Prepare the geometries within the elements of the stencil
            auto fvGeometry = localView(globalFvGeometry);
            fvGeometry.bind(element);

            auto elemVolVars = localView(problem.model().curGlobalVolVars());
            elemVolVars.bind(element, fvGeometry, problem.model().curSol());

            // prepare all the caches of the scvfs inside the corresponding interaction volume
            for (auto&& scvf : scvfs(fvGeometry))
                if (!fluxVarsCache_[scvf.index()].isUpdated())
                    filler.fill(*this, fluxVarsCache_[scvf.index()], element, fvGeometry, elemVolVars, scvf);
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

    void clear_()
    {
        fluxVarsCache_.clear();
        interactionVolumes_.clear();
        boundaryInteractionVolumes_.clear();
        ivDataHandles_.clear();
        boundaryIvDataHandles_.clear();
    }

    // TODO Hopefully we can get rid of this prop tag in the non-coupled framework
    template<class T = TypeTag>
    typename std::enable_if<GET_PROP_VALUE(T, MpfaFacetCoupling)>::type
    prepareContext_(Problem& problem, const Element& element)
    { problem.couplingManager().setCouplingContext(element); }

    template<class T = TypeTag>
    typename std::enable_if<!GET_PROP_VALUE(T, MpfaFacetCoupling)>::type
    prepareContext_(Problem& problem, const Element& element) {}

    // access operators in the case of caching
    const FluxVariablesCache& operator [](const SubControlVolumeFace& scvf) const
    { return fluxVarsCache_[scvf.index()]; }

    FluxVariablesCache& operator [](const SubControlVolumeFace& scvf)
    { return fluxVarsCache_[scvf.index()]; }

    const FluxVariablesCache& operator [](IndexType scvfIdx) const
    { return fluxVarsCache_[scvfIdx]; }

    FluxVariablesCache& operator [](IndexType scvfIdx)
    { return fluxVarsCache_[scvfIdx]; }

    const Problem& problem_() const
    { return *problemPtr_; }

    const Problem* problemPtr_;

    std::vector<FluxVariablesCache> fluxVarsCache_;

    // store the interaction volumes and handles
    std::vector<InteractionVolume> interactionVolumes_;
    std::vector<BoundaryInteractionVolume> boundaryInteractionVolumes_;
    std::vector<DataHandle> ivDataHandles_;
    std::vector<DataHandle> boundaryIvDataHandles_;
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
