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
#ifndef DUMUX_DISCRETIZATION_CCMPFA_ELEMENT_FLUXVARSCACHE_HH
#define DUMUX_DISCRETIZATION_CCMPFA_ELEMENT_FLUXVARSCACHE_HH

#include <dumux/implicit/properties.hh>

#include "fluxvariablescachefiller.hh"

namespace Dumux
{

/*!
 * \ingroup ImplicitModel
 * \brief Base class for the stencil local flux variables cache
 */
template<class TypeTag, bool EnableGlobalFluxVariablesCache>
class CCMpfaElementFluxVariablesCache;

/*!
 * \ingroup ImplicitModel
 * \brief Spezialization when caching globally
 */
template<class TypeTag>
class CCMpfaElementFluxVariablesCache<TypeTag, true>
{
    using Problem = typename GET_PROP_TYPE(TypeTag, Problem);
    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
    using IndexType = typename GridView::IndexSet::IndexType;
    using Element = typename GridView::template Codim<0>::Entity;
    using FVElementGeometry = typename GET_PROP_TYPE(TypeTag, FVElementGeometry);
    using ElementVolumeVariables = typename GET_PROP_TYPE(TypeTag, ElementVolumeVariables);
    using FluxVariablesCache = typename GET_PROP_TYPE(TypeTag, FluxVariablesCache);
    using GlobalFluxVariablesCache = typename GET_PROP_TYPE(TypeTag, GlobalFluxVariablesCache);
    using SubControlVolumeFace = typename GET_PROP_TYPE(TypeTag, SubControlVolumeFace);

public:
    CCMpfaElementFluxVariablesCache(const GlobalFluxVariablesCache& global)
    : globalFluxVarsCachePtr_(&global) {}

    // Specialization for the global caching being enabled - do nothing here
    void bindElement(const Element& element,
                     const FVElementGeometry& fvGeometry,
                     const ElementVolumeVariables& elemVolVars) {}

    // Specialization for the global caching being enabled - do nothing here
    void bind(const Element& element,
              const FVElementGeometry& fvGeometry,
              const ElementVolumeVariables& elemVolVars) {}

    // Specialization for the global caching being enabled - do nothing here
    void bindScvf(const Element& element,
                  const FVElementGeometry& fvGeometry,
                  const ElementVolumeVariables& elemVolVars,
                  const SubControlVolumeFace& scvf) {}

    // access operators in the case of caching
    const FluxVariablesCache& operator [](const SubControlVolumeFace& scvf) const
    { return (*globalFluxVarsCachePtr_)[scvf.index()]; }

    //! The global object we are a restriction of
    const GlobalFluxVariablesCache& globalFluxVarsCache() const
    {  return *globalFluxVarsCachePtr_; }

private:
    const GlobalFluxVariablesCache* globalFluxVarsCachePtr_;
};

/*!
 * \ingroup ImplicitModel
 * \brief Spezialization when not using global caching
 */
template<class TypeTag>
class CCMpfaElementFluxVariablesCache<TypeTag, false>
{
    using Problem = typename GET_PROP_TYPE(TypeTag, Problem);
    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
    using IndexType = typename GridView::IndexSet::IndexType;
    using Element = typename GridView::template Codim<0>::Entity;
    using FVElementGeometry = typename GET_PROP_TYPE(TypeTag, FVElementGeometry);
    using ElementVolumeVariables = typename GET_PROP_TYPE(TypeTag, ElementVolumeVariables);
    using FluxVariablesCache = typename GET_PROP_TYPE(TypeTag, FluxVariablesCache);
    using GlobalFluxVariablesCache = typename GET_PROP_TYPE(TypeTag, GlobalFluxVariablesCache);
    using SubControlVolumeFace = typename GET_PROP_TYPE(TypeTag, SubControlVolumeFace);
    using FluxVariablesCacheFiller = CCMpfaFluxVariablesCacheFiller<TypeTag>;

public:
    CCMpfaElementFluxVariablesCache(const GlobalFluxVariablesCache& global)
    : globalFluxVarsCachePtr_(&global) {}

    // This function has to be called prior to flux calculations on the element.
    // Prepares the transmissibilities of the scv faces in an element. The FvGeometry is assumed to be bound.
    void bindElement(const Element& element,
                     const FVElementGeometry& fvGeometry,
                     const ElementVolumeVariables& elemVolVars)
    {
        DUNE_THROW(Dune::InvalidStateException, "Does local flux var cache binding make sense in general?");
    }

    // This function is called by the CCLocalResidual before flux calculations during assembly.
    // Prepares the transmissibilities of the scv faces in the stencil. The FvGeometries are assumed to be bound.
    void bind(const Element& element,
              const FVElementGeometry& fvGeometry,
              const ElementVolumeVariables& elemVolVars)
    {
        clear();
        const auto& problem = globalFluxVarsCache().problem_();
        const auto& globalFvGeometry = problem.model().globalFvGeometry();

        // // first prepare all the global indices
        for (auto&& scvf : scvfs(fvGeometry))
        {
            const bool boundary = globalFvGeometry.scvfTouchesBoundary(scvf);
            const auto& ivSeed = boundary ? globalFvGeometry.boundaryInteractionVolumeSeed(scvf) : globalFvGeometry.interactionVolumeSeed(scvf);

            // loop over all the scvfs in the interaction region
            for (auto scvfIdx : ivSeed.globalScvfIndices())
                globalScvfIndices_.push_back(scvfIdx);
        }
        // make global indices unique
        std::sort(globalScvfIndices_.begin(), globalScvfIndices_.end());
        globalScvfIndices_.erase(std::unique(globalScvfIndices_.begin(), globalScvfIndices_.end()), globalScvfIndices_.end());

        // prepare all the caches of the scvfs inside the corresponding interaction volume using helper class
        fluxVarsCache_.resize(globalScvfIndices_.size());
        for (auto&& scvf : scvfs(fvGeometry))
        {
            if (!(*this)[scvf].isUpdated())
                FluxVariablesCacheFiller::fillFluxVarCache(problem, element, fvGeometry, elemVolVars, scvf, *this);
        }
    }

    void bindScvf(const Element& element,
                  const FVElementGeometry& fvGeometry,
                  const ElementVolumeVariables& elemVolVars,
                  const SubControlVolumeFace& scvf)
    {
        DUNE_THROW(Dune::InvalidStateException, "Does binding of one scvf make sense in general?");
    }

    // access operators in the case of no caching
    const FluxVariablesCache& operator [](const SubControlVolumeFace& scvf) const
    { return fluxVarsCache_[getLocalScvfIdx_(scvf.index())]; }

    const FluxVariablesCache& operator [](const IndexType scvfIdx) const
    { return fluxVarsCache_[getLocalScvfIdx_(scvfIdx)]; }

    FluxVariablesCache& operator [](const SubControlVolumeFace& scvf)
    { return fluxVarsCache_[getLocalScvfIdx_(scvf.index())]; }

    FluxVariablesCache& operator [](const IndexType scvfIdx)
    { return fluxVarsCache_[getLocalScvfIdx_(scvfIdx)]; }

    //! The global object we are a restriction of
    const GlobalFluxVariablesCache& globalFluxVarsCache() const
    {  return *globalFluxVarsCachePtr_; }

private:
    const GlobalFluxVariablesCache* globalFluxVarsCachePtr_;

    void clear()
    {
        fluxVarsCache_.clear();
        globalScvfIndices_.clear();
    }

    // get index of scvf in the local container
    int getLocalScvfIdx_(const int scvfIdx) const
    {
        auto it = std::lower_bound(globalScvfIndices_.begin(), globalScvfIndices_.end(), scvfIdx);
        assert(it != globalScvfIndices_.end() && "Could not find the flux vars cache for scvfIdx");
        assert(globalScvfIndices_[std::distance(globalScvfIndices_.begin(), it)] == scvfIdx && "Could not find the flux vars cache for scvfIdx");
        return std::distance(globalScvfIndices_.begin(), it);
    }

    std::vector<FluxVariablesCache> fluxVarsCache_;
    std::vector<IndexType> globalScvfIndices_;
};

} // end namespace

#endif
