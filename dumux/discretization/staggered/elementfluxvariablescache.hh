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
 * \ingroup StaggeredDiscretization
 * \copydoc Dumux::StaggeredElementFluxVariablesCache
 */
#ifndef DUMUX_DISCRETIZATION_STAGGERED_ELEMENT_FLUXVARSCACHE_HH
#define DUMUX_DISCRETIZATION_STAGGERED_ELEMENT_FLUXVARSCACHE_HH

#include <dumux/common/properties.hh>

namespace Dumux
{

/*!
 * \ingroup StaggeredDiscretization
 * \brief Base class for the stencil local flux variables cache for the staggered model
 */
template<class TypeTag, bool EnableGridFluxVariablesCache>
class StaggeredElementFluxVariablesCache;

/*!
 * \ingroup StaggeredDiscretization
 * \brief Class for the stencil local flux variables cache for the staggered model.
          Specialization for the case of storing the fluxvars cache globally.
 */
template<class TypeTag>
class StaggeredElementFluxVariablesCache<TypeTag, true>
{
    using Problem = typename GET_PROP_TYPE(TypeTag, Problem);
    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
    using IndexType = typename GridView::IndexSet::IndexType;
    using Element = typename GridView::template Codim<0>::Entity;
    using FVElementGeometry = typename GET_PROP_TYPE(TypeTag, FVElementGeometry);
    using ElementVolumeVariables = typename GET_PROP_TYPE(TypeTag, ElementVolumeVariables);
    using FluxVariablesCache = typename GET_PROP_TYPE(TypeTag, FluxVariablesCache);
    using GridFluxVariablesCache = typename GET_PROP_TYPE(TypeTag, GridFluxVariablesCache);
    using SubControlVolumeFace = typename FVElementGeometry::SubControlVolumeFace;

public:
    StaggeredElementFluxVariablesCache(const GridFluxVariablesCache& global)
    : gridFluxVarsCachePtr_(&global) {}

    //! Specialization for the global caching being enabled - do nothing here
    void bindElement(const Element& element,
                     const FVElementGeometry& fvGeometry,
                     const ElementVolumeVariables& elemVolVars) {}

    //! Specialization for the global caching being enabled - do nothing here
    void bind(const Element& element,
              const FVElementGeometry& fvGeometry,
              const ElementVolumeVariables& elemVolVars) {}

    //! Specialization for the global caching being enabled - do nothing here
    void bindScvf(const Element& element,
                  const FVElementGeometry& fvGeometry,
                  const ElementVolumeVariables& elemVolVars,
                  const SubControlVolumeFace& scvf) {}

    //! operators in the case of caching
    const FluxVariablesCache& operator [](const SubControlVolumeFace& scvf) const
    { return (*gridFluxVarsCachePtr_)[scvf.index()]; }

    //! The global object we are a restriction of
    const GridFluxVariablesCache& gridFluxVarsCache() const
    {  return *gridFluxVarsCachePtr_; }

private:
    const GridFluxVariablesCache* gridFluxVarsCachePtr_;
};

/*!
 * \ingroup StaggeredDiscretization
 * \brief Class for the stencil local flux variables cache for the staggered model.
          Specialization for the case of not storing the fluxvars cache globally.
 */
template<class TypeTag>
class StaggeredElementFluxVariablesCache<TypeTag, false>
{
    using Problem = typename GET_PROP_TYPE(TypeTag, Problem);
    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
    using IndexType = typename GridView::IndexSet::IndexType;
    using Element = typename GridView::template Codim<0>::Entity;
    using FVElementGeometry = typename GET_PROP_TYPE(TypeTag, FVElementGeometry);
    using ElementVolumeVariables = typename GET_PROP_TYPE(TypeTag, ElementVolumeVariables);
    using FluxVariablesCache = typename GET_PROP_TYPE(TypeTag, FluxVariablesCache);
    using GridFluxVariablesCache = typename GET_PROP_TYPE(TypeTag, GridFluxVariablesCache);
    using SubControlVolumeFace = typename FVElementGeometry::SubControlVolumeFace;

public:
    StaggeredElementFluxVariablesCache(const GridFluxVariablesCache& global)
    : gridFluxVarsCachePtr_(&global) {}

    //! This function has to be called prior to flux calculations on the element.
    //! Prepares the transmissibilities of the scv faces in an element. The FvGeometry is assumed to be bound.
    void bindElement(const Element& element,
                     const FVElementGeometry& fvGeometry,
                     const ElementVolumeVariables& elemVolVars)
    {
        // resizing of the cache
        const auto numScvf = fvGeometry.numScvf();
        fluxVarsCache_.resize(numScvf);
        globalScvfIndices_.resize(numScvf);
        IndexType localScvfIdx = 0;

        // fill the containers
        for (auto&& scvf : scvfs(fvGeometry))
        {
            fluxVarsCache_[localScvfIdx].update(gridFluxVarsCache().problem_(), element, fvGeometry, elemVolVars, scvf);
            globalScvfIndices_[localScvfIdx] = scvf.index();
            localScvfIdx++;
        }
    }

    //! This function is called by the StaggeredLocalResidual before flux calculations during assembly.
    //! Prepares the transmissibilities of the scv faces in the stencil. The FvGeometries are assumed to be bound.
    void bind(const Element& element,
              const FVElementGeometry& fvGeometry,
              const ElementVolumeVariables& elemVolVars)
    {
        const auto globalI = gridFluxVarsCache().problem_().elementMapper().index(element);
        const auto& neighborStencil = gridFluxVarsCache().problem_().model().stencils(element).neighborStencil();
        const auto numNeighbors = neighborStencil.size();

        // find the number of scv faces that need to be prepared
        auto numScvf = fvGeometry.numScvf();
        for (IndexType localIdxJ = 0; localIdxJ < numNeighbors; ++localIdxJ)
        {
            const auto& fluxVarIndicesJ = gridFluxVarsCache().problem_().model().localJacobian().assemblyMap()[globalI][localIdxJ];
            numScvf += fluxVarIndicesJ.size();
        }

        // fill the containers with the data on the scv faces inside the actual element
        fluxVarsCache_.resize(numScvf);
        globalScvfIndices_.resize(numScvf);
        IndexType localScvfIdx = 0;
        for (auto&& scvf : scvfs(fvGeometry))
        {
            fluxVarsCache_[localScvfIdx].update(gridFluxVarsCache().problem_(), element, fvGeometry, elemVolVars, scvf);
            globalScvfIndices_[localScvfIdx] = scvf.index();
            localScvfIdx++;
        }

        // add required data on the scv faces in the neighboring elements
        for (IndexType localIdxJ = 0; localIdxJ < numNeighbors; ++localIdxJ)
        {
            const auto& fluxVarIndicesJ = gridFluxVarsCache().problem_().model().localJacobian().assemblyMap()[globalI][localIdxJ];
            const auto elementJ = fvGeometry.fvGridGeometry().element(neighborStencil[localIdxJ]);

            for (auto fluxVarIdx : fluxVarIndicesJ)
            {
                auto&& scvfJ = fvGeometry.scvf(fluxVarIdx);
                fluxVarsCache_[localScvfIdx].update(gridFluxVarsCache().problem_(), elementJ, fvGeometry, elemVolVars, scvfJ);
                globalScvfIndices_[localScvfIdx] = scvfJ.index();
                localScvfIdx++;
            }
        }
    }

    void bindScvf(const Element& element,
                  const FVElementGeometry& fvGeometry,
                  const ElementVolumeVariables& elemVolVars,
                  const SubControlVolumeFace& scvf)
    {
        fluxVarsCache_.resize(1);
        globalScvfIndices_.resize(1);

        fluxVarsCache_[0].update(gridFluxVarsCache().problem_(), element, fvGeometry, elemVolVars, scvf);
        globalScvfIndices_[0] = scvf.index();
    }

    //! access operators in the case of no caching
    const FluxVariablesCache& operator [](const SubControlVolumeFace& scvf) const
    { return fluxVarsCache_[getLocalScvfIdx_(scvf.index())]; }

    FluxVariablesCache& operator [](const SubControlVolumeFace& scvf)
    { return fluxVarsCache_[getLocalScvfIdx_(scvf.index())]; }

    //! The global object we are a restriction of
    const GridFluxVariablesCache& gridFluxVarsCache() const
    {  return *gridFluxVarsCachePtr_; }

private:
    const GridFluxVariablesCache* gridFluxVarsCachePtr_;

    // get index of scvf in the local container
    int getLocalScvfIdx_(const int scvfIdx) const
    {
        auto it = std::find(globalScvfIndices_.begin(), globalScvfIndices_.end(), scvfIdx);
        assert(it != globalScvfIndices_.end() && "Could not find the flux vars cache for scvfIdx");
        return std::distance(globalScvfIndices_.begin(), it);
    }

    std::vector<FluxVariablesCache> fluxVarsCache_;
    std::vector<IndexType> globalScvfIndices_;
};

} // end namespace

#endif
