// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup CCTpfaDiscretization
 * \brief The flux variables caches for an element
 */
#ifndef DUMUX_DISCRETIZATION_CCTPFA_ELEMENT_FLUXVARSCACHE_HH
#define DUMUX_DISCRETIZATION_CCTPFA_ELEMENT_FLUXVARSCACHE_HH

#include <algorithm>
#include <cassert>
#include <vector>
#include <utility>

#include <dune/common/exceptions.hh>

namespace Dumux {

/*!
 * \ingroup CCTpfaDiscretization
 * \brief The flux variables caches for an element
 * \note The class is specialized for a version with and without caching
 * If grid caching is enabled the flux caches are stored for the whole gridview in the corresponding
 * GridFluxVariablesCache which is memory intensive but faster. For caching disabled the
 * flux caches are locally computed for each element whenever needed.
 */
template<class GFVC, bool cachingEnabled>
class CCTpfaElementFluxVariablesCache;

/*!
 * \ingroup CCTpfaDiscretization
 * \brief The flux variables caches for an element with caching enabled
 */
template<class GFVC>
class CCTpfaElementFluxVariablesCache<GFVC, true>
{
    //! the type of the flux variables cache filler
    using FluxVariablesCacheFiller = typename GFVC::Traits::FluxVariablesCacheFiller;

public:
    //! export the type of the grid flux variables cache
    using GridFluxVariablesCache = GFVC;

    //! export the type of the flux variables cache
    using FluxVariablesCache = typename GFVC::FluxVariablesCache;

    CCTpfaElementFluxVariablesCache(const GridFluxVariablesCache& global)
    : gridFluxVarsCachePtr_(&global) {}

    /*!
    * \brief bind the local view (r-value overload)
    * This overload is called when an instance of this class is a temporary in the usage context
    * This allows a usage like this: `const auto view = localView(...).bind(element);`
    */
    template<class FVElementGeometry, class ElementVolumeVariables>
    CCTpfaElementFluxVariablesCache bindElement(const typename FVElementGeometry::GridGeometry::GridView::template Codim<0>::Entity& element,
                                                const FVElementGeometry& fvGeometry,
                                                const ElementVolumeVariables& elemVolVars) &&
    {
        this->bindElement(element, fvGeometry, elemVolVars);
        return std::move(*this);
    }

    //! Specialization for the global caching being enabled - do nothing here
    template<class FVElementGeometry, class ElementVolumeVariables>
    void bindElement(const typename FVElementGeometry::GridGeometry::GridView::template Codim<0>::Entity& element,
                     const FVElementGeometry& fvGeometry,
                     const ElementVolumeVariables& elemVolVars) & {}

    /*!
    * \brief bind the local view (r-value overload)
    * This overload is called when an instance of this class is a temporary in the usage context
    * This allows a usage like this: `const auto view = localView(...).bind(element);`
    */
    template<class FVElementGeometry, class ElementVolumeVariables>
    CCTpfaElementFluxVariablesCache bind(const typename FVElementGeometry::GridGeometry::GridView::template Codim<0>::Entity& element,
                                         const FVElementGeometry& fvGeometry,
                                         const ElementVolumeVariables& elemVolVars) &&
    {
        this->bind(element, fvGeometry, elemVolVars);
        return std::move(*this);
    }

    //! Specialization for the global caching being enabled - do nothing here
    template<class FVElementGeometry, class ElementVolumeVariables>
    void bind(const typename FVElementGeometry::GridGeometry::GridView::template Codim<0>::Entity& element,
              const FVElementGeometry& fvGeometry,
              const ElementVolumeVariables& elemVolVars) & {}

    /*!
    * \brief bind the local view (r-value overload)
    * This overload is called when an instance of this class is a temporary in the usage context
    * This allows a usage like this: `const auto view = localView(...).bind(element);`
    */
    template<class FVElementGeometry, class ElementVolumeVariables>
    CCTpfaElementFluxVariablesCache bindScvf(const typename FVElementGeometry::GridGeometry::GridView::template Codim<0>::Entity& element,
                                             const FVElementGeometry& fvGeometry,
                                             const ElementVolumeVariables& elemVolVars,
                                             const typename FVElementGeometry::SubControlVolumeFace& scvf) &&
    {
        this->bindScvf(element, fvGeometry, elemVolVars, scvf);
        return std::move(*this);
    }

    //! Specialization for the global caching being enabled - do nothing here
    template<class FVElementGeometry, class ElementVolumeVariables>
    void bindScvf(const typename FVElementGeometry::GridGeometry::GridView::template Codim<0>::Entity& element,
                  const FVElementGeometry& fvGeometry,
                  const ElementVolumeVariables& elemVolVars,
                  const typename FVElementGeometry::SubControlVolumeFace& scvf) & {}

    //! Specialization for the global caching being enabled - do nothing here
    template<class FVElementGeometry, class ElementVolumeVariables>
    void update(const typename FVElementGeometry::GridGeometry::GridView::template Codim<0>::Entity& element,
                const FVElementGeometry& fvGeometry,
                const ElementVolumeVariables& elemVolVars) {}

    //! access operators in the case of caching
    template<class SubControlVolumeFace>
    const FluxVariablesCache& operator [](const SubControlVolumeFace& scvf) const
    { return gridFluxVarsCache()[scvf]; }

    //! The global object we are a restriction of
    const GridFluxVariablesCache& gridFluxVarsCache() const
    {  return *gridFluxVarsCachePtr_; }

private:
    const GridFluxVariablesCache* gridFluxVarsCachePtr_;
};

/*!
 * \ingroup CCTpfaDiscretization
 * \brief The flux variables caches for an element with caching disabled
 */
template<class GFVC>
class CCTpfaElementFluxVariablesCache<GFVC, false>
{
    //! the type of the flux variables cache filler
    using FluxVariablesCacheFiller = typename GFVC::Traits::FluxVariablesCacheFiller;

public:
    //! export the type of the grid flux variables cache
    using GridFluxVariablesCache = GFVC;

    //! export the type of the flux variables cache
    using FluxVariablesCache = typename GFVC::FluxVariablesCache;

    CCTpfaElementFluxVariablesCache(const GridFluxVariablesCache& global)
    : gridFluxVarsCachePtr_(&global) {}

    /*!
    * \brief bind the local view (r-value overload)
    * This overload is called when an instance of this class is a temporary in the usage context
    * This allows a usage like this: `const auto view = localView(...).bind(element);`
    */
    template<class FVElementGeometry, class ElementVolumeVariables>
    CCTpfaElementFluxVariablesCache bindElement(const typename FVElementGeometry::GridGeometry::GridView::template Codim<0>::Entity& element,
                                                const FVElementGeometry& fvGeometry,
                                                const ElementVolumeVariables& elemVolVars) &&
    {
        this->bindElement_(element, fvGeometry, elemVolVars);
        return std::move(*this);
    }

    //! Specialization for the global caching being enabled - do nothing here
    template<class FVElementGeometry, class ElementVolumeVariables>
    void bindElement(const typename FVElementGeometry::GridGeometry::GridView::template Codim<0>::Entity& element,
                     const FVElementGeometry& fvGeometry,
                     const ElementVolumeVariables& elemVolVars) &
    { this->bindElement_(element, fvGeometry, elemVolVars); }

    /*!
    * \brief bind the local view (r-value overload)
    * This overload is called when an instance of this class is a temporary in the usage context
    * This allows a usage like this: `const auto view = localView(...).bind(element);`
    */
    template<class FVElementGeometry, class ElementVolumeVariables>
    CCTpfaElementFluxVariablesCache bind(const typename FVElementGeometry::GridGeometry::GridView::template Codim<0>::Entity& element,
                                         const FVElementGeometry& fvGeometry,
                                         const ElementVolumeVariables& elemVolVars) &&
    {
        this->bind_(element, fvGeometry, elemVolVars);
        return std::move(*this);
    }

    //! Specialization for the global caching being enabled - do nothing here
    template<class FVElementGeometry, class ElementVolumeVariables>
    void bind(const typename FVElementGeometry::GridGeometry::GridView::template Codim<0>::Entity& element,
              const FVElementGeometry& fvGeometry,
              const ElementVolumeVariables& elemVolVars) &
    { this->bind_(element, fvGeometry, elemVolVars); }

    /*!
    * \brief bind the local view (r-value overload)
    * This overload is called when an instance of this class is a temporary in the usage context
    * This allows a usage like this: `const auto view = localView(...).bind(element);`
    */
    template<class FVElementGeometry, class ElementVolumeVariables>
    CCTpfaElementFluxVariablesCache bindScvf(const typename FVElementGeometry::GridGeometry::GridView::template Codim<0>::Entity& element,
                                             const FVElementGeometry& fvGeometry,
                                             const ElementVolumeVariables& elemVolVars,
                                             const typename FVElementGeometry::SubControlVolumeFace& scvf) &&
    {
        this->bindScvf_(element, fvGeometry, elemVolVars, scvf);
        return std::move(*this);
    }

    //! Specialization for the global caching being enabled - do nothing here
    template<class FVElementGeometry, class ElementVolumeVariables>
    void bindScvf(const typename FVElementGeometry::GridGeometry::GridView::template Codim<0>::Entity& element,
                  const FVElementGeometry& fvGeometry,
                  const ElementVolumeVariables& elemVolVars,
                  const typename FVElementGeometry::SubControlVolumeFace& scvf) &
    { this->bindScvf_(element, fvGeometry, elemVolVars, scvf); }

    /*!
     * \brief Update the transmissibilities if the volume variables have changed
     * \note Results in undefined behaviour if called before bind() or with a different element
     */
    template<class FVElementGeometry, class ElementVolumeVariables>
    void update(const typename FVElementGeometry::GridGeometry::GridView::template Codim<0>::Entity& element,
                const FVElementGeometry& fvGeometry,
                const ElementVolumeVariables& elemVolVars)
    {
        if (FluxVariablesCacheFiller::isSolDependent)
        {
            const auto& problem = gridFluxVarsCache().problem();
            const auto globalI = fvGeometry.gridGeometry().elementMapper().index(element);

            // instantiate filler class
            FluxVariablesCacheFiller filler(problem);

            // let the filler class update the caches
            for (unsigned int localScvfIdx = 0; localScvfIdx < fluxVarsCache_.size(); ++localScvfIdx)
            {
                const auto& scvf = fvGeometry.scvf(globalScvfIndices_[localScvfIdx]);

                const auto scvfInsideScvIdx = scvf.insideScvIdx();
                const auto& insideElement = scvfInsideScvIdx == globalI ?
                                            element :
                                            fvGeometry.gridGeometry().element(scvfInsideScvIdx);

                filler.fill(*this, fluxVarsCache_[localScvfIdx], insideElement, fvGeometry, elemVolVars, scvf);
            }
        }
    }

    //! access operators in the case of no caching
    template<class SubControlVolumeFace>
    const FluxVariablesCache& operator [](const SubControlVolumeFace& scvf) const
    { return fluxVarsCache_[getLocalScvfIdx_(scvf.index())]; }

    //! access operators in the case of no caching
    template<class SubControlVolumeFace>
    FluxVariablesCache& operator [](const SubControlVolumeFace& scvf)
    { return fluxVarsCache_[getLocalScvfIdx_(scvf.index())]; }

    //! The global object we are a restriction of
    const GridFluxVariablesCache& gridFluxVarsCache() const
    {  return *gridFluxVarsCachePtr_; }

private:

    /*!
     * \brief Prepares the transmissibilities of the scv faces in an element
     * \note the fvGeometry is assumed to be bound to the same element
     * \note this function has to be called prior to flux calculations on the element.
     */
    template<class FVElementGeometry, class ElementVolumeVariables>
    void bindElement_(const typename FVElementGeometry::GridGeometry::GridView::template Codim<0>::Entity& element,
                      const FVElementGeometry& fvGeometry,
                      const ElementVolumeVariables& elemVolVars)
    {
        // resizing of the cache
        const auto numScvf = fvGeometry.numScvf();
        fluxVarsCache_.resize(numScvf);
        globalScvfIndices_.resize(numScvf);

        // instantiate helper class to fill the caches
        FluxVariablesCacheFiller filler(gridFluxVarsCache().problem());

        std::size_t localScvfIdx = 0;
        // fill the containers
        for (auto&& scvf : scvfs(fvGeometry))
        {
            filler.fill(*this, fluxVarsCache_[localScvfIdx], element, fvGeometry, elemVolVars, scvf, true);
            globalScvfIndices_[localScvfIdx] = scvf.index();
            localScvfIdx++;
        }
    }

    /*!
     * \brief Prepares the transmissibilities of the scv faces in the stencil of an element
     * \note the fvGeometry is assumed to be bound to the same element
     * \note this function has to be called prior to flux calculations on the element.
     */
    template<class FVElementGeometry, class ElementVolumeVariables>
    void bind_(const typename FVElementGeometry::GridGeometry::GridView::template Codim<0>::Entity& element,
               const FVElementGeometry& fvGeometry,
               const ElementVolumeVariables& elemVolVars)
    {
        const auto& problem = gridFluxVarsCache().problem();
        const auto& gridGeometry = fvGeometry.gridGeometry();
        const auto globalI = gridGeometry.elementMapper().index(element);
        const auto& connectivityMapI = gridGeometry.connectivityMap()[globalI];
        const auto numNeighbors = connectivityMapI.size();

        // instantiate helper class to fill the caches
        FluxVariablesCacheFiller filler(problem);

        // find the number of scv faces that need to be prepared
        auto numScvf = fvGeometry.numScvf();
        for (unsigned int localIdxJ = 0; localIdxJ < numNeighbors; ++localIdxJ)
            numScvf += connectivityMapI[localIdxJ].scvfsJ.size();

        // fill the containers with the data on the scv faces inside the actual element
        fluxVarsCache_.resize(numScvf);
        globalScvfIndices_.resize(numScvf);
        unsigned int localScvfIdx = 0;
        for (auto&& scvf : scvfs(fvGeometry))
        {
            filler.fill(*this, fluxVarsCache_[localScvfIdx], element, fvGeometry, elemVolVars, scvf, true);
            globalScvfIndices_[localScvfIdx] = scvf.index();
            localScvfIdx++;
        }

        // add required data on the scv faces in the neighboring elements
        for (unsigned int localIdxJ = 0; localIdxJ < numNeighbors; ++localIdxJ)
        {
            const auto elementJ = gridGeometry.element(connectivityMapI[localIdxJ].globalJ);
            for (auto scvfIdx : connectivityMapI[localIdxJ].scvfsJ)
            {
                auto&& scvfJ = fvGeometry.scvf(scvfIdx);
                filler.fill(*this, fluxVarsCache_[localScvfIdx], elementJ, fvGeometry, elemVolVars, scvfJ, true);
                globalScvfIndices_[localScvfIdx] = scvfJ.index();
                localScvfIdx++;
            }
        }
    }

    /*!
     * \brief Prepares the transmissibilities of a single scv face
     * \note the fvGeometry is assumed to be bound to the same element
     * \note this function has to be called prior to flux calculations on the element.
     */
    template<class FVElementGeometry, class ElementVolumeVariables>
    void bindScvf_(const typename FVElementGeometry::GridGeometry::GridView::template Codim<0>::Entity& element,
                   const FVElementGeometry& fvGeometry,
                   const ElementVolumeVariables& elemVolVars,
                   const typename FVElementGeometry::SubControlVolumeFace& scvf)
    {
        fluxVarsCache_.resize(1);
        globalScvfIndices_.resize(1);

        // instantiate helper class to fill the caches
        FluxVariablesCacheFiller filler(gridFluxVarsCache().problem());

        filler.fill(*this, fluxVarsCache_[0], element, fvGeometry, elemVolVars, scvf, true);
        globalScvfIndices_[0] = scvf.index();
    }

    const GridFluxVariablesCache* gridFluxVarsCachePtr_;

    //! get index of scvf in the local container
    int getLocalScvfIdx_(const int scvfIdx) const
    {
        auto it = std::find(globalScvfIndices_.begin(), globalScvfIndices_.end(), scvfIdx);
        assert(it != globalScvfIndices_.end() && "Could not find the flux vars cache for scvfIdx");
        return std::distance(globalScvfIndices_.begin(), it);
    }

    std::vector<FluxVariablesCache> fluxVarsCache_;
    std::vector<std::size_t> globalScvfIndices_;
};

} // end namespace Dumux

#endif
