// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup PNMTwoPModel
 * \brief Element flux variable cache
 */
#ifndef DUMUX_PNM_2P_ELEMNT_FLUXVARSCACHE_HH
#define DUMUX_PNM_2P_ELEMNT_FLUXVARSCACHE_HH

#include <cstddef>
#include <vector>
#include <utility>

#include <dumux/discretization/cvfe/elementfluxvariablescache.hh>

namespace Dumux::PoreNetwork {

/*!
 * \ingroup PNMTwoPModel
 * \brief The flux variables caches for an element
 */
template<class GFVC, bool cachingEnabled>
class PNMTwoPElementFluxVariablesCache;

/*!
 * \ingroup PNMTwoPModel
 * \brief The flux variables caches for an element with caching enabled
 */
template<class GFVC>
class PNMTwoPElementFluxVariablesCache<GFVC, true> : public CVFEElementFluxVariablesCache<GFVC, true>
{
    using ParentType = CVFEElementFluxVariablesCache<GFVC, true>;
public:
    using ParentType::ParentType;
};

/*!
 * \ingroup PNMTwoPModel
 * \brief The flux variables caches for an element with caching disabled
 */
template<class GFVC>
class PNMTwoPElementFluxVariablesCache<GFVC, false>
{
public:
    //! export the type of the grid flux variables cache
    using GridFluxVariablesCache = GFVC;

    //! export the type of the flux variables cache
    using FluxVariablesCache = typename GFVC::FluxVariablesCache;

    PNMTwoPElementFluxVariablesCache(const GridFluxVariablesCache& global)
    : gridFluxVarsCachePtr_(&global) {}

    // Function is called by the BoxLocalJacobian prior to flux calculations on the element.
    // We assume the FVGeometries to be bound at this point
    template<class FVElementGeometry, class ElementVolumeVariables>
    void bind(const typename FVElementGeometry::GridGeometry::GridView::template Codim<0>::Entity& element,
              const FVElementGeometry& fvGeometry,
              const ElementVolumeVariables& elemVolVars) &
    { bindElement(element, fvGeometry, elemVolVars); }

    /*!
    * \brief bind the local view (r-value overload)
    * This overload is called when an instance of this class is a temporary in the usage context
    * This allows a usage like this: `const auto view = localView(...).bind(element);`
    */
    template<class FVElementGeometry, class ElementVolumeVariables>
    PNMTwoPElementFluxVariablesCache bind(const typename FVElementGeometry::GridGeometry::GridView::template Codim<0>::Entity& element,
                                          const FVElementGeometry& fvGeometry,
                                          const ElementVolumeVariables& elemVolVars) &&
    {
        this->bind(element, fvGeometry, elemVolVars);
        return std::move(*this);
    }

    template<class FVElementGeometry, class ElementVolumeVariables>
    void bindElement(const typename FVElementGeometry::GridGeometry::GridView::template Codim<0>::Entity& element,
                     const FVElementGeometry& fvGeometry,
                     const ElementVolumeVariables& elemVolVars) &
    {
        for (auto&& scvf : scvfs(fvGeometry))
            fluxVarsCache_.update(gridFluxVarsCache().problem(), element, fvGeometry, elemVolVars, scvf, gridFluxVarsCache().invasionState().invaded(element));
    }

    /*!
    * \brief bind the local view (r-value overload)
    * This overload is called when an instance of this class is a temporary in the usage context
    * This allows a usage like this: `const auto view = localView(...).bind(element);`
    */
    template<class FVElementGeometry, class ElementVolumeVariables>
    PNMTwoPElementFluxVariablesCache bindElement(const typename FVElementGeometry::GridGeometry::GridView::template Codim<0>::Entity& element,
                                                 const FVElementGeometry& fvGeometry,
                                                 const ElementVolumeVariables& elemVolVars) &&
    {
        this->bindElement(element, fvGeometry, elemVolVars);
        return std::move(*this);
    }

    template<class FVElementGeometry, class ElementVolumeVariables>
    void bindScvf(const typename FVElementGeometry::GridGeometry::GridView::template Codim<0>::Entity& element,
                  const FVElementGeometry& fvGeometry,
                  const ElementVolumeVariables& elemVolVars,
                  const typename FVElementGeometry::SubControlVolumeFace& scvf) &
    {
        fluxVarsCache_.update(gridFluxVarsCache().problem(), element, fvGeometry, elemVolVars, scvf, gridFluxVarsCache().invasionState().invaded(element));
    }

    /*!
    * \brief bind the local view (r-value overload)
    * This overload is called when an instance of this class is a temporary in the usage context
    * This allows a usage like this: `const auto view = localView(...).bind(element);`
    */
    template<class FVElementGeometry, class ElementVolumeVariables>
    PNMTwoPElementFluxVariablesCache bindScvf(const typename FVElementGeometry::GridGeometry::GridView::template Codim<0>::Entity& element,
                                              const FVElementGeometry& fvGeometry,
                                              const ElementVolumeVariables& elemVolVars,
                                              const typename FVElementGeometry::SubControlVolumeFace& scvf) &&
    {
        this->bindScvf(element, fvGeometry, elemVolVars, scvf);
        return std::move(*this);
    }

    /*!
     * \brief Update the caches if the volume variables have changed and the cache is solution-dependent
     * \note Results in undefined behaviour if called before bind() or with a different element
     */
    template<class FVElementGeometry, class ElementVolumeVariables>
    void update(const typename FVElementGeometry::Element& element,
                const FVElementGeometry& fvGeometry,
                const ElementVolumeVariables& elemVolVars)
    {
        if constexpr (FluxVariablesCache::isSolDependent)
        {
            for (const auto& scvf : scvfs(fvGeometry))
                fluxVarsCache_.update(
                    gridFluxVarsCache().problem(), element, fvGeometry, elemVolVars, scvf,
                    gridFluxVarsCache().invasionState().invaded(element)
                );
        }
    }

    //! access operator
    template<class SubControlVolumeFace>
    const FluxVariablesCache& operator [](const SubControlVolumeFace& scvf) const
    { return fluxVarsCache_; }

    //! access operator
    template<class SubControlVolumeFace>
    FluxVariablesCache& operator [](const SubControlVolumeFace& scvf)
    { return fluxVarsCache_; }

    //! The global object we are a restriction of
    const GridFluxVariablesCache& gridFluxVarsCache() const
    {  return *gridFluxVarsCachePtr_; }

private:
    const GridFluxVariablesCache* gridFluxVarsCachePtr_;
    FluxVariablesCache fluxVarsCache_;
};

} // end namespace Dumux::PoreNetwork

#endif
