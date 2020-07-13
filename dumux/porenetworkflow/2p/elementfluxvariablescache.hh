// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   See the file COPYING for full copying permissions.                      *
 *                                                                           *
 *   This program is free software: you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation, either version 3 of the License, or       *
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
 * \ingroup BoxDiscretization
 * \brief Global flux variable cache
 */
#ifndef DUMUX_PNM_2P_ELEMNT_FLUXVARSCACHE_HH
#define DUMUX_PNM_2P_ELEMNT_FLUXVARSCACHE_HH

#include <cstddef>
#include <vector>
#include <dumux/discretization/box/elementfluxvariablescache.hh>

namespace Dumux {

/*!
 * \ingroup BoxDiscretization
 * \brief The flux variables caches for an element
 * \note The class is specialized for a version with and without caching
 * If grid caching is enabled the flux caches are stored for the whole gridview in the corresponding
 * GridFluxVariablesCache which is memory intensive but faster. For caching disabled the
 * flux caches are locally computed for each element whenever needed.
 */
template<class GFVC, bool cachingEnabled>
class PNMTwoPElementFluxVariablesCache;

/*!
 * \ingroup BoxDiscretization
 * \brief The flux variables caches for an element with caching enabled
 */
template<class GFVC>
class PNMTwoPElementFluxVariablesCache<GFVC, true> : public BoxElementFluxVariablesCache<GFVC, true>
{
    using ParentType = BoxElementFluxVariablesCache<GFVC, true>;
public:
    using ParentType::ParentType;
};

/*!
 * \ingroup BoxDiscretization
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
              const ElementVolumeVariables& elemVolVars)
    {
        bindElement(element, fvGeometry, elemVolVars);
    }

    template<class FVElementGeometry, class ElementVolumeVariables>
    void bindElement(const typename FVElementGeometry::GridGeometry::GridView::template Codim<0>::Entity& element,
                     const FVElementGeometry& fvGeometry,
                     const ElementVolumeVariables& elemVolVars)
    {
        // temporary resizing of the cache
        fluxVarsCache_.resize(fvGeometry.numScvf());
        for (auto&& scvf : scvfs(fvGeometry))
            (*this)[scvf].update(gridFluxVarsCache().problem(), element, fvGeometry, elemVolVars, scvf, gridFluxVarsCache().invasionState().invaded(element));
    }

    template<class FVElementGeometry, class ElementVolumeVariables>
    void bindScvf(const typename FVElementGeometry::GridGeometry::GridView::template Codim<0>::Entity& element,
                  const FVElementGeometry& fvGeometry,
                  const ElementVolumeVariables& elemVolVars,
                  const typename FVElementGeometry::SubControlVolumeFace& scvf)
    {
        fluxVarsCache_.resize(fvGeometry.numScvf());
        (*this)[scvf].update(gridFluxVarsCache().problem(), element, fvGeometry, elemVolVars, scvf, gridFluxVarsCache().invasionState().invaded(element));
    }

    // access operator
    template<class SubControlVolumeFace>
    const FluxVariablesCache& operator [](const SubControlVolumeFace& scvf) const
    { return fluxVarsCache_[scvf.index()]; }

    // access operator
    template<class SubControlVolumeFace>
    FluxVariablesCache& operator [](const SubControlVolumeFace& scvf)
    { return fluxVarsCache_[scvf.index()]; }

    //! The global object we are a restriction of
    const GridFluxVariablesCache& gridFluxVarsCache() const
    {  return *gridFluxVarsCachePtr_; }

private:
    const GridFluxVariablesCache* gridFluxVarsCachePtr_;
    std::vector<FluxVariablesCache> fluxVarsCache_;
};

} // end namespace Dumux

#endif
