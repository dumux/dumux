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
 * \ingroup CCWMpfaDiscretization
 * \brief The flux variables caches for an element
 */
#ifndef DUMUX_DISCRETIZATION_CCWMPFA_ELEMENT_FLUXVARSCACHE_HH
#define DUMUX_DISCRETIZATION_CCWMPFA_ELEMENT_FLUXVARSCACHE_HH

#include <algorithm>
#include <cassert>
#include <vector>
#include <dune/common/exceptions.hh>

namespace Dumux {

template<class InterpolationOperator, class FaceDataHandle>
struct DataStorage
{
    std::vector<InterpolationOperator> interpolationOperators;
    std::vector<FaceDataHandle> facesDataHandles;
};

/*!
 * \ingroup CCWMpfaDiscretization
 * \brief The flux variables caches for an element
 * \note The class is specialized for a version with and without caching
 * If grid caching is enabled the flux caches are stored for the whole gridview in the corresponding
 * GridFluxVariablesCache which is memory intensive but faster. For caching disabled the
 * flux caches are locally computed for each element whenever needed.
 */
template<class GFVC, bool cachingEnabled>
class CCWMpfaElementFluxVariablesCache;

/*!
 * \ingroup CCWMpfaDiscretization
 * \brief The flux variables caches for an element with caching enabled
 */
template<class GFVC>
class CCWMpfaElementFluxVariablesCache<GFVC, true>
{
    //! the type of the flux variables cache filler
    using FluxVariablesCacheFiller = typename GFVC::Traits::FluxVariablesCacheFiller;

public:
    //! export the type of the grid flux variables cache
    using GridFluxVariablesCache = GFVC;

    //! export the type of the flux variables cache
    using FluxVariablesCache = typename GFVC::FluxVariablesCache;

        //! export the data handle types used
    using InterpolationOperator = typename GFVC::InterpolationOperator;
    using FaceDataHandle = typename GFVC::FaceDataHandle;

    CCWMpfaElementFluxVariablesCache(const GridFluxVariablesCache& global)
    : gridFluxVarsCachePtr_(&global) {}

    //! Specialization for the global caching being enabled - do nothing here
    template<class FVElementGeometry, class ElementVolumeVariables>
    void bindElement(const typename FVElementGeometry::GridGeometry::GridView::template Codim<0>::Entity& element,
                     const FVElementGeometry& fvGeometry,
                     const ElementVolumeVariables& elemVolVars) {}

    //! Specialization for the global caching being enabled - do nothing here
    template<class FVElementGeometry, class ElementVolumeVariables>
    void bind(const typename FVElementGeometry::GridGeometry::GridView::template Codim<0>::Entity& element,
              const FVElementGeometry& fvGeometry,
              const ElementVolumeVariables& elemVolVars) {}

    //! Specialization for the global caching being enabled - do nothing here
    template<class FVElementGeometry, class ElementVolumeVariables>
    void bindScvf(const typename FVElementGeometry::GridGeometry::GridView::template Codim<0>::Entity& element,
                  const FVElementGeometry& fvGeometry,
                  const ElementVolumeVariables& elemVolVars,
                  const typename FVElementGeometry::SubControlVolumeFace& scvf) {}

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

} // end namespace Dumux

#endif
