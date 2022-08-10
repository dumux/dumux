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
#ifndef DUMUX_DISCRETIZATION_BOX_GRID_FLUXVARSCACHE_HH
#define DUMUX_DISCRETIZATION_BOX_GRID_FLUXVARSCACHE_HH

#warning "This header is deprecated and will be removed after 3.6"

#include <dumux/discretization/cvfe/gridfluxvariablescache.hh>

namespace Dumux {

/*!
 * \ingroup BoxDiscretization
 * \brief Flux variable caches traits
 */
template<class P, class FVC>
struct BoxDefaultGridFVCTraits
{
    using Problem = P;
    using FluxVariablesCache = FVC;

    template<class GridFluxVariablesCache, bool cachingEnabled>
    using LocalView = CVFEElementFluxVariablesCache<GridFluxVariablesCache, cachingEnabled>;
};

/*!
 * \ingroup BoxDiscretization
 * \brief Flux variable caches on a gridview
 * \note The class is specialized for a version with and without grid caching
 */
template<class Problem,
         class FluxVariablesCache,
         bool cachingEnabled = false,
         class Traits = CVFEDefaultGridFVCTraits<Problem, FluxVariablesCache> >
using BoxGridFluxVariablesCache [[deprecated("Will be removed after 3.6")]]
    = CVFEGridFluxVariablesCache<Problem, FluxVariablesCache, cachingEnabled, Traits>;

} // end namespace Dumux

#endif
