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
 * \ingroup CvfeDiscretization
 * \brief Global flux variable cache
 */
#ifndef DUMUX_DISCRETIZATION_CVFE_GRID_FLUXVARSCACHE_HH
#define DUMUX_DISCRETIZATION_CVFE_GRID_FLUXVARSCACHE_HH

// make the local view function available whenever we use this class
#include <dumux/discretization/localview.hh>
#include <dumux/discretization/cvfe/elementfluxvariablescache.hh>
#include <dumux/discretization/box/gridfluxvariablescache.hh>

namespace Dumux {

/*!
 * \ingroup CvfeDiscretization
 * \brief Flux variable caches traits
 */
template<class P, class FVC>
struct CvfeDefaultGridFVCTraits
{
    using Problem = P;
    using FluxVariablesCache = FVC;

    template<class GridFluxVariablesCache, bool cachingEnabled>
    using LocalView = CvfeElementFluxVariablesCache<GridFluxVariablesCache, cachingEnabled>;
};

/*!
 * \ingroup CvfeDiscretization
 * \brief Flux variable caches on a gridview
 * \note The class is specialized for a version with and without grid caching
 */
template<class Problem,
         class FluxVariablesCache,
         bool cachingEnabled = false,
         class Traits = CvfeDefaultGridFVCTraits<Problem, FluxVariablesCache> >
using CvfeGridFluxVariablesCache = BoxGridFluxVariablesCache<Problem, FluxVariablesCache, cachingEnabled, Traits>;

} // end namespace Dumux

#endif
