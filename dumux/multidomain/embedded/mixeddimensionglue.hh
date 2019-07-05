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
 * \ingroup MixedDimension
 * \brief A class glueing two grids of different dimension geometrically
 *        Intersections are computed using axis-aligned bounding box trees
 */
#ifndef DUMUX_MULTIDOMAIN_MIXEDDIMENSION_GLUE_HH
#define DUMUX_MULTIDOMAIN_MIXEDDIMENSION_GLUE_HH

#warning "This header is deprecated and will be removed after release 3.1. Use multidomain/glue.hh"

#include <dumux/multidomain/glue.hh>
#include <dune/grid/common/mcmgmapper.hh>

namespace Dumux {

// forward declaration
template<class BulkGridView, class LowDimGridView,
         class BulkMapper = Dune::MultipleCodimMultipleGeomTypeMapper<BulkGridView>,
         class LowDimMapper = Dune::MultipleCodimMultipleGeomTypeMapper<LowDimGridView>>
using MixedDimensionGlue [[deprecated("Use MultiDomainGlue instead. Will be removed after 3.1!")]]
    = MultiDomainGlue<BulkGridView, LowDimGridView, BulkMapper, LowDimMapper>;

} // end namespace Dumux

#endif
