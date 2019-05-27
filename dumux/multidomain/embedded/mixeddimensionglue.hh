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

#include <dumux/multidomain/glue.hh>
#include <dune/grid/common/mcmgmapper.hh>

namespace Dumux {

// forward declaration
template<class BulkGridView, class LowDimGridView,
         class BulkMapper = Dune::MultipleCodimMultipleGeomTypeMapper<BulkGridView>,
         class LowDimMapper = Dune::MultipleCodimMultipleGeomTypeMapper<LowDimGridView>>
class DUNE_DEPRECATED_MSG("Use MultiDomainGlue instead. Will be removed after 3.1!")
MixedDimensionGlue
: public MultiDomainGlue<LowDimGridView, BulkGridView, LowDimMapper, BulkMapper>
{
    using ParentType = MultiDomainGlue<LowDimGridView, BulkGridView, LowDimMapper, BulkMapper>;
public:
    //! Default constructor
    MixedDimensionGlue() = default;

    /*!
     * \brief constructor from bounding box trees
     * \note The definition of Glue::Intersection::inside() has changed!
     *       inside() used to return the lower-dimensional entities
     *       (i.e. from LowDimGridView). Now, inside() returns the entities
     *       of the template argument DomainGridView, i.e. the Dune::GridView
     *       given as the first template argument. In order to be backwards
     *       compatible, we flip the order here (and in the definition of ParentType).
     */
    template<class BulkTree, class LowDimTree>
    MixedDimensionGlue(const BulkTree& bulkTree, const LowDimTree& lowDimTree)
    : ParentType(lowDimTree, bulkTree)
    {}

    /*!
     * \brief construction from bounding box trees
     * \note The definition of Glue::Intersection::inside() has changed!
     *       inside() used to return the lower-dimensional entities
     *       (i.e. from LowDimGridView). Now, inside() returns the entities
     *       of the template argument DomainGridView, i.e. the Dune::GridView
     *       given as the first template argument. In order to be backwards
     *       compatible, we flip the order here (and in the definition of ParentType).
     */
    template<class BulkTree, class LowDimTree>
    void build(const BulkTree& bulkTree, const LowDimTree& lowDimTree)
    {
        ParentType::build(lowDimTree, bulkTree);
    }
};

} // end namespace Dumux

#endif
