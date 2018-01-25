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
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the           *
 *   GNU General Public License for more details.                            *
 *                                                                           *
 *   You should have received a copy of the GNU General Public License       *
 *   along with this program.  If not, see <http://www.gnu.org/licenses/>.   *
 *****************************************************************************/
/*!
 * \file
 * \ingroup Common
 * \brief A class storing intersections from intersecting two bounding box trees
 */
#ifndef DUMUX_BOUNDING_BOX_TREE_INTERSECTION_HH
#define DUMUX_BOUNDING_BOX_TREE_INTERSECTION_HH

#include <dune/common/fvector.hh>
#include <dune/common/promotiontraits.hh>

namespace Dumux
{

/*!
 * \ingroup Common
 * \brief An intersection object resulting from the intersection of two bounding box tree primitives
 *
 * After is has been found that two leaf bounding boxes intersect a primitive test has to be
 * performed to see if the actual entities inside the bounding box intersect too. The result
 * if such an intersection is found as an object of this class containing the indices of the
 * intersecting entities and the corners of the intersection object.
 */
template<class EntitySet0, class EntitySet1>
class BoundingBoxTreeIntersection
{
    enum { dimworld = EntitySet0::dimensionworld };
    using ctype = typename Dune::PromotionTraits<typename EntitySet0::ctype, typename EntitySet1::ctype>::PromotedType;
    using GlobalPosition = Dune::FieldVector<ctype, dimworld>;

public:
    template<class Corners>
    explicit BoundingBoxTreeIntersection(std::size_t a,
                                         std::size_t b,
                                         Corners&& c)
    : a_(a)
    , b_(b)
    , corners_(c.begin(), c.end())
    {
        static_assert(int(EntitySet0::dimensionworld) == int(EntitySet1::dimensionworld),
          "Can only store intersections of entity sets with the same world dimension");
    }

    //! Get the index of the intersecting entity belonging to this grid
    std::size_t first() const
    { return a_; }

    //! Get the index of the intersecting entity belonging to the other grid
    std::size_t second() const
    { return b_; }

    //! Get the corners of the intersection geometry
    std::vector<GlobalPosition> corners() const
    { return corners_; }

    /*!
     * \brief Check if the corners of this intersection match with the given corners
     * \note This is useful to check if the intersection geometry of two intersections coincide.
     */
    bool cornersMatch(const std::vector<GlobalPosition>& otherCorners) const
    {
        if (otherCorners.size() != corners_.size())
            return false;

        const auto eps = 1.5e-7*(corners_[1] - corners_[0]).two_norm();
        for (int i = 0; i < corners_.size(); ++i)
            if ((corners_[i] - otherCorners[i]).two_norm() > eps)
                return false;

        return true;
    }

private:
    std::size_t a_, b_; //!< Indices of the intersection elements
    std::vector<GlobalPosition> corners_; //!< the corner points of the intersection geometry
};

} // end namespace Dumux

#endif
