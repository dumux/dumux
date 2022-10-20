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
 * \ingroup Discretization
 * \copydoc Dumux::IndexedIntersections
 */
#ifndef DUMUX_DISCRETIZATION_INDEXED_INTERSECTIONS_HH
#define DUMUX_DISCRETIZATION_INDEXED_INTERSECTIONS_HH

#include <memory>
#include <concepts>
#include <algorithm>
#include <ranges>
#include <vector>
#include <tuple>
#include <optional>

#include <dune/common/float_cmp.hh>
#include <dune/common/iteratorfacades.hh>
#include <dumux/geometry/diameter.hh>

#include <dumux/experimental/new_assembly/dumux/common/storage.hh>
#include <dumux/experimental/new_assembly/dumux/common/size.hh>
#include <dumux/experimental/new_assembly/dumux/discretization/concepts.hh>

namespace Dumux {

#ifndef DOXYGEN
namespace Detail {

template<std::forward_iterator IndexIterator,
         std::forward_iterator IntersectionIterator>
using IndexedIntersectionIteratorValueType = std::tuple<
    std::iter_reference_t<IndexIterator>,
    std::iter_reference_t<IntersectionIterator>
>;

template<std::forward_iterator FaceIndexIterator,
         std::forward_iterator IntersectionIterator>
class IndexedIntersectionIterator
: public Dune::ForwardIteratorFacade<
    IndexedIntersectionIterator<FaceIndexIterator, IntersectionIterator>,
    IndexedIntersectionIteratorValueType<FaceIndexIterator, IntersectionIterator>,
    IndexedIntersectionIteratorValueType<FaceIndexIterator, IntersectionIterator>>
{
    using VT = IndexedIntersectionIteratorValueType<FaceIndexIterator, IntersectionIterator>;

public:
    IndexedIntersectionIterator(FaceIndexIterator indexIterator,
                                IntersectionIterator intersectionIterator)
    : indexIterator_(indexIterator)
    , intersectionIterator_(intersectionIterator)
    {}

    void increment()
    {
        ++indexIterator_;
        ++intersectionIterator_;
    }

    bool equals(const IndexedIntersectionIterator& other) const
    { return intersectionIterator_ == other.intersectionIterator_; }

    VT dereference() const
    { return {*indexIterator_, *intersectionIterator_}; }

private:
    FaceIndexIterator indexIterator_;
    IntersectionIterator intersectionIterator_;
};

} // namespace Detail
#endif // DOXYGEN

/*!
 * \ingroup Discretization
 * \brief Iterator range over tuples of (element-local index, intersection).
 *        On network grids, multiple intersections may be defined at the same geometric
 *        location (with different neighbors). In this case, geometrically identical
 *        intersections map to the same element-local index.
 * \note On network grids, in index map is set up internally upon construction with
 *       an element. Since this may involve dynamic memory allocation (depending on
 *       the template argument `maxNumElementIntersections`), it can be more efficient
 *       to construct this class only once and then set new elements with the provided
 *       `set` function.
 * \note The second template argument can be used to specify a maximum number of
 *       intersections per element. This template argument only has an effect for
 *       the specialization on network grids, where providing an upper bound for
 *       the number of intersections per element at compile-time can lead to a
 *       speedup by avoiding dynamic memory allocation.
 * \note This identifies geometrically identical intersections by checking their center
 *       points for fuzzy equality.
 */
template<Concepts::GridView GV,
         Concepts::Size auto maxNumElementIntersections = Size::dynamic>
class IndexedIntersections
{
    using Element = typename GV::template Codim<0>::Entity;
    using Indices = std::ranges::iota_view<int>;
    using Intersections = decltype(intersections(
        std::declval<const GV&>(),
        std::declval<const Element&>()
    ));

    static_assert(std::ranges::forward_range<Indices>);
    static_assert(std::ranges::forward_range<Intersections>);

public:
    using GridView = GV;
    using Iterator = Detail::IndexedIntersectionIterator<
        std::ranges::iterator_t<Indices>,
        std::ranges::iterator_t<Intersections>
    >;

    IndexedIntersections(const GridView& gv) : gridView_(gv) {}
    IndexedIntersections(const GridView& gv, const Element& element)
    : gridView_(gv)
    { set(element); }

    //! Prepare the intersection range for the given element
    void set(const Element& element)
    {
        indices_ = std::make_unique<Indices>(std::views::iota(0));
        intersections_ = std::make_unique<Intersections>(intersections(gridView_, element));
    }

    Iterator begin() const { return {indices_->begin(), intersections_->begin()}; }
    Iterator end() const { return {indices_->begin(), intersections_->end()}; }

private:
    const GridView& gridView_;
    std::unique_ptr<Indices> indices_{nullptr};
    std::unique_ptr<Intersections> intersections_{nullptr};
};

/*!
 * \ingroup Discretization
 * \brief Specialization for network grids
 * \note While not in line with the specialization for non-network grids,
 *       this specialization allows to (optionally) pass in a tolerance
 *       to be used when identifying geometrically identical intersections.
 *       If none is passed, a default tolerance is computed from the size
 *       of the given element. The main reason for this constructor argument
 *       is to speed up computations, but be careful to choose an adequate
 *       tolerance for the entire grid.
 */
template<Concepts::NetworkGridView GV,
         Concepts::Size auto maxNumElementIntersections>
class IndexedIntersections<GV, maxNumElementIntersections>
{
    using ctype = typename GV::ctype;
    using Element = typename GV::template Codim<0>::Entity;
    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;

    using Indices = DefaultStorage<int, maxNumElementIntersections>;
    using Intersections = decltype(intersections(
        std::declval<const GV&>(),
        std::declval<const Element&>()
    ));

    static_assert(std::ranges::forward_range<Indices>);
    static_assert(std::ranges::forward_range<Intersections>);

public:
    using GridView = GV;
    using Iterator = Detail::IndexedIntersectionIterator<
        typename Indices::const_iterator,
        std::ranges::iterator_t<Intersections>
    >;

    IndexedIntersections(const GridView& gv,
                         std::optional<ctype> tolerance = {})
    : gridView_(gv)
    , tolerance_(tolerance)
    {}

    IndexedIntersections(const GridView& gv,
                         const Element& element,
                         std::optional<ctype> tolerance = {})
    : gridView_(gv)
    , tolerance_(tolerance)
    { set(element); }

    //! Prepare the intersection range for the given element
    void set(const Element& element)
    {
        resetIndices_(element);
        setIndices_(element);
        intersections_ = std::make_unique<Intersections>(intersections(gridView_, element));
    }

    Iterator begin() const { return {indices_.begin(), intersections_->begin()}; }
    Iterator end() const { return {indices_.begin(), intersections_->end()}; }

private:
    void resetIndices_(const Element& element)
    {
        indices_.clear();
        centers_.clear();
        indices_.reserve(element.subEntities(1));
        centers_.reserve(element.subEntities(1));
    }

    void setIndices_(const Element& element)
    {
        const ctype tolerance = tolerance_.has_value() ? *tolerance_ : 1e-6*diameter(element.geometry());
        for (const auto& is : intersections(gridView_, element))
        {
            const auto& center = is.geometry().center();
            if (auto fIdx = getIndex_(center, tolerance); fIdx)
                indices_.push_back(*fIdx);
            else
            {
                indices_.push_back(centers_.size());
                centers_.push_back(center);
            }
        }
    }

    std::optional<int> getIndex_(const GlobalPosition& center, ctype tolerance) const
    {
        auto it = std::ranges::find_if(centers_, [&] (const GlobalPosition& c) {
            return Dune::FloatCmp::eq(center, c, tolerance);
        });
        if (it != std::ranges::end(centers_))
            return std::distance(std::ranges::begin(centers_), it);
        return std::nullopt;
    }

    const GridView& gridView_;
    std::optional<ctype> tolerance_;

    Indices indices_;
    std::vector<GlobalPosition> centers_;
    std::unique_ptr<Intersections> intersections_{nullptr};
};

/*!
 * \ingroup Discretization
 * \brief Convenience function to return a range over tuples of (element-local-index, intersection).
 */
template<Concepts::GridView GV>
auto indexedIntersections(const GV& gridView, const typename GV::template Codim<0>::Entity& element)
{ return IndexedIntersections{gridView, element}; }

} // namespace Dumux

#endif
