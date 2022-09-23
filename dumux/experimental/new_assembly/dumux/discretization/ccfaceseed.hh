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
 * \ingroup CCDiscretization
 * \copydoc Dumux::CCFaceSeed
 */
#ifndef DUMUX_DISCRETIZATION_CC_FACE_SEED_HH
#define DUMUX_DISCRETIZATION_CC_FACE_SEED_HH

#include <ranges>
#include <cassert>
#include <cstdint>
#include <concepts>
#include <type_traits>

#include <dumux/experimental/new_assembly/dumux/common/size.hh>
#include <dumux/experimental/new_assembly/dumux/common/storage.hh>
#include <dumux/experimental/new_assembly/dumux/discretization/facet.hh>

namespace Dumux {

/*!
 * \ingroup CCDiscretization
 * \brief Stores connectivity information for a grid face.
 */
template<Concepts::Size auto maxNumFaceNeighbors = Size::dynamic,
         std::integral GridIndex = std::size_t,
         std::integral LocalIndex = std::uint_least8_t>
class CCFaceSeed
{
public:
    using Facet = Dumux::Facet<GridIndex, LocalIndex>;

    //! Constructor from an element facet
    explicit CCFaceSeed(Facet facet, bool boundary) noexcept(std::is_nothrow_move_constructible_v<Facet>)
    : facets_({std::move(facet)})
    , boundary_(boundary)
    {}

    void addOutsideFacet(Facet&& facet)
    { facets_.push_back(std::move(facet)); }

    bool onBoundary() const { return boundary_; }
    std::size_t numNeighbors() const { return facets_.size(); }
    std::size_t numOutsideNeighbors() const { return numNeighbors() - 1; }

    const auto& facets() const  { return facets_; }
    const Facet& facet(unsigned int i) const { return facets_[i]; }
    const Facet& outsideFacet(unsigned int i = 0) const { return facet(i+1); }
    const Facet& insideFacet() const { return facets_[0]; }

    std::ranges::range auto outsideFacets() const
    {
        return std::views::transform(
            std::views::iota(std::size_t{0}, numOutsideNeighbors()),
            [&] (std::size_t i) { return outsideFacet(i); }
        );
    }

private:
    DefaultStorage<Facet, maxNumFaceNeighbors> facets_;
    bool boundary_;
};

} // end namespace Dumux

#endif
