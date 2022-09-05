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
 * \ingroup CCTpfaDiscretization
 * \brief TODO: Doc me
 */
#ifndef DUMUX_DISCRETIZATION_CCTPFA_SCVF_SEED_HH
#define DUMUX_DISCRETIZATION_CCTPFA_SCVF_SEED_HH

#include <type_traits>
#include <concepts>
#include <cstdint>
#include <ranges>

#include <dumux/experimental/new_assembly/dumux/common/concepts.hh>
#include <dumux/experimental/new_assembly/dumux/common/storage.hh>

namespace Dumux::CCTpfa {

template<std::integral Index = std::size_t,
         std::integral LocalIndex = std::uint_least8_t>
struct Facet
{
    Index elementIndex;
    LocalIndex facetIndex;
};

template<auto maxNumFaceNeighbors = Size::dynamic,
         std::integral Index = std::size_t,
         std::integral LocalIndex = std::uint_least8_t>
class ScvfSeed
{
public:
    using Facet = CCTpfa::Facet<Index, LocalIndex>;

    //! Constructor from a primary face seed
    explicit ScvfSeed(Facet facet) noexcept(std::is_nothrow_move_constructible_v<Facet>)
    : facets_({std::move(facet)})
    {}

    const Facet& facet() const
    { return facets_[0]; }

    std::size_t numNeighbors() const
    { return facets_.size(); }

    const auto& neighboringFacets() const
    { return facets_; }

    std::ranges::range auto outsideFacets() const
    {
        return std::views::transform(
            std::views::iota(std::size_t{1}, facets_.size()),
            [&] (std::size_t i) { return facets_[i]; }
        );
    }

    void addOutsideFacet(Facet&& facet)
    { facets_.push_back(std::move(facet)); }

private:
    DefaultStorage<Facet, maxNumFaceNeighbors> facets_;
};

} // end namespace Dumux::CCTpfa

#endif
