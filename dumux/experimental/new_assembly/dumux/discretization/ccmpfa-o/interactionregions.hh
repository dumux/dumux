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
 * \ingroup CCMpfaDiscretization
 * \copydoc Dumux::CCMpfaOInteractionRegions
 */
#ifndef DUMUX_DISCRETIZATION_CCMPFA_O_INTERACTION_REGIONS_HH
#define DUMUX_DISCRETIZATION_CCMPFA_O_INTERACTION_REGIONS_HH

#include <vector>
#include <cstdint>
#include <utility>
#include <cassert>

namespace Dumux::CCMpfaO {

struct SubCell
{
    std::size_t cellIndex;
    std::uint_least8_t localCornerIndex;
};

class InteractionRegion
{
public:
    void add(SubCell&& subCell)
    { subCells_.emplace_back(subCell); }

    auto size() const { return subCells_.size(); }

    auto begin() { return subCells_.begin(); }
    auto begin() const { return subCells_.begin(); }

    auto end() { return subCells_.end(); }
    auto end() const { return subCells_.end(); }

    const SubCell& operator[](std::size_t i) const
    {
        assert(i < subCells_.size());
        return subCells_[i];
    }

private:
    std::vector<SubCell> subCells_;
};

/*!
 * \ingroup CCMpfaDiscretization
 * \todo doc me!
 */
class InteractionRegions
{
public:
    template<typename GridView,
             typename ElementMapper,
             typename VertexMapper>
    InteractionRegions(const GridView& gridView,
                       const ElementMapper& elementMapper,
                       const VertexMapper& vertexMapper)
    { construct_(gridView, elementMapper, vertexMapper); }

    std::size_t size() const
    { return interactionRegions_.size(); }

    const InteractionRegion& operator[](std::size_t vertexIndex) const
    { return interactionRegions_[vertexIndex]; }

private:
    template<typename GridView,
             typename ElementMapper,
             typename VertexMapper>
    void construct_(const GridView& gridView,
                    const ElementMapper& elementMapper,
                    const VertexMapper& vertexMapper)
    {
        static constexpr int dim = GridView::dimension;
        interactionRegions_.resize(gridView.size(dim));
        for (const auto& element : elements(gridView))
            for (int corner = 0; corner < element.subEntities(dim); ++corner)
                interactionRegions_[vertexMapper.subIndex(element, corner, dim)].add(
                    SubCell{
                        elementMapper.index(element),
                        static_cast<std::uint_least8_t>(corner)
                    }
                );
    }

    std::vector<InteractionRegion> interactionRegions_;
};

} // end namespace Dumux::CCMpfaO

#endif
