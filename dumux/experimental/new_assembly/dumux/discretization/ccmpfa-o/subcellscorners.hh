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
 * \todo Doc me!
 */
#ifndef DUMUX_DISCRETIZATION_CCMPFA_O_SUB_CELLS_CORNERS_HH
#define DUMUX_DISCRETIZATION_CCMPFA_O_SUB_CELLS_CORNERS_HH

#include <array>
#include <cstdint>
#include <cassert>

#include <dune/geometry/type.hh>

namespace Dumux::CCMpfaO {

struct SubEntity
{
    std::uint_least8_t codim;
    std::uint_least8_t index;
};

template<Dune::GeometryType::Id Id>
struct SubCellsLocalCorners;

template<>
struct SubCellsLocalCorners<Dune::GeometryTypes::quadrilateral>
{
    static constexpr int numSubCells = 4;
    static constexpr int numSubCellFaces = 2;
    static constexpr int numSubFacesPerSubCell = 2;
    static constexpr int numCornersPerSubFace = 2;

    using FacetIndices = std::array<std::uint_least8_t, numSubCellFaces>;
    using SubFaceCorners = std::array<SubEntity, numCornersPerSubFace>;
    using SubCellSubFaceCorners = std::array<SubFaceCorners, numSubFacesPerSubCell>;

    static constexpr const FacetIndices& facetIndices(unsigned int subCellIndex)
    {
        assert(subCellIndex < numSubCells);
        return facetIndexMap_[subCellIndex];
    }

    static constexpr const SubCellSubFaceCorners& subFaceCorners(unsigned int subCellIndex)
    {
        assert(subCellIndex < numSubCells);
        return subFaceCorners_[subCellIndex];
    }

private:
    static constexpr std::array<FacetIndices, numSubCells> facetIndexMap_{{
        {{0, 2}},
        {{2, 1}},
        {{3, 0}},
        {{1, 3}}
    }};

    static constexpr std::array<SubCellSubFaceCorners, numSubCells> subFaceCorners_{{
        {{
            {{{2, 0}, {1, 0}}}, // sub-cell 0 -> sub-face 0 (cell facet 0)
            {{{1, 2}, {2, 0}}}, // sub-cell 0 -> sub-face 1 (cell facet 2)
        }},
        {{
            {{{2, 1}, {1, 2}}}, // sub-cell 1 -> sub-face 0 (cell facet 2)
            {{{1, 1}, {2, 1}}}, // sub-cell 1 -> sub-face 1 (cell facet 1)
        }},
        {{
            {{{2, 2}, {1, 3}}}, // sub-cell 2 -> sub-face 0 (cell facet 3)
            {{{1, 0}, {2, 2}}}, // sub-cell 2 -> sub-face 1 (cell facet 0)
        }},
        {{
            {{{2, 3}, {1, 1}}}, // sub-cell 3 -> sub-face 0 (cell facet 1)
            {{{1, 3}, {2, 3}}}, // sub-cell 3 -> sub-face 1 (cell facet 3)
        }}
    }};
};

} // end namespace Dumux::CCMpfaO

#endif
