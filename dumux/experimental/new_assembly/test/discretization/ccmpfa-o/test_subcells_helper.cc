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
 * \brief Tests for the `CCMpfaO` subcells helper functions.
 */
#include <iostream>
#include <cmath>
#include <array>

#include <dune/common/exceptions.hh>
#include <dune/common/fvector.hh>

#include <dune/geometry/type.hh>
#include <dune/geometry/multilineargeometry.hh>

#include <dumux/common/math.hh>
#include <dumux/experimental/new_assembly/dumux/discretization/ccmpfa-o/subcellshelper.hh>


template<typename Position, std::size_t N>
bool isRightHandSystem(const std::array<Position, N>& basis)
{
    static_assert(N == 2);

    std::cout << "Checking basis:" << std::endl;
    std::ranges::for_each(basis, [] (const auto& b) { std::cout << b << "\n"; });
    std::cout << std::endl;

    if constexpr (N == 2)
        return !std::signbit(Dumux::crossProduct(basis[0], basis[1]));
}

template<typename Position, std::size_t N>
auto computeSubFaceArea(const std::array<Position, N>& corners)
{
    static_assert(N == 2);
    if constexpr (N == 2)
        return (corners[1] - corners[0]).two_norm();
}


void testQuadrilateral()
{
    Dune::MultiLinearGeometry<double, 2, 2> geo{
            Dune::GeometryTypes::quadrilateral,
            std::vector<Dune::FieldVector<double, 2>>{{
                {0.0, 0.0}, {1.0, 0.0},
                {0.0, 1.0}, {1.0, 1.0}
            }}
        };

    Dumux::CCMpfaO::SubCellsHelper subCellsHelper{geo};
    for (int i = 0; i < subCellsHelper.numSubCells(); ++i)
    {
        const auto numSubCellFaces = std::ranges::size(subCellsHelper.facetIndices(i));
        if (numSubCellFaces != subCellsHelper.numSubFacesPerSubCell)
            DUNE_THROW(Dune::InvalidStateException, "Unexpected number of faces");

        const auto basis = subCellsHelper.makeBasis(i);
        const auto area = Dumux::crossProduct(basis[0], basis[1]);
        if (!isRightHandSystem(basis))
            DUNE_THROW(Dune::InvalidStateException, "Local basis not a right hand system");
        if (std::abs(area - 0.5*0.5) > 1e-7)
            DUNE_THROW(Dune::InvalidStateException, "Unexpected sub cell area");

        for (int fIdx = 0; fIdx < subCellsHelper.numSubFacesPerSubCell; ++fIdx)
        {
            const auto& subFaceCornerIndices = subCellsHelper.subFaceCornerIndices(i, fIdx);
            if (std::ranges::size(subFaceCornerIndices) != subCellsHelper.numCornersPerSubFace)
                DUNE_THROW(Dune::InvalidStateException, "Unexpected number of corners");

            const auto& subFaceCorners = subCellsHelper.subFaceCorners(i, fIdx);
            const auto subFaceArea = computeSubFaceArea(subFaceCorners);
            if (std::abs(subFaceArea - 0.5) > 1e-7)
                DUNE_THROW(Dune::InvalidStateException, "Unexpected sub-face area");
        }
    }
}

int main (int argc, char *argv[])
{
    testQuadrilateral();

    std::cout << "\nAll tests passed" << std::endl;
    return 0;
}
