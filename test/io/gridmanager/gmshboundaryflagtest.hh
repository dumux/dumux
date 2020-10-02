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
 *
 * \brief Test whether GmshBoundaryFlag works as expected
 * \note GmshBoundaryFlag only exists for ALUGrid
 *
 * This tests whether the boundary marker IDs accessed using boundary flags are
 * correct when using ALUGrid, a Gmsh mesh file, and the class GmshBoundaryFlag.
 */
#ifndef DUMUX_TEST_IO_GMSHBOUNDARYFLAG_TEST_HH
#define DUMUX_TEST_IO_GMSHBOUNDARYFLAG_TEST_HH

#include <dumux/io/grid/gridmanager.hh>
#include <dumux/discretization/method.hh>

#include <iostream>

namespace Dumux {

template<class Grid>
class GmshBoundaryFlagTest
{
    using GridView = typename Grid::LeafGridView;
    using Scalar = double;
    static const int dim = Grid::dimension;
    using GridManager = typename Dumux::GridManager<Grid>;

public:

    template<class GridGeometry>
    static void testGmshBoundaryFlag(const GridView& leafGridView,
                                     std::shared_ptr<const GridGeometry> gridGeometry,
                                     std::shared_ptr<const GridData<Grid>> gridData)
    {

        for(const auto& element : elements(leafGridView))
        {
            auto fvGeometry = localView(*gridGeometry);
            fvGeometry.bind(element);

            for (auto&& scvf : scvfs(fvGeometry))
            {
                if (scvf.boundary())
                {
                    const auto boundaryMarkerId = gridData->getBoundaryDomainMarker(scvf.boundaryFlag());
                    const auto& pos = scvf.center();
                    std::cout << "z-coordinate: " << pos[dim-1] << ", actual ID = " << boundaryMarkerId << ", ";

                    /* According to unitcube.geo:
                     * top = 1, bottom = 2, side = 3
                     *
                     * According to unitcube.msh:
                     * top = 1, bottom = 0, side = 2
                     */
                    const int topId = 1;
                    const int bottomId = 0;
                    const int sideId = 2;

                    const bool isTop = pos[dim-1] > 1.0 - eps_;
                    const bool isBottom = pos[dim-1] < eps_;
                    const bool isSide = !isTop && !isBottom;

                    if (isTop)
                    {
                        std::cout << "correct ID = " << topId << " (is top surface)" << std::endl;
                        if (boundaryMarkerId != topId)
                            DUNE_THROW(Dune::Exception, "BoundaryMarkerId for top is wrong!");
                    }
                    else if (isBottom)
                    {
                        std::cout << "correct ID = " << bottomId << " (is bottom surface)" << std::endl;
                        if (boundaryMarkerId != bottomId)
                            DUNE_THROW(Dune::Exception, "BoundaryMarkerId for bottom is wrong!");
                    }
                    else if (isSide)
                    {
                        std::cout << "correct ID = " << sideId << " (is side surface)" << std::endl;
                        if (boundaryMarkerId != sideId)
                            DUNE_THROW(Dune::Exception, "BoundaryMarkerId for side is wrong!");
                    }
                } // end if boundary
            } // end scvf loop
        } // end element loop
    } // end testGmshBoundaryFlag

private:
    static constexpr Scalar eps_ = 1e-4;

};


} // end namespace Dumux


#endif
