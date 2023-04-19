//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
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

#include <dumux/io/grid/gridmanager_base.hh>
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
        auto fvGeometry = localView(*gridGeometry);
        for(const auto& element : elements(leafGridView))
        {
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
