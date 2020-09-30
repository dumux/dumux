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
 * \brief Test for the grid finite volume element geometry for non-conforming grids.
 *        A square grid with 9 elements is created of which the central element is refined.
 *        Subsequently, all directions & connectivity of the scvfs are checked for correctness.
 */
#include <config.h>

#include <cmath>
#include <iostream>
#include <utility>

#include <dune/common/float_cmp.hh>

#include <dune/grid/utility/structuredgridfactory.hh>
#include <dune/alugrid/grid.hh>

#include <dumux/adaptive/markelements.hh>
#include <dumux/discretization/cellcentered/tpfa/fvgridgeometry.hh>

#ifndef DOXGEN
namespace Dumux {
namespace Test {

//! epsilon for checking direction of scvf normals
constexpr double eps = 1e-6;

//! returns true if a given position is inside the central element (before refinement)
bool isInCentralElement(const Dune::FieldVector<double, 2>& pos)
{
    return pos[0] < 2.0 && pos[0] > 1.0 && pos[1] > 1.0 && pos[1] < 2.0;
}

//! returns if the element with given center position is corner element
bool isCornerElement(const Dune::FieldVector<double, 2>& center)
{
    const auto distVec = center - Dune::FieldVector<double, 2>({1.5, 1.5});
    return (std::abs(distVec[0]) > 1.0 - eps && std::abs(distVec[1]) > 1.0 - eps);
}

//! returns if the element is middle left element
bool isMiddleLeftElement(const Dune::FieldVector<double, 2>& center)
{
    const auto distVec = center - Dune::FieldVector<double, 2>({1.5, 1.5});
    return distVec[0] < -1.0 + eps && std::abs(distVec[1]) < eps;
}

//! returns if the element is middle right element
bool isMiddleRightElement(const Dune::FieldVector<double, 2>& center)
{
    const auto distVec = center - Dune::FieldVector<double, 2>({1.5, 1.5});
    return distVec[0] > 1.0 - eps && std::abs(distVec[1]) < eps;
}

//! returns if the element is middle upper element
bool isMiddleUpperElement(const Dune::FieldVector<double, 2>& center)
{
    const auto distVec = center - Dune::FieldVector<double, 2>({1.5, 1.5});
    return distVec[1] > 1.0 - eps && std::abs(distVec[0]) < eps;
}

//! returns if the element is middle lower element
bool isMiddleLowerElement(const Dune::FieldVector<double, 2>& center)
{
    const auto distVec = center - Dune::FieldVector<double, 2>({1.5, 1.5});
    return distVec[1] < -1.0 + eps && std::abs(distVec[0]) < eps;
}

//! returns for a given element center the element type name used in this test
std::string elementTypeName(const Dune::FieldVector<double, 2>& center)
{
    if (isInCentralElement(center))
        return "Central element";
    if (isCornerElement(center))
        return "Corner element";
    if (isMiddleLeftElement(center))
        return "Middle left";
    if (isMiddleRightElement(center))
        return "Middle right";
    if (isMiddleUpperElement(center))
        return "Middle upper";
    if (isMiddleLowerElement(center))
        return "Middle lower";

    DUNE_THROW(Dune::InvalidStateException, "Element center position could not be interpreted.");
}
} // end namespace Test
} // end namespace Dumux
#endif

int main (int argc, char *argv[])
{
    using namespace Dumux;
    using namespace Dumux::Test;

    // maybe initialize mpi
    Dune::MPIHelper::instance(argc, argv);

    std::cout << "Checking the FVGeometries, SCVs and SCV faces on a non-conforming grid" << std::endl;

    using Grid = Dune::ALUGrid<2, 2, Dune::cube, Dune::nonconforming>;

    constexpr int dim = Grid::dimension;

    using GridGeometry = CCTpfaFVGridGeometry<typename Grid::LeafGridView, ENABLE_CACHING>;
    using Element = typename GridGeometry::GridView::template Codim<0>::Entity;
    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;

    //! make a grid
    GlobalPosition lower(0.0);
    GlobalPosition upper(3.0);
    std::array<unsigned int, dim> els{{3, 3}};
    std::shared_ptr<Grid> grid = Dune::StructuredGridFactory<Grid>::createCubeGrid(lower, upper, els);

    //! refine the central element once
    auto leafGridView = grid->leafGridView();
    markElements(*grid, [&](const auto& e){ return isInCentralElement(e.geometry().center()) ? 1 : 0; });
    grid->preAdapt();
    grid->adapt();
    grid->postAdapt();

    //! if the leaf now does not have 12 elements, something went wrong
    if (leafGridView.size(0) != 12)
        DUNE_THROW(Dune::InvalidStateException, "Refined grid does not have exactly 12 elements!");

    //! instantiate and update gridGeometry
    GridGeometry gridGeometry(leafGridView);
    gridGeometry.update();

    //! We should have constructed 12 scvfs
    if (gridGeometry.numScv() != 12)
        DUNE_THROW(Dune::InvalidStateException, "FvGridGeometry does not have exactly 12 scvs!");

    //! We should have constructed 52 scvfs
    if (gridGeometry.numScvf() != 52)
        DUNE_THROW(Dune::InvalidStateException, "FvGridGeometry does not have exactly 52 scvfs!");

    //! iterate over elements and check for each element the number of scvfs
    for (const auto& element : elements(leafGridView))
    {
        auto fvGeometry = localView(gridGeometry);
        fvGeometry.bind(element);

        //! For the tpfa scheme there is always one scv per element
        if (fvGeometry.numScv() != 1)
            DUNE_THROW(Dune::InvalidStateException, "More than one scv found in an element!");

        //! make sure the scv has the same center as the element
        const auto elementCenter = element.geometry().center();
        for (const auto& scv : scvs(fvGeometry))
        {
            const auto scvCenter = scv.center();
            for (unsigned int dir = 0; dir < dim; ++dir)
                if ( !Dune::FloatCmp::eq(elementCenter[dir], scvCenter[dir], eps) )
                    DUNE_THROW(Dune::InvalidStateException, "Element center - " << elementCenter << " - and scv center - " << scvCenter << " - do not coincide");
        }

        //! check if the right number of scvfs point in the right direction
        std::size_t numScvfsPosX = 0;
        std::size_t numScvfsNegX = 0;
        std::size_t numScvfsPosY = 0;
        std::size_t numScvfsNegY = 0;
        std::size_t numScvfsTotal = 0;

        //! Also, scheck how many neighbors on the two levels the element has
        //! and store the neighboring centers for output
        std::size_t numL0Neighbors = 0;
        std::size_t numL1Neighbors = 0;
        std::vector<GlobalPosition> neighborCenters;
        for (const auto& scvf : scvfs(fvGeometry))
        {
            const auto n = scvf.unitOuterNormal();

            //! the normal vector should point either in x- or y-direction with length 1
            using std::abs;
            if (std::abs(n[0]) > eps && std::abs(n[1]) > eps)
                DUNE_THROW(Dune::InvalidStateException, "Wrong unit outer normal vector");

            //! Outer normal must be pointing outwards
            if ( n*(scvf.center() - elementCenter) < 0 )
                DUNE_THROW(Dune::InvalidStateException, "Normal vector does not point outwards");

            //! the face should never have more the one neighbor
            if (scvf.numOutsideScvs() > 1)
                DUNE_THROW(Dune::InvalidStateException, "Scvf has more than one neighbor");

            //! center must always be between the corners
            const auto d1 = scvf.corner(0) - scvf.center();
            const auto d2 = scvf.corner(1) - scvf.center();
            if ( d1 * d2 >= 0 )
                DUNE_THROW(Dune::InvalidStateException, "Center is not between the two corners");

            //! count up faces depending on direction of n
            if ( Dune::FloatCmp::eq(n[0],  1.0, eps) )
                numScvfsPosX++;
            if ( Dune::FloatCmp::eq(n[0], -1.0, eps) )
                numScvfsNegX++;
            if ( Dune::FloatCmp::eq(n[1],  1.0, eps) )
                numScvfsPosY++;
            if ( Dune::FloatCmp::eq(n[1], -1.0, eps) )
                numScvfsNegY++;

            //! keep track of total number of scvfs
            numScvfsTotal++;

            //! check levels of neighbors
            if (!scvf.boundary())
            {
                const auto outsideElement = gridGeometry.element(scvf.outsideScvIdx());
                const auto outsideCenter = outsideElement.geometry().center();

                if (isInCentralElement(outsideCenter))
                    numL1Neighbors++;
                else
                    numL0Neighbors++;

                //! check if ipGlobal & area make sense
                const auto distVec = scvf.ipGlobal() - elementCenter;
                if ( isInCentralElement(outsideCenter) && !isInCentralElement(elementCenter) )
                {
                    if ( Dune::FloatCmp::eq(distVec[0], 0.0, eps) || Dune::FloatCmp::eq(distVec[1], 0.0, eps) )
                        DUNE_THROW(Dune::InvalidStateException, "Vector from element center to ipGlobal must NOT be axis-parallel " <<
                                                                "when going from lower to higher level");
                    if ( !Dune::FloatCmp::eq(scvf.area(), 0.5, eps) )
                        DUNE_THROW(Dune::InvalidStateException, "Area of scvf towards level 1 element is not 0.5!");
                }
                else
                {
                    if ( !Dune::FloatCmp::eq(distVec[0], 0.0, eps) && !Dune::FloatCmp::eq(distVec[1], 0.0, eps) )
                        DUNE_THROW(Dune::InvalidStateException, "Vector from element center to ipGlobal must be axis-parallel " <<
                                                                "when neighboring levels are identical or outside level is lower");

                    if ( !isInCentralElement(outsideCenter) && !isInCentralElement(elementCenter) )
                    {
                        if ( !Dune::FloatCmp::eq(scvf.area(), 1.0, eps) )
                            DUNE_THROW(Dune::InvalidStateException, "Area of scvf between level 0 element is not 1.0!");
                    }
                    else
                    {
                        if ( !Dune::FloatCmp::eq(scvf.area(), 0.5, eps) )
                            DUNE_THROW(Dune::InvalidStateException, "Area of scvf between different levels is not 0.5!");
                    }
                }

                //! store outside center for output
                neighborCenters.push_back(outsideCenter);
            }
        }

        //! make sure the number of found scvfs makes sense
        const std::size_t sumScvfs = numScvfsPosX+numScvfsNegX+numScvfsPosY+numScvfsNegY;
        if (numScvfsTotal != sumScvfs)
            DUNE_THROW(Dune::InvalidStateException, "Number of total scvfs is not equal to number of faces in individual directions.\n"
                                                        << "Total number: " << numScvfsTotal << ", sum individual faces: " << sumScvfs);

        //! print found combination in element to terminal
        std::cout << elementTypeName(elementCenter) << " with element center " << elementCenter << " has " << numScvfsTotal << " faces:\n"
                  << "    - " << numScvfsPosX << " in positive x-direction\n"
                  << "    - " << numScvfsNegX << " in negative x-direction\n"
                  << "    - " << numScvfsPosY << " in positive y-direction\n"
                  << "    - " << numScvfsNegY << " in negative x-direction\n"
                  << "The neighboring element centers are:" << std::endl;
        for (const auto& p : neighborCenters) std::cout << "\t\t - " << p << std::endl;

        //! in elements within the original central element or corner elements we should have only 4 scvfs (1 in each direction)
        using std::abs;
        if ( isInCentralElement(elementCenter) || isCornerElement(elementCenter) )
        {
            if (numScvfsTotal != 4)
                DUNE_THROW(Dune::InvalidStateException, "Did not find exactly 4 scvfs in a central/corner element");
            if (numScvfsPosX != 1)
                DUNE_THROW(Dune::InvalidStateException, "Did not find exactly 1 scvfs in a central/corner element in positive x-direction");
            if (numScvfsPosY != 1)
                DUNE_THROW(Dune::InvalidStateException, "Did not find exactly 1 scvfs in a central/corner element in positive y-direction");
            if (numScvfsNegX != 1)
                DUNE_THROW(Dune::InvalidStateException, "Did not find exactly 1 scvfs in a central/corner element in negative x-direction");
            if (numScvfsNegY != 1)
                DUNE_THROW(Dune::InvalidStateException, "Did not find exactly 1 scvf1 in a central/corner element in negative y-direction");

            //! In the corners we should find two level 0 neighbor elements
            if ( isCornerElement(elementCenter) && numL1Neighbors != 0 && numL0Neighbors != 2 )
                DUNE_THROW(Dune::InvalidStateException, "Did not find exactly two level 0 neighbors and zero level 1 neighbors in corner element");

            //! In the center elements we should find two level 0 and 2 level 1 neighbor elements
            if ( isInCentralElement(elementCenter) && numL1Neighbors != 2 && numL0Neighbors != 2 )
                DUNE_THROW(Dune::InvalidStateException, "Did not find exactly two level 0 neighbors and two level 1 neighbors in central element");
        }
        //! In the element on the left to the center we should have 2 scvfs in positive x-direction
        else if ( isMiddleLeftElement(elementCenter) ) //! to the left of the center
        {
            if (numScvfsTotal != 5)
                DUNE_THROW(Dune::InvalidStateException, "Did not find exactly 5 scvfs in middle left element");
            if (numScvfsPosX != 2)
                DUNE_THROW(Dune::InvalidStateException, "Did not find exactly 2 scvfs in middle left element in positive x-direction");
            if (numScvfsPosY != 1)
                DUNE_THROW(Dune::InvalidStateException, "Did not find exactly 1 scvfs in middle left element in positive y-direction");
            if (numScvfsNegX != 1)
                DUNE_THROW(Dune::InvalidStateException, "Did not find exactly 1 scvfs in middle left element in negative x-direction");
            if (numScvfsNegY != 1)
                DUNE_THROW(Dune::InvalidStateException, "Did not find exactly 1 scvfs in middle left element in negative y-direction");

            //! we should find two level 0 and 2 level 1 neighbor elements
            if (numL1Neighbors != 2 && numL0Neighbors != 2)
                DUNE_THROW(Dune::InvalidStateException, "Did not find exactly two level 0 neighbors and two level 1 neighbors in middle left element");
        }
        //! In the element on the right to the center we should have 2 scvfs in negative x-direction
        else if ( isMiddleRightElement(elementCenter) ) //! to the right of the center
        {
            if (numScvfsTotal != 5)
                DUNE_THROW(Dune::InvalidStateException, "Did not find exactly 5 scvfs in middle right element");
            if (numScvfsPosX != 1)
                DUNE_THROW(Dune::InvalidStateException, "Did not find exactly 1 scvfs in middle right element in positive x-direction");
            if (numScvfsPosY != 1)
                DUNE_THROW(Dune::InvalidStateException, "Did not find exactly 1 scvfs in middle right element in positive y-direction");
            if (numScvfsNegX != 2)
                DUNE_THROW(Dune::InvalidStateException, "Did not find exactly 2 scvfs in middle right element in negative x-direction");
            if (numScvfsNegY != 1)
                DUNE_THROW(Dune::InvalidStateException, "Did not find exactly 1 scvfs in middle right element in negative y-direction");

            //! we should find two level 0 and 2 level 1 neighbor elements
            if ( numL1Neighbors != 2 && numL0Neighbors != 2 )
                DUNE_THROW(Dune::InvalidStateException, "Did not find exactly two level 0 neighbors and two level 1 neighbors in middle right element");
        }
        //! In the element above the center we should have 2 scvfs in negative y-direction
        else if ( isMiddleUpperElement(elementCenter) ) //! above the center
        {
            if (numScvfsTotal != 5)
                DUNE_THROW(Dune::InvalidStateException, "Did not find exactly 5 scvfs in middle upper element");
            if (numScvfsPosX != 1)
                DUNE_THROW(Dune::InvalidStateException, "Did not find exactly 1 scvfs in middle upper element in positive x-direction");
            if (numScvfsPosY != 1)
                DUNE_THROW(Dune::InvalidStateException, "Did not find exactly 1 scvfs in middle upper element in positive y-direction");
            if (numScvfsNegX != 1)
                DUNE_THROW(Dune::InvalidStateException, "Did not find exactly 1 scvfs in middle upper element in negative x-direction");
            if (numScvfsNegY != 2)
                DUNE_THROW(Dune::InvalidStateException, "Did not find exactly 2 scvfs in middle upper element in negative y-direction");

            //! we should find two level 0 and 2 level 1 neighbor elements
            if ( numL1Neighbors != 2 && numL0Neighbors != 2 )
                DUNE_THROW(Dune::InvalidStateException, "Did not find exactly two level 0 neighbors and two level 1 neighbors in middle upper element");
        }
        //! In the element below the center we should have 2 scvfs in positive y-direction
        else if ( isMiddleLowerElement(elementCenter) ) //! below the center
        {
            if (numScvfsTotal != 5)
                DUNE_THROW(Dune::InvalidStateException, "Did not find exactly 5 scvfs in middle lower element");
            if (numScvfsPosX != 1)
                DUNE_THROW(Dune::InvalidStateException, "Did not find exactly 1 scvfs in middle lower element in positive x-direction");
            if (numScvfsPosY != 2)
                DUNE_THROW(Dune::InvalidStateException, "Did not find exactly 1 scvfs in middle lower element in positive y-direction");
            if (numScvfsNegX != 1)
                DUNE_THROW(Dune::InvalidStateException, "Did not find exactly 1 scvfs in middle lower element in negative x-direction");
            if (numScvfsNegY != 1)
                DUNE_THROW(Dune::InvalidStateException, "Did not find exactly 1 scvfs in middle lower element in negative y-direction");

            //! we should find two level 0 and 2 level 1 neighbor elements
            if ( numL1Neighbors != 2 && numL0Neighbors != 2 )
                DUNE_THROW(Dune::InvalidStateException, "Did not find exactly two level 0 neighbors and two level 2 neighbors in middle lower element");
        }
        else
            DUNE_THROW(Dune::InvalidStateException, "Element center position could not be interpreted.");
    }
}
