//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \brief Test for the cake grid manager
 */
#include <config.h>
#include <iostream>
#include <cmath>
#include <string>

#include <dune/common/fvector.hh>
#include <dune/common/timer.hh>
#include <dune/grid/io/file/vtk.hh>
#include <dune/grid/yaspgrid.hh>
#if HAVE_DUNE_ALUGRID
#include <dune/alugrid/grid.hh>
#endif
#if HAVE_DUNE_SPGRID
#include <dune/grid/spgrid.hh>
#endif

#include <dumux/io/grid/gridmanager_yasp.hh>
#include <dumux/io/grid/gridmanager_alu.hh>
#include <dumux/io/grid/gridmanager_sub.hh>
#include <dumux/common/initialize.hh>
#include <dumux/common/parameters.hh>


/*!
 * \brief A method providing an () operator in order to select elements for a subgrid.
 */
template<class GlobalPosition>
class CircleSelector
{
public:
    CircleSelector(const GlobalPosition& center) : center_(center) {}

    //! Select all elements within a circle around a center point.
    template<class Element>
    int operator() (const Element& element) const
    {
        const auto x = element.geometry().center()[0];
        const auto y = element.geometry().center()[1];
        const double radius = 0.3;
        return std::hypot(x-center_[0], y-center_[1]) < radius;
    }
private:
    const GlobalPosition center_;
};

template<int dim, class HostGrid>
void testSubGrid(const std::string& hostGridName)
{
    using SubGrid = Dune::SubGrid<dim, HostGrid>;

    using HostGridManager = Dumux::GridManager<HostGrid>;
    HostGridManager externalHostGridManager;
    externalHostGridManager.init("External");
    auto& hostGrid = externalHostGridManager.grid();

    // Calculate the bounding box of the host grid view.
    using GlobalPosition = Dune::FieldVector<double, dim>;
    GlobalPosition bBoxMin(std::numeric_limits<double>::max());
    GlobalPosition bBoxMax(std::numeric_limits<double>::min());
    for (const auto& vertex : vertices(hostGrid.leafGridView()))
    {
        for (int i=0; i<dim; i++)
        {
            using std::min;
            using std::max;
            bBoxMin[i] = min(bBoxMin[i], vertex.geometry().corner(0)[i]);
            bBoxMax[i] = max(bBoxMax[i], vertex.geometry().corner(0)[i]);
        }
    }

    // Get the center of the hostgrid's domain.
    const GlobalPosition center{bBoxMin[0]+0.5*bBoxMax[0], bBoxMin[1]+0.5*bBoxMax[1]};

    // Write out the host grid and the subgrids.
    {
        Dune::VTKWriter<typename HostGrid::LeafGridView> vtkWriter(hostGrid.leafGridView());
        vtkWriter.write("hostgrid");
    }

    // Create different subgrids from the same hostgrid
    {
        std::cout << "Constructing SubGrid with host grid and lambda element selector" << std::endl;

        // Select all elements right of the center.
        auto elementSelector = [&center](const auto& element)
        {
            return element.geometry().center()[0] > center[0];
        };

        Dumux::GridManager<SubGrid> subgridManager;
        subgridManager.init(hostGrid, elementSelector);
        Dune::VTKWriter<typename SubGrid::LeafGridView> vtkWriter(subgridManager.grid().leafGridView());
        vtkWriter.write("subgrid_right");
    }

    // create without constructing host grid first
    {
        std::cout << "Constructing SubGrid from lambda without specifying host grid" << std::endl;

        // Select all elements left of the center.
        auto elementSelector = [&center](const auto& element)
        {
            return element.geometry().center()[0] < center[0];
        };

        Dumux::GridManager<Dune::SubGrid<2, HostGrid>> subgridManager;
        subgridManager.init(elementSelector, "Internal");
        Dune::VTKWriter<typename  SubGrid::LeafGridView> vtkWriter(subgridManager.grid().leafGridView());
        vtkWriter.write("subgrid_left");
    }

    {
        std::cout << "Constructing SubGrid with host grid and functor selector" << std::endl;

        // Select all elements within a circle around the center.
        // Instead of a lambda, we use a class providing an () operator.
        // Of course, a lambda would be possible here, too.
        CircleSelector<GlobalPosition> elementSelector(center);

        Dumux::GridManager<SubGrid> subgridManager;
        subgridManager.init(hostGrid, elementSelector);
        Dune::VTKWriter<typename SubGrid::LeafGridView> vtkWriter(subgridManager.grid().leafGridView());
        vtkWriter.write("subgrid_circle_" + hostGridName);
    }
}

int main(int argc, char** argv)
{
    using namespace Dumux;

    // maybe initialize MPI and/or multithreading backend
    initialize(argc, argv);

    // First read parameters from input file.
    Dumux::Parameters::init(argc, argv);

    constexpr int dim = 2;

    {
        Dune::Timer timer;
        using HostGrid = Dune::YaspGrid<dim, Dune::TensorProductCoordinates<double, dim> >;
        testSubGrid<dim, HostGrid>("yasp");
        std::cout << "Constructing a yasp host grid and three subgrids took "  << timer.elapsed() << " seconds.\n";
    }

    {
#if HAVE_DUNE_ALUGRID
        Dune::Timer timer;
        using HostGrid = Dune::ALUGrid<dim, dim, Dune::cube, Dune::nonconforming>;
        testSubGrid<dim, HostGrid>("alu");
        std::cout << "Constructing a alu host grid and three subgrids took "  << timer.elapsed() << " seconds.\n";
#else
        std::cout << "Skipped test with ALUGrid as host grid.\n";
#endif
    }

    {
        std::cout << "Constructing SubGrid from binary image" << std::endl;
        using HostGrid = Dune::YaspGrid<2, Dune::EquidistantOffsetCoordinates<double, 2>>;
        using GridManager = Dumux::GridManager<Dune::SubGrid<2, HostGrid>>;
        GridManager subgridManager; subgridManager.init("Image");
        Dune::VTKWriter<GridManager::Grid::LeafGridView> vtkWriter(subgridManager.grid().leafGridView());
        vtkWriter.write("subgrid_binary_image");
    }

    {
        std::cout << "Constructing SubGrid from binary image" << std::endl;
        using HostGrid = Dune::YaspGrid<dim, Dune::EquidistantOffsetCoordinates<double, dim>>;
        using GridManager = Dumux::GridManager<Dune::SubGrid<dim, HostGrid>>;
        GridManager subgridManager;
        subgridManager.init("RepeatedImage");
        Dune::VTKWriter<GridManager::Grid::LeafGridView> vtkWriter(subgridManager.grid().leafGridView());
        vtkWriter.write("repeatedsubgrid_binary_image");
    }

#if HAVE_DUNE_SPGRID
        Dune::Timer timer;
        using HostGrid = Dune::SPGrid<double, dim>;
        using HostGridManager = GridManager<HostGrid>;
        HostGridManager hostGridManager;
        std::string gridName = "SPGrid";
        hostGridManager.init(gridName);
        auto& hostGrid = hostGridManager.grid();

        using SubGridManager = Dumux::GridManager<Dune::SubGrid<dim, HostGrid>>;
        SubGridManager subGridManager;

        // The subgrid is created from a 3x3 periodic host grid.
        // Two cells are removed (one internal, one at a boundary) such that all possible intersection cases are seen.

        std::vector<bool> img = {1, 0, 0, 0, 1, 0, 0, 0, 0};

        auto gridSelector = [&hostGrid, &img](const auto& element)
        {
            auto eIdx = hostGrid.leafGridView().indexSet().index(element);
            return img[eIdx] == 0;
        };

        subGridManager.init(hostGrid, gridSelector, gridName);
        std::cout << "Constructing a sp host grid and one basic image grid in "  << timer.elapsed() << " seconds.\n";
        Dune::VTKWriter<SubGridManager::Grid::LeafGridView> vtkWriter(subGridManager.grid().leafGridView());
        vtkWriter.write("periodic_subgrid_binary_image");

        // The intersections for the periodic subgrid should match the following properties.
        // Each entry contains a pair with bool values that should match (it.Boundary(), it.Neighbor())
        std::vector<std::pair<bool, bool>> intersectionProperties {{true, false}, {false, true}, {true, true}, {true, false},
                                                                   {false, true}, {true, false}, {true, true}, {false, true},
                                                                   {true, false}, {true, false}, {true, false}, {false, true},
                                                                   {true, false}, {true, false}, {false, true}, {false, true},
                                                                   {true, false}, {false, true}, {false, true}, {true, false},
                                                                   {false, true}, {false, true}, {true, false}, {true, true},
                                                                   {false, true}, {true, false}, {false, true}, {true, true}};
        int countSub = 0;
        int intersectionIdx = 0;
        for (const auto& element : elements(subGridManager.grid().leafGridView()))
        {
            for (auto& intersection : intersections(subGridManager.grid().leafGridView(), element))
            {
                if ((intersection.boundary() != intersectionProperties[intersectionIdx].first) ||
                    (intersection.neighbor() != intersectionProperties[intersectionIdx].second) )
                    DUNE_THROW(Dune::Exception, "The Subgrid intersection's (" << intersectionIdx << ") properties are incorrect. \n" <<
                                                "Intersection Center is : " << intersection.geometry().center() <<
                                                "\n Boundary is " <<  intersection.boundary() << " and should be " <<
                                                intersectionProperties[intersectionIdx].first <<
                                                "\n Neighbor is " <<  intersection.neighbor() << " and should be " <<
                                                intersectionProperties[intersectionIdx].second );
                if (intersection.boundary())
                    if (intersection.neighbor())
                        countSub++;
            intersectionIdx++;
            }
        }
        if (countSub == 0)
            DUNE_THROW(Dune::Exception, "subGrid is not periodic!");

#else
        std::cout << "Skipped test with SPGRID as host grid.\n";
#endif
    }

    return 0;
}
