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
 * \brief Test for the cake grid manager
 */
#include <config.h>
#include <iostream>
#include <cmath>
#include <string>

#include <dune/common/parallel/mpihelper.hh>
#include <dune/common/fvector.hh>
#include <dune/common/timer.hh>
#include <dune/grid/io/file/vtk.hh>
#include <dune/grid/yaspgrid.hh>
#if HAVE_DUNE_ALUGRID
#include <dune/alugrid/grid.hh>
#endif

#include <dumux/io/grid/gridmanager_sub.hh>
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

    // create without contructing host grid first
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

    // Initialize MPI, finalize is done automatically on exit.
    Dune::MPIHelper::instance(argc, argv);

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

    return 0;
}
