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

#include <dune/common/parallel/mpihelper.hh>
#include <dune/common/fvector.hh>
#include <dune/common/timer.hh>
#include <dune/grid/io/file/vtk.hh>
#include <dune/grid/yaspgrid.hh>

#include <dumux/common/parameters.hh>
#include <dumux/io/grid/gridmanager.hh>
#include <dumux/io/grid/subgridmanager.hh>
#include <dumux/discretization/method.hh>

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
        return std::sqrt((x - center_[0])*(x - center_[0]) + (y - center_[1])*(y - center_[1])) < radius;
    }
private:
    const GlobalPosition center_;
};

int main(int argc, char** argv) try
{
    using namespace Dumux;

    // Initialize MPI, finalize is done automatically on exit.
    Dune::MPIHelper::instance(argc, argv);

    // First read parameters from input file.
    Dumux::Parameters::init(argc, argv);

    constexpr int dim = 2;
    using GlobalPosition = Dune::FieldVector<double, dim>;

    Dune::Timer timer;
    using HostGrid = Dune::YaspGrid<dim, Dune::TensorProductCoordinates<double, dim> >;

    using HostGridManager = Dumux::GridManager<HostGrid>;
    HostGridManager hostGridManager;
    hostGridManager.init();
    auto& hostGrid = hostGridManager.grid();

    // Calculate the bounding box of the host grid view.
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

    // Select all elements right of the center.
    auto elementSelectorOne = [&center](const auto& element)
    {
        return element.geometry().center()[0] > center[0];
    };

    // Select all elements left of the center.
    auto elementSelectorTwo = [&center](const auto& element)
    {
        return element.geometry().center()[0] < center[0];
    };

    // Select all elements within a circle around the center.
    // Instead of a lambda, we use a class providing an () operator.
    // Of course, a lambda would be possible here, too.
    CircleSelector<GlobalPosition> elementSelectorThree(center);

    // Create three different subgrids from the same hostgrid.
    auto subgridPtrOne = SubgridManager<HostGrid>::makeGrid(hostGrid, elementSelectorOne);
    auto subgridPtrTwo = SubgridManager<HostGrid>::makeGrid(hostGrid, elementSelectorTwo);
    auto subgridPtrThree = SubgridManager<HostGrid>::makeGrid(hostGrid, elementSelectorThree);

    std::cout << "Constructing a host grid and three subgrids took "  << timer.elapsed() << " seconds.\n";

    // Write out the host grid and the subgrids.
    {
        Dune::VTKWriter<HostGrid::LeafGridView> vtkWriter(hostGrid.leafGridView());
        vtkWriter.write("hostgrid");
    }
    {
        Dune::VTKWriter<SubgridManager<HostGrid>::Grid::LeafGridView> vtkWriter(subgridPtrOne->leafGridView());
        vtkWriter.write("subgrid_one");
    }
    {
        Dune::VTKWriter<SubgridManager<HostGrid>::Grid::LeafGridView> vtkWriter(subgridPtrTwo->leafGridView());
        vtkWriter.write("subgrid_two");
    }
    {
        Dune::VTKWriter<SubgridManager<HostGrid>::Grid::LeafGridView> vtkWriter(subgridPtrThree->leafGridView());
        vtkWriter.write("subgrid_three");
    }

    return 0;
}
///////////////////////////////////////
//////// Error handler ////////////////
///////////////////////////////////////
catch (Dumux::ParameterException &e) {
    std::cerr << e << ". Abort!\n";
    return 1;
}
catch (Dune::Exception &e) {
    std::cerr << "Dune reported error: " << e << std::endl;
    return 3;
}
catch (...) {
    std::cerr << "Unknown exception thrown!\n";
    return 4;
}
