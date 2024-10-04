//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
#include <config.h>

#include <dune/common/float_cmp.hh>
#include <dune/common/fvector.hh>
#include <dune/common/exceptions.hh>

#include <dumux/common/parameters.hh>
#include <dumux/common/initialize.hh>
#include <dumux/discretization/basegridgeometry.hh>
#include <dumux/io/grid/gridmanager_yasp.hh>
#include <dumux/io/grid/gridmanager_sub.hh>

#include <dumux/geometry/boundingboxtree.hh>
#include <dumux/geometry/geometricentityset.hh>
#include <dumux/geometry/intersectionentityset.hh>
#include <dumux/common/quadrature.hh>

int main (int argc, char *argv[])
{
    // maybe initialize MPI and/or multithreading backend
    Dumux::initialize(argc, argv);

    // initialize parameters
    Dumux::Parameters::init([](auto& p){
        p["Grid.UpperRight"] = "1 1";
        p["Grid.Cells"] = "10 10";
    });

    // initialize a grid
    static constexpr int dim = 2;
    using HostGrid = Dune::YaspGrid<dim>;
    using Grid = Dune::SubGrid<dim, HostGrid>;

    // Cut rectangle with size [0.2,0.8] x [0.2,0.8] out of domain
    auto elementSelector = [](const auto& element)
    {
        auto center = element.geometry().center();
        center[0] -= 0.5; center[1] -= 0.5;
        return center.infinity_norm() > 0.3;
    };

    Dumux::GridManager<Grid> gridManager;
    gridManager.init(elementSelector);
    const auto leafGridView = gridManager.grid().leafGridView();

    using Scalar = typename Grid::ctype;
    using GlobalPosition = Dune::FieldVector<Scalar, 2>;

    // a linear function
    auto f = [&] (const GlobalPosition& pos) { return 0.85*pos[0] + 2.3*pos[1]; };

    // perform volume averaging on a coarser grid
    using CoarseGrid = Dune::YaspGrid<dim>;
    const auto coarseGrid = Dune::StructuredGridFactory<CoarseGrid>::createCubeGrid(
        {0.0, 0.0}, {1.0, 1.0}, {6, 6}
    );

    Scalar totalVolumeFraction(0.0);

    using GridEntitySet = Dumux::GridViewGeometricEntitySet<typename Grid::LeafGridView, 0>;
    using GeometryEntitySet = Dumux::SingleGeometryEntitySet<typename CoarseGrid::template Codim<0>::Entity::Geometry>;
    Dumux::IntersectionEntitySet<GridEntitySet, GeometryEntitySet> iset;
    const auto gridEntitySet = std::make_shared<GridEntitySet>(leafGridView);

    for(const auto& coarseElement : elements(coarseGrid->leafGridView()))
    {
        Scalar fraction = 0.0;
        Scalar f_avg = 0.0;

        iset.build(gridEntitySet, std::make_unique<GeometryEntitySet>(coarseElement.geometry()));
        for (const auto& is : intersections(iset))
        {
            const auto geometry = is.geometry();
            const auto volume = geometry.volume();
            const auto center = geometry.center();

            fraction += volume;
            f_avg += volume*f(center);
        }

        totalVolumeFraction += fraction;
        Scalar coarseElementVolume = coarseElement.geometry().volume();
        fraction /= coarseElementVolume; f_avg /= coarseElementVolume;

        // If the whole coarse element is embedded within the fine grid
        // the average is exact and given by evalution at the centroid of the averaging volume
        if(Dune::FloatCmp::eq(fraction, 1.0))
        {
            Scalar exactAverage = f(coarseElement.geometry().center());
            if (Dune::FloatCmp::ne(f_avg, exactAverage))
                DUNE_THROW(Dune::InvalidStateException, "Unexpected volume average" << f_avg << " vs. " << exactAverage);
        }
    }

    // Expected fraction is domain area minus the interior rectangle (porosity)
    Scalar expectedVolumeFraction = 1.0-0.6*0.6;
    if (Dune::FloatCmp::ne(totalVolumeFraction, expectedVolumeFraction))
        DUNE_THROW(Dune::InvalidStateException, "Unexpected volume fraction" << totalVolumeFraction << " vs. " << expectedVolumeFraction);

    std::cout << "All tests passed!" << std::endl;
    return 0;
}
