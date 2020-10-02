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
 * \brief Test for the projection of discrete solutions from
 *        a two-dimensional discrete solution to a one-dimensional.
 */
#include <config.h>

#include <iostream>
#include <fstream>
#include <cmath>

#include <dune/common/fvector.hh>
#include <dune/istl/bvector.hh>

#include <dune/grid/utility/structuredgridfactory.hh>
#include <dune/grid/io/file/vtk/vtkwriter.hh>
#include <dune/grid/yaspgrid.hh>
#include <dune/foamgrid/foamgrid.hh>

#include <dune/functions/functionspacebases/lagrangebasis.hh>
#include <dune/functions/functionspacebases/interpolate.hh>
#include <dune/functions/gridfunctions/discreteglobalbasisfunction.hh>

#include <dumux/io/container.hh>
#include <dumux/common/parameters.hh>
#include <dumux/common/integrate.hh>

#include <dumux/discretization/box/fvgridgeometry.hh>
#include <dumux/discretization/fem/fegridgeometry.hh>
#include <dumux/discretization/cellcentered/tpfa/fvgridgeometry.hh>
#include <dumux/discretization/projection/projector.hh>
#include <dumux/multidomain/glue.hh>

int main (int argc, char *argv[])
{
    // maybe initialize mpi
    Dune::MPIHelper::instance(argc, argv);

    // initialize parameter tree
    Dumux::Parameters::init(argc, argv);

    using Grid1 = Dune::YaspGrid<2>;
    using Grid2 = Dune::FoamGrid<1, 2>;

    using GridView1 = typename Grid1::LeafGridView;
    using GridView2 = typename Grid2::LeafGridView;
    using FEBasis1 = Dune::Functions::LagrangeBasis<GridView1, 2>;
    using FEBasis2 = Dune::Functions::LagrangeBasis<GridView2, 2>;

    using Scalar = double;
    using FEGridGeometry1 = Dumux::FEGridGeometry<FEBasis1>;
    using FEGridGeometry2 = Dumux::FEGridGeometry<FEBasis2>;

    // make the 2d grid
    using GlobalPosition = Dune::FieldVector<Scalar, 2>;
    GlobalPosition lower1(0.0);
    GlobalPosition upper1(1.0);
    std::array<unsigned int, 2> els1{{10, 10}};
    std::shared_ptr<Grid1> grid1 = Dune::StructuredGridFactory<Grid1>::createCubeGrid(lower1, upper1, els1);

    // make the 1d grid
    GlobalPosition lower2({0., 0.5});
    GlobalPosition upper2({1., 0.5});
    std::array<unsigned int, 1> els2{{13}};
    std::shared_ptr<Grid2> grid2 = Dune::StructuredGridFactory<Grid2>::createSimplexGrid(lower2, upper2, els2);

    // create grid geometries
    FEGridGeometry1 feGridGeometry1(std::make_shared<FEBasis1>(grid1->leafGridView()));
    FEGridGeometry2 feGridGeometry2(std::make_shared<FEBasis2>(grid2->leafGridView()));

    // create solution vectors
    using BlockType = Dune::FieldVector<Scalar, 1>;
    using SolutionVector = Dune::BlockVector< BlockType >;
    SolutionVector sol1;
    SolutionVector sol2;

    // function to distribute on grids
    auto f = [] (const auto& pos) -> BlockType { return {pos.two_norm2()}; };

    Dune::Functions::interpolate(feGridGeometry1.feBasis(), sol1, f);
    Dune::Functions::interpolate(feGridGeometry2.feBasis(), sol2, f);

    // make projections and check them
    const auto projectors = Dumux::makeProjectorPair( feGridGeometry1.feBasis(),
                                                      feGridGeometry2.feBasis(),
                                                      makeGlue(feGridGeometry1, feGridGeometry2) );
    const auto& forwardProjector = projectors.first;
    const auto& backwardProjector = projectors.second;

    auto params = forwardProjector.defaultParams();
    params.residualReduction = 1e-16;

    const auto p12 = forwardProjector.project(sol1, params);
    const auto p21 = backwardProjector.project(sol2, params);

    using namespace Dune::Functions;

    auto af1 = makeAnalyticGridViewFunction(f, feGridGeometry1.gridView());
    auto af2 = makeAnalyticGridViewFunction(f, feGridGeometry2.gridView());

    auto gf1 = makeDiscreteGlobalBasisFunction<BlockType>(feGridGeometry1.feBasis(), p21);
    auto gf2 = makeDiscreteGlobalBasisFunction<BlockType>(feGridGeometry2.feBasis(), p12);

    if (Dumux::getParam<bool>("Vtk.EnableOutput", false))
    {
        Dune::VTKWriter<typename FEGridGeometry1::GridView> writer1(feGridGeometry1.gridView());
        Dune::VTKWriter<typename FEGridGeometry2::GridView> writer2(feGridGeometry2.gridView());

        const auto fieldInfoAnalytic = Dune::VTK::FieldInfo({"u", Dune::VTK::FieldInfo::Type::scalar, 1});
        const auto fieldInfoProjected = Dune::VTK::FieldInfo({"u_p", Dune::VTK::FieldInfo::Type::scalar, 1});

        writer1.addVertexData(af1, fieldInfoAnalytic);
        writer1.addVertexData(gf1, fieldInfoProjected);

        writer2.addVertexData(af2, fieldInfoAnalytic);
        writer2.addVertexData(gf2, fieldInfoProjected);

        writer1.write("sol_2d");
        writer2.write("sol_1d");
    }

    // make sure that solution projected on 2d grid is "exact" on the intersecting edges
    for (const auto& element : elements(feGridGeometry1.gridView()))
    {
        auto lf1 = localFunction(gf1);
        lf1.bind(element);

        for (const auto& is : intersections(feGridGeometry1.gridView(), element))
        {
            const auto isGeom = is.geometry();
            const auto p = isGeom.center();

            if (p[1] > 0.5 - 1e-6 && p[1] < 0.5 + 1e-6)
            {
                const auto error = lf1(element.geometry().local(p)) - f(p);

                using std::sqrt;
                if ( sqrt(error*error) > 1e-14 )
                    DUNE_THROW(Dune::MathError, "Error norm on 2d grid too high (" << sqrt(error*error) << ")!");
            }
        }
    }

    // one-dimensional solution should still be "exact"
    const auto l2Norm = Dumux::integrateL2Error( feGridGeometry2.gridView(), af2, gf2, 3);
    if ( l2Norm > 1e-14 )
        DUNE_THROW(Dune::MathError, "Error norm on 1d grid too high!");

    return 0;
}
