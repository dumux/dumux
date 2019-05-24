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

#include <dumux/discretization/box/fvgridgeometry.hh>
#include <dumux/discretization/fem/fegridgeometry.hh>
#include <dumux/discretization/cellcentered/tpfa/fvgridgeometry.hh>
#include <dumux/discretization/projection/projector.hh>
#include <dumux/multidomain/glue.hh>

template< class Scalar, class SolutionVector, class GridGeometry1, class GridGeometry2 >
void writeProjection(const GridGeometry1& gg1, const GridGeometry2& gg2,
                     const SolutionVector& sol1, const SolutionVector& sol2,
                     const std::string& acro1, const std::string& acro2)
{
    const auto projectors = Dumux::makeProjectorPair( getFunctionSpaceBasis(gg1),
                                                      getFunctionSpaceBasis(gg2),
                                                      makeGlue(gg1, gg2) );
    const auto& forwardProjector = projectors.first;
    const auto& backwardProjector = projectors.second;

    SolutionVector p2_1, p1_2;
    forwardProjector.project(sol1, p2_1);
    backwardProjector.project(sol2, p1_2);

    const std::string filename1 = "2d_" + acro1 + "_" + acro2;
    const std::string filename2 = "1d_" + acro1 + "_" + acro2;

    Dumux::writeContainerToFile(p1_2, filename1 + ".csv");
    Dumux::writeContainerToFile(p2_1, filename2 + ".csv");

    static const bool writeVtk = Dumux::getParam<bool>("Vtk.EnableOutput", false);
    if (writeVtk)
    {
        using namespace Dune;
        VTKWriter<typename GridGeometry1::GridView> writer1(gg1.gridView());
        VTKWriter<typename GridGeometry2::GridView> writer2(gg2.gridView());

        const auto& basis1 = Dumux::getFunctionSpaceBasis(gg1);
        const auto& basis2 = Dumux::getFunctionSpaceBasis(gg2);

        const auto uInfo = VTK::FieldInfo("u", VTK::FieldInfo::Type::scalar, 1);
        const auto upInfo = VTK::FieldInfo("u_p", VTK::FieldInfo::Type::scalar, 1);

        const auto uFunction1 = Functions::makeDiscreteGlobalBasisFunction<Scalar>(basis1, sol1);
        const auto uFunction2 = Functions::makeDiscreteGlobalBasisFunction<Scalar>(basis2, sol2);
        const auto upFunction1 = Functions::makeDiscreteGlobalBasisFunction<Scalar>(basis1, p1_2);
        const auto upFunction2 = Functions::makeDiscreteGlobalBasisFunction<Scalar>(basis2, p2_1);

        if (GridGeometry1::discMethod == Dumux::DiscretizationMethod::cctpfa)
        {
            writer1.addCellData(uFunction1, uInfo);
            writer1.addCellData(upFunction1, upInfo);
        }
        else
        {
            writer1.addVertexData(uFunction1, uInfo);
            writer1.addVertexData(upFunction1, upInfo);
        }

        if (GridGeometry2::discMethod == Dumux::DiscretizationMethod::cctpfa)
        {
            writer2.addCellData(uFunction2, uInfo);
            writer2.addCellData(upFunction2, upInfo);
        }
        else
        {
            writer2.addVertexData(uFunction2, uInfo);
            writer2.addVertexData(upFunction2, upInfo);
        }

        writer1.write(filename1);
        writer2.write(filename2);
    }
}

int main (int argc, char *argv[]) try
{
    // maybe initialize mpi
    Dune::MPIHelper::instance(argc, argv);

    // initialize parameter tree
    Dumux::Parameters::init(argc, argv);

    using Grid1 = Dune::YaspGrid<2>;
    using Grid2 = Dune::FoamGrid<1, 2>;

    using GridView1 = typename Grid1::LeafGridView;
    using GridView2 = typename Grid2::LeafGridView;
    using FEBasis1 = Dune::Functions::LagrangeBasis<GridView1, 1>;
    using FEBasis2 = Dune::Functions::LagrangeBasis<GridView2, 1>;

    using Scalar = double;
    using TpfaGridGeometry1 = Dumux::CCTpfaFVGridGeometry<GridView1>;
    using BoxGridGeometry1 = Dumux::BoxFVGridGeometry<Scalar, GridView1>;
    using FEGridGeometry1 = Dumux::FEGridGeometry<FEBasis1>;

    using TpfaGridGeometry2 = Dumux::CCTpfaFVGridGeometry<GridView2>;
    using BoxGridGeometry2 = Dumux::BoxFVGridGeometry<Scalar, GridView2>;
    using FEGridGeometry2 = Dumux::FEGridGeometry<FEBasis2>;

    // make the 2d grid
    using GlobalPosition = Dune::FieldVector<Scalar, 2>;
    GlobalPosition lower1(0.0);
    GlobalPosition upper1(1.0);
    std::array<unsigned int, 2> els1{{10, 10}};
    std::shared_ptr<Grid1> grid1 = Dune::StructuredGridFactory<Grid1>::createCubeGrid(lower1, upper1, els1);

    // make the 1d grid
    GlobalPosition lower2({0., 0.5});
    GlobalPosition upper2({1. ,0.5});
    std::array<unsigned int, 1> els2{{13}};
    std::shared_ptr<Grid2> grid2 = Dune::StructuredGridFactory<Grid2>::createSimplexGrid(lower2, upper2, els2);

    Dune::Functions::LagrangeBasis<GridView1, 0> q0Basis1(grid1->leafGridView());
    Dune::Functions::LagrangeBasis<GridView2, 0> q0Basis2(grid2->leafGridView());
    Dune::Functions::LagrangeBasis<GridView1, 1> q1Basis1(grid1->leafGridView());
    Dune::Functions::LagrangeBasis<GridView2, 1> q1Basis2(grid2->leafGridView());

    // create all grid geometries
    TpfaGridGeometry1 tpfaGridGeometry1(grid1->leafGridView());
    BoxGridGeometry1 boxGridGeometry1(grid1->leafGridView());
    FEGridGeometry1 feGridGeometry1(std::make_shared<FEBasis1>(grid1->leafGridView()));

    TpfaGridGeometry2 tpfaGridGeometry2(grid2->leafGridView());
    BoxGridGeometry2 boxGridGeometry2(grid2->leafGridView());
    FEGridGeometry2 feGridGeometry2(std::make_shared<FEBasis2>(grid2->leafGridView()));

    // function to distribute on grids
    auto evalFunction1 = [] (const auto& pos)
    { return std::sin(2*M_PI*pos[0])*std::cos(4*M_PI*pos[1]); };
    auto evalFunction2 = [] (const auto& pos)
    { return std::cos(4*M_PI*(pos[0] - M_PI))*std::sin(2*M_PI*(pos[1]-M_PI/2.0)); };

    // create solution vectors
    using SolutionVector = Dune::BlockVector< Dune::FieldVector<Scalar, 1> >;
    SolutionVector sol2d_tpfa, sol2d_box, sol2d_fem;
    SolutionVector sol1d_tpfa, sol1d_box, sol1d_fem;

    Dune::Functions::interpolate(q0Basis1, sol2d_tpfa, evalFunction1);
    Dune::Functions::interpolate(q0Basis2, sol1d_tpfa, evalFunction2);

    Dune::Functions::interpolate(q1Basis1, sol2d_box, evalFunction1);
    Dune::Functions::interpolate(q1Basis2, sol1d_box, evalFunction2);

    Dune::Functions::interpolate(feGridGeometry1.feBasis(), sol2d_fem, evalFunction1);
    Dune::Functions::interpolate(feGridGeometry2.feBasis(), sol1d_fem, evalFunction2);

    // tpfa - tpfa
    writeProjection<Scalar>(tpfaGridGeometry1, tpfaGridGeometry2, sol2d_tpfa, sol1d_tpfa, "tpfa", "tpfa");
    // box - box
    writeProjection<Scalar>(boxGridGeometry1, boxGridGeometry2, sol2d_box, sol1d_box, "box", "box");
    // fem - fem
    writeProjection<Scalar>(feGridGeometry1, feGridGeometry2, sol2d_fem, sol1d_fem, "fem", "fem");
    // box - tpfa
    writeProjection<Scalar>(boxGridGeometry1, tpfaGridGeometry2, sol2d_box, sol1d_tpfa, "box", "tpfa");
    // tpfa - box
    writeProjection<Scalar>(tpfaGridGeometry1, boxGridGeometry2, sol2d_tpfa, sol1d_box, "tpfa", "box");
    // fem - tpfa
    writeProjection<Scalar>(feGridGeometry1, tpfaGridGeometry2, sol2d_fem, sol1d_tpfa, "fem", "tpfa");
    // fem - box
    writeProjection<Scalar>(feGridGeometry1, boxGridGeometry2, sol2d_fem, sol1d_box, "fem", "box");
    // tpfa - fem
    writeProjection<Scalar>(tpfaGridGeometry1, feGridGeometry2, sol2d_tpfa, sol1d_fem, "tpfa", "fem");

    return 0;
}
// //////////////////////////////////
//   Error handler
// /////////////////////////////////
catch (Dune::Exception &e) {

    std::cout << e << std::endl;
    return 1;
}
