// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \brief Test for the L2 projection of analytic solutions
 */
#include <config.h>

#include <iostream>
#include <cmath>
#include <memory>
#include <array>

#include <dune/common/exceptions.hh>

#include <dune/grid/utility/structuredgridfactory.hh>
#include <dune/grid/io/file/vtk/subsamplingvtkwriter.hh>
#include <dune/grid/yaspgrid.hh>

#include <dune/functions/functionspacebases/lagrangebasis.hh>
#include <dune/functions/functionspacebases/interpolate.hh>
#include <dune/functions/gridfunctions/discreteglobalbasisfunction.hh>

#include <dumux/common/initialize.hh>

#include <dumux/common/integrate.hh>
#include <dumux/discretization/projection/l2_projection.hh>

int main (int argc, char *argv[])
{
    // maybe initialize MPI and/or multithreading backend
    Dumux::initialize(argc, argv);

    // make a 2d grid
    using Grid = Dune::YaspGrid<2>;
    using GlobalPosition = Dune::FieldVector<double, 2>;
    GlobalPosition lower(0.0);
    GlobalPosition upper(1.0);
    std::array<unsigned int, 2> els{{10, 10}};
    std::shared_ptr<Grid> grid = Dune::StructuredGridFactory<Grid>::createCubeGrid(lower, upper, els);

    using GridView = typename Grid::LeafGridView;
    using FEBasis = Dune::Functions::LagrangeBasis<GridView, 2>;
    const auto feBasis = std::make_shared<FEBasis>(grid->leafGridView());

    auto source = [] (const auto& x) { return std::sin(2*M_PI*x[0])*std::sin(2*M_PI*x[1]); };
    Dumux::L2Projection projection(*feBasis);
    const auto projectedCoefficients = projection.project(source);

    auto interpolatedCoefficients = projectedCoefficients;
    Dune::Functions::interpolate(*feBasis, interpolatedCoefficients, source);

    auto analytic = Dune::Functions::makeAnalyticGridViewFunction(source, grid->leafGridView());
    auto projected = Dune::Functions::makeDiscreteGlobalBasisFunction<double>(*feBasis, projectedCoefficients);
    auto interpolated = Dune::Functions::makeDiscreteGlobalBasisFunction<double>(*feBasis, interpolatedCoefficients);

    Dune::SubsamplingVTKWriter<GridView> vtkWriter(grid->leafGridView(), Dune::refinementLevels(2));
    vtkWriter.addVertexData(analytic, Dune::VTK::FieldInfo("u", Dune::VTK::FieldInfo::Type::scalar, 1));
    vtkWriter.addVertexData(projected, Dune::VTK::FieldInfo("u_p", Dune::VTK::FieldInfo::Type::scalar, 1));
    vtkWriter.addVertexData(interpolated, Dune::VTK::FieldInfo("u_i", Dune::VTK::FieldInfo::Type::scalar, 1));
    vtkWriter.write("test_l2_projection");

    auto l2ErrorProjected = Dumux::integrateL2Error(grid->leafGridView(), analytic, projected, 5);
    auto l2ErrorInterpolated = Dumux::integrateL2Error(grid->leafGridView(), analytic, interpolated, 5);

    std::cout << "L2-Error of L2-projected discrete solution: " << l2ErrorProjected << std::endl;
    std::cout << "L2-Error of interpolated discrete solution: " << l2ErrorInterpolated << std::endl;

    if (l2ErrorProjected > l2ErrorInterpolated)
        DUNE_THROW(Dune::Exception, "L2-error of L2-projected solution is larger than the interpolated solution");

    return 0;
}
