// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \brief Test for writing parallel files with the generic grid file format writer.
 */
#include <config.h>

#include <dune/grid/yaspgrid.hh>
#if HAVE_DUNE_FUNCTIONS
#include <dune/functions/gridfunctions/analyticgridviewfunction.hh>
#endif

#include <dumux/common/initialize.hh>
#include <dumux/io/gridwriter.hh>

// The rank is added as cell field, although the writer
// writes it out as meta data per default (which is much
// more space-efficient). However, VTK seems to have issues
// displaying field data in parallel, and I could not find
// any helpful info in the code or documentation.

using namespace Dumux::IO;

template<class GlobalPosition>
double testFunction(const GlobalPosition& x)
{ return x[0]*x[1]; }

template<class Writer>
void addFields(Writer& writer, int rank)
{
    writer.setPointField("pfield", [=] (const auto& p) {
        return testFunction(p.geometry().center());
    });
    writer.setCellField("cfield", [=] (const auto& e) {
        return testFunction(e.geometry().center());
    });
    writer.setCellField("rank", [=] (const auto& e) {
        return rank;
    });
}

template<class GridView>
void testFirstOrder(const GridView& gridView)
{
    GridWriter writer{Format::vti, gridView};
    addFields(writer, gridView.comm().rank());
    std::cout << "Wrote " << writer.write("test_grid_writer_parallel") << std::endl;
}

template<class GridView>
void testSecondOrder(const GridView& gridView)
{
#if HAVE_DUNE_FUNCTIONS
    GridWriter writer{Format::vtu, gridView, order<2>};
    auto f = Dune::Functions::makeAnalyticGridViewFunction([] (const auto& x) {
        return testFunction(x);
    }, gridView);
    auto rank = Dune::Functions::makeAnalyticGridViewFunction([rank=gridView.comm().rank()] (const auto& x) {
        return rank;
    }, gridView);
    writer.setPointField("pfunc", f);
    writer.setCellField("cfunc", f);
    writer.setCellField("rank", rank);
    std::cout << "Wrote " << writer.write("test_grid_writer_parallel_second_order") << std::endl;
#endif
}

int main(int argc, char** argv)
{
    Dumux::initialize(argc, argv);

    // write a parallel vti file
    Dune::YaspGrid<2> grid{
        {1.0, 1.0},
        {20, 20},
        std::bitset<2>(), // periodicity
        0                 // overlap
    };
    grid.loadBalance();
    const auto& gridView = grid.leafGridView();

    testFirstOrder(gridView);
    testSecondOrder(gridView);

    return 0;
}
