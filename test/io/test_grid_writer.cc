// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \brief Test for reading/writing from/to different formats with the generic grid writers & readers.
 */
#include <config.h>

#include <dune/grid/yaspgrid.hh>
#if HAVE_DUNE_FUNCTIONS
#include <dune/functions/gridfunctions/analyticgridviewfunction.hh>
#endif

#include <dumux/common/initialize.hh>
#include <dumux/io/gridwriter.hh>

using namespace Dumux::IO;

template<class Writer>
void setFields(Writer& writer)
{
    writer.setPointField("pcoords", [] (const auto& p) {
        return p.geometry().center();
    });
    writer.setCellField("ccoords", [] (const auto& e) {
        return e.geometry().center();
    });
}

template<typename GridView>
auto writeFirstOrder(const GridView& gridView)
{
    GridWriter writer{Format::vtu, gridView, order<1>};
    setFields(writer);
    std::cout << "Wrote " << writer.write("test_grid_writer_first_order") << std::endl;
}

template<typename GridView>
auto writeSecondOrder(const GridView& gridView)
{
#if HAVE_DUNE_FUNCTIONS
    GridWriter writer{Format::vtu, gridView, order<2>};
    auto f = Dune::Functions::makeAnalyticGridViewFunction([&] (const auto& x) {
        return x[0]*x[1];
    }, gridView);
    writer.setPointField("point_function", f);
    writer.setCellField("cell_function", f);
    std::cout << "Wrote " << writer.write("test_grid_writer_second_order") << std::endl;
#endif
}

template<typename GridView>
auto writeSpecifiedFormat(const GridView& gridView)
{
    GridWriter writer{Format::vtu.with({.encoder = Encoding::raw}), gridView};
    setFields(writer);
    std::cout << "Wrote " << writer.write("test_grid_writer_format_selection") << std::endl;;
}

int main(int argc, char** argv)
{
    Dumux::initialize(argc, argv);

    Dune::YaspGrid<2> grid{
        {1.0, 1.0},
        {20, 20},
        std::bitset<2>(), // periodicity
        0                // overlap
    };

    writeFirstOrder(grid.leafGridView());
    writeSecondOrder(grid.leafGridView());
    writeSpecifiedFormat(grid.leafGridView());

    return 0;
}
