// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \brief Output of the errors of a run for the convergence test driver bin/testing/dumux_convergencetest.py
 */
#ifndef DUMUX_TEST_CONVERGENCETEST_HH
#define DUMUX_TEST_CONVERGENCETEST_HH

#include <algorithm>
#include <fstream>
#include <string>
#include <utility>
#include <vector>

#include <dune/grid/common/partitionset.hh>
#include <dune/grid/common/rangegenerators.hh>

#include <dumux/geometry/diameter.hh>
#include <dumux/io/format.hh>

namespace Dumux::Test {

/*!
 * \brief Write the mesh size and the named errors of a run to <name>-errors.csv
 *
 * The file has a header row and one data row, with the mesh size in the first column.
 * It is overwritten by every run.
 */
template<class Scalar>
void writeConvergenceTestErrors(const std::string& name, Scalar h,
                                const std::vector<std::pair<std::string, Scalar>>& errors)
{
    std::ofstream file(name + "-errors.csv");
    file << "h";
    for (const auto& error : errors)
        file << "," << error.first;
    file << "\n" << Fmt::format("{:.16e}", h);
    for (const auto& error : errors)
        file << Fmt::format(",{:.16e}", error.second);
    file << "\n";
}

/*!
 * \brief The largest diameter of the elements of a grid view
 */
template<class GridView>
typename GridView::ctype maxElementDiameter(const GridView& gridView)
{
    typename GridView::ctype h = 0.0;
    for (const auto& element : elements(gridView, Dune::Partitions::interior))
        h = std::max(h, diameter(element.geometry()));
    return gridView.comm().max(h);
}

/*!
 * \brief The largest extent of the elements of a grid view along a coordinate direction
 * \note The mesh size of a series that refines in one direction only
 */
template<class GridView>
typename GridView::ctype maxElementExtent(const GridView& gridView, int direction)
{
    typename GridView::ctype h = 0.0;
    for (const auto& element : elements(gridView, Dune::Partitions::interior))
    {
        const auto geometry = element.geometry();
        auto min = geometry.corner(0)[direction];
        auto max = min;
        for (int i = 1; i < geometry.corners(); ++i)
        {
            min = std::min(min, geometry.corner(i)[direction]);
            max = std::max(max, geometry.corner(i)[direction]);
        }
        h = std::max(h, max - min);
    }
    return gridView.comm().max(h);
}

} // end namespace Dumux::Test

#endif
