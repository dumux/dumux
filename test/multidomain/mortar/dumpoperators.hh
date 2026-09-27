// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup MultiDomainTests
 * \brief Dense dumps of the operators on a mortar space, for an eigenvalue analysis outside
 *        the simulation.
 */
#ifndef DUMUX_MORTAR_TEST_DUMP_OPERATORS_HH
#define DUMUX_MORTAR_TEST_DUMP_OPERATORS_HH

#include <cstddef>
#include <fstream>
#include <string>

namespace Dumux::Mortar::Test {

//! Write a vector on the mortar space as one line of entries
template<typename SolutionVector>
void dumpMortarVector(const std::string& fileName, const SolutionVector& v)
{
    std::ofstream out(fileName);
    out.precision(17);
    for (std::size_t i = 0; i < v.size(); ++i)
        out << v[i][0] << (i + 1 < v.size() ? " " : "\n");
}

//! Write the dense matrix of an operator on the mortar space, one line per column
template<typename SolutionVector, typename Apply>
void dumpMortarOperator(const std::string& fileName, std::size_t numDofs, const Apply& apply)
{
    std::ofstream out(fileName);
    out.precision(17);
    SolutionVector unit(numDofs), column(numDofs);
    for (std::size_t j = 0; j < numDofs; ++j)
    {
        unit = 0.0;
        unit[j] = 1.0;
        apply(unit, column);
        for (std::size_t i = 0; i < numDofs; ++i)
            out << column[i][0] << (i + 1 < numDofs ? " " : "\n");
    }
}

} // end namespace Dumux::Mortar::Test

#endif
