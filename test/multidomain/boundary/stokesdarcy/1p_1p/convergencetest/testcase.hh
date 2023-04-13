// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup BoundaryTests
 * \brief The different test cases
 */

#ifndef DUMUX_CONVERGENCE_TEST_TESTCASE_HH
#define DUMUX_CONVERGENCE_TEST_TESTCASE_HH

namespace Dumux {

enum class TestCase
{
    ShiueExampleOne, ShiueExampleTwo, Rybak, Schneider
};

} // end namespace Dumux

#endif
