// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup Discretization
 * \brief Element solution classes and factory functions
 */
#ifndef DUMUX_DISCRETIZATION_ELEMENT_SOLUTION_HH
#define DUMUX_DISCRETIZATION_ELEMENT_SOLUTION_HH

#include <dumux/discretization/cellcentered/elementsolution.hh>
#include <dumux/discretization/cvfe/elementsolution.hh>
#include <dumux/discretization/staggered/elementsolution.hh>

namespace Dumux {

struct EmptyElementSolution {};

} // end namespace Dumux

#endif
