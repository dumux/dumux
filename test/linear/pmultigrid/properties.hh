// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \brief The scheme (SCHEME_MODEL) and grid (GRIDTYPE) of the Laplace test of the p-multigrid preconditioner
 */
#ifndef DUMUX_TEST_LINEAR_PMULTIGRID_PROPERTIES_HH
#define DUMUX_TEST_LINEAR_PMULTIGRID_PROPERTIES_HH

#include <dune/grid/yaspgrid.hh>
#if HAVE_DUNE_ALUGRID
#include <dune/alugrid/grid.hh>
#endif

#include <dumux/discretization/pq2.hh>
#include <dumux/discretization/pq3.hh>

#include "laplace.hh"

namespace Dumux::Properties {

namespace TTag {
struct LaplacePMultigrid { using InheritsFrom = std::tuple<LaplaceTest, SCHEME_MODEL>; };
} // end namespace TTag

template<class TypeTag>
struct Grid<TypeTag, TTag::LaplacePMultigrid> { using type = GRIDTYPE; };

} // end namespace Dumux::Properties

#endif
