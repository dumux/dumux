// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup Discretization
 * \brief Check the overlap size for different discretization methods
 */
#ifndef DUMUX_DISCRETIZATION_CHECK_OVERLAP_SIZE_HH
#define DUMUX_DISCRETIZATION_CHECK_OVERLAP_SIZE_HH

#include <dumux/discretization/method.hh>

namespace Dumux {

/*!
 * \ingroup Discretization
 * \brief Check if the overlap size is valid for a given discretization method
 * \note the default checks if the grid has at least an overlap of one if there are no ghosts
 * \note for sequential grids every overlap is fine
 * \note specialize this for your discretization method if the default doesn't apply
 */
template<class DiscretizationMethod>
struct CheckOverlapSize
{
    template<class GridView>
    static bool isValid(const GridView& gridView) noexcept
    { return gridView.comm().size() <= 1 || gridView.overlapSize(0) + gridView.ghostSize(0) > 0; }
};

//! specialization for the box method which requires an overlap size of 0
template<>
struct CheckOverlapSize<DiscretizationMethods::Box>
{
    template<class GridView>
    static bool isValid(const GridView& gridView) noexcept
    { return gridView.comm().size() <= 1 || gridView.overlapSize(0) == 0; }
};

//! specialization for the finite element method which requires an overlap size of 0
//! \note Overloads for bases that require overlap regions can be defined in the future
template<>
struct CheckOverlapSize<DiscretizationMethods::FEM>
{
    template<class FEBasis>
    static bool isValid(const FEBasis& feBasis) noexcept
    { return feBasis.gridView().comm().size() <= 1 || feBasis.gridView().overlapSize(0) == 0; }
};

// fc staggered requires an overlap of exactly 1
template<>
struct CheckOverlapSize<DiscretizationMethods::FCStaggered>
{
    template<class GridView>
    static bool isValid(const GridView& gridView) noexcept
    { return gridView.comm().size() <= 1 || gridView.overlapSize(0) == 1; }
};

} // end namespace Dumux

#endif
