// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup StaggeredDiscretization
 * \brief The local element solution class for staggered methods
 */
#ifndef DUMUX_STAGGERED_ELEMENT_SOLUTION_HH
#define DUMUX_STAGGERED_ELEMENT_SOLUTION_HH

#include <type_traits>
#include <dune/istl/bvector.hh>
#include <dumux/discretization/method.hh>

namespace Dumux {

/*!
 * \ingroup StaggeredDiscretization
 * \brief  Helper function to create a PrimaryVariables object from CellCenterPrimaryVariables
 * \tparam PrimaryVariables The type of the desired primary variables object
 * \tparam CellCenterPrimaryVariables The type of the cell center (input) primary variables object
 * \param cellCenterPriVars The cell center (input) primary variables object
 */
template<class PrimaryVariables, class CellCenterPrimaryVariables>
PrimaryVariables makePriVarsFromCellCenterPriVars(const CellCenterPrimaryVariables& cellCenterPriVars)
{
    static_assert(int(PrimaryVariables::dimension) > int(CellCenterPrimaryVariables::dimension),
                  "PrimaryVariables' size must be greater than the one of CellCenterPrimaryVariables");

    PrimaryVariables priVars(0.0);
    constexpr auto offset = PrimaryVariables::dimension - CellCenterPrimaryVariables::dimension;
    for (std::size_t i = 0; i < cellCenterPriVars.size(); ++i)
        priVars[i + offset] = cellCenterPriVars[i];
    return priVars;
}

template<class PrimaryVariables>
using StaggeredElementSolution = Dune::BlockVector<PrimaryVariables>;

/*!
 * \ingroup Discretization
 * \brief  Make an element solution for staggered schemes
 * \note This is e.g. used to construct an element solution at Dirichlet boundaries
 */
template<class FVElementGeometry, class PrimaryVariables>
auto elementSolution(PrimaryVariables&& priVars)
-> std::enable_if_t<FVElementGeometry::GridGeometry::discMethod == DiscretizationMethods::staggered,
                    StaggeredElementSolution<PrimaryVariables>>
{
    return StaggeredElementSolution<PrimaryVariables>({std::move(priVars)});
}

/*!
 * \ingroup Discretization
 * \brief  Helper function to create an elementSolution from cell center primary variables
 * \tparam PrimaryVariables The type of the desired primary variables object
 * \tparam CellCenterPrimaryVariables The type of the cell center (input) primary variables object
 * \param cellCenterPriVars The cell center (input) primary variables object
 */
template<class PrimaryVariables, class CellCenterPrimaryVariables>
StaggeredElementSolution<PrimaryVariables> makeElementSolutionFromCellCenterPrivars(const CellCenterPrimaryVariables& cellCenterPriVars)
{
    return StaggeredElementSolution<PrimaryVariables>({makePriVarsFromCellCenterPriVars<PrimaryVariables>(cellCenterPriVars)});
}

} // end namespace Dumux

#endif
