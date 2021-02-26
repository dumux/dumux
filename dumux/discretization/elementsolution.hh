// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   See the file COPYING for full copying permissions.                      *
 *                                                                           *
 *   This program is free software: you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation, either version 3 of the License, or       *
 *   (at your option) any later version.                                     *
 *                                                                           *
 *   This program is distributed in the hope that it will be useful,         *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of          *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the            *
 *   GNU General Public License for more details.                            *
 *                                                                           *
 *   You should have received a copy of the GNU General Public License       *
 *   along with this program.  If not, see <http://www.gnu.org/licenses/>.   *
 *****************************************************************************/
/*!
 * \file
 * \ingroup Discretization
 * \brief Element solution classes and factory functions
 */
#ifndef DUMUX_DISCRETIZATION_ELEMENT_SOLUTION_HH
#define DUMUX_DISCRETIZATION_ELEMENT_SOLUTION_HH

#include <dumux/discretization/method.hh>
#include <dumux/discretization/cellcentered/elementsolution.hh>
#include <dumux/discretization/box/elementsolution.hh>
#include <dumux/discretization/staggered/elementsolution.hh>

namespace Dumux {

namespace Detail {
    template<class GGLV, class PV, DiscretizationMethod dm>
    struct ElementSolutionHelper;
} // end namespace Detail

struct EmptyElementSolution {};

/*!
 * \ingroup Discretization
 * \brief Convenience alias for discretization-specific element solution.
 * \tparam GGLocalView Local view on the grid geometry
 * \tparam PrimaryVariables The primary variables type
 */
template<class GGLocalView, class PrimaryVariables>
using ElementSolution
    = typename Detail::ElementSolutionHelper<GGLocalView,
                                             PrimaryVariables,
                                             GGLocalView::GridGeometry::discMethod>::Type;

namespace Detail {

    // TODO: We could also require a generic template interface for elemsols
    template<class GGLV, class PV>
    struct ElementSolutionHelper<GGLV, PV, DiscretizationMethod::box>
    { using Type = BoxElementSolution<GGLV, PV>; };

    template<class GGLV, class PV>
    struct ElementSolutionHelper<GGLV, PV, DiscretizationMethod::cctpfa>
    { using Type = CCElementSolution<GGLV, PV>; };

    template<class GGLV, class PV>
    struct ElementSolutionHelper<GGLV, PV, DiscretizationMethod::ccmpfa>
    { using Type = CCElementSolution<GGLV, PV>; };

    template<class GGLV, class PV>
    struct ElementSolutionHelper<GGLV, PV, DiscretizationMethod::staggered>
    { using Type = StaggeredElementSolution<PV>; };

} // end namespace Detail
} // end namespace Dumux

#endif
