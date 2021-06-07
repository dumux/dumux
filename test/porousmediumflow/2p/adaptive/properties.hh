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
 * \ingroup TwoPTests
 * \brief Soil contamination problem where DNAPL infiltrates a fully
 *        water saturated medium.
 */

#ifndef DUMUX_LENSPROBLEM_ADAPTIVE_PROPERTIES_HH
#define DUMUX_LENSPROBLEM_ADAPTIVE_PROPERTIES_HH

#if HAVE_DUNE_UGGRID
#include <dune/grid/uggrid.hh>
#endif
#if HAVE_DUNE_ALUGRID
#include <dune/alugrid/grid.hh>
#endif
#include <test/porousmediumflow/2p/incompressible/problem.hh>
#include "../incompressible/properties.hh"
#include "problem.hh"
#include "pointsourceproblem.hh"

// Type tags for the adaptive versions of the two-phase incompressible problem
namespace Dumux::Properties {

//! Type Tags for the adaptive tests
// Create new type tags
namespace TTag {
struct TwoPIncompressibleAdaptiveTpfa { using InheritsFrom = std::tuple<TwoPIncompressibleTpfa>; };
struct TwoPIncompressibleAdaptiveMpfa { using InheritsFrom = std::tuple<TwoPIncompressibleMpfa>; };
struct TwoPIncompressibleAdaptiveBox { using InheritsFrom = std::tuple<TwoPIncompressibleBox>; };
struct TwoPAdaptivePointSource { using InheritsFrom = std::tuple<TwoPIncompressibleAdaptiveTpfa>; };
} // end namespace TTag

//! Use non-conforming refinement in the cell-centered tests, conforming for box
#if HAVE_DUNE_ALUGRID
template<class TypeTag>
struct Grid<TypeTag, TTag::TwoPIncompressibleAdaptiveTpfa> { using type = Dune::ALUGrid<2, 2, Dune::cube, Dune::nonconforming>; };
template<class TypeTag>
struct Grid<TypeTag, TTag::TwoPIncompressibleAdaptiveMpfa> { using type = Dune::ALUGrid<2, 2, Dune::cube, Dune::nonconforming>; };
#endif
#if HAVE_DUNE_UGGRID
template<class TypeTag>
struct Grid<TypeTag, TTag::TwoPIncompressibleAdaptiveBox> { using type = Dune::UGGrid<2>; };
#endif
template<class TypeTag>
struct Problem<TypeTag, TTag::TwoPAdaptivePointSource> { using type = PointSourceTestProblem<TypeTag>; };
template<class TypeTag>
struct Problem<TypeTag, TTag::TwoPIncompressibleAdaptiveTpfa> { using type = TwoPTestProblemAdaptive<TypeTag>; };
template<class TypeTag>
struct Problem<TypeTag, TTag::TwoPIncompressibleAdaptiveMpfa> { using type = TwoPTestProblemAdaptive<TypeTag>; };
template<class TypeTag>
struct Problem<TypeTag, TTag::TwoPIncompressibleAdaptiveBox> { using type = TwoPTestProblemAdaptive<TypeTag>; };

} // end namespace Dumux::Properties

#endif
