// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
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
