// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup Discretization
 * \brief Provides helper aliases and functionality to obtain the types
 *        and instances of Dune::Functions function space bases that
 *        underlie different discretization schemes
 */
#ifndef DUMUX_DISCRETIZATION_FUNCTION_SPACE_BASIS_HH
#define DUMUX_DISCRETIZATION_FUNCTION_SPACE_BASIS_HH

#if HAVE_DUNE_FUNCTIONS

#include <dumux/discretization/method.hh>
#include <dune/functions/functionspacebases/lagrangebasis.hh>

namespace Dumux {

/*!
 * \ingroup Discretization
 * \brief Traits class that contains information on the
 *        function space basis used by the discretization
 *        scheme underlying a grid geometry class. Specializations
 *        of the traits for different schemes are provided below.
 */
template< class GridGeometry, class DiscretizationMethod = typename GridGeometry::DiscretizationMethod >
struct FunctionSpaceBasisTraits;

/*!
 * \ingroup Discretization
 * \brief Creates a Dune::Functions object of the underlying basis
 *        of a discretization scheme that is not finite elements.
 */
template<class GridGeometry, std::enable_if_t<GridGeometry::discMethod != DiscretizationMethods::fem, int> = 0>
typename FunctionSpaceBasisTraits<GridGeometry>::GlobalBasis
getFunctionSpaceBasis(const GridGeometry& gridGeometry)
{ return {gridGeometry.gridView()}; }

/*!
 * \ingroup Discretization
 * \brief Returns the Dune::Functions object for the basis of a finite element scheme.
 */
template<class GridGeometry, std::enable_if_t<GridGeometry::discMethod == DiscretizationMethods::fem, int> = 0>
const typename FunctionSpaceBasisTraits<GridGeometry>::GlobalBasis&
getFunctionSpaceBasis(const GridGeometry& gridGeometry)
{ return gridGeometry.feBasis(); }


///////////////////////////////////////////////////////////
// Specializations of the FunctionSpaceBasisTraits class //
///////////////////////////////////////////////////////////

//! Traits specialization: box scheme uses lagrange basis of order 1
template< class GridGeometry >
struct FunctionSpaceBasisTraits<GridGeometry, DiscretizationMethods::Box>
{ using GlobalBasis = Dune::Functions::LagrangeBasis<typename GridGeometry::GridView, /*order*/1>; };

//! Traits specialization: cc schemes use lagrange bases of order 0
template< class GridGeometry >
struct FunctionSpaceBasisTraits<GridGeometry, DiscretizationMethods::CCTpfa>
{ using GlobalBasis = Dune::Functions::LagrangeBasis<typename GridGeometry::GridView, /*order*/0>; };

//! Traits specialization: cc schemes use lagrange bases of order 0
template< class GridGeometry >
struct FunctionSpaceBasisTraits<GridGeometry, DiscretizationMethods::CCMpfa>
{ using GlobalBasis = Dune::Functions::LagrangeBasis<typename GridGeometry::GridView, /*order*/0>; };

//! Traits specialization: fem defines its basis
template< class GridGeometry >
struct FunctionSpaceBasisTraits<GridGeometry, DiscretizationMethods::FEM>
{ using GlobalBasis = typename GridGeometry::FEBasis; };

} // end namespace Dumux

#endif // HAVE_DUNE_FUNCTIONS
#endif
