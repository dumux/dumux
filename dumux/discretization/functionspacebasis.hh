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
template< class GridGeometry,
          DiscretizationMethod dm = GridGeometry::discMethod >
struct FunctionSpaceBasisTraits;

/*!
 * \ingroup Discretization
 * \brief Creates a Dune::Functions object of the underlying basis
 *        of a discretization scheme that is not finite elements.
 */
template<class GridGeometry, std::enable_if_t<GridGeometry::discMethod != DiscretizationMethod::fem, int> = 0>
typename FunctionSpaceBasisTraits<GridGeometry>::GlobalBasis
getFunctionSpaceBasis(const GridGeometry& gridGeometry)
{ return {gridGeometry.gridView()}; }

/*!
 * \ingroup Discretization
 * \brief Returns the Dune::Functions object for the basis of a finite element scheme.
 */
template<class GridGeometry, std::enable_if_t<GridGeometry::discMethod == DiscretizationMethod::fem, int> = 0>
const typename FunctionSpaceBasisTraits<GridGeometry>::GlobalBasis&
getFunctionSpaceBasis(const GridGeometry& gridGeometry)
{ return gridGeometry.feBasis(); }


///////////////////////////////////////////////////////////
// Specializations of the FunctionSpaceBasisTraits class //
///////////////////////////////////////////////////////////

//! Traits specialization: box scheme uses lagrange basis of order 1
template< class GridGeometry >
struct FunctionSpaceBasisTraits<GridGeometry, DiscretizationMethod::box>
{ using GlobalBasis = Dune::Functions::LagrangeBasis<typename GridGeometry::GridView, /*order*/1>; };

//! Traits specialization: cc schemes use lagrange bases of order 0
template< class GridGeometry >
struct FunctionSpaceBasisTraits<GridGeometry, DiscretizationMethod::cctpfa>
{ using GlobalBasis = Dune::Functions::LagrangeBasis<typename GridGeometry::GridView, /*order*/0>; };

//! Traits specialization: cc schemes use lagrange bases of order 0
template< class GridGeometry >
struct FunctionSpaceBasisTraits<GridGeometry, DiscretizationMethod::ccmpfa>
{ using GlobalBasis = Dune::Functions::LagrangeBasis<typename GridGeometry::GridView, /*order*/0>; };

//! Traits specialization: fem defines its basis
template< class GridGeometry >
struct FunctionSpaceBasisTraits<GridGeometry, DiscretizationMethod::fem>
{ using GlobalBasis = typename GridGeometry::FEBasis; };

} // end namespace Dumux

#endif // HAVE_DUNE_FUNCTIONS
#endif
