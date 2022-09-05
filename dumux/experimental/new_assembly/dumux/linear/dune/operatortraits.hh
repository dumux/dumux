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
 * \ingroup Linear
 * \brief Traits related to dune linear operators
 */
#ifndef DUMUX_LINEAR_DUNE_OPERATOR_TRAITS_HH
#define DUMUX_LINEAR_DUNE_OPERATOR_TRAITS_HH

#include <type_traits>

#include <dune/istl/operators.hh>

#include <dumux/experimental/new_assembly/dumux/linear/operator.hh>
#include <dumux/experimental/new_assembly/dumux/linear/dune/detail.hh>

namespace Dumux::Linear::Traits {

template<LinearSystem::Detail::DuneMatrix M,
         LinearSystem::Detail::DuneVector V>
struct MatrixOperator<M, V>
: public std::type_identity<Dune::MatrixAdapter<M, V, V>>
{};

} // namespace Dumux::Linear::Traits

#endif
