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
 * \ingroup Common
 * \brief Free functions to get element-local solutions.
 */
#ifndef DUMUX_DISCRETIZATION_ELEMENT_SOLUTION_HH
#define DUMUX_DISCRETIZATION_ELEMENT_SOLUTION_HH

#include <concepts>
#include <type_traits>

#include <dune/common/fvector.hh>
#include <dune/istl/bvector.hh>

#include <dumux/common/enumerate.hh>
#include <dumux/common/reservedblockvector.hh>
#include <dumux/discretization/method.hh>

#include <dumux/experimental/new_assembly/dumux/common/concepts.hh>
#include <dumux/experimental/new_assembly/dumux/common/typetraits.hh>
#include <dumux/experimental/new_assembly/dumux/common/indexstrategies.hh>
#include <dumux/experimental/new_assembly/dumux/discretization/concepts.hh>
#include <dumux/experimental/new_assembly/dumux/linear/system.hh>

namespace Dumux {

#ifndef DOXYGEN
namespace Detail {

template<typename PrimaryVariables, auto size>
struct ElementSolution;

template<typename PrimaryVariables>
struct ElementSolution<PrimaryVariables, Size::dynamic>
: public std::type_identity<Dune::BlockVector<PrimaryVariables>>
{};

template<typename PrimaryVariables, std::integral auto size>
struct ElementSolution<PrimaryVariables, size>
: public std::type_identity<ReservedBlockVector<PrimaryVariables, size>>
{};

} // namespace Detail
#endif // DOXYGEN

template<typename PrimaryVariables, Concepts::Size auto maxSize>
using ElementSolution = typename Detail::ElementSolution<PrimaryVariables, maxSize>::type;


#ifndef DOXYGEN
namespace Detail { using ElemSolMI = Dumux::MultiIndex<std::size_t, std::size_t>; }
#endif // DOXYGEN

template<std::size_t numEq,
         Concepts::CCTpfaGridGeometryLocalView LocalView,
         Concepts::IndexStrategy<Detail::ElemSolMI> IndexStrategy,
         LinearSystem::Interoperable Dofs>
auto elementSolution(const LocalView& localGridGeom,
                     const IndexStrategy& indexStrategy,
                     const Dofs& dofs)
{
    using Scalar = LinearSystem::ScalarType<Dofs>;
    using PrimaryVariables = Dune::FieldVector<Scalar, numEq>;
    ElementSolution<PrimaryVariables, 1> result(localGridGeom.numScv());
    for (const auto& [i, scv] : enumerate(scvs(localGridGeom)))
        for (std::size_t eqIdx = 0; eqIdx < numEq; ++eqIdx)
            result[i][eqIdx] = LinearSystem::get(dofs,
                indexStrategy[MultiIndex{localGridGeom.dofIndex(scv), eqIdx}]
            );
    return result;
}

template<Concepts::CCTpfaGridGeometryLocalView LocalView,
         Concepts::IndexStrategy<Detail::ElemSolMI> IndexStrategy,
         LinearSystem::Interoperable Dofs,
         std::integral auto numEq>
auto elementSolution(const LocalView& localGridGeom,
                     const IndexStrategy& indexStrategy,
                     const Dofs& dofs,
                     const std::integral_constant<decltype(numEq), numEq>&)
{ return elementSolution<numEq>(localGridGeom, indexStrategy, dofs); }

} // namespace Dumux

#endif
