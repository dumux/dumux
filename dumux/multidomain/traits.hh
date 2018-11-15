// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   See the file COPYING for full copying permissions.                      *
 *                                                                           *
 *   This program is free software: you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation, either version 2 of the License, or       *
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
 * \ingroup MultiDomain
 * \brief Traits for multidomain problems
 */

#ifndef DUMUX_MULTIDOMAIN_TRAITS_HH
#define DUMUX_MULTIDOMAIN_TRAITS_HH

#include <type_traits>
#include <tuple>
#include <utility>

#include <dune/common/fmatrix.hh>
#include <dune/common/indices.hh>

#include <dune/istl/bcrsmatrix.hh>
#include <dune/istl/multitypeblockvector.hh>
#include <dune/istl/multitypeblockmatrix.hh>

#include <dumux/common/properties.hh>
#include <dumux/common/typetraits/matrix.hh>
#include <dumux/common/typetraits/utility.hh>

namespace Dumux {

namespace Detail {

//! a helper class to create a multitype matrix given the diagonal matrix blocks
template<class Scalar, class... JacobianBlocks>
class createMultiTypeBlockMatrixType
{
    //! TODO: replace by std::conjuction in C++17
    template<bool...> struct boolPack;
    template<bool... bools>
    using all_true = std::is_same<boolPack<bools..., true>, boolPack<true, bools...>>;

    static_assert(all_true<isBCRSMatrix<JacobianBlocks>::value...>::value, "Jacobian blocks have to be BCRSMatrices!");

    template<std::size_t id>
    using JacobianDiagBlock = typename std::tuple_element_t<id, std::tuple<JacobianBlocks...>>;

    template<std::size_t id>
    static constexpr decltype(auto) numEq()
    { return JacobianDiagBlock<id>::block_type::rows; }

    template <std::size_t id, class I> struct makeRow;

    template <std::size_t id, std::size_t... Is>
    struct makeRow<id, std::index_sequence<Is...>>
    {
        using type = Dune::MultiTypeBlockVector<Dune::BCRSMatrix<Dune::FieldMatrix<Scalar, numEq<id>(), numEq<Is>()>>...>;
    };

    template <class I> struct makeMatrix;

    template <std::size_t... Is>
    struct makeMatrix<std::index_sequence<Is...>>
    {
        using type = Dune::MultiTypeBlockMatrix<typename makeRow<Is, std::index_sequence<Is...>>::type...>;
    };

    using Indices = std::index_sequence_for<JacobianBlocks...>;
public:
    using type = typename makeMatrix<Indices>::type;
};

} // end namespace Detail

/*
 * \ingroup MultiDomain
 * \brief A traits class every multidomain model has to provide
 * \tparam SubDomainTypeTags the TypeTags of the sub domain problems
 * \note should export the types
 * \code
 *       //! the type tag of the sub domain problem with id
 *       template<std::size_t id>
 *       using SubDomainTypeTag = ...
 *
 *       //! the index to access sub domain matrices and vectors
 *       //! to use with multitype matrices and vectors
 *       template<std::size_t id>
 *       using DomainIdx = ...
 *
 *       //! the scalar type
 *       using Scalar = ...
 *
 *       //! the solution vector type
 *       using SolutionVector = ...
 *
 *       //! the jacobian type
 *       using JacobianMatrix = ...
 * \endcode
 */
template<typename... SubDomainTypeTags>
struct MultiDomainTraits
{
    static constexpr std::size_t numSubDomains = sizeof...(SubDomainTypeTags);

    //! the type tag of a sub domain problem
    template<std::size_t id>
    using SubDomainTypeTag = typename std::tuple_element_t<id, std::tuple<SubDomainTypeTags...>>;

    //! the static domain indices
    template<std::size_t id>
    using DomainIdx = Dune::index_constant<id>;

    //! the sub domain geometry
    template<std::size_t id>
    using SubDomainFVGridGeometry = GetPropType<SubDomainTypeTag<id>, Properties::FVGridGeometry>;

private:
    using Indices = std::make_index_sequence<numSubDomains>;

    template<std::size_t id>
    using SolutionSubVector = GetPropType<SubDomainTypeTag<id>, Properties::SolutionVector>;

    template<std::size_t id>
    using SubDomainScalar = GetPropType<SubDomainTypeTag<id>, Properties::Scalar>;

    template<std::size_t id>
    using SubDomainProblem = std::shared_ptr<const GetPropType<SubDomainTypeTag<id>, Properties::Problem>>;

    template<std::size_t id>
    using SubDomainFVGridGeometryPtr = std::shared_ptr<const GetPropType<SubDomainTypeTag<id>, Properties::FVGridGeometry>>;

    template<std::size_t id>
    using SubDomainGridVariables = std::shared_ptr<GetPropType<SubDomainTypeTag<id>, Properties::GridVariables>>;

    template<std::size_t id>
    using JacobianDiagBlock = GetPropType<SubDomainTypeTag<id>, Properties::JacobianMatrix>;

public:

    //! the scalar type
    using Scalar = typename makeFromIndexedType<std::common_type_t, SubDomainScalar, Indices>::type;

    //! the solution vector type
    using SolutionVector = typename makeFromIndexedType<Dune::MultiTypeBlockVector, SolutionSubVector, Indices>::type;

    template<typename... MatrixBlocks>
    using createMatrixType = typename Detail::createMultiTypeBlockMatrixType<Scalar, MatrixBlocks...>::type::type;

    //! the jacobian type
    using JacobianMatrix = typename makeFromIndexedType<createMatrixType, JacobianDiagBlock, Indices>::type;

    //! the tuple of problems
    using ProblemTuple = typename makeFromIndexedType<std::tuple, SubDomainProblem, Indices>::type;

    //! the tuple of fv grid geometries
    using FVGridGeometryTuple = typename makeFromIndexedType<std::tuple, SubDomainFVGridGeometryPtr, Indices>::type;

    //! the tuple of grid variables
    using GridVariablesTuple = typename makeFromIndexedType<std::tuple, SubDomainGridVariables, Indices>::type;

    //! convenience alias to create tuple from type
    template<template<std::size_t> class T>
    using MakeTuple = typename makeFromIndexedType<std::tuple, T, Indices>::type;
};

} //end namespace Dumux

#endif
