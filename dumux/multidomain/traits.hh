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
 * \ingroup MultiDomain
 * \brief Traits for multidomain problems
 */

#ifndef DUMUX_MULTIDOMAIN_TRAITS_HH
#define DUMUX_MULTIDOMAIN_TRAITS_HH

#include <type_traits>
#include <tuple>
#include <utility>
#include <memory>

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
    static_assert(std::conjunction_v<isBCRSMatrix<JacobianBlocks>...>, "Jacobian blocks have to be BCRSMatrices!");

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

//! helper alias to create a tuple of shared_ptr<...> from an indexed type
template<template<std::size_t> class T, class Indices>
struct MultiDomainTupleSharedPtr
{
    template<std::size_t i>
    using PtrType = std::shared_ptr<T<i>>;

    using type = typename makeFromIndexedType<std::tuple, PtrType, Indices>::type;
};

//! helper alias to create a tuple of shared_ptr<const ...> from an indexed type
template<template<std::size_t> class T, class Indices>
struct MultiDomainTupleSharedPtrConst
{
    template<std::size_t i>
    using PtrType = std::shared_ptr<const T<i>>;

    using type = typename makeFromIndexedType<std::tuple, PtrType, Indices>::type;
};

//! helper alias to create the JacobianMatrix type
template<template<std::size_t> class SubDomainDiagBlocks, class Indices, class Scalar>
struct MultiDomainMatrixType
{
    template<typename... MatrixBlocks>
    using M = typename createMultiTypeBlockMatrixType<Scalar, MatrixBlocks...>::type::type;

    using type = typename makeFromIndexedType<M, SubDomainDiagBlocks, Indices>::type;
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
    //! the number of subdomains
    static constexpr std::size_t numSubDomains = sizeof...(SubDomainTypeTags);

private:

    //! the type tag of a sub domain problem
    template<std::size_t id>
    using SubDomainTypeTag = typename std::tuple_element_t<id, std::tuple<SubDomainTypeTags...>>;

    //! helper alias to construct derived multidomain types like tuples
    using Indices = std::make_index_sequence<numSubDomains>;

    //! the scalar type of each sub domain
    template<std::size_t id>
    using SubDomainScalar = GetPropType<SubDomainTypeTag<id>, Properties::Scalar>;

    //! the jacobian type of each sub domain
    template<std::size_t id>
    using SubDomainJacobianMatrix = GetPropType<SubDomainTypeTag<id>, Properties::JacobianMatrix>;

    //! the solution type of each sub domain
    template<std::size_t id>
    using SubDomainSolutionVector = GetPropType<SubDomainTypeTag<id>, Properties::SolutionVector>;

public:

    /*
     * \brief sub domain types
     */
    //\{

    template<std::size_t id>
    struct SubDomain
    {
        using Index = Dune::index_constant<id>;
        using TypeTag = SubDomainTypeTag<id>;
        using Grid = GetPropType<SubDomainTypeTag<id>, Properties::Grid>;
        using GridGeometry = GetPropType<SubDomainTypeTag<id>, Properties::GridGeometry>;
        using Problem = GetPropType<SubDomainTypeTag<id>, Properties::Problem>;
        using GridVariables =GetPropType<SubDomainTypeTag<id>, Properties::GridVariables>;
        using IOFields = GetPropType<SubDomainTypeTag<id>, Properties::IOFields>;
        using SolutionVector = GetPropType<SubDomainTypeTag<id>, Properties::SolutionVector>;
    };

    //\}

    /*
     * \brief multi domain types
     */
    //\{

    //! the scalar type
    using Scalar = typename makeFromIndexedType<std::common_type_t, SubDomainScalar, Indices>::type;

    //! the solution vector type
    using SolutionVector = typename makeFromIndexedType<Dune::MultiTypeBlockVector, SubDomainSolutionVector, Indices>::type;

    //! the jacobian type
    using JacobianMatrix = typename Detail::MultiDomainMatrixType<SubDomainJacobianMatrix, Indices, Scalar>::type;

    //\}

    /*
     * \brief helper aliases to contruct derived tuple types
     */
    //\{

    //! helper alias to create tuple<...> from indexed type
    template<template<std::size_t> class T>
    using Tuple = typename makeFromIndexedType<std::tuple, T, Indices>::type;

    //! helper alias to create tuple<std::shared_ptr<...>> from indexed type
    template<template<std::size_t> class T>
    using TupleOfSharedPtr = typename Detail::MultiDomainTupleSharedPtr<T, Indices>::type;

    //! helper alias to create tuple<std::shared_ptr<const ...>> from indexed type
    template<template<std::size_t> class T>
    using TupleOfSharedPtrConst = typename Detail::MultiDomainTupleSharedPtrConst<T, Indices>::type;

    //\}
};

} //end namespace Dumux

#endif
