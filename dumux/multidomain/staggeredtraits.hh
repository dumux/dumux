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
 * \ingroup StaggeredDiscretization
 * \brief Linear algebra traits for mixeddimension problems
 */

#ifndef DUMUX_STAGGERED_MULTIDOMAIN_TRAITS_HH
#define DUMUX_STAGGERED_MULTIDOMAIN_TRAITS_HH

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

#include "traits.hh"

namespace Dumux {
namespace Detail {
namespace Staggered {

//////////////////////////////////////////////////////////
template<template<std::size_t> class SubDomainTypeTag, std::size_t i>
struct SubDomainFVGridGeometryImpl
{ using type = GetPropType<SubDomainTypeTag<i>, Properties::GridGeometry>; };

template<template<std::size_t> class SubDomainTypeTag>
struct SubDomainFVGridGeometryImpl<SubDomainTypeTag, 0>
{ using type = typename GetPropType<SubDomainTypeTag<0>, Properties::GridGeometry>::FaceFVGridGeometryType; };

template<template<std::size_t> class SubDomainTypeTag>
struct SubDomainFVGridGeometryImpl<SubDomainTypeTag, 1>
{ using type = typename GetPropType<SubDomainTypeTag<0>, Properties::GridGeometry>::CellCenterFVGridGeometryType; };
//////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////
template<template<std::size_t> class SubDomainTypeTag, std::size_t i>
struct SubDomainGridVariablesImpl
{ using type = GetPropType<SubDomainTypeTag<i>, Properties::GridVariables>; };

template<template<std::size_t> class SubDomainTypeTag>
struct SubDomainGridVariablesImpl<SubDomainTypeTag, 0>
{ using type = typename GetPropType<SubDomainTypeTag<0>, Properties::GridVariables>::FaceGridVariablesType; };

template<template<std::size_t> class SubDomainTypeTag>
struct SubDomainGridVariablesImpl<SubDomainTypeTag, 1>
{ using type = typename GetPropType<SubDomainTypeTag<0>, Properties::GridVariables>::CellCenterGridVariablesType; };
//////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////
template<template<std::size_t> class SubDomainTypeTag, std::size_t i>
struct SubDomainPrimaryVariablesImpl
{ using type = GetPropType<SubDomainTypeTag<i>, Properties::PrimaryVariables>; };

template<template<std::size_t> class SubDomainTypeTag>
struct SubDomainPrimaryVariablesImpl<SubDomainTypeTag, 0>
{ using type = GetPropType<SubDomainTypeTag<0>, Properties::FacePrimaryVariables>; };

template<template<std::size_t> class SubDomainTypeTag>
struct SubDomainPrimaryVariablesImpl<SubDomainTypeTag, 1>
{ using type = GetPropType<SubDomainTypeTag<0>, Properties::CellCenterPrimaryVariables>; };
//////////////////////////////////////////////////////////

template<class Scalar, int numEq>
struct JacobianTypeImpl
{
    private:
        using MatrixBlock = typename Dune::FieldMatrix<Scalar, numEq, numEq>;
    public:
        using type = typename Dune::BCRSMatrix<MatrixBlock>;
};

//////////////////////////////////////////////////////////
template<template<std::size_t> class SubDomainTypeTag, std::size_t i>
struct SubDomainJacobianMatrixImpl
{ using type = GetPropType<SubDomainTypeTag<i>, Properties::JacobianMatrix>; };

template<template<std::size_t> class SubDomainTypeTag>
struct SubDomainJacobianMatrixImpl<SubDomainTypeTag, 0>
{ using type = typename JacobianTypeImpl<GetPropType<SubDomainTypeTag<0>, Properties::Scalar>,
                                         getPropValue<SubDomainTypeTag<0>, Properties::NumEqFace>()>::type; };

template<template<std::size_t> class SubDomainTypeTag>
struct SubDomainJacobianMatrixImpl<SubDomainTypeTag, 1>
{ using type = typename JacobianTypeImpl<GetPropType<SubDomainTypeTag<1>, Properties::Scalar>,
                                         getPropValue<SubDomainTypeTag<0>, Properties::NumEqCellCenter>()>::type; };
//////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////
template<template<std::size_t> class SubDomainTypeTag, std::size_t i>
struct SubDomainSolutionVectorImpl
{ using type = GetPropType<SubDomainTypeTag<i>, Properties::SolutionVector>; };

template<template<std::size_t> class SubDomainTypeTag>
struct SubDomainSolutionVectorImpl<SubDomainTypeTag, 0>
{ using type = GetPropType<SubDomainTypeTag<0>, Properties::FaceSolutionVector>; };

template<template<std::size_t> class SubDomainTypeTag>
struct SubDomainSolutionVectorImpl<SubDomainTypeTag, 1>
{ using type = GetPropType<SubDomainTypeTag<0>, Properties::CellCenterSolutionVector>; };
//////////////////////////////////////////////////////////

} // end namespace Staggered
} // end namespace Detail

/*
 * \ingroup MultiDomain
 * \ingroup StaggeredDiscretization
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
struct StaggeredMultiDomainTraits
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

    template<std::size_t id>
    using SubDomainJacobianMatrix = typename Detail::Staggered::SubDomainJacobianMatrixImpl<SubDomainTypeTag, id>::type;

    template<std::size_t id>
    using SubDomainSolutionVector = typename Detail::Staggered::SubDomainSolutionVectorImpl<SubDomainTypeTag, id>::type;

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
        using Problem = GetPropType<SubDomainTypeTag<id>, Properties::Problem>;
        using GridGeometry = typename Detail::Staggered::SubDomainFVGridGeometryImpl<SubDomainTypeTag, id>::type;
        using GridVariables = typename Detail::Staggered::SubDomainGridVariablesImpl<SubDomainTypeTag, id>::type;
        using SolutionVector = typename Detail::Staggered::SubDomainSolutionVectorImpl<SubDomainTypeTag, id>::type;
        using PrimaryVariables = typename Detail::Staggered::SubDomainPrimaryVariablesImpl<SubDomainTypeTag, id>::type;
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
