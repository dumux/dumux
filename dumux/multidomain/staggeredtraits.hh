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
    static constexpr std::size_t numSubDomains = sizeof...(SubDomainTypeTags);

    //! the type tag of a sub domain problem
    template<std::size_t id>
    using SubDomainTypeTag = typename std::tuple_element_t<id, std::tuple<SubDomainTypeTags...>>;

    //! the static domain indices
    template<std::size_t id>
    using DomainIdx = Dune::index_constant<id>;

private:
    using Indices = std::make_index_sequence<numSubDomains>;

    template<std::size_t id>
    using SolutionSubVector = std::conditional_t<(id < 2),
                                                 std::conditional_t<(id == 0),
                                                                    typename GET_PROP_TYPE(SubDomainTypeTag<0>, CellCenterSolutionVector),
                                                                    typename GET_PROP_TYPE(SubDomainTypeTag<0>, FaceSolutionVector)>,
                                                 typename GET_PROP_TYPE(SubDomainTypeTag<id>, SolutionVector)>;

    template<std::size_t id>
    using SubDomainScalar = typename GET_PROP_TYPE(SubDomainTypeTag<id>, Scalar);

    template<std::size_t id>
    using SubDomainProblem = std::shared_ptr<const typename GET_PROP_TYPE(SubDomainTypeTag<id>, Problem)>;

    template<std::size_t id>
    using SubDomainFVGridGeometry = std::shared_ptr<std::conditional_t<(id < 2),
                                                                       std::conditional_t<(id == 0),
                                                                                          typename GET_PROP_TYPE(SubDomainTypeTag<0>, FVGridGeometry)::CellCenterFVGridGeometryType,
                                                                                          typename GET_PROP_TYPE(SubDomainTypeTag<0>, FVGridGeometry)::FaceFVGridGeometryType>,
                                                                       typename GET_PROP_TYPE(SubDomainTypeTag<id>, FVGridGeometry)>>;

    template<std::size_t id>
    using SubDomainGridVariables = std::shared_ptr<std::conditional_t<(id < 2),
                                                                      std::conditional_t<(id == 0),
                                                                                         typename GET_PROP_TYPE(SubDomainTypeTag<0>, GridVariables)::CellCenterGridVariablesType,
                                                                                         typename GET_PROP_TYPE(SubDomainTypeTag<0>, GridVariables)::FaceGridVariablesType>,
                                                                      typename GET_PROP_TYPE(SubDomainTypeTag<id>, GridVariables)>>;

    template<class Scalar, int numEq>
    struct JacobianType
    {
        private:
            using MatrixBlock = typename Dune::FieldMatrix<Scalar, numEq, numEq>;
        public:
            using type = typename Dune::BCRSMatrix<MatrixBlock>;
    };


    template<std::size_t id>
    using JacobianDiagBlock = std::conditional_t<(id < 2),
                                                 std::conditional_t<(id == 0),
                                                                    typename JacobianType<typename GET_PROP_TYPE(SubDomainTypeTag<0>, Scalar), GET_PROP_VALUE(SubDomainTypeTag<0>, NumEqCellCenter)>::type,
                                                                    typename JacobianType<typename GET_PROP_TYPE(SubDomainTypeTag<0>, Scalar), GET_PROP_VALUE(SubDomainTypeTag<0>, NumEqFace)>::type>,
                                                 typename GET_PROP_TYPE(SubDomainTypeTag<id>, JacobianMatrix)>;

public:

    //! the scalar type
    using Scalar = typename makeFromIndexedType<std::common_type_t, SubDomainScalar, Indices>::type;

    //! the solution vector type
    using SolutionVector = typename makeFromIndexedType<Dune::MultiTypeBlockVector, SolutionSubVector, Indices>::type;

    template<std::size_t id>
    using PrimaryVariables = std::conditional_t<(id < 2),
                                                 std::conditional_t<(id == 0), typename GET_PROP_TYPE(SubDomainTypeTag<0>, CellCenterPrimaryVariables),
                                                                               typename GET_PROP_TYPE(SubDomainTypeTag<0>, FacePrimaryVariables)>,
                                                 typename GET_PROP_TYPE(SubDomainTypeTag<id>, PrimaryVariables)>;

    template<typename... MatrixBlocks>
    using createMatrixType = typename Detail::createMultiTypeBlockMatrixType<Scalar, MatrixBlocks...>::type::type;

    //! the jacobian type
    using JacobianMatrix = typename makeFromIndexedType<createMatrixType, JacobianDiagBlock, Indices>::type;

    //! the tuple of problems
    using ProblemTuple = typename makeFromIndexedType<std::tuple, SubDomainProblem, Indices>::type;

    //! the tuple of fv grid geometries
    using FVGridGeometryTuple = typename makeFromIndexedType<std::tuple, SubDomainFVGridGeometry, Indices>::type;

    //! the tuple of grid variables
    using GridVariablesTuple = typename makeFromIndexedType<std::tuple, SubDomainGridVariables, Indices>::type;

    //! convenience alias to create tuple from type
    template<template<std::size_t> class T>
    using MakeTuple = typename makeFromIndexedType<std::tuple, T, Indices>::type;
};

} //end namespace Dumux

#endif
