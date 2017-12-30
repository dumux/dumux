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
 * \ingroup MixedDimension
 * \brief Linear algebra traits for mixeddimension problems
 */

#ifndef DUMUX_MIXEDDIMENSION_TRAITS_HH
#define DUMUX_MIXEDDIMENSION_TRAITS_HH

#include <dune/common/fvector.hh>
#include <dune/common/fmatrix.hh>
#include <dune/common/indices.hh>
#include <dune/common/promotiontraits.hh>
#include <dune/istl/bcrsmatrix.hh>
#include <dune/istl/multitypeblockvector.hh>
#include <dune/istl/multitypeblockmatrix.hh>

#include <dumux/common/properties.hh>

namespace Dumux {

template<class BTypeTag, class LTypeTag>
struct MixedDimensionTraits
{
    using BulkTypeTag = BTypeTag;
    using LowDimTypeTag = LTypeTag;

    enum {
        numEqBulk = GET_PROP_VALUE(BulkTypeTag, NumEq),
        numEqLowDim = GET_PROP_VALUE(LowDimTypeTag, NumEq)
    };

    using Scalar = typename Dune::PromotionTraits<
                                typename GET_PROP_TYPE(BulkTypeTag, Scalar),
                                typename GET_PROP_TYPE(LowDimTypeTag, Scalar)
                            >::PromotedType;

    using SolutionVectorBulk = typename GET_PROP_TYPE(BulkTypeTag, SolutionVector);
    using SolutionVectorLowDim = typename GET_PROP_TYPE(LowDimTypeTag, SolutionVector);
    using SolutionVector = Dune::MultiTypeBlockVector<SolutionVectorBulk, SolutionVectorLowDim>;

    // the sub-blocks
    using MatrixLittleBlockBulk = Dune::FieldMatrix<Scalar, numEqBulk, numEqBulk>;
    using MatrixLittleBlockBulkCoupling = Dune::FieldMatrix<Scalar, numEqBulk, numEqLowDim>;
    using MatrixLittleBlockLowDim = Dune::FieldMatrix<Scalar, numEqLowDim, numEqLowDim>;
    using MatrixLittleBlockLowDimCoupling = Dune::FieldMatrix<Scalar, numEqLowDim, numEqBulk>;

    // the BCRS matrices of the subproblems as big blocks
    using MatrixBlockBulk = Dune::BCRSMatrix<MatrixLittleBlockBulk>;
    using MatrixBlockBulkCoupling = Dune::BCRSMatrix<MatrixLittleBlockBulkCoupling>;
    using MatrixBlockLowDim = Dune::BCRSMatrix<MatrixLittleBlockLowDim>;
    using MatrixBlockLowDimCoupling = Dune::BCRSMatrix<MatrixLittleBlockLowDimCoupling>;

    // the row types
    using RowBulk = Dune::MultiTypeBlockVector<MatrixBlockBulk, MatrixBlockBulkCoupling>;
    using RowLowDim = Dune::MultiTypeBlockVector<MatrixBlockLowDimCoupling, MatrixBlockLowDim>;

    // the jacobian matrix
    using JacobianMatrix = Dune::MultiTypeBlockMatrix<RowBulk, RowLowDim>;

    // Definition of the indices of the subproblems in the global solution vector
    using bulkIdx = Dune::index_constant<0>;
    using lowDimIdx = Dune::index_constant<1>;
};

} //end namespace Dumux

#endif
