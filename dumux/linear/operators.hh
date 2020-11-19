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
 * \brief Linear operators for iterative solvers
 */
#ifndef DUMUX_LINEAR_OPERATORS_HH
#define DUMUX_LINEAR_OPERATORS_HH

#include <tuple>

#include <dune/common/hybridutilities.hh>
#include <dune/istl/solvercategory.hh>
#include <dune/istl/operators.hh>

namespace Dumux {

/*!
 * \ingroup Linear
 * \brief A linear operator based on a tuple of linear operators
 *
 * \tparam Vector usually a MultiTypeBlockVector
 * \tparam Matrix usually a MultiTypeBlockMatrix
 * \tparam LinearOperatorTuple tuple of linear operators
 *
 * Vector, Matrix and LinearOperatorTuple have to be consistent in the sense
 * that a tuple element implements the linear operator for a corresponding
 * diagonal block of a Matrix and block of a Vector.
 */
template<class Vector, class Matrix, class LinearOperatorTuple>
class TupleLinearOperator: public Dune::LinearOperator<Vector, Vector> {
public:
    //! The type of the domain of the operator.
    typedef Vector domain_type;
    //! The type of the range of the operator.
    typedef Vector range_type;
    //! The field type of the operator.
    typedef typename Vector::field_type field_type;

    /*! \brief Construct the linear operator
     *
     *  \param lops tuple of linear operators
     *  \param m matrix
     *
     *  The linear operator tuple and the matrix are both needed.
     *  The first realizes the action of the diagonal blocks and additionally
     *  carries the parallel information. With the second, the action of the
     *  off-diagonal blocks is implemented.
     */
    TupleLinearOperator (const LinearOperatorTuple& lops, const Matrix& m)
    : linearOperators_(lops), matrix_(m)
    {}

    /*! \brief Apply operator to x:  \f$ y = A(x) \f$
     *
     *  The input vector is consistent and the output must also be
     *  consistent on the interior+border partition.
     */
    void apply (const Vector& x, Vector& y) const override
    {
        using namespace Dune::Hybrid;
        forEach(integralRange(Dune::Hybrid::size(x)), [&](const auto i)
        {
            std::get<i>(linearOperators_)->apply(x[i], y[i]);

            forEach(integralRange(Dune::Hybrid::size(x)), [&](const auto j)
            {
                if (i != j)
                    matrix_[i][j].umv(x[j], y[i]);
            });
        });
    }

    //! Apply operator to x, scale and add:  \f$ y = y + \alpha A(x) \f$
    void applyscaleadd (field_type alpha, const Vector& x, Vector& y) const override
    {
        using namespace Dune::Hybrid;
        forEach(integralRange(Dune::Hybrid::size(x)), [&](const auto i)
        {
            std::get<i>(linearOperators_)->applyscaleadd(alpha, x[i], y[i]);

            forEach(integralRange(Dune::Hybrid::size(x)), [&](const auto j)
            {
                if (i != j)
                    matrix_[i][j].usmv(alpha, x[j], y[i]);
            });
        });
    }

    /*! \brief Category of the linear operator
     *
     *  While each component may be of a different category,
     *  overlapping is selected in parallel for the overall
     *  operator because no adequate value exists. Has to be
     *  consistent with the categories for the scalar product
     *  and the preconditioner.
     */
    Dune::SolverCategory::Category category() const
    {
        if (std::get<0>(linearOperators_)->category() == Dune::SolverCategory::sequential)
            return Dune::SolverCategory::sequential;

        return Dune::SolverCategory::overlapping;
    }

private:
    const LinearOperatorTuple& linearOperators_;
    const Matrix& matrix_;
};

} // end namespace Dumux

#endif
