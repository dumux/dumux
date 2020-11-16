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
 * \brief Dumux linear operators for iterative solvers
 */
#ifndef DUMUX_LINEAR_OPERATORS_HH
#define DUMUX_LINEAR_OPERATORS_HH

#include <tuple>

#include <dune/common/hybridutilities.hh>
#include <dune/istl/solvercategory.hh>
#include <dune/istl/operators.hh>

namespace Dumux {

template<class Vector, class Matrix, class LinearOperatorTuple>
class TupleLinearOperator: public Dune::LinearOperator<Vector, Vector> {
public:
    //! The type of the domain of the operator.
    typedef Vector domain_type;
    //! The type of the range of the operator.
    typedef Vector range_type;
    //! The field type of the operator.
    typedef typename Vector::field_type field_type;

    TupleLinearOperator (const LinearOperatorTuple& lops, const Matrix& m)
    : lops_(lops), m_(m)
    {}

    /*! \brief apply operator to x:  \f$ y = A(x) \f$
     *          The input vector is consistent and the output must also be
     *       consistent on the interior+border partition.
     */
    void apply (const Vector& x, Vector& y) const override
    {
        using namespace Dune::Hybrid;
        forEach(integralRange(Dune::Hybrid::size(x)), [&](const auto i)
        {
            std::get<i>(lops_)->apply(x[i], y[i]);

            forEach(integralRange(Dune::Hybrid::size(x)), [&](const auto j)
            {
                if (i != j)
                    m_[i][j].umv(x[j], y[i]);
            });
        });
    }

    //! apply operator to x, scale and add:  \f$ y = y + \alpha A(x) \f$
    void applyscaleadd (field_type alpha, const Vector& x, Vector& y) const override
    {
        using namespace Dune::Hybrid;
        forEach(integralRange(Dune::Hybrid::size(x)), [&](const auto i)
        {
            std::get<i>(lops_)->applyscaleadd(alpha, x[i], y[i]);

            forEach(integralRange(Dune::Hybrid::size(x)), [&](const auto j)
            {
                if (i != j)
                    m_[i][j].usmv(alpha, x[j], y[i]);
            });
        });
    }

    //! Category of the linear operator (see SolverCategory::Category)
    Dune::SolverCategory::Category category() const
    {
        return Dune::SolverCategory::overlapping;
    }

private:
    const LinearOperatorTuple& lops_;
    const Matrix& m_;
};

} // end namespace Dumux

#endif
