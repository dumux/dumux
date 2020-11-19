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
 * \brief Scalar products for iterative solvers
 */
#ifndef DUMUX_LINEAR_SCALARPRODUCTS_HH
#define DUMUX_LINEAR_SCALARPRODUCTS_HH

#include <tuple>

#include <dune/common/hybridutilities.hh>
#include <dune/istl/scalarproducts.hh>
#include <dune/istl/operators.hh>

namespace Dumux {

/*!
 * \ingroup Linear
 * \brief A scalar product based on a tuple of scalar products
 *
 * \tparam Vector usually a MultiTypeBlockVector
 * \tparam ScalarProductTuple tuple of scalar products
 *
 * Vector and ScalarProductTuple have to be consistent in the sense
 * that a tuple element implements the scalar product for a corresponding
 * subtype of Vector.
 */
template<class Vector, class ScalarProductTuple>
class TupleScalarProduct : public Dune::ScalarProduct<Vector>
{
    using field_type = typename Dune::ScalarProduct<Vector>::field_type;
    using real_type = typename Dune::ScalarProduct<Vector>::real_type;

public:
    /*! \brief Construct the scalar product
     *
     *  \param sps tuple of scalar products
     */
    TupleScalarProduct (const ScalarProductTuple& sps)
    : scalarProducts_(sps)
    {}

    /*! \brief Dot product of two vectors
     *
     *  Each vector component must be consistent on the interior+border partition.
     */
    field_type dot (const Vector& x, const Vector& y) const override
    {
        field_type result(0);

        using namespace Dune::Hybrid;
        forEach(integralRange(Dune::Hybrid::size(x)), [&](const auto i)
        {
            result += std::get<i>(scalarProducts_)->dot(x[i], y[i]);
        });

        return result;
    }

    /*! \brief Norm of a right-hand side vector
     *
     *  Each vector component must be consistent on the interior+border partition.
     */
    real_type norm (const Vector& x) const override
    {
        using std::sqrt;
        return sqrt(dot(x, x));
    }

    /*! \brief Category of the scalar product
     *
     *  While each component may be of a different category,
     *  overlapping is selected in parallel for the overall
     *  product because no adequate value exists. Has to be
     *  consistent with the categories for the linear operator
     *  and the preconditioner.
     */
    Dune::SolverCategory::Category category() const
    {
        if (std::get<0>(scalarProducts_)->category() == Dune::SolverCategory::sequential)
            return Dune::SolverCategory::sequential;

        return Dune::SolverCategory::overlapping;
    }

private:
    const ScalarProductTuple& scalarProducts_;
};

} // end namespace Dumux

#endif
