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
#ifndef DUMUX_LINEAR_SCALARPRODUCTS_HH
#define DUMUX_LINEAR_SCALARPRODUCTS_HH

#include <tuple>

#include <dune/common/hybridutilities.hh>
#include <dune/istl/scalarproducts.hh>
#include <dune/istl/operators.hh>

namespace Dumux {

template<class Vector, class ScalarProductTuple>
class TupleScalarProduct : public Dune::ScalarProduct<Vector>
{
    using field_type = typename Dune::ScalarProduct<Vector>::field_type;
    using real_type = typename Dune::ScalarProduct<Vector>::real_type;

public:
    /*!
     * \param comm The communication object for syncing overlap and copy
     * data points.
     * \param cat parallel solver category (nonoverlapping or overlapping)
     */
    TupleScalarProduct (const ScalarProductTuple& sps)
    : sps_(sps)
    {}

    /*! \brief Dot product of two vectors.
     *       It is assumed that the vectors are consistent on the interior+border
     *       partition.
     */
    field_type dot (const Vector& x, const Vector& y) const override
    {
        field_type result(0);

        using namespace Dune::Hybrid;
        forEach(integralRange(Dune::Hybrid::size(x)), [&](const auto i)
        {
            result += std::get<i>(sps_)->dot(x[i], y[i]);
        });

        return result;
    }

    /*! \brief Norm of a right-hand side vector.
     *       The vector must be consistent on the interior+border partition
     */
    real_type norm (const Vector& x) const override
    {
        using std::sqrt;
        return sqrt(dot(x, x));
    }

    //! Category of the scalar product (see SolverCategory::Category)
    Dune::SolverCategory::Category category() const
    {
        return Dune::SolverCategory::overlapping;
    }

private:
    const ScalarProductTuple& sps_;
};

} // end namespace Dumux

#endif
