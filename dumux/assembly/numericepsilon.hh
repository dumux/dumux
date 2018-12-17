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
 * \ingroup Assembly
 * \brief An adapter class for local assemblers using numeric differentiation
 */
#ifndef DUMUX_ASSEMBLY_NUMERIC_EPSILON_HH
#define DUMUX_ASSEMBLY_NUMERIC_EPSILON_HH

#include <dune/common/fvector.hh>
#include <dumux/common/parameters.hh>
#include <dumux/common/numericdifferentiation.hh>

namespace Dumux {

/*!
 * \ingroup Assembly
 * \brief A helper class for local assemblers using numeric differentiation to determine the epsilon
 * \tparam Scalar the scalar type
 * \tparam numEq the number of primary variables / equations
 */
template<class Scalar, int numEq>
class NumericEpsilon
{
    using NumEqVector = Dune::FieldVector<Scalar, numEq>;

public:
    explicit NumericEpsilon(const std::string& paramGroup = "")
    {
        // get epsilons from input file with invalid default
        baseEps_ = getParamFromGroup<Scalar>(paramGroup, "Assembly.NumericDifference.BaseEpsilon", 1e-10);
        magnitude_ = getParamFromGroup<NumEqVector>(paramGroup, "Assembly.NumericDifference.PriVarMagnitude", NumEqVector(-1));
    }

    /*!
     * \brief get the epsilon
     * \note If no user input was specified -> try to estimate magnitude from primary variable value
     *       else -> use given magnitude for the primary variable times the base epsilon (default 1e-10)
     */
    Scalar operator() (Scalar priVar, int priVarIdx) const noexcept
    {
        return magnitude_[priVarIdx] > 0.0 ? baseEps_*magnitude_[priVarIdx]
                                           : NumericDifferentiation::epsilon(priVar, baseEps_);
    }

private:
    Scalar baseEps_;
    NumEqVector magnitude_;
};

} // end namespace Dumux

#endif
