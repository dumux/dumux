// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
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
