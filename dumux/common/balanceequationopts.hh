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
 * \ingroup Common
 * \brief Traits class to set options used by the local residual when
 *        when evaluating the balance equations.
 */
#ifndef BALANCE_EQUATION_OPTIONS_HH
#define BALANCE_EQUATION_OPTIONS_HH

namespace Dumux {

/*!
 * \ingroup Common
 * \brief Traits class to set options used by the local residual when
 *        when evaluating the balance equations.
 * \todo include useMoles here
 * \todo include replaceCompIdx here
 */
template <class TypeTag>
class BalanceEquationOptions
{
public:

    //! If a certain component is balanced in this model
    // per default all phases are balanced. See e.g. Richards for an example where
    // the air component exists but is not balanced. Or the tracer model where the
    // carrier phase main component exists but is not balanced.
    static constexpr bool mainComponentIsBalanced(int phaseIdx)
    { return true; }
};

} // end namespace Dumux

#endif
