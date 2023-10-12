// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup Core
 * \brief Traits class to set options used by the local residual when
 *        when evaluating the balance equations.
 */
#ifndef BALANCE_EQUATION_OPTIONS_HH
#define BALANCE_EQUATION_OPTIONS_HH

namespace Dumux {

/*!
 * \ingroup Core
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
