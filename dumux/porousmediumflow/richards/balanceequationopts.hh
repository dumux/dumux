// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup RichardsModel
 * \brief Traits class to set options used by the local residual when
 *        when evaluating the balance equations.
 */
#ifndef DUMUX_POROUSMEDIUMFLOW_RICHARDS_BALANCE_EQUATION_OPTIONS_HH
#define DUMUX_POROUSMEDIUMFLOW_RICHARDS_BALANCE_EQUATION_OPTIONS_HH

namespace Dumux {

/*!
 * \ingroup RichardsModel
 * \brief Traits class to set options used by the local residual when
 *        when evaluating the balance equations.
 */
template <class FluidSystem>
struct RichardsBalanceEquationOptions
{
    /*
     * The main component in the liquid phase (phase 0) is always balanced
     * and the main component in the gas phase (phase 1) is never balanced
     * For the index convention see RichardsVolumeVariables.
     */
    static constexpr bool mainComponentIsBalanced(int phaseIdx)
    { return phaseIdx == FluidSystem::phase0Idx; }
};

} // end namespace Dumux

#endif
