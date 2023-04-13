// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup NIModel
 * \brief Defines the indices used by the non-isothermal two-phase two-component model.
 */

#ifndef DUMUX_ENERGY_INDICES_HH
#define DUMUX_ENERGY_INDICES_HH

namespace Dumux {

/*!
 * \ingroup NIModel
 * \brief Indices for the non-isothermal two-phase two-component model.
 *
 * \tparam IsothermalIndices The indices of the isothermal model
 * \tparam numEq the number of equations of the non-isothermal model
 */
template <class IsothermalIndices, int numEq>
struct EnergyIndices : public IsothermalIndices
{
    static const int temperatureIdx = numEq - 1; //!< The index for temperature in primary variable vectors.
    static const int energyEqIdx = numEq - 1;    //!< The index for energy in equation vectors.
};

} // end namespace Dumux

#endif
