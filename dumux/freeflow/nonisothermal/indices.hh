// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup FreeflowNIModel
 * \copydoc Dumux::FreeflowNonIsothermalIndices
 */
#ifndef DUMUX_FREEFLOW_NI_INDICES_HH
#define DUMUX_FREEFLOW_NI_INDICES_HH

namespace Dumux {

/*!
 * \ingroup FreeflowNIModel
 * \brief Indices for the non-isothermal Navier-Stokes model.
 *
 * \tparam IsothermalIndices The isothermal indices class
 * \tparam numEq the number of equations of the non-isothermal model
 */
template <class IsothermalIndices, int numEq>
class FreeflowNonIsothermalIndices : public IsothermalIndices
{
public:
    static constexpr int energyEqIdx = numEq - 1;
    static constexpr int temperatureIdx = numEq - 1;
};

} // end namespace Dumux

#endif
