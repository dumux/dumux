// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup AllenCahnModel
 * \brief Variables for the Allen-Cahn model.
 */
#ifndef DUMUX_PHASEFIELD_ALLENCAHN_VARIABLES_HH
#define DUMUX_PHASEFIELD_ALLENCAHN_VARIABLES_HH

#include <dumux/phasefield/common/variables.hh>

namespace Dumux {

/*!
 * \ingroup AllenCahnModel
 * \brief Variables for the Allen-Cahn model. Allen-Cahn has a single primary
 *        variable (the phase field itself), so this is simply the shared
 *        `Dumux::PhaseField::ScalarVariables` for the model's own `Traits`.
 */
template<class Traits>
using AllenCahnVariables = PhaseField::ScalarVariables<Traits>;

} // end namespace Dumux

#endif
