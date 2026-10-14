// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup Flux
 * \brief The density in the gravity term of Darcy's law
 */
#ifndef DUMUX_FLUX_GRAVITY_TERM_DENSITY_HH
#define DUMUX_FLUX_GRAVITY_TERM_DENSITY_HH

namespace Dumux {

/*!
 * \ingroup Flux
 * \brief The phase density is the density in the gravity term of Darcy's law
 */
struct PhaseDensityInGravityTerm
{
    template<class Problem, class Variables>
    static auto density(const Problem&, const Variables& vars, int phaseIdx)
    { return vars.density(phaseIdx); }
};

} // end namespace Dumux

#endif
