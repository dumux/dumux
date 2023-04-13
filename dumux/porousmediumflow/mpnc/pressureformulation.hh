// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup MPNCModel
 * \brief Enumeration of the formulations accepted by the MpNc model.
 */

#ifndef DUMUX_MPNC_PRESSUREFORMULATION_HH
#define DUMUX_MPNC_PRESSUREFORMULATION_HH

namespace Dumux
{

/*!
 * \ingroup MPNCModel
 * \brief Enumerates the formulations which the MpNc model accepts.
 */
enum class MpNcPressureFormulation
{
    mostWettingFirst, leastWettingFirst
};

}

#endif
