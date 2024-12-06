// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup OnePTests
 * \ingroup MultiDomainTests
 * \brief Additional property declarations for the single-phase darcy-darcy mortar-coupling test
 */
#ifndef DUMUX_MORTAR_DARCY_ONEP_TEST_PROPERTY_DECLARATIONS_HH
#define DUMUX_MORTAR_DARCY_ONEP_TEST_PROPERTY_DECLARATIONS_HH

#include <dumux/common/properties.hh>

namespace Dumux::Properties {

DUMUX_DEFINE_PROPERTY(MortarGrid);
DUMUX_DEFINE_PROPERTY(MortarSolutionVector);

}  // namespace Dumux::Properties

#endif
