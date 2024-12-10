// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup MultiDomain
 * \ingroup MortarCoupling
 * \brief Additional property declarations for mortar-coupling models.
 */
#ifndef DUMUX_MULTIDOMAIN_MORTAR_PROPERTIES_HH
#define DUMUX_MULTIDOMAIN_MORTAR_PROPERTIES_HH

#include <dumux/common/properties.hh>

namespace Dumux::Properties {

DUMUX_DEFINE_PROPERTY(MortarGrid);
DUMUX_DEFINE_PROPERTY(MortarSolutionVector);

}  // namespace Dumux::Properties

#endif
