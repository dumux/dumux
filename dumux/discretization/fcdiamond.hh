// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup Discretization
 * \brief Defines a type tag and some properties for models using the diamond scheme.
          This scheme features degrees of freedom at the elements' centers and intersections (faces).
 */

#ifndef DUMUX_DISCRETIZATION_FACECENTERED_DIAMOND_HH
#define DUMUX_DISCRETIZATION_FACECENTERED_DIAMOND_HH

#include <dumux/discretization/pq1nonconforming.hh>

namespace Dumux::Properties::TTag {

//! Alias for backwards compatibility
using FaceCenteredDiamondModel = PQ1NonconformingFVModel;

} // end namespace Dumux::Properties::TTag

#endif
