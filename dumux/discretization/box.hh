// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup Discretization
 * \brief Defines a type tag and some properties for models using the box scheme.
 */

#ifndef DUMUX_DISCRETIZTAION_BOX_HH
#define DUMUX_DISCRETIZTAION_BOX_HH

#include <dumux/discretization/pq1.hh>

namespace Dumux::Properties::TTag {

//! Alias for Box model
using BoxModel = PQ1FVModel;

} // end namespace Dumux::Properties::TTag

namespace Dumux::Detail {

template<class TypeTag>
concept BoxModel = PQ1FVModel<TypeTag>;

} // end namespace Dumux::Detail

#endif
