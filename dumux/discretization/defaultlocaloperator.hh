// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup Discretization
 * \brief The default local operator than can be specialized for each discretization scheme
 */

#ifndef DUMUX_DISCRETIZATION_DEFAULT_LOCAL_OPERTAOR_HH
#define DUMUX_DISCRETIZATION_DEFAULT_LOCAL_OPERTAOR_HH

namespace Dumux::Detail {

template<class TypeTag>
struct DiscretizationDefaultLocalOperator;

} // end namespace Dumux::Detail

namespace Dumux {

template<class TypeTag>
using DiscretizationDefaultLocalOperator
    = typename Detail::DiscretizationDefaultLocalOperator<TypeTag>::type;

} // end namespace Dumux

#endif
