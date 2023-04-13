// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup MultiDomain
 * \brief Some useful type traits
 */
#ifndef DUMUX_MULTIDOMAIN_FREEFLOW_TYPETRAITS_HH
#define DUMUX_MULTIDOMAIN_FREEFLOW_TYPETRAITS_HH

namespace Dumux::Detail {

template<class Traits>
struct MomentumDiscretizationMethod
{ using type = typename Traits::template SubDomain<0>::GridGeometry::DiscretizationMethod; };

} // end namespace Dumux::Detail

#endif
