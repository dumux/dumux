// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup EmbeddedCoupling
 * \brief Coupling manager for embedded fractures
 */

#ifndef DUMUX_MULTIDOMAIN_EMBEDDED_COUPLINGMANAGER_2D3D_HH
#define DUMUX_MULTIDOMAIN_EMBEDDED_COUPLINGMANAGER_2D3D_HH

#include <dumux/multidomain/embedded/couplingmanagerbase.hh>

namespace Dumux {

/*!
 * \ingroup EmbeddedCoupling
 * \brief Coupling manager for embedded fractures
 * \note we just use the default coupling manager
 */
template<class MDTraits>
class EmbeddedCouplingManager2d3d
: public EmbeddedCouplingManagerBase<MDTraits,
                                     EmbeddedCouplingManager2d3d<MDTraits>>
{
    using ThisType = EmbeddedCouplingManager2d3d<MDTraits>;
    using ParentType = EmbeddedCouplingManagerBase<MDTraits, ThisType>;
public:
    using ParentType::ParentType;
};

//! we support multithreaded assembly
template<class MDTraits>
struct CouplingManagerSupportsMultithreadedAssembly<EmbeddedCouplingManager2d3d<MDTraits>>
: public std::true_type {};

} // end namespace Dumux

#endif
