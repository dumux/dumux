// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup CCMpfaDiscretization
 * \brief Stores the face indices corresponding to the neighbors of an element
 *        that contribute to the derivative calculation. Depending on if an
 *        mpfa scheme leads to a symmetric/unsymmetric sparsity pattern, the
 *        adequate implementation of the connectivity map is chosen.
 */
#ifndef DUMUX_CC_MPFA_CONNECTIVITY_MAP_HH
#define DUMUX_CC_MPFA_CONNECTIVITY_MAP_HH

#include <dumux/discretization/cellcentered/mpfa/methods.hh>
#include <dumux/discretization/cellcentered/connectivitymap.hh>

namespace Dumux {

//! Forward declaration of method specific implementation of the assembly map
template<class GridGeometry, MpfaMethods method>
class CCMpfaConnectivityMap;

//! The o-method can use the simple (symmetric) assembly map
template<class GridGeometry>
class CCMpfaConnectivityMap<GridGeometry, MpfaMethods::oMethod> : public CCSimpleConnectivityMap<GridGeometry> {};
} // end namespace Dumux

#endif
