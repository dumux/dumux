// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup Assembly
 * \brief An enum class to define various differentiation methods available in order to compute
          the derivatives of the residual i.e. the entries in the jacobian matrix.
 */
#ifndef DUMUX_JACOBIAN_DIFFERENTIATION_METHODS_HH
#define DUMUX_JACOBIAN_DIFFERENTIATION_METHODS_HH

namespace Dumux {

/*!
 * \ingroup Assembly
 * \brief Differentiation methods in order to compute the derivatives
 *        of the residual i.e. the entries in the jacobian matrix.
 * \todo automatic differentation is not yet implemented
 */
enum class DiffMethod
{
    numeric, analytic, automatic
};

} // end namespace Dumux

#endif
