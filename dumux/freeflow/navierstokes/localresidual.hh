// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup NavierStokesModel
 * \copydoc Dumux::NavierStokesResidual
 */
#ifndef DUMUX_NAVIERSTOKES_LOCAL_RESIDUAL_HH
#define DUMUX_NAVIERSTOKES_LOCAL_RESIDUAL_HH

#include <dumux/common/properties.hh>
#include <dumux/discretization/method.hh>
#include <dumux/freeflow/navierstokes/staggered/localresidual.hh>

namespace Dumux {

// forward declaration
template<class TypeTag, class DiscretizationMethod, class Implementation>
class NavierStokesResidualImpl;

/*!
 * \ingroup NavierStokesModel
 * \brief The local residual class for the Navier-Stokes model (balance equations).
          This is a convenience alias for the actual,
          discretization-specific local residual.
 * \note Not all specializations are currently implemented
 * \tparam TypeTag The model type tag
 * \tparam Implementation optional implementation using CRTP (use void if NavierStokesResidual is itself the implementation)
 */
template<class TypeTag, class Implementation = void>
using NavierStokesResidual = NavierStokesResidualImpl<
    TypeTag,
    typename GetPropType<TypeTag, Properties::GridGeometry>::DiscretizationMethod,
    Implementation
>;

} // end namespace Dumux

#endif
