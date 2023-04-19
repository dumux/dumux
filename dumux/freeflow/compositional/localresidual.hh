// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
 /*!
  * \file
  * \ingroup FreeflowNCModel
  * \copydoc Dumux::FreeflowNCResidual
  */
#ifndef DUMUX_FREEFLOW_NC_LOCAL_RESIDUAL_HH
#define DUMUX_FREEFLOW_NC_LOCAL_RESIDUAL_HH

#include <dumux/common/properties.hh>
#include <dumux/discretization/method.hh>
#include <dumux/freeflow/navierstokes/localresidual.hh>
#include <dumux/freeflow/compositional/staggered/localresidual.hh>

namespace Dumux {

// forward declaration
template<class TypeTag, class DiscretizationMethod>
class FreeflowNCResidualImpl;

/*!
 * \ingroup FreeflowNCModel
 * \brief The local residual class for the multi-component free-flow model (balance equations).
          This is a convenience alias for the actual,
          discretization-specific local residual.
 * \note  Not all specializations are currently implemented
 */
template<class TypeTag>
using FreeflowNCResidual = FreeflowNCResidualImpl<TypeTag, typename GetPropType<TypeTag, Properties::GridGeometry>::DiscretizationMethod>;

}

#endif
