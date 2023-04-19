// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
 /*!
  * \file
  * \ingroup FreeflowNCModel
  * \copydoc Dumux::FreeflowNCFluxVariables
  */
#ifndef DUMUX_FREELOW_NC_FLUXVARIABLES_HH
#define DUMUX_FREELOW_NC_FLUXVARIABLES_HH

#include <dumux/common/properties.hh>
#include <dumux/freeflow/compositional/staggered/fluxvariables.hh>

namespace Dumux {

// forward declaration
template<class TypeTag, class DiscretizationMethod>
class FreeflowNCFluxVariablesImpl;

/*!
 * \ingroup FreeflowNCModel
 * \brief The flux variables class for the multi-component free-flow model.
          This is a convenience alias for that actual,
          discretization-specific flux variables.
 * \note  Not all specializations are currently implemented
 */
template<class TypeTag>
using FreeflowNCFluxVariables = FreeflowNCFluxVariablesImpl<TypeTag, typename GetPropType<TypeTag, Properties::GridGeometry>::DiscretizationMethod>;


} // end namespace

#endif
