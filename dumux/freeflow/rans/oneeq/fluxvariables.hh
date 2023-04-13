// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup OneEqModel
 * \copydoc Dumux::OneEqFluxVariables
 */
#ifndef DUMUX_ONEEQFLUXVARIABLES_HH
#define DUMUX_ONEEQFLUXVARIABLES_HH

#include <dumux/common/properties.hh>
#include <dumux/freeflow/navierstokes/fluxvariables.hh>
#include <dumux/freeflow/rans/oneeq/staggered/fluxvariables.hh>

namespace Dumux {

// forward declaration
template<class TypeTag, class BaseFluxVariables, class DiscretizationMethod>
class OneEqFluxVariablesImpl;

/*!
 * \ingroup OneEqModel
 * \brief The flux variables class for the one-equation turbulence model by Spalart-Allmaras.
          This is a convenience alias for that actual,
          discretization-specific flux variables.
 * \note  Not all specializations are currently implemented
 */
template<class TypeTag, class BaseFluxVariables>
using OneEqFluxVariables = OneEqFluxVariablesImpl<TypeTag, BaseFluxVariables, typename GetPropType<TypeTag, Properties::GridGeometry>::DiscretizationMethod>;

} // end namespace

#endif
