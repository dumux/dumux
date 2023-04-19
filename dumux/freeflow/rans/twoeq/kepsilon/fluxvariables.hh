// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup KEpsilonModel
 * \copydoc Dumux::KEpsilonFluxVariables
 */
#ifndef DUMUX_KEPSILON_FLUXVARIABLES_HH
#define DUMUX_KEPSILON_FLUXVARIABLES_HH

#include <dumux/common/properties.hh>
#include <dumux/freeflow/rans/twoeq/kepsilon/staggered/fluxvariables.hh>

namespace Dumux {

// forward declaration
template<class TypeTag, class BaseFluxVariables, class DiscretizationMethod>
class KEpsilonFluxVariablesImpl;

/*!
 * \ingroup KEpsilonModel
 * \brief The flux variables class for the k-epsilon model.
          This is a convenience alias for that actual,
          discretization-specific flux variables.
 * \note  Not all specializations are currently implemented
 */
template<class TypeTag, class BaseFluxVariables>
using KEpsilonFluxVariables = KEpsilonFluxVariablesImpl<TypeTag, BaseFluxVariables, typename GetPropType<TypeTag, Properties::GridGeometry>::DiscretizationMethod>;

} // end namespace

#endif
