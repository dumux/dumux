// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
 /*!
  * \file
  * \ingroup KEpsilonModel
  * \copydoc Dumux::KEpsilonResidual
  */
#ifndef DUMUX_KEPSILON_LOCAL_RESIDUAL_HH
#define DUMUX_KEPSILON_LOCAL_RESIDUAL_HH

#include <dumux/common/properties.hh>
#include <dumux/discretization/method.hh>
#include <dumux/freeflow/navierstokes/localresidual.hh>
#include <dumux/freeflow/rans/twoeq/kepsilon/staggered/localresidual.hh>

namespace Dumux {

// forward declaration
template<class TypeTag, class BaseLocalResidual, class DiscretizationMethod>
class KEpsilonResidualImpl;

/*!
 * \ingroup KEpsilonModel
 * \brief The local residual class for the k-epsilon model.
          This is a convenience alias for the actual,
          discretization-specific local residual.
 * \note  Not all specializations are currently implemented
 */
template<class TypeTag, class BaseLocalResidual>
using KEpsilonResidual = KEpsilonResidualImpl<TypeTag, BaseLocalResidual, typename GetPropType<TypeTag, Properties::GridGeometry>::DiscretizationMethod>;

}

#endif
