// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup OneEqModel
 * \copydoc Dumux::OneEqResidual
 */
#ifndef DUMUX_ONEEQ_LOCAL_RESIDUAL_HH
#define DUMUX_ONEEQ_LOCAL_RESIDUAL_HH

#include <dumux/common/properties.hh>
#include <dumux/discretization/method.hh>
#include <dumux/freeflow/navierstokes/localresidual.hh>
#include <dumux/freeflow/rans/oneeq/staggered/localresidual.hh>

namespace Dumux {

// forward declaration
template<class TypeTag, class BaseLocalResidual, class DiscretizationMethod>
class OneEqResidualImpl;

/*!
 * \ingroup OneEqModel
 * \brief The local residual class for the one-equation turbulence model by Spalart-Allmaras
          This is a convenience alias for the actual,
          discretization-specific local residual.
 * \note  Not all specializations are currently implemented
 */
template<class TypeTag, class BaseLocalResidual>
using OneEqResidual = OneEqResidualImpl<TypeTag, BaseLocalResidual, typename GetPropType<TypeTag, Properties::GridGeometry>::DiscretizationMethod>;

} // end namespace Dumux

#endif
