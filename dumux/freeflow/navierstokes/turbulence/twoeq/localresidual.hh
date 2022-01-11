// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   See the file COPYING for full copying permissions.                      *
 *                                                                           *
 *   This program is free software: you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation, either version 3 of the License, or       *
 *   (at your option) any later version.                                     *
 *                                                                           *
 *   This program is distributed in the hope that it will be useful,         *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of          *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the            *
 *   GNU General Public License for more details.                            *
 *                                                                           *
 *   You should have received a copy of the GNU General Public License       *
 *   along with this program.  If not, see <http://www.gnu.org/licenses/>.   *
 *****************************************************************************/
/*!
 * \file
 * \ingroup TwoEqModel
 * \copydoc Dumux::TwoEqResidual
 */
#ifndef DUMUX_TWOEQ_LOCAL_RESIDUAL_HH
#define DUMUX_TWOEQ_LOCAL_RESIDUAL_HH

#include <dumux/common/properties.hh>
#include <dumux/discretization/method.hh>
#include <dumux/freeflow/navierstokes/turbulence/twoeq/cellcentered/localresidual.hh>

namespace Dumux {

// forward declaration
template<class TypeTag, class BaseLocalResidual, class DiscretizationMethod>
class TwoEqResidualImpl;

/*!
 * \ingroup TwoEqModel
 * \brief The local residual class for the two-eq turbulence models
          This is a convenience alias for the actual,
          discretization-specific local residual.
 * \note  At the moment, only the Cell Centered version is currently implemented
 */
template<class TypeTag, class BaseLocalResidual>
using TwoEqResidual = TwoEqResidualImpl<TypeTag, BaseLocalResidual, typename GetPropType<TypeTag, Properties::GridGeometry>::DiscretizationMethod>;

}

#endif
