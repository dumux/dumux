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
 * \ingroup Assembly
 * \copydoc Dumux::OperatorWeights.
 */
#ifndef DUMUX_ASSEMBLY_OPERATOR_WEIGHTS_HH
#define DUMUX_ASSEMBLY_OPERATOR_WEIGHTS_HH

#include <concepts>
#include <optional>

namespace Dumux {

/*!
 * \file
 * \ingroup Assembly
 * \brief Small class to store weights associated with spatial and temporal
 *        operators of a time-dependent PDE. This is used e.g. in the context
 *        of multi-stage time integration methods.
 */
template<std::floating_point Scalar>
struct OperatorWeights
{
    std::optional<Scalar> spatialWeight;
    std::optional<Scalar> temporalWeight;
};

} // end namespace Dumux

#endif
