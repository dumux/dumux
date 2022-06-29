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
 * \ingroup SpatialParameters
 * \brief The base class for spatial parameters of one-phase problems
 * using a fully implicit discretization method.
 */
#ifndef DUMUX_FV_SPATIAL_PARAMS_ONE_P_HH
#define DUMUX_FV_SPATIAL_PARAMS_ONE_P_HH

#warning "This file is deprecated, use FVPorousMediumFlowSpatialParamsOneP from dumux/porousmediumflow/fvspatialparams1p.hh instead!"

#include <dune/common/exceptions.hh>
#include <dune/common/fmatrix.hh>

#include <dumux/common/parameters.hh>
#include <dumux/common/math.hh>
#include <dumux/common/typetraits/isvalid.hh>

namespace Dumux {

#ifndef DOXYGEN
namespace Detail {
// helper struct detecting if the user-defined spatial params class has a permeabilityAtPos function
template<class GlobalPosition>
struct hasPermeabilityAtPos
{
    template<class SpatialParams>
    auto operator()(const SpatialParams& a)
    -> decltype(a.permeabilityAtPos(std::declval<GlobalPosition>()))
    {}
};

template<class GlobalPosition, class SolidSystem>
struct hasInertVolumeFractionAtPos
{
    template<class SpatialParams>
    auto operator()(const SpatialParams& a)
    -> decltype(a.template inertVolumeFractionAtPos<SolidSystem>(std::declval<GlobalPosition>(), 0))
    {}
};

template<class GlobalPosition>
struct hasPorosityAtPos
{
    template<class SpatialParams>
    auto operator()(const SpatialParams& a)
    -> decltype(a.porosityAtPos(std::declval<GlobalPosition>()))
    {}
};
} // end namespace Detail
#endif

} // namespace Dumux

#endif
