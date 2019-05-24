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
 * \ingroup Common
 * \brief Helpers for deprecation
 */

#ifndef DUMUX_COMMON_DEPRECATED_HH
#define DUMUX_COMMON_DEPRECATED_HH

#include <type_traits>

#include <dune/common/deprecated.hh>

#include <dumux/common/typetraits/isvalid.hh>

namespace Dumux {

#ifndef DOXYGEN // hide from doxygen
// Helper classes/functions for deprecation
// Each implementation has to state after which release
// it will be removed. Implementations in the Deprecated
// namespace will be removed without
// deprecation after their usage in the code exprired,
// so most likely you don't want to use this in your code
namespace Deprecated {

////////////////////////////////////////////////////////
///// REMOVE THIS AFTER RELEASE 3.1
////////////////////////////////////////////////////////


// support old interface of the effective thermal conductivity laws
template<class VV>
struct HasNewEffThermCondIF
{
    template<class ETC>
    auto operator()(ETC&& e) -> decltype(e.effectiveThermalConductivity(std::declval<const VV&>())) {}
};

template<class ETC, class VV, class SpatialParams, class Element, class FVGeometry,
         typename std::enable_if_t<!decltype(isValid(HasNewEffThermCondIF<VV>()).template check<ETC>())::value, int> = 0>
auto effectiveThermalConductivity(const VV& volVars,
                                  const SpatialParams& spatialParams,
                                  const Element& element,
                                  const FVGeometry& fvGeometry,
                                  const typename FVGeometry::SubControlVolume& scv)
{
    return ETC::effectiveThermalConductivity(volVars, spatialParams, element, fvGeometry, scv);
}

template<class ETC, class VV, class SpatialParams, class Element, class FVGeometry,
         typename std::enable_if_t<decltype(isValid(HasNewEffThermCondIF<VV>()).template check<ETC>())::value, int> = 0>
auto effectiveThermalConductivity(const VV& volVars,
                                  const SpatialParams& spatialParams,
                                  const Element& element,
                                  const FVGeometry& fvGeometry,
                                  const typename FVGeometry::SubControlVolume& scv)
{
    return ETC::effectiveThermalConductivity(volVars);
}
////////////////////////////////////////////////////////
////////////////////////////////////////////////////////

} // end namespace Deprecated
#endif

} // end namespace Dumux

#endif
