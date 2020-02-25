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

// support old interface of the neumann() function on problems
template<class E, class FVEG, class EVV, class EFVC>
class HasNewNeumannIF
{
    using SCVF = typename FVEG::SubControlVolumeFace;

public:
    template<class P>
    auto operator()(P&& p) -> decltype(p.neumann(std::declval<const E&>(),
                                                 std::declval<const FVEG&>(),
                                                 std::declval<const EVV&>(),
                                                 std::declval<const EFVC&>(),
                                                 std::declval<const SCVF&>()))
    {}
};

template<class P, class E, class FVEG, class EVV, class EFVC,
         typename std::enable_if_t<!decltype(isValid(HasNewNeumannIF<E, FVEG, EVV, EFVC>()).template check<P>())::value, int> = 0>
auto DUNE_DEPRECATED_MSG("Use new neumann() interface (see common/fvproblem.hh) that additionally receives the element flux variables cache in your problem!. Will be removed after 3.1 release")
neumann(const P& problem,
        const E& element,
        const FVEG& fvGeometry,
        const EVV& elemVolVars,
        const EFVC& elemFluxVarsCache,
        const typename FVEG::SubControlVolumeFace& scvf)
{
    return problem.neumann(element, fvGeometry, elemVolVars, scvf);
}

template<class P, class E, class FVEG, class EVV, class EFVC,
         typename std::enable_if_t<decltype(isValid(HasNewNeumannIF<E, FVEG, EVV, EFVC>()).template check<P>())::value, int> = 0>
auto neumann(const P& problem,
             const E& element,
             const FVEG& fvGeometry,
             const EVV& elemVolVars,
             const EFVC& elemFluxVarsCache,
             const typename FVEG::SubControlVolumeFace& scvf)
{
    return problem.neumann(element, fvGeometry, elemVolVars, elemFluxVarsCache, scvf);
}

constexpr auto hasEffTherCondImpl = Dumux::isValid([](auto&& v) -> decltype(v.effectiveThermalConductivity()){return 0;});
template<class VolumeVariables> constexpr bool hasEffTherCond = decltype(hasEffTherCondImpl(std::declval<VolumeVariables>())){};

constexpr auto hasEffDiffCoeffImpl = Dumux::isValid([](auto&& v) -> decltype(v.effectiveDiffusionCoefficient(0,0,0)){return 0;});
template<class VolumeVariables> constexpr bool hasEffDiffCoeff = decltype(hasEffDiffCoeffImpl(std::declval<VolumeVariables>())){};

////////////////////////////////////////////////////////
////////////////////////////////////////////////////////

} // end namespace Deprecated
#endif

} // end namespace Dumux

#endif
