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

#include <dune/common/version.hh>
#include <dune/common/exceptions.hh>
#include <dune/common/std/type_traits.hh>

namespace Dumux {

#ifndef DOXYGEN // hide from doxygen
// Helper classes/functions for deprecation
// Each implementation has to state after which release
// it will be removed. Implementations in the Deprecated
// namespace will be removed without
// deprecation after their usage in the code expired,
// so most likely you don't want to use this in your code
namespace Deprecated {

template <class Mapper, class GridView>
using GridViewDetector = decltype(std::declval<Mapper>().update(std::declval<GridView>()));

template<class Mapper, class GridView>
static constexpr bool hasUpdateGridView()
{ return Dune::Std::is_detected<GridViewDetector, Mapper, GridView>::value; }

// helper class to print deprecated message
template <class Mapper>
#if DUNE_VERSION_GTE(DUNE_GRID,2,8)
[[deprecated("The interface mapper.update() is deprecated. All mappers now have to implement `update(gridView)` instead (with a gridView as argument). Only mappers with the new interface will be support for dune-grid 2.7 is dropped.")]]
#endif
void update(Mapper& mapper)
{ mapper.update(); };

template <typename Problem, typename GlobalPosition>
using HasIsOnWallDetector = decltype(std::declval<Problem>().isOnWallAtPos(std::declval<GlobalPosition>()));

template<class Problem, typename GlobalPosition>
static constexpr bool hasIsOnWall()
{ return Dune::Std::is_detected<HasIsOnWallDetector, Problem, GlobalPosition>::value; }

template <typename BcTypes>
using HasWallBCDetector = decltype(std::declval<BcTypes>().hasWall());

template<class BcTypes>
static constexpr bool hasHasWallBC()
{ return Dune::Std::is_detected<HasWallBCDetector, BcTypes>::value; }

template <typename ModelTraits>
using HasEnableCompositionalDispersionDetector = decltype(ModelTraits::enableCompositionalDispersion());

template<class ModelTraits>
static constexpr bool hasEnableCompositionalDispersion()
{ return Dune::Std::is_detected<HasEnableCompositionalDispersionDetector, ModelTraits>::value; }

template <typename ModelTraits>
using HasEnableThermalDispersionDetector = decltype(ModelTraits::enableThermalDispersion());

template<class ModelTraits>
static constexpr bool hasEnableThermalDispersion()
{ return Dune::Std::is_detected<HasEnableThermalDispersionDetector, ModelTraits>::value; }

template<class SpatialParams, class Element, class Scv, class ElemSol>
using HasTemperatureDetector = decltype(std::declval<SpatialParams>().temperature(
    std::declval<Element>(),
    std::declval<Scv>(),
    std::declval<ElemSol>()
));

template<typename Problem, typename Element, typename Scv, typename ElemSol>
decltype(auto) temperature(const Problem& problem, const Element& element, const Scv& scv, const ElemSol& elemSol)
{
    using SpatialParams = std::decay_t<decltype(problem.spatialParams())>;
    if constexpr (Dune::Std::is_detected<HasTemperatureDetector, SpatialParams, Element, Scv, ElemSol>::value)
        return problem.spatialParams().temperature(element, scv, elemSol);
    else
        return problem.temperatureAtPos(scv.dofPosition());
}

template<class SpatialParams, class Element, class SubControlVolume, class ElementSolution>
using HasExtrusionFactorDetector = decltype(std::declval<SpatialParams>().extrusionFactor(
        std::declval<Element>(),
        std::declval<SubControlVolume>(),
        std::declval<ElementSolution>()
    ));

template<class Problem, class Element, class SubControlVolume, class ElementSolution>
using HasBaseProblemExtrusionFactorDetector = decltype(std::declval<Problem>().extrusionFactor(
        std::declval<Element>(),
        std::declval<SubControlVolume>(),
        std::declval<ElementSolution>(),
        double{}
    ));

template<class Problem, class GlobalPosition>
using HasBaseProblemExtrusionFactorAtPosDetector = decltype(std::declval<Problem>().extrusionFactorAtPos(
        std::declval<GlobalPosition>(),
        double{}
    ));

template<class Problem, class Element, class SubControlVolume, class ElementSolution>
decltype(auto) extrusionFactor(const Problem& problem,
                               const Element& element,
                               const SubControlVolume& scv,
                               const ElementSolution& elemSol)
{
    using SpatialParams = std::decay_t<decltype(problem.spatialParams())>;
    using GlobalPosition = std::decay_t<decltype(scv.center())>;

    static constexpr bool hasNewSpatialParamsInterface = Dune::Std::is_detected<
        HasExtrusionFactorDetector, SpatialParams, Element, SubControlVolume, ElementSolution
    >::value;

    static constexpr bool hasBaseProblemInterface = Dune::Std::is_detected<
        HasBaseProblemExtrusionFactorDetector, Problem, Element, SubControlVolume, ElementSolution
    >::value;

    static constexpr bool hasBaseProblemAtPosInterface = Dune::Std::is_detected<
        HasBaseProblemExtrusionFactorAtPosDetector, Problem, GlobalPosition
    >::value;

    static constexpr bool hasUserDefinedProblemExtrusionFactor = !hasBaseProblemInterface || !hasBaseProblemAtPosInterface;

    if constexpr (hasNewSpatialParamsInterface && hasUserDefinedProblemExtrusionFactor)
        DUNE_THROW(Dune::InvalidStateException,
                   "Extrusion factor defined both in problem implementation (deprecated interface) and spatial params (new interface). "
                   "Please move the overload in your problem implementation to your spatial parameters.");

    if constexpr (hasNewSpatialParamsInterface)
        return problem.spatialParams().extrusionFactor(element, scv, elemSol);
    else
        return problem.extrusionFactor(element, scv, elemSol);
}

} // end namespace Deprecated
#endif

} // end namespace Dumux
#endif
