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

////////////////////////////////////////////
// Remove the following after Release 3.2 //
////////////////////////////////////////////

///////////////////////////////////////////////////////////////
// Deprecation warnings for the new material law //
///////////////////////////////////////////////////////////////

// support old interface of the effective thermal conductivity laws
template<class E, class SCV, class Sol>
struct HasNewFIAIF
{
    template<class S>
    auto operator()(S&& sp)
    -> decltype(sp.fluidMatrixInteraction(std::declval<const E&>(),
                                          std::declval<const SCV&>(),
                                          std::declval<const Sol&>())) {}
};

template<class Pos>
struct HasNewFIAIFAtPos
{
    template<class S>
    auto operator()(S&& sp)
    -> decltype(sp.fluidMatrixInteractionAtPos(std::declval<const Pos&>())) {}
};

// support new and old twop material law interface
template<class Scalar, class SpatialParams, class Element, class Scv, class ElemSol>
class TwoPMaterialLawWrapper
{
    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;

    template<class S>
    static constexpr bool hasNew()
    { return decltype(isValid(HasNewFIAIF<Element, Scv, ElemSol>()).template check<S>())::value; }

    template<class S>
    static constexpr bool hasNewAtPos()
    { return decltype(isValid(HasNewFIAIFAtPos<GlobalPosition>()).template check<S>())::value; }

public:
    // pass scalar so template arguments can all be deduced
    TwoPMaterialLawWrapper(const Scalar& scalar,
                           const SpatialParams& sp,
                           const Element& element,
                           const Scv& scv,
                           const ElemSol& elemSol)
    : spatialParams_(sp), element_(element), scv_(scv), elemSol_(elemSol)
    {
        // maybe throw deprecation warning
        checkUsingOldMaterialLaw();
    }

    template<class S = SpatialParams, typename std::enable_if_t<!hasNew<S>() && !hasNewAtPos<S>(), int> = 0>
    [[deprecated("The material laws have been overhauled! Your spatial params implement the old interface. Use the new style material laws. Old material laws will no longer be supported after release 3.3")]]
    void checkUsingOldMaterialLaw() {}

    template<class S = SpatialParams, typename std::enable_if_t<hasNew<S>() || hasNewAtPos<S>(), int> = 0>
    void checkUsingOldMaterialLaw() {}

    /*!
     * \brief The capillary pressure-saturation curve
     */
    template<class S = SpatialParams, typename std::enable_if_t<!hasNew<S>() && !hasNewAtPos<S>(), int> = 0>
    Scalar pc(const Scalar sw) const
    {
        const auto& params = spatialParams_.materialLawParamsDeprecated(element_, scv_, elemSol_);
        return S::MaterialLaw::pc(params, sw);
    }

    template<class S = SpatialParams, typename std::enable_if_t<hasNew<S>(), int> = 0>
    Scalar pc(const Scalar sw) const
    {
        return spatialParams_.fluidMatrixInteraction(element_, scv_, elemSol_).pcSw().pc(sw);
    }

    template<class S = SpatialParams, typename std::enable_if_t<hasNewAtPos<S>() && !hasNew<S>(), int> = 0>
    Scalar pc(const Scalar sw) const
    {
        return spatialParams_.fluidMatrixInteractionAtPos(scv_.center()).pc(sw);
    }

    template<class S = SpatialParams, typename std::enable_if_t<!hasNew<S>() && !hasNewAtPos<S>(), int> = 0>
    Scalar dpc_dsw(const Scalar sw) const
    {
        const auto& params = spatialParams_.materialLawParamsDeprecated(element_, scv_, elemSol_);
        return S::MaterialLaw::dpc_dsw(params, sw);
    }

    template<class S = SpatialParams, typename std::enable_if_t<hasNew<S>(), int> = 0>
    Scalar dpc_dsw(const Scalar sw) const
    {
        return spatialParams_.fluidMatrixInteraction(element_, scv_, elemSol_).dpc_dsw(sw);
    }

    template<class S = SpatialParams, typename std::enable_if_t<hasNewAtPos<S>() && !hasNew<S>(), int> = 0>
    Scalar dpc_dsw(const Scalar sw) const
    {
        return spatialParams_.fluidMatrixInteractionAtPos(scv_.center()).dpc_dsw(sw);
    }

    template<class S = SpatialParams, typename std::enable_if_t<!hasNew<S>() && !hasNewAtPos<S>(), int> = 0>
    Scalar endPointPc() const
    {
        const auto& params = spatialParams_.materialLawParamsDeprecated(element_, scv_, elemSol_);
        return S::MaterialLaw::endPointPc(params);
    }

    template<class S = SpatialParams, typename std::enable_if_t<hasNew<S>(), int> = 0>
    Scalar endPointPc() const
    {
        return spatialParams_.fluidMatrixInteraction(element_, scv_, elemSol_).pcSw().endPointPc();
    }

    template<class S = SpatialParams, typename std::enable_if_t<hasNewAtPos<S>() && !hasNew<S>(), int> = 0>
    Scalar endPointPc() const
    {
        return spatialParams_.fluidMatrixInteractionAtPos(scv_.center()).endPointPc();
    }

    template<class S = SpatialParams, typename std::enable_if_t<!hasNew<S>() && !hasNewAtPos<S>(), int> = 0>
    Scalar sw(const Scalar pc) const
    {
        const auto& params = spatialParams_.materialLawParamsDeprecated(element_, scv_, elemSol_);
        return S::MaterialLaw::sw(params, pc);
    }

    template<class S = SpatialParams, typename std::enable_if_t<hasNew<S>(), int> = 0>
    Scalar sw(const Scalar pc) const
    {
        return spatialParams_.fluidMatrixInteraction(element_, scv_, elemSol_).pcSw().sw(pc);
    }

    template<class S = SpatialParams, typename std::enable_if_t<hasNewAtPos<S>() && !hasNew<S>(), int> = 0>
    Scalar sw(const Scalar pc) const
    {
        return spatialParams_.fluidMatrixInteractionAtPos(scv_.center()).sw(pc);
    }

    template<class S = SpatialParams, typename std::enable_if_t<!hasNew<S>() && !hasNewAtPos<S>(), int> = 0>
    Scalar dsw_dpc(const Scalar pc) const
    {
        const auto& params = spatialParams_.materialLawParamsDeprecated(element_, scv_, elemSol_);
        return S::MaterialLaw::dsw_dpc(params, pc);
    }

    template<class S = SpatialParams, typename std::enable_if_t<hasNew<S>(), int> = 0>
    Scalar dsw_dpc(const Scalar pc) const
    {
        return spatialParams_.fluidMatrixInteraction(element_, scv_, elemSol_).dsw_dpc(pc);
    }

    template<class S = SpatialParams, typename std::enable_if_t<hasNewAtPos<S>() && !hasNew<S>(), int> = 0>
    Scalar dsw_dpc(const Scalar pc) const
    {
        return spatialParams_.fluidMatrixInteractionAtPos(scv_.center()).dsw_dpc(pc);
    }



    /*!
     * \brief The relative permeability for the wetting phase
     */
    template<class S = SpatialParams, typename std::enable_if_t<!hasNew<S>() && !hasNewAtPos<S>(), int> = 0>
    Scalar krw(const Scalar sw) const
    {
        const auto& params = spatialParams_.materialLawParamsDeprecated(element_, scv_, elemSol_);
        return S::MaterialLaw::krw(params, sw);
    }

    template<class S = SpatialParams, typename std::enable_if_t<hasNew<S>(), int> = 0>
    Scalar krw(const Scalar sw) const
    {
        return spatialParams_.fluidMatrixInteraction(element_, scv_, elemSol_).krSw().krw(sw);
    }

    template<class S = SpatialParams, typename std::enable_if_t<hasNewAtPos<S>() && !hasNew<S>(), int> = 0>
    Scalar krw(const Scalar sw) const
    {
        return spatialParams_.fluidMatrixInteractionAtPos(scv_.center()).krw(sw);
    }

    /*!
     * \brief The derivative of the relative permeability for the wetting phase w.r.t. saturation
     */
    template<class S = SpatialParams, typename std::enable_if_t<!hasNew<S>() && !hasNewAtPos<S>(), int> = 0>
    Scalar dkrw_dsw(const Scalar sw) const
    {
        const auto& params = spatialParams_.materialLawParamsDeprecated(element_, scv_, elemSol_);
        return S::MaterialLaw::dkrw_dsw(params, sw);
    }

    template<class S = SpatialParams, typename std::enable_if_t<hasNew<S>(), int> = 0>
    Scalar dkrw_dsw(const Scalar sw) const
    {
        return spatialParams_.fluidMatrixInteraction(element_, scv_, elemSol_).dkrw_dsw(sw);
    }

    template<class S = SpatialParams, typename std::enable_if_t<hasNewAtPos<S>() && !hasNew<S>(), int> = 0>
    Scalar dkrw_dsw(const Scalar sw) const
    {
        return spatialParams_.fluidMatrixInteractionAtPos(scv_.center()).dkrw_dsw(sw);
    }

    /*!
     * \brief The relative permeability for the non-wetting phase
     */
    template<class S = SpatialParams, typename std::enable_if_t<!hasNew<S>() && !hasNewAtPos<S>(), int> = 0>
    Scalar krn(const Scalar sw) const
    {
        const auto& params = spatialParams_.materialLawParamsDeprecated(element_, scv_, elemSol_);
        return S::MaterialLaw::krn(params, sw);
    }

    template<class S = SpatialParams, typename std::enable_if_t<hasNew<S>(), int> = 0>
    Scalar krn(const Scalar sw) const
    {
        return spatialParams_.fluidMatrixInteraction(element_, scv_, elemSol_).krn(sw);
    }

    template<class S = SpatialParams, typename std::enable_if_t<hasNewAtPos<S>() && !hasNew<S>(), int> = 0>
    Scalar krn(const Scalar sw) const
    {
        return spatialParams_.fluidMatrixInteractionAtPos(scv_.center()).krn(sw);
    }

    /*!
     * \brief The derivative of the relative permeability for the non-wetting phase w.r.t. saturation
     */
    template<class S = SpatialParams, typename std::enable_if_t<!hasNew<S>() && !hasNewAtPos<S>(), int> = 0>
    Scalar dkrn_dsw(const Scalar sw) const
    {
        const auto& params = spatialParams_.materialLawParamsDeprecated(element_, scv_, elemSol_);
        return S::MaterialLaw::dkrn_dsw(params, sw);
    }

    template<class S = SpatialParams, typename std::enable_if_t<hasNew<S>(), int> = 0>
    Scalar dkrn_dsw(const Scalar sw) const
    {
        return spatialParams_.fluidMatrixInteraction(element_, scv_, elemSol_).dkrn_dsw(sw);
    }

    template<class S = SpatialParams, typename std::enable_if_t<hasNewAtPos<S>() && !hasNew<S>(), int> = 0>
    Scalar dkrn_dsw(const Scalar sw) const
    {
        return spatialParams_.fluidMatrixInteractionAtPos(scv_.center()).dkrn_dsw(sw);
    }

private:
    const SpatialParams& spatialParams_;
    const Element& element_;
    const Scv& scv_;
    const ElemSol& elemSol_;
};

///////////////////////////////////////////////////////////////
// Deprecation warnings for the kinetic surface areas //
///////////////////////////////////////////////////////////////

// support old interface of surface area
template<class E, class SCV, class Sol>
struct HasNewAns
{
    template<class S>
    auto operator()(S&& sp)
    -> decltype(sp.nonwettingSolidInterfacialArea(std::declval<const E&>(),
                                                  std::declval<const SCV&>(),
                                                  std::declval<const Sol&>())) {}
};

template<class Pos>
struct HasNewAnsAtPos
{
    template<class S>
    auto operator()(S&& sp)
    -> decltype(sp.nonwettingSolidInterfacialAreaAtPos(std::declval<const Pos&>())) {}
};

// support old interface of surface area
template<class E, class SCV, class Sol>
struct HasNewAnw
{
    template<class S>
    auto operator()(S&& sp)
    -> decltype(sp.wettingNonwettingInterfacialArea(std::declval<const E&>(),
                                                    std::declval<const SCV&>(),
                                                    std::declval<const Sol&>())) {}
};

template<class Pos>
struct HasNewAnwAtPos
{
    template<class S>
    auto operator()(S&& sp)
    -> decltype(sp.wettingNonwettingInterfacialAreaAtPos(std::declval<const Pos&>())) {}
};

// support old interface of surface area
template<class E, class SCV, class Sol>
struct HasNewAws
{
    template<class S>
    auto operator()(S&& sp)
    -> decltype(sp.wettingSolidInterfacialArea(std::declval<const E&>(),
                                               std::declval<const SCV&>(),
                                               std::declval<const Sol&>())) {}
};

template<class Pos>
struct HasNewAwsAtPos
{
    template<class S>
    auto operator()(S&& sp)
    -> decltype(sp.wettingSolidInterfacialAreaAtPos(std::declval<const Pos&>())) {}
};

// support new and old twop material law interface
template<class Scalar, class SpatialParams, class Element, class Scv, class ElemSol>
class WettingNonwettingInterfacialAreaWrapper
{
    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;

    template<class S>
    static constexpr bool hasNewAnw()
    { return decltype(isValid(HasNewAnw<Element, Scv, ElemSol>()).template check<S>())::value; }

    template<class S>
    static constexpr bool hasNewAnwAtPos()
    { return decltype(isValid(HasNewAnwAtPos<GlobalPosition>()).template check<S>())::value; }


public:
    // pass scalar so template arguments can all be deduced
    WettingNonwettingInterfacialAreaWrapper(const Scalar& scalar,
                                            const SpatialParams& sp,
                                            const Element& element,
                                            const Scv& scv,
                                            const ElemSol& elemSol)
    : spatialParams_(sp), element_(element), scv_(scv), elemSol_(elemSol)
    {
        maybePrintWarning();
    }

    template<class S = SpatialParams>
    void maybePrintWarning() const
    {
      if constexpr (!hasNewAnw<S>() && !hasNewAnwAtPos<S>())
        printWarning();
    }

    [[deprecated("The material laws have been overhauled! Your spatial params implement the old interface. Use the new style material laws. Old material laws will no longer be supported after release 3.3")]]
    void printWarning() const {}

    template<class S = SpatialParams, typename std::enable_if_t<!hasNewAnw<S>() && !hasNewAnwAtPos<S>(), int> = 0>
    const auto& basicParams() const
    {
      printWarning();
      return spatialParams_.aWettingNonWettingSurfaceParams(element_, scv_, elemSol_);
    }

    template<class S = SpatialParams, typename std::enable_if_t<hasNewAnw<S>(), int> = 0>
    const auto& basicParams() const
    {
      return spatialParams_.wettingNonwettingInterfacialArea(element_, scv_, elemSol_).basicParams();
    }

    template<class S = SpatialParams, typename std::enable_if_t<hasNewAnwAtPos<S>() && !hasNewAnw<S>(), int> = 0>
    const auto& basicParams() const
    {
      return spatialParams_.wettingNonwettingInterfacialAreaAtPos(scv_.center()).basicParams();
    }

    /*!
     * \brief The capillary pressure-saturation curve
     */
    template<class S = SpatialParams, typename std::enable_if_t<!hasNewAnw<S>() && !hasNewAnwAtPos<S>(), int> = 0>
    Scalar area(const Scalar sw, const Scalar pc) const
    {
        printWarning();
        const auto& surfaceParams = spatialParams_.aWettingNonWettingSurfaceParams(element_, scv_, elemSol_);
        const auto& materialParams = spatialParams_.materialLawParams(element_, scv_, elemSol_);
        using AwnSurface = typename SpatialParams::AwnSurface;
        return AwnSurface::interfacialArea(surfaceParams, materialParams, sw, pc);
    }

    template<class S = SpatialParams, typename std::enable_if_t<hasNewAnw<S>(), int> = 0>
    Scalar area(const Scalar sw, const Scalar pc) const
    {
        return spatialParams_.wettingNonwettingInterfacialArea(element_, scv_, elemSol_).area(sw, pc);
    }

    template<class S = SpatialParams, typename std::enable_if_t<hasNewAnwAtPos<S>() && !hasNewAnw<S>(), int> = 0>
    Scalar area(const Scalar sw, const Scalar pc) const
    {
        return spatialParams_.wettingNonwettingInterfacialAreaAtPos(scv_.center()).area(sw, pc);
    }

    template<class S = SpatialParams, typename std::enable_if_t<!hasNewAnw<S>() && !hasNewAnwAtPos<S>(), int> = 0>
    Scalar darea_dpc(const Scalar sw, const Scalar pc)
    {
      printWarning();
      const auto& surfaceParams = spatialParams_.aWettingNonWettingSurfaceParams(element_, scv_, elemSol_);
      using AwnSurface = typename SpatialParams::AwnSurface;
      return AwnSurface::dawn_dpc(surfaceParams, sw, pc);
    }

    template<class S = SpatialParams, typename std::enable_if_t<hasNewAnw<S>(), int> = 0>
    Scalar darea_dpc(const Scalar sw, const Scalar pc)
    {
      return spatialParams_.wettingNonwettingInterfacialArea(element_, scv_, elemSol_).darea_dpc(sw, pc);
    }

    template<class S = SpatialParams, typename std::enable_if_t<hasNewAnwAtPos<S>() && !hasNewAnw<S>(), int> = 0>
    Scalar darea_dpc(const Scalar sw, const Scalar pc)
    {
      return spatialParams_.wettingNonwettingInterfacialAreaAtPos(scv_.center()).darea_dpc(sw, pc);
    }

    template<class S = SpatialParams, typename std::enable_if_t<!hasNewAnw<S>() && !hasNewAnwAtPos<S>(), int> = 0>
    Scalar darea_dsw(const Scalar sw, const Scalar pc)
    {
      printWarning();
      const auto& surfaceParams = spatialParams_.aWettingNonWettingSurfaceParams(element_, scv_, elemSol_);
      using AwnSurface = typename SpatialParams::AwnSurface;
      return AwnSurface::dawn_dsw(surfaceParams, sw, pc);
    }

    template<class S = SpatialParams, typename std::enable_if_t<hasNewAnw<S>(), int> = 0>
    Scalar darea_dsw(const Scalar sw, const Scalar pc)
    {
      return spatialParams_.wettingNonwettingInterfacialArea(element_, scv_, elemSol_).darea_dsw(sw, pc);
    }

    template<class S = SpatialParams, typename std::enable_if_t<hasNewAnwAtPos<S>() && !hasNewAnw<S>(), int> = 0>
    Scalar darea_dsw(const Scalar sw, const Scalar pc)
    {
      return spatialParams_.wettingNonwettingInterfacialAreaAtPos(scv_.center()).darea_dsw(sw, pc);
    }

private:
    const SpatialParams& spatialParams_;
    const Element& element_;
    const Scv& scv_;
    const ElemSol& elemSol_;
};

// support new and old twop material law interface
template<class Scalar, class SpatialParams, class Element, class Scv, class ElemSol>
class NonwettingSolidInterfacialAreaWrapper
{
    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;

    template<class S>
    static constexpr bool hasNewAns()
    { return decltype(isValid(HasNewAns<Element, Scv, ElemSol>()).template check<S>())::value; }

    template<class S>
    static constexpr bool hasNewAnsAtPos()
    { return decltype(isValid(HasNewAnsAtPos<GlobalPosition>()).template check<S>())::value; }


public:
    // pass scalar so template arguments can all be deduced
    NonwettingSolidInterfacialAreaWrapper(const Scalar& scalar,
                                          const SpatialParams& sp,
                                          const Element& element,
                                          const Scv& scv,
                                          const ElemSol& elemSol)
    : spatialParams_(sp), element_(element), scv_(scv), elemSol_(elemSol)
    {
        maybePrintWarning();
    }

    template<class S = SpatialParams>
    void maybePrintWarning() const
    {
      if constexpr (!hasNewAns<S>() && !hasNewAnsAtPos<S>())
        printWarning();
    }

    [[deprecated("The material laws have been overhauled! Your spatial params implement the old interface. Use the new style material laws. Old material laws will no longer be supported after release 3.3")]]
    void printWarning() const {}

    template<class S = SpatialParams, typename std::enable_if_t<!hasNewAns<S>() && !hasNewAnsAtPos<S>(), int> = 0>
    const auto& basicParams() const
    {
      printWarning();
      return spatialParams_.aNonWettingSolidSurfaceParams(element_, scv_, elemSol_);
    }

    template<class S = SpatialParams, typename std::enable_if_t<hasNewAns<S>(), int> = 0>
    const auto& basicParams() const
    {
      return spatialParams_.nonwettingSolidInterfacialArea(element_, scv_, elemSol_).basicParams();
    }

    template<class S = SpatialParams, typename std::enable_if_t<hasNewAnsAtPos<S>() && !hasNewAns<S>(), int> = 0>
    const auto& basicParams() const
    {
      return spatialParams_.nonwettingSolidInterfacialAreaAtPos(scv_.center()).basicParams();
    }

    /*!
     * \brief The capillary pressure-saturation curve
     */
    template<class S = SpatialParams, typename std::enable_if_t<!hasNewAns<S>() && !hasNewAnsAtPos<S>(), int> = 0>
    Scalar area(const Scalar sw, const Scalar pc) const
    {
        printWarning();
        const auto& surfaceParams = spatialParams_.aNonWettingSolidSurfaceParams(element_, scv_, elemSol_);
        const auto& materialParams = spatialParams_.materialLawParams(element_, scv_, elemSol_);
        using AnsSurface = typename SpatialParams::AnsSurface;
        return AnsSurface::interfacialArea(surfaceParams, materialParams, sw, pc);
    }

    template<class S = SpatialParams, typename std::enable_if_t<hasNewAns<S>(), int> = 0>
    Scalar area(const Scalar sw, const Scalar pc) const
    {
        return spatialParams_.nonwettingSolidInterfacialArea(element_, scv_, elemSol_).area(sw, pc);
    }

    template<class S = SpatialParams, typename std::enable_if_t<hasNewAnsAtPos<S>() && !hasNewAns<S>(), int> = 0>
    Scalar area(const Scalar sw, const Scalar pc) const
    {
        return spatialParams_.nonwettingSolidInterfacialAreaAtPos(scv_.center()).area(sw, pc);
    }

    template<class S = SpatialParams, typename std::enable_if_t<!hasNewAns<S>() && !hasNewAnsAtPos<S>(), int> = 0>
    Scalar darea_dpc(const Scalar sw, const Scalar pc)
    {
      printWarning();
      const auto& surfaceParams = spatialParams_.aNonWettingSolidSurfaceParams(element_, scv_, elemSol_);
      using AnsSurface = typename SpatialParams::AnsSurface;
      return AnsSurface::dawn_dpc(surfaceParams, sw, pc);
    }

    template<class S = SpatialParams, typename std::enable_if_t<hasNewAns<S>(), int> = 0>
    Scalar darea_dpc(const Scalar sw, const Scalar pc)
    {
      return spatialParams_.nonwettingSolidInterfacialArea(element_, scv_, elemSol_).darea_dpc(sw, pc);
    }

    template<class S = SpatialParams, typename std::enable_if_t<hasNewAnsAtPos<S>() && !hasNewAns<S>(), int> = 0>
    Scalar darea_dpc(const Scalar sw, const Scalar pc)
    {
      return spatialParams_.nonwettingSolidInterfacialAreaAtPos(scv_.center()).darea_dpc(sw, pc);
    }

    template<class S = SpatialParams, typename std::enable_if_t<!hasNewAns<S>() && !hasNewAnsAtPos<S>(), int> = 0>
    Scalar darea_dsw(const Scalar sw, const Scalar pc)
    {
      printWarning();
      const auto& surfaceParams = spatialParams_.aNonWettingSolidSurfaceParams(element_, scv_, elemSol_);
      using AnsSurface = typename SpatialParams::AnsSurface;
      return AnsSurface::dawn_dsw(surfaceParams, sw, pc);
    }

    template<class S = SpatialParams, typename std::enable_if_t<hasNewAns<S>(), int> = 0>
    Scalar darea_dsw(const Scalar sw, const Scalar pc)
    {
      return spatialParams_.nonwettingSolidInterfacialArea(element_, scv_, elemSol_).darea_dsw(sw, pc);
    }

    template<class S = SpatialParams, typename std::enable_if_t<hasNewAnsAtPos<S>() && !hasNewAns<S>(), int> = 0>
    Scalar darea_dsw(const Scalar sw, const Scalar pc)
    {
      return spatialParams_.nonwettingSolidInterfacialAreaAtPos(scv_.center()).darea_dsw(sw, pc);
    }

private:
    const SpatialParams& spatialParams_;
    const Element& element_;
    const Scv& scv_;
    const ElemSol& elemSol_;
};

// support new and old twop material law interface
template<class Scalar, class SpatialParams, class Element, class Scv, class ElemSol>
class WettingSolidInterfacialAreaWrapper
{
    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;

    template<class S>
    static constexpr bool hasNewAws()
    { return decltype(isValid(HasNewAws<Element, Scv, ElemSol>()).template check<S>())::value; }

    template<class S>
    static constexpr bool hasNewAwsAtPos()
    { return decltype(isValid(HasNewAwsAtPos<GlobalPosition>()).template check<S>())::value; }


public:
    // pass scalar so template arguments can all be deduced
    WettingSolidInterfacialAreaWrapper(const Scalar& scalar,
                                       const SpatialParams& sp,
                                       const Element& element,
                                       const Scv& scv,
                                       const ElemSol& elemSol)
    : spatialParams_(sp), element_(element), scv_(scv), elemSol_(elemSol)
    {
        maybePrintWarning();
    }

    template<class S = SpatialParams>
    void maybePrintWarning() const
    {
      if constexpr (!hasNewAws<S>() && !hasNewAwsAtPos<S>())
        printWarning();
    }

    [[deprecated("The material laws have been overhauled! Your spatial params implement the old interface. Use the new style material laws. Old material laws will no longer be supported after release 3.3")]]
    void printWarning() const {}

    template<class S = SpatialParams, typename std::enable_if_t<!hasNewAws<S>() && !hasNewAwsAtPos<S>(), int> = 0>
    const auto& basicParams() const
    {
      printWarning();
      return spatialParams_.aWettingSolidSurfaceParams(element_, scv_, elemSol_);
    }

    template<class S = SpatialParams, typename std::enable_if_t<hasNewAws<S>(), int> = 0>
    const auto& basicParams() const
    {
      return spatialParams_.wettingSolidInterfacialArea(element_, scv_, elemSol_).basicParams();
    }

    template<class S = SpatialParams, typename std::enable_if_t<hasNewAwsAtPos<S>() && !hasNewAws<S>(), int> = 0>
    const auto& basicParams() const
    {
      return spatialParams_.wettingSolidInterfacialAreaAtPos(scv_.center()).basicParams();
    }

    /*!
     * \brief The capillary pressure-saturation curve
     */
    template<class S = SpatialParams, typename std::enable_if_t<!hasNewAws<S>() && !hasNewAwsAtPos<S>(), int> = 0>
    Scalar area(const Scalar sw, const Scalar pc) const
    {
        printWarning();
        const auto& surfaceParams = spatialParams_.aWettingSolidSurfaceParams(element_, scv_, elemSol_);
        const auto& materialParams = spatialParams_.materialLawParams(element_, scv_, elemSol_);
        using AwsSurface = typename SpatialParams::AwsSurface;
        return AwsSurface::interfacialArea(surfaceParams, materialParams, sw, pc);
    }

    template<class S = SpatialParams, typename std::enable_if_t<hasNewAws<S>(), int> = 0>
    Scalar area(const Scalar sw, const Scalar pc) const
    {
        return spatialParams_.wettingSolidInterfacialArea(element_, scv_, elemSol_).area(sw, pc);
    }

    template<class S = SpatialParams, typename std::enable_if_t<hasNewAwsAtPos<S>() && !hasNewAws<S>(), int> = 0>
    Scalar area(const Scalar sw, const Scalar pc) const
    {
        return spatialParams_.wettingSolidInterfacialAreaAtPos(scv_.center()).area(sw, pc);
    }

    template<class S = SpatialParams, typename std::enable_if_t<!hasNewAws<S>() && !hasNewAwsAtPos<S>(), int> = 0>
    Scalar darea_dpc(const Scalar sw, const Scalar pc)
    {
      printWarning();
      const auto& surfaceParams = spatialParams_.aWettingSolidSurfaceParams(element_, scv_, elemSol_);
      using AwsSurface = typename SpatialParams::AwsSurface;
      return AwsSurface::dawn_dpc(surfaceParams, sw, pc);
    }

    template<class S = SpatialParams, typename std::enable_if_t<hasNewAws<S>(), int> = 0>
    Scalar darea_dpc(const Scalar sw, const Scalar pc)
    {
      return spatialParams_.wettingSolidInterfacialArea(element_, scv_, elemSol_).darea_dpc(sw, pc);
    }

    template<class S = SpatialParams, typename std::enable_if_t<hasNewAwsAtPos<S>() && !hasNewAws<S>(), int> = 0>
    Scalar darea_dpc(const Scalar sw, const Scalar pc)
    {
      return spatialParams_.wettingSolidInterfacialAreaAtPos(scv_.center()).darea_dpc(sw, pc);
    }

    template<class S = SpatialParams, typename std::enable_if_t<!hasNewAws<S>() && !hasNewAwsAtPos<S>(), int> = 0>
    Scalar darea_dsw(const Scalar sw, const Scalar pc)
    {
      printWarning();
      const auto& surfaceParams = spatialParams_.aWettingSolidSurfaceParams(element_, scv_, elemSol_);
      using AwsSurface = typename SpatialParams::AwsSurface;
      return AwsSurface::dawn_dsw(surfaceParams, sw, pc);
    }

    template<class S = SpatialParams, typename std::enable_if_t<hasNewAws<S>(), int> = 0>
    Scalar darea_dsw(const Scalar sw, const Scalar pc)
    {
      return spatialParams_.wettingSolidInterfacialArea(element_, scv_, elemSol_).darea_dsw(sw, pc);
    }

    template<class S = SpatialParams, typename std::enable_if_t<hasNewAwsAtPos<S>() && !hasNewAws<S>(), int> = 0>
    Scalar darea_dsw(const Scalar sw, const Scalar pc)
    {
      return spatialParams_.wettingSolidInterfacialAreaAtPos(scv_.center()).darea_dsw(sw, pc);
    }

private:
    const SpatialParams& spatialParams_;
    const Element& element_;
    const Scv& scv_;
    const ElemSol& elemSol_;
};
} // end namespace Deprecated
#endif

} // end namespace Dumux
#endif
