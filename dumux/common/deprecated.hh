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
#include <dumux/material/fluidmatrixinteractions/fluidmatrixinteraction.hh>
#include <dumux/material/fluidmatrixinteractions/mp/mpadapter.hh>

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


template<class ScalarT, class SpatialParams, class Element, class Scv, class ElemSol>
class PcKrSwHelper : public FluidMatrix::Adapter<PcKrSwHelper<ScalarT, SpatialParams, Element, Scv, ElemSol>, FluidMatrix::PcKrSw>
{
public:
    using Scalar = ScalarT;

    // pass scalar so template arguments can all be deduced
    PcKrSwHelper(const Scalar& scalar,
                 const SpatialParams& sp,
                 const Element& element,
                 const Scv& scv,
                 const ElemSol& elemSol)
    : spatialParams_(sp), element_(element), scv_(scv), elemSol_(elemSol)
    {}

    Scalar krw(const Scalar sw) const
    {
        const auto& params = spatialParams_.materialLawParamsDeprecated(element_, scv_, elemSol_);
        return SpatialParams::MaterialLaw::krw(params, sw);
    }

    Scalar krn(const Scalar sw) const
    {
        const auto& params = spatialParams_.materialLawParamsDeprecated(element_, scv_, elemSol_);
        return SpatialParams::MaterialLaw::krn(params, sw);
    }

    Scalar pc(const Scalar sw) const
    {
        const auto& params = spatialParams_.materialLawParamsDeprecated(element_, scv_, elemSol_);
        return SpatialParams::MaterialLaw::pc(params, sw);
    }

    Scalar dpc_dsw(const Scalar sw) const
    {
        const auto& params = spatialParams_.materialLawParamsDeprecated(element_, scv_, elemSol_);
        return SpatialParams::MaterialLaw::dpc_dsw(params, sw);
    }

    Scalar endPointPc() const
    {
        const auto& params = spatialParams_.materialLawParamsDeprecated(element_, scv_, elemSol_);
        return SpatialParams::MaterialLaw::endPointPc(params);
    }

    Scalar sw(const Scalar pc) const
    {
        const auto& params = spatialParams_.materialLawParamsDeprecated(element_, scv_, elemSol_);
        return SpatialParams::MaterialLaw::sw(params, pc);
    }

    Scalar dsw_dpc(const Scalar pc) const
    {
        const auto& params = spatialParams_.materialLawParamsDeprecated(element_, scv_, elemSol_);
        return SpatialParams::MaterialLaw::dsw_dpc(params, pc);
    }

    Scalar dkrw_dsw(const Scalar sw) const
    {
        const auto& params = spatialParams_.materialLawParamsDeprecated(element_, scv_, elemSol_);
        return SpatialParams::MaterialLaw::dkrw_dsw(params, sw);
    }

    Scalar dkrn_dsw(const Scalar sw) const
    {
        const auto& params = spatialParams_.materialLawParamsDeprecated(element_, scv_, elemSol_);
        return SpatialParams::MaterialLaw::dkrn_dsw(params, sw);
    }

    const auto& basicParams() const
    { return spatialParams_.materialLawParamsDeprecated(element_, scv_, elemSol_); }

    const auto& effToAbsParams() const
    { return spatialParams_.materialLawParamsDeprecated(element_, scv_, elemSol_); }

private:
    const SpatialParams& spatialParams_;
    const Element& element_;
    const Scv& scv_;
    const ElemSol& elemSol_;
};

template<class ScalarT, class SpatialParams, class Element, class Scv, class ElemSol>
class PcKrSwThreePHelper : public FluidMatrix::Adapter<PcKrSwThreePHelper<ScalarT, SpatialParams, Element, Scv, ElemSol>, FluidMatrix::ThreePhasePcKrSw>
{
public:
    using Scalar = ScalarT;

    // pass scalar so template arguments can all be deduced
    PcKrSwThreePHelper(const Scalar& scalar,
                       const SpatialParams& sp,
                       const Element& element,
                       const Scv& scv,
                       const ElemSol& elemSol)
    : spatialParams_(sp), element_(element), scv_(scv), elemSol_(elemSol)
    {}

    Scalar pcgw(const Scalar sw, const Scalar /*dummySn*/) const
    {
        const auto& params = spatialParams_.materialLawParamsDeprecated(element_, scv_, elemSol_);
        return SpatialParams::MaterialLaw::pcgw(params, sw);
    }

    Scalar pcnw(const Scalar sw, const Scalar /*dummySn*/) const
    {
        const auto& params = spatialParams_.materialLawParamsDeprecated(element_, scv_, elemSol_);
        return SpatialParams::MaterialLaw::pcnw(params, sw);
    }

    Scalar pcgn(const Scalar sw, const Scalar sn) const
    {
        const auto& params = spatialParams_.materialLawParamsDeprecated(element_, scv_, elemSol_);
        return SpatialParams::MaterialLaw::pcgn(params, sw + sn);
    }

    Scalar pcAlpha(const Scalar /*dummySw*/, const Scalar sn) const
    {
        const auto& params = spatialParams_.materialLawParamsDeprecated(element_, scv_, elemSol_);
        return SpatialParams::MaterialLaw::pcAlpha(params, sn);
    }

    Scalar krw(const Scalar sw, const Scalar sn) const
    {
        const auto& params = spatialParams_.materialLawParamsDeprecated(element_, scv_, elemSol_);
        return SpatialParams::MaterialLaw::krw(params, sw, sn);
    }

    Scalar krn(const Scalar sw, const Scalar sn) const
    {
        const auto& params = spatialParams_.materialLawParamsDeprecated(element_, scv_, elemSol_);
        return SpatialParams::MaterialLaw::krw(params, sw, sn);
    }

    Scalar krg(const Scalar sw, const Scalar sn) const
    {
        const auto& params = spatialParams_.materialLawParamsDeprecated(element_, scv_, elemSol_);
        return SpatialParams::MaterialLaw::krg(params, sw, sn);
    }

    Scalar kr(const int phaseIdx, const Scalar sw, const Scalar sn) const
    {
        const auto& params = spatialParams_.materialLawParamsDeprecated(element_, scv_, elemSol_);
        return SpatialParams::MaterialLaw::kr(params, phaseIdx, sw, sn, 1 - sw - sn);
    }

    const auto& basicParams() const
    { return spatialParams_.materialLawParamsDeprecated(element_, scv_, elemSol_); }

    const auto& effToAbsParams() const
    { return spatialParams_.materialLawParamsDeprecated(element_, scv_, elemSol_); }

private:
    const SpatialParams& spatialParams_;
    const Element& element_;
    const Scv& scv_;
    const ElemSol& elemSol_;
};

// for implicit models
template<int numPhases = 2, class Scalar, class SpatialParams, class Element, class Scv, class ElemSol>
auto makePcKrSw(const Scalar& scalar,
                const SpatialParams& sp,
                const Element& element,
                const Scv& scv,
                const ElemSol& elemSol)
{
    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;
    constexpr bool hasNew = decltype(isValid(HasNewFIAIF<Element, Scv, ElemSol>()).template check<SpatialParams>())::value;
    [[maybe_unused]] constexpr bool hasNewAtPos = decltype(isValid(HasNewFIAIFAtPos<GlobalPosition>()).template check<SpatialParams>())::value;
    if constexpr (hasNew)
        return sp.fluidMatrixInteraction(element, scv, elemSol);
    else if constexpr (hasNewAtPos)
        return sp.fluidMatrixInteractionAtPos(scv.center());
    else
    {
        if constexpr (numPhases == 2)
            return makeFluidMatrixInteraction(PcKrSwHelper(scalar, sp, element, scv, elemSol));
        else
            return makeFluidMatrixInteraction(PcKrSwThreePHelper(scalar, sp, element, scv, elemSol));
    }
}

// for sequential models
template<class Scalar, class SpatialParams, class Element>
auto makePcKrSw(const Scalar& scalar,
                const SpatialParams& sp,
                const Element& element)
{
    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;
    constexpr bool hasNewAtPos = decltype(isValid(HasNewFIAIFAtPos<GlobalPosition>()).template check<SpatialParams>())::value;
    if constexpr (hasNewAtPos)
        return sp.fluidMatrixInteractionAtPos(element.geometry().center());
    else
    {
        using DummyScv = int;
        using DummyElemSol = int;
        return makeFluidMatrixInteraction(PcKrSwHelper(scalar, sp, element, DummyScv(), DummyElemSol()));
    }
}


///////////////////////////////////////////////////////////////
// Deprecation warnings for the mp material law stuff //
///////////////////////////////////////////////////////////////

template<class ScalarT, class SpatialParams, class Element, class Scv, class ElemSol, class NumPhases>
class PcKrSwMPHelper : public FluidMatrix::Adapter<PcKrSwMPHelper<ScalarT, SpatialParams, Element, Scv, ElemSol, NumPhases>, FluidMatrix::MultiPhasePcKrSw>
{
    using MaterialLaw = typename SpatialParams::MaterialLaw;
    using MPAdapter = Dumux::MPAdapter<MaterialLaw, NumPhases{}()>;
public:
    using Scalar = ScalarT;

    // pass scalar so template arguments can all be deduced
    PcKrSwMPHelper(const Scalar& scalar,
                   const SpatialParams& sp,
                   const Element& element,
                   const Scv& scv,
                   const ElemSol& elemSol,
                   const NumPhases& np)
    : spatialParams_(sp), element_(element), scv_(scv), elemSol_(elemSol)
    {}

    template<class FluidState>
    auto capillaryPressures(const FluidState& fs, int wPhaseIdx) const
    {
        const auto& params = spatialParams_.materialLawParamsDeprecated(element_, scv_, elemSol_);
        Dune::FieldVector<Scalar, NumPhases{}()> pc;
        MPAdapter::capillaryPressures(pc, params, fs, wPhaseIdx);
        return pc;
    }

    template<class FluidState>
    auto relativePermeabilities(const FluidState& fs, int wPhaseIdx) const
    {
        const auto& params = spatialParams_.materialLawParamsDeprecated(element_, scv_, elemSol_);
        Dune::FieldVector<Scalar, NumPhases{}()> kr;
        MPAdapter::capillaryPressures(kr, params, fs, wPhaseIdx);
        return kr;
    }

private:
    const SpatialParams& spatialParams_;
    const Element& element_;
    const Scv& scv_;
    const ElemSol& elemSol_;
};

template<class Scalar, class SpatialParams, class Element, class Scv, class ElemSol, class NumPhases>
auto makeMPPcKrSw(const Scalar& scalar,
                  const SpatialParams& sp,
                  const Element& element,
                  const Scv& scv,
                  const ElemSol& elemSol,
                  const NumPhases& np)
{
    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;
    constexpr bool hasNew = decltype(isValid(HasNewFIAIF<Element, Scv, ElemSol>()).template check<SpatialParams>())::value;
    constexpr bool hasNewAtPos = decltype(isValid(HasNewFIAIFAtPos<GlobalPosition>()).template check<SpatialParams>())::value;
    if constexpr (hasNew)
        return sp.fluidMatrixInteraction(element, scv, elemSol);
    else if constexpr (hasNewAtPos)
        return sp.fluidMatrixInteractionAtPos(scv.center());
    else
        return makeFluidMatrixInteraction(PcKrSwMPHelper(scalar, sp, element, scv, elemSol, np));
}

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

template<class ScalarT, class SpatialParams, class Element, class Scv, class ElemSol>
class WettingNonwettingInterfacialArea : public FluidMatrix::Adapter<WettingNonwettingInterfacialArea<ScalarT, SpatialParams, Element, Scv, ElemSol>, FluidMatrix::WettingNonwettingInterfacialAreaPcSw>
{
public:
    using Scalar = ScalarT;

    WettingNonwettingInterfacialArea(const Scalar& scalar,
                                     const SpatialParams& sp,
                                     const Element& element,
                                     const Scv& scv,
                                     const ElemSol& elemSol)
    : spatialParams_(sp), element_(element), scv_(scv), elemSol_(elemSol)
    {}

    const auto& basicParams() const
    {
      return spatialParams_.aWettingNonWettingSurfaceParams(element_, scv_, elemSol_);
    }

    Scalar area(const Scalar sw, const Scalar pc) const
    {
        const auto& surfaceParams = spatialParams_.aWettingNonWettingSurfaceParams(element_, scv_, elemSol_);
        const auto& materialParams = spatialParams_.materialLawParams(element_, scv_, elemSol_);
        using AwnSurface = typename SpatialParams::AwnSurface;
        return AwnSurface::interfacialArea(surfaceParams, materialParams, sw, pc);
    }

    Scalar darea_dpc(const Scalar sw, const Scalar pc)
    {
      const auto& surfaceParams = spatialParams_.aWettingNonWettingSurfaceParams(element_, scv_, elemSol_);
      using AwnSurface = typename SpatialParams::AwnSurface;
      return AwnSurface::dawn_dpc(surfaceParams, sw, pc);
    }

    Scalar darea_dsw(const Scalar sw, const Scalar pc)
    {
      const auto& surfaceParams = spatialParams_.aWettingNonWettingSurfaceParams(element_, scv_, elemSol_);
      using AwnSurface = typename SpatialParams::AwnSurface;
      return AwnSurface::dawn_dsw(surfaceParams, sw, pc);
    }

private:
    const SpatialParams& spatialParams_;
    const Element& element_;
    const Scv& scv_;
    const ElemSol& elemSol_;
};

template<class ScalarT, class SpatialParams, class Element, class Scv, class ElemSol>
class NonwettingSolidInterfacialArea : public FluidMatrix::Adapter<NonwettingSolidInterfacialArea<ScalarT, SpatialParams, Element, Scv, ElemSol>, FluidMatrix::NonwettingSolidInterfacialAreaPcSw>
{
public:
    using Scalar = ScalarT;

    NonwettingSolidInterfacialArea(const Scalar& scalar,
                                   const SpatialParams& sp,
                                   const Element& element,
                                   const Scv& scv,
                                   const ElemSol& elemSol)
    : spatialParams_(sp), element_(element), scv_(scv), elemSol_(elemSol)
    {}

    const auto& basicParams() const
    {
      return spatialParams_.aNonWettingSolidSurfaceParams(element_, scv_, elemSol_);
    }

    Scalar area(const Scalar sw, const Scalar pc) const
    {
        const auto& surfaceParams = spatialParams_.aNonWettingSolidSurfaceParams(element_, scv_, elemSol_);
        const auto& materialParams = spatialParams_.materialLawParams(element_, scv_, elemSol_);
        using AnsSurface = typename SpatialParams::AnsSurface;
        return AnsSurface::interfacialArea(surfaceParams, materialParams, sw, pc);
    }

    Scalar darea_dpc(const Scalar sw, const Scalar pc)
    {
      const auto& surfaceParams = spatialParams_.aNonWettingSolidSurfaceParams(element_, scv_, elemSol_);
      using AnsSurface = typename SpatialParams::AnsSurface;
      return AnsSurface::dawn_dpc(surfaceParams, sw, pc);
    }

    Scalar darea_dsw(const Scalar sw, const Scalar pc)
    {
      const auto& surfaceParams = spatialParams_.aNonWettingSolidSurfaceParams(element_, scv_, elemSol_);
      using AnsSurface = typename SpatialParams::AnsSurface;
      return AnsSurface::dawn_dsw(surfaceParams, sw, pc);
    }

private:
    const SpatialParams& spatialParams_;
    const Element& element_;
    const Scv& scv_;
    const ElemSol& elemSol_;
};

template<class ScalarT, class SpatialParams, class Element, class Scv, class ElemSol>
class WettingSolidInterfacialArea : public FluidMatrix::Adapter<WettingSolidInterfacialArea<ScalarT, SpatialParams, Element, Scv, ElemSol>, FluidMatrix::WettingSolidInterfacialAreaPcSw>
{
public:
    using Scalar = ScalarT;

    WettingSolidInterfacialArea(const Scalar& scalar,
                                const SpatialParams& sp,
                                const Element& element,
                                const Scv& scv,
                                const ElemSol& elemSol)
    : spatialParams_(sp), element_(element), scv_(scv), elemSol_(elemSol)
    {}

    const auto& basicParams() const
    {
      return spatialParams_.aWettingSolidSurfaceParams(element_, scv_, elemSol_);
    }

    Scalar area(const Scalar sw, const Scalar pc) const
    {
        const auto& surfaceParams = spatialParams_.aWettingSolidSurfaceParams(element_, scv_, elemSol_);
        const auto& materialParams = spatialParams_.materialLawParams(element_, scv_, elemSol_);
        using AwsSurface = typename SpatialParams::AwsSurface;
        return AwsSurface::interfacialArea(surfaceParams, materialParams, sw, pc);
    }

    Scalar darea_dpc(const Scalar sw, const Scalar pc)
    {
      const auto& surfaceParams = spatialParams_.aWettingSolidSurfaceParams(element_, scv_, elemSol_);
      using AwsSurface = typename SpatialParams::AwsSurface;
      return AwsSurface::dawn_dpc(surfaceParams, sw, pc);
    }

    Scalar darea_dsw(const Scalar sw, const Scalar pc)
    {
      const auto& surfaceParams = spatialParams_.aWettingSolidSurfaceParams(element_, scv_, elemSol_);
      using AwsSurface = typename SpatialParams::AwsSurface;
      return AwsSurface::dawn_dsw(surfaceParams, sw, pc);
    }

private:
    const SpatialParams& spatialParams_;
    const Element& element_;
    const Scv& scv_;
    const ElemSol& elemSol_;
};

template<class Scalar, class SpatialParams, class Element, class Scv, class ElemSol>
auto makeInterfacialArea(const Scalar& scalar,
                               const SpatialParams& sp,
                               const Element& element,
                               const Scv& scv,
                               const ElemSol& elemSol)
{
    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;
    constexpr bool hasNew = decltype(isValid(HasNewFIAIF<Element, Scv, ElemSol>()).template check<SpatialParams>())::value;
    constexpr bool hasNewAtPos = decltype(isValid(HasNewFIAIFAtPos<GlobalPosition>()).template check<SpatialParams>())::value;
    if constexpr (hasNew)
        return sp.fluidMatrixInteraction(element, scv, elemSol);
    else if constexpr (hasNewAtPos)
        return sp.fluidMatrixInteractionAtPos(scv.center());
    else
        return makeFluidMatrixInteraction(WettingNonwettingInterfacialArea(scalar, sp, element, scv, elemSol),
                                          NonwettingSolidInterfacialArea(scalar, sp, element, scv, elemSol),
                                          WettingSolidInterfacialArea(scalar, sp, element, scv, elemSol));
}

} // end namespace Deprecated
#endif

} // end namespace Dumux
#endif
