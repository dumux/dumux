// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup OnePTests
 * \brief A test for marking some faces as internal boundaries
 */

#ifndef DUMUX_INCOMPRESSIBLE_ONEP_TEST_PROBLEM_INTERNAL_BARRIER_PROPERTIES_HH
#define DUMUX_INCOMPRESSIBLE_ONEP_TEST_PROBLEM_INTERNAL_BARRIER_PROPERTIES_HH

#include <memory>
#include <test/porousmediumflow/1p/incompressible/properties.hh>
#include <dumux/discretization/cellcentered/tpfa/fvgridgeometry.hh>
#include <dumux/common/parameters.hh>

namespace Dumux {

// this is useful to mark some faces as internal boundaries
// it is also useful for network grids where maybe some faces should be considered boundaries
// although internal neighbors exist (e.g. in a network grid with a kink on the boundary)
template<class GridView, class MapperTraits = DefaultMapperTraits<GridView>>
struct CCTpfaInternalBarrierGridGeometryTraits
: public CCTpfaDefaultGridGeometryTraits<GridView, MapperTraits>
{
    using Element = typename GridView::template Codim<0>::Entity;
    using Intersection = typename GridView::Intersection;

    // usually we would precompute some flags here but for this tests we just check the position
    // otherwise consider passing a mapper (to get the element index) and a map to the constructor.
    explicit CCTpfaInternalBarrierGridGeometryTraits(const bool enableInternalBarrier)
    : enableInternalBarrier_(enableInternalBarrier) {}

    //! Whether an intersection has a neighbor
    bool neighbor(const Element& element, const Intersection& intersection) const
    {
        if (intersection.neighbor())
            return false;

        using std::abs;
        const auto center = intersection.geometry().center();
        if (enableInternalBarrier_ && abs(intersection.centerUnitOuterNormal()[1]) > 0.5 && center[0] < 0.8 && center[0] > 0.2 && center[1] < 0.55 && center[1] > 0.45)
            return false;

        return true;
    }

    //! Whether an intersection is at a boundary
    bool boundary(const Element& element, const Intersection& intersection) const
    {
        if (intersection.boundary())
            return true;

        using std::abs;
        const auto center = intersection.geometry().center();
        if (enableInternalBarrier_ && abs(intersection.centerUnitOuterNormal()[1]) > 0.5 && center[0] < 0.8 && center[0] > 0.2 && center[1] < 0.55 && center[1] > 0.45)
        {
            std::cout << "internal barrier at " << center << std::endl;
            return true;
        }

        return false;
    }

private:
    bool enableInternalBarrier_;
};

} // end namespace Dumux

namespace Dumux::Properties {
// Create new type tags
namespace TTag {
struct OnePInternalBarrier {};
struct OnePInternalBarrierTpfa
{
    using InheritsFrom = std::tuple<OnePInternalBarrier, OnePIncompressibleTpfa>;
    using EnableGridGeometryCache = std::true_type;
};
} // end namespace TTag

template<class TypeTag>
struct GridGeometry<TypeTag, TTag::OnePInternalBarrierTpfa>
{
private:
    static constexpr bool enableCache = getPropValue<TypeTag, Properties::EnableGridGeometryCache>();
    using GridView = typename GetPropType<TypeTag, Properties::Grid>::LeafGridView;
    using Traits = CCTpfaInternalBarrierGridGeometryTraits<GridView>;
public:
    using type = CCTpfaFVGridGeometry<GridView, enableCache, Traits>;
};

} // end namespace Dumux::Properties

namespace Dumux {

template<class T> struct GridGeometryFactory;

template<>
struct GridGeometryFactory<Properties::TTag::OnePInternalBarrierTpfa>
{
    using TypeTag = Properties::TTag::OnePInternalBarrierTpfa;
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using GridView = typename GridGeometry::GridView;

    std::unique_ptr<GridGeometry> create(const GridView& gridView) const
    {
        const bool enableInternalBarrier = getParam<bool>("Problem.InternalBarrier", true);
        CCTpfaInternalBarrierGridGeometryTraits<GridView> strategy(enableInternalBarrier);
        return std::make_unique<GridGeometry>(gridView, strategy);
    }
};

} // end namespace Dumux

#endif
