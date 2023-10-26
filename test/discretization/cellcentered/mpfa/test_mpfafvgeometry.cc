// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \brief Test for the mpfa finite volume element geometry.
 */
#include <config.h>

#include <iostream>
#include <vector>
#include <utility>
#include <algorithm>

#include <dune/common/fvector.hh>
#include <dune/common/float_cmp.hh>
#include <dune/grid/yaspgrid.hh>

#include <dumux/common/initialize.hh>
#include <dumux/common/properties.hh>
#include <dumux/common/parameters.hh>
#include <dumux/discretization/ccmpfa.hh>

#ifndef ENABLECACHING
#define ENABLECACHING 0
#endif

namespace Dumux::Properties {
namespace TTag { struct TestTag { using InheritsFrom = std::tuple<CCMpfaModel>; }; }

template<typename TypeTag>
struct Grid<TypeTag, TTag::TestTag> { using type = Dune::YaspGrid<2>; };

template<typename TypeTag>
struct Scalar<TypeTag, TTag::TestTag> { using type = double; };

template<typename TypeTag>
struct EnableGridGeometryCache<TypeTag, TTag::TestTag> {
    static constexpr bool value = ENABLECACHING;
};

} // namespace Dumux::Properties

template<class Range>
auto size(const Range& r)
{ return std::distance(std::begin(r), std::end(r)); }

template<class FVElementGeometry>
auto numBoundaryScvfs(const FVElementGeometry& fvGeometry)
{
    unsigned int count = 0;
    for (const auto& scvf : scvfs(fvGeometry))
        if (scvf.boundary())
            count++;
    return count;
}

template<typename Geometry>
auto popped(std::vector<typename Geometry::GlobalCoordinate> expected, const Geometry& geo)
{
    using T = std::vector<typename Geometry::GlobalCoordinate>;
    for (int c = 0; c < geo.corners(); ++c)
    {
        auto it = std::find_if(
            expected.begin(),
            expected.end(),
            [corner=geo.corner(c)] (const auto& exp) {
                return Dune::FloatCmp::eq(corner, exp);
        });
        if (it == expected.end())
        {
            std::cout << "Could not find corner " << geo.corner(c) << std::endl;
            return std::optional<T>{};
        }
        expected.erase(it);
    }
    return std::optional<T>{std::move(expected)};
}

int main (int argc, char *argv[])
{
    using namespace Dumux;

    initialize(argc, argv);
    Parameters::init(argc, argv);

    std::cout << "Checking the FVGeometries, SCVs and SCV faces" << std::endl;

    using TypeTag = Properties::TTag::TestTag;
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using FVElementGeometry = typename GridGeometry::LocalView;
    using SubControlVolumeFace = typename FVElementGeometry::SubControlVolumeFace;

    GetPropType<TypeTag, Properties::Grid> grid{{1.0, 1.0}, {1, 1}};
    auto leafGridView = grid.leafGridView();
    GridGeometry gridGeometry(leafGridView);

    for (const auto& element : elements(leafGridView))
    {
        auto fvGeometry = localView(gridGeometry);

        if (fvGeometry.isBound())
            DUNE_THROW(Dune::Exception, "Local view should not be bound at this point");

        fvGeometry.bind(element);
        if (!fvGeometry.isBound())
            DUNE_THROW(Dune::Exception, "Local view should be bound at this point");

        const auto eIdx = gridGeometry.elementMapper().index(element);
        const auto eIdxBound = gridGeometry.elementMapper().index(fvGeometry.element());
        if (eIdx != eIdxBound)
            DUNE_THROW(Dune::Exception, "Bound element index does not match");

        if (size(scvs(fvGeometry)) != 1) DUNE_THROW(Dune::InvalidStateException, "Unexpected number of scvs");
        if (size(scvfs(fvGeometry)) != 8) DUNE_THROW(Dune::InvalidStateException, "Unexpected number of scvfs");
        if (numBoundaryScvfs(fvGeometry) != 8) DUNE_THROW(Dune::InvalidStateException, "Unexpected number of boundary scvfs");

        for (const auto& scv : scvs(fvGeometry))
        {
            const auto geo = fvGeometry.geometry(scv);
            if (!Dune::FloatCmp::eq(geo.volume(), 1.0))
                DUNE_THROW(Dune::InvalidStateException, "Unexpected scv volume");
            const auto result = popped({
                {0.0, 0.0},
                {1.0, 0.0},
                {1.0, 1.0},
                {0.0, 1.0}
            }, geo);
            if (!result.has_value() || !result.value().empty())
                DUNE_THROW(Dune::InvalidStateException, "Unexpected scv corners");
        }

        std::vector<typename SubControlVolumeFace::GlobalPosition> expected{
            {0.0, 0.0}, {0.5, 0.0},
            {0.5, 0.0}, {1.0, 0.0},
            {1.0, 0.0}, {1.0, 0.5},
            {1.0, 0.5}, {1.0, 1.0},
            {1.0, 1.0}, {0.5, 1.0},
            {0.5, 1.0}, {0.0, 1.0},
            {0.0, 1.0}, {0.0, 0.5},
            {0.0, 0.5}, {0.0, 0.0}
        };
        for (const auto& scvf : scvfs(fvGeometry))
        {
            const auto geo = fvGeometry.geometry(scvf);
            if (!Dune::FloatCmp::eq(geo.volume(), 0.5))
                DUNE_THROW(Dune::InvalidStateException, "Unexpected scvf area " << geo.volume());
            const auto shrunk = popped(expected, geo);
            if (!shrunk.has_value())
                DUNE_THROW(Dune::InvalidStateException, "Unexpected scvf corners");
            expected = shrunk.value();
        }
        if (!expected.empty())
            DUNE_THROW(Dune::InvalidStateException, "Could not find all expected scvf corners");
    }
}
