// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup Discretization
 * \brief Test the interpolation point data of element variables for all combinations of storing
 *        the variables and the interpolation point data for the whole grid
 *
 * Data that are not stored for the whole grid are computed on the first access after binding.
 * All combinations have to give the data stored for the whole grid, and shape values and
 * gradients have to reproduce linear functions exactly.
 */
#include <config.h>

#include <array>
#include <atomic>
#include <cmath>
#include <cstddef>
#include <iostream>
#include <memory>
#include <string>
#include <type_traits>
#include <vector>

#include <dune/common/exceptions.hh>
#include <dune/common/fvector.hh>
#include <dune/grid/yaspgrid.hh>
#include <dune/istl/bvector.hh>

#include <dumux/common/initialize.hh>
#include <dumux/common/typetraits/localdofs_.hh>
#include <dumux/discretization/box/fvgridgeometry.hh>
#include <dumux/discretization/pq2/fvgridgeometry.hh>
#include <dumux/discretization/pq1/fegriddiscretization.hh>
#include <dumux/discretization/pq2/fegriddiscretization.hh>
#include <dumux/discretization/cvfe/quadraturerules.hh>
#include <dumux/discretization/cvfe/interpolationpointdata.hh>
#include <dumux/discretization/cvfe/gridvariablescache_.hh>
#include <dumux/discretization/cvfe/hybrid/gridvariablescache.hh>
#include <dumux/discretization/fem/gridvariablescache.hh>

namespace Dumux {

struct TestProblem {};

struct TestVariables
{
    template<class... Args>
    void update(const Args&...) {}
};

template<class GridDiscretization>
class CountingInterpolationPointData : public CVFE::LocalBasisInterpolationPointData<GridDiscretization>
{
    using ParentType = CVFE::LocalBasisInterpolationPointData<GridDiscretization>;
public:
    template<class... Args>
    void update(const Args&... args)
    {
        ++numUpdates;
        ParentType::update(args...);
    }

    static inline std::atomic<std::size_t> numUpdates = 0;
};

//! The grid variables cache of control-volume schemes and where it provides interpolation point data
struct ControlVolume
{
    template<class InterpolationPointData, bool enableCache, bool enableIpDataCache>
    using Cache = Experimental::CVFE::CVFEGridVariablesCache<
        Experimental::CVFE::CVFEDefaultGridVariablesCacheTraits<TestProblem, TestVariables, InterpolationPointData>, enableCache, enableIpDataCache>;
    static constexpr bool hasScvfData = true;
    static constexpr bool hasElementData = false;
};

//! The grid variables cache of hybrid schemes and where it provides interpolation point data
struct Hybrid
{
    template<class InterpolationPointData, bool enableCache, bool enableIpDataCache>
    using Cache = Experimental::CVFE::HybridCVFEGridVariablesCache<
        Experimental::CVFE::HybridCVFEDefaultGridVariablesCacheTraits<TestProblem, TestVariables, InterpolationPointData>, enableCache, enableIpDataCache>;
    static constexpr bool hasScvfData = true;
    static constexpr bool hasElementData = true;
};

//! The grid variables cache of finite element schemes and where it provides interpolation point data
struct FiniteElement
{
    template<class InterpolationPointData, bool enableCache, bool enableIpDataCache>
    using Cache = Experimental::FEGridVariablesCache<
        Experimental::FEDefaultGridVariablesCacheTraits<TestProblem, TestVariables, InterpolationPointData>, enableCache, enableIpDataCache>;
    static constexpr bool hasScvfData = false;
    static constexpr bool hasElementData = true;
};

//! Call f for each interpolation point at which element variables of the given scheme provide data
template<class Scheme, class ElementDiscretization, class F>
void forEachCachedInterpolationPoint(const ElementDiscretization& elemDisc, F&& f)
{
    if constexpr (Scheme::hasScvfData)
        for (const auto& scvf : scvfs(elemDisc))
            for (const auto& qpData : CVFE::quadratureRule(elemDisc, scvf))
                f(qpData.ipData());

    if constexpr (Scheme::hasElementData)
    {
        for (const auto& qpData : CVFE::quadratureRule(elemDisc, elemDisc.element()))
            f(qpData.ipData());

        for (const auto& boundaryFace : boundaryFaces(elemDisc))
            for (const auto& qpData : CVFE::quadratureRule(elemDisc, boundaryFace))
                f(qpData.ipData());
    }
}

template<class Data, class ReferenceData>
void checkIdentical(const Data& data, const ReferenceData& reference, std::size_t numLocalDofs)
{
    for (std::size_t i = 0; i < numLocalDofs; ++i)
        if (data.shapeValues()[i] != reference.shapeValues()[i] || data.gradN(i) != reference.gradN(i))
            DUNE_THROW(Dune::Exception, "Interpolation point data differ from the data cached on the grid");
}

//! Shape values and gradients interpolate the dof positions to the point and its gradient to the identity
template<class ElementDiscretization, class Data, class GlobalPosition>
void checkLinearReproduction(const ElementDiscretization& elemDisc, const Data& data, const GlobalPosition& ipPosition)
{
    constexpr int dimWorld = GlobalPosition::dimension;
    constexpr double eps = 1e-12;

    double sumShapeValues = 0.0;
    GlobalPosition position(0.0);
    GlobalPosition sumGradients(0.0);
    std::array<GlobalPosition, dimWorld> gradPosition;
    gradPosition.fill(GlobalPosition(0.0));

    for (const auto& localDof : localDofs(elemDisc))
    {
        const auto i = localDof.index();
        const auto dofPosition = ipData(elemDisc, localDof).global();
        const auto shapeValue = data.shapeValues()[i][0];
        sumShapeValues += shapeValue;
        position.axpy(shapeValue, dofPosition);
        sumGradients += data.gradN(i);
        for (int k = 0; k < dimWorld; ++k)
            gradPosition[k].axpy(dofPosition[k], data.gradN(i));
    }

    if (std::abs(sumShapeValues - 1.0) > eps)
        DUNE_THROW(Dune::Exception, "Shape values sum to " << sumShapeValues);
    if ((position - ipPosition).two_norm() > eps)
        DUNE_THROW(Dune::Exception, "Interpolated position " << position << " instead of " << ipPosition);
    if (sumGradients.two_norm() > eps)
        DUNE_THROW(Dune::Exception, "Shape function gradients sum to " << sumGradients);
    for (int k = 0; k < dimWorld; ++k)
        for (int l = 0; l < dimWorld; ++l)
            if (std::abs(gradPosition[k][l] - (k == l ? 1.0 : 0.0)) > eps)
                DUNE_THROW(Dune::Exception, "Gradient of the interpolated position is not the identity");
}

template<class Scheme, bool enableCache, bool enableIpDataCache, class GridDiscretization>
void checkInterpolationPointDataCache(const GridDiscretization& gridDiscretization, const std::string& name)
{
    using IpData = CountingInterpolationPointData<GridDiscretization>;
    using Cache = typename Scheme::template Cache<IpData, enableCache, enableIpDataCache>;
    using ReferenceCache = typename Scheme::template Cache<CVFE::LocalBasisInterpolationPointData<GridDiscretization>, true, true>;
    static_assert(Cache::interpolationPointDataCachingEnabled == enableIpDataCache);

    TestProblem problem;
    Dune::BlockVector<Dune::FieldVector<double, 1>> sol(gridDiscretization.numDofs());
    sol = 0.0;

    IpData::numUpdates = 0;
    Cache gridCache(problem);
    gridCache.init(gridDiscretization, sol);
    const std::size_t numUpdatesInit = IpData::numUpdates;
    ReferenceCache referenceCache(problem);
    referenceCache.init(gridDiscretization, sol);

    auto elemDisc = localView(gridDiscretization);
    auto elemVars = localView(gridCache);
    auto referenceElemVars = localView(referenceCache);

    std::size_t numChecked = 0;
    for (const auto& element : elements(gridDiscretization.gridView()))
    {
        elemDisc.bind(element);
        referenceElemVars.bind(element, elemDisc, sol);

        // a view moved after binding keeps working
        auto movedElemVars = localView(gridCache).bind(element, elemDisc, sol);

        IpData::numUpdates = 0;
        elemVars.bind(element, elemDisc, sol);
        if (IpData::numUpdates != 0)
            DUNE_THROW(Dune::Exception, name << ": binding computed interpolation point data");

        const auto numLocalDofs = Detail::LocalDofs::numLocalDofs(elemDisc);
        std::size_t numPoints = 0;
        forEachCachedInterpolationPoint<Scheme>(elemDisc, [&](const auto& ipData)
        {
            const auto& data = cache(elemVars, ipData);
            checkIdentical(data, cache(referenceElemVars, ipData), numLocalDofs);
            checkIdentical(cache(movedElemVars, ipData), cache(referenceElemVars, ipData), numLocalDofs);
            checkLinearReproduction(elemDisc, data, ipData.global());
            ++numPoints;
        });

        if (numPoints == 0)
            DUNE_THROW(Dune::Exception, name << ": no interpolation points checked");

        // data stored for the whole grid are computed once, the others once per accessed view
        const std::size_t expectedUpdates = enableIpDataCache ? 0 : 2*numPoints;
        if (IpData::numUpdates != expectedUpdates)
            DUNE_THROW(Dune::Exception, name << ": " << IpData::numUpdates << " updates for "
                                              << numPoints << " points in two views, expected " << expectedUpdates);
        numChecked += numPoints;
    }

    if (numUpdatesInit != (enableIpDataCache ? numChecked : 0))
        DUNE_THROW(Dune::Exception, name << ": " << numUpdatesInit << " updates at initialization for "
                                          << numChecked << " points");

    std::cout << name << " (variables " << (enableCache ? "stored" : "local")
              << ", interpolation point data " << (enableIpDataCache ? "stored" : "local")
              << "): checked " << numChecked << " interpolation points" << std::endl;
}

template<class Scheme, class GridDiscretization>
void checkInterpolationPointDataCache(const GridDiscretization& gridDiscretization, const std::string& name)
{
    checkInterpolationPointDataCache<Scheme, false, false>(gridDiscretization, name);
    checkInterpolationPointDataCache<Scheme, false, true>(gridDiscretization, name);
    checkInterpolationPointDataCache<Scheme, true, false>(gridDiscretization, name);
    checkInterpolationPointDataCache<Scheme, true, true>(gridDiscretization, name);
}

} // end namespace Dumux

int main(int argc, char** argv)
{
    using namespace Dumux;
    initialize(argc, argv);

    // anisotropic cells, so that the shape function gradients depend on the inverse Jacobian
    auto makeGrid = []<int dim>(std::integral_constant<int, dim>)
    {
        using Grid = Dune::YaspGrid<dim, Dune::TensorProductCoordinates<double, dim>>;
        std::array<std::vector<double>, dim> coordinates;
        for (int k = 0; k < dim; ++k)
            coordinates[k] = {0.0, 0.1*(k + 1), 0.5, 1.0 + k};
        return std::make_unique<Grid>(coordinates);
    };

    const auto grid2d = makeGrid(std::integral_constant<int, 2>{});
    const auto grid3d = makeGrid(std::integral_constant<int, 3>{});
    using GridView2d = std::decay_t<decltype(grid2d->leafGridView())>;
    using GridView3d = std::decay_t<decltype(grid3d->leafGridView())>;

    checkInterpolationPointDataCache<ControlVolume>(
        BoxFVGridGeometry<double, GridView2d, false>(grid2d->leafGridView()), "box 2d");
    checkInterpolationPointDataCache<ControlVolume>(
        BoxFVGridGeometry<double, GridView3d, false>(grid3d->leafGridView()), "box 3d");
    checkInterpolationPointDataCache<Hybrid>(
        PQ2FVGridGeometry<double, GridView2d, false>(grid2d->leafGridView()), "pq2 hybrid 2d");
    checkInterpolationPointDataCache<FiniteElement>(
        Experimental::PQ1FEGridDiscretization<double, GridView2d, false>(grid2d->leafGridView()), "pq1 fe 2d");
    checkInterpolationPointDataCache<FiniteElement>(
        Experimental::PQ2FEGridDiscretization<double, GridView3d, false>(grid3d->leafGridView()), "pq2 fe 3d");

    return 0;
}
