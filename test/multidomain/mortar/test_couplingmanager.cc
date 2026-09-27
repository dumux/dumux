// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \brief Tests the evaluation of subdomain trace couplings: the pointwise lookup of trace
 *        data on elements coupled to several trace cells, for data given per trace cell and
 *        per trace vertex, and the refusal of data in the wrong layout.
 */
#include <config.h>

#include <cmath>
#include <memory>
#include <iostream>
#include <string>
#include <type_traits>

#include <dune/common/exceptions.hh>
#include <dune/common/fvector.hh>
#include <dune/grid/yaspgrid.hh>
#include <dune/foamgrid/foamgrid.hh>
#include <dune/istl/bvector.hh>

#include <dumux/common/initialize.hh>
#include <dumux/common/concepts/field_.hh>
#include <dumux/common/concepts/mortarcouplingmanager_.hh>
#include <dumux/discretization/box/fvgridgeometry.hh>
#include <dumux/multidomain/mortar/trace.hh>
#include <dumux/multidomain/mortar/couplingmanager.hh>

namespace {

//! a constant of seven, integrated over each face
struct ConstantField
{
    template<typename Element>
    auto bind(const Element&) const
    { return [] (const auto& scvf) { return Dune::FieldVector<double, 1>{7.0*scvf.area()}; }; }
};

} // end anonymous namespace

int main(int argc, char** argv)
{
    using namespace Dumux;
    initialize(argc, argv);

    using Grid = Dune::YaspGrid<2>;
    using MortarGrid = Dune::FoamGrid<1, 2>;
    using GridGeometry = BoxFVGridGeometry<double, typename Grid::LeafGridView>;
    using MortarSolution = Dune::BlockVector<Dune::FieldVector<double, 1>>;
    using CouplingManager = Mortar::SubDomainCouplingManager<GridGeometry, MortarGrid, MortarSolution>;
    using GlobalPosition = Dune::FieldVector<double, 2>;

    // a bulk subdomain has a facet trace and composes the manager for it, which models the
    // interface a subdomain problem is written against; a grid geometry does not
    static_assert(std::is_same_v<CouplingManager, Mortar::FacetTraceCouplingManager<GridGeometry, MortarGrid, MortarSolution>>);
    static_assert(Concept::MortarSubDomainCouplingManager<CouplingManager>);
    static_assert(!Concept::MortarSubDomainCouplingManager<GridGeometry>);

    // what a trace is assembled from has to be a field over the faces of the subdomain
    using Element = typename GridGeometry::GridView::template Codim<0>::Entity;
    static_assert(Concept::FaceField<ConstantField, Element, typename GridGeometry::SubControlVolumeFace>);
    static_assert(!Concept::FaceField<GridGeometry, Element, typename GridGeometry::SubControlVolumeFace>);

    // an L-shaped trace along the two boundaries through the origin: the corner element
    // couples to trace cells on both of its boundary facets, so a pointwise lookup must
    // pick the cell containing the queried position, not merely a nearby one
    Grid grid{{1.0, 1.0}, {4, 4}};
    auto gridGeometry = std::make_shared<GridGeometry>(grid.leafGridView());

    auto trace = std::make_shared<Mortar::FacetTrace<Grid, MortarGrid>>(grid, [] (const auto&, const auto& is) {
        const auto c = is.geometry().center();
        return c[1] < 1e-10 || c[0] < 1e-10;
    });
    const auto& traceGridView = trace->gridView();
    if (trace->size() != 8)
        DUNE_THROW(Dune::InvalidStateException, "Expected the eight boundary facets of the L-shaped trace");

    CouplingManager couplingManager{gridGeometry};
    couplingManager.registerTrace(trace, 0);

    // a box subdomain takes its trace data per trace vertex by default; this test starts
    // with data per trace cell
    if (couplingManager.traceDataOrder() != 1)
        DUNE_THROW(Dune::InvalidStateException, "A box subdomain should take trace data per vertex by default");
    couplingManager.setTraceDataOrder(0);
    if (couplingManager.numTraceDofs(0) != traceGridView.size(0))
        DUNE_THROW(Dune::InvalidStateException, "Data per trace cell has one entry per cell");

    // data in the wrong layout is refused
    {
        bool threw = false;
        try { MortarSolution wrong(traceGridView.size(1)); wrong = 0.0; couplingManager.setTraceVariables(0, wrong); }
        catch (const Dune::InvalidStateException&) { threw = true; }
        if (!threw)
            DUNE_THROW(Dune::InvalidStateException, "Data per vertex was accepted for a manager reading data per cell");
    }

    // distinct value per trace cell, distinguishable between the two legs of the L
    MortarSolution data(traceGridView.size(0));
    Dune::MultipleCodimMultipleGeomTypeMapper<typename MortarGrid::LeafGridView> traceMapper(
        traceGridView, Dune::mcmgElementLayout());
    for (const auto& element : elements(traceGridView))
    {
        const auto c = element.geometry().center();
        data[traceMapper.index(element)] = c[1] < 1e-10 ? 100.0 + std::floor(c[0]*4.0)
                                                        : 200.0 + std::floor(c[1]*4.0);
    }
    couplingManager.setTraceVariables(0, data);

    const auto elementAt = [&] (const GlobalPosition& pos) {
        for (const auto& element : elements(gridGeometry->gridView()))
        {
            const auto geo = element.geometry();
            bool inside = true;
            for (int d = 0; d < 2; ++d)
                if (pos[d] < geo.corner(0)[d] - 1e-12 || pos[d] > geo.corner(3)[d] + 1e-12)
                    inside = false;
            if (inside)
                return element;
        }
        DUNE_THROW(Dune::InvalidStateException, "No element at " << pos);
    };

    // pointwise lookup: values of the containing trace cell, including positions on the
    // corner element where the two legs of the L would confuse a proximity-only match
    const auto check = [&] (const GlobalPosition& pos, double expected) {
        const auto value = Mortar::traceVariablesAt(couplingManager, elementAt(pos), pos)[0];
        if (std::abs(value - expected) > 1e-12)
            DUNE_THROW(Dune::MathError, "Pointwise lookup at " << pos << ": "
                       << value << ", expected " << expected);
    };
    for (const double x : {0.05, 0.2, 0.45, 0.62, 0.9})
        check({x, 0.0}, 100.0 + std::floor(x*4.0));
    for (const double y : {0.05, 0.2, 0.7})
        check({0.0, y}, 200.0 + std::floor(y*4.0));
    check({0.01, 0.0}, 100.0);
    check({0.0, 0.01}, 200.0);
    std::cout << "Pointwise containment lookups match the containing trace cell" << std::endl;

    // the position query: true on both legs of the L, false on the other boundaries and inside
    for (const double t : {0.05, 0.5, 0.95})
        if (!couplingManager.isCoupledAtPos({t, 0.0}) || !couplingManager.isCoupledAtPos({0.0, t}))
            DUNE_THROW(Dune::InvalidStateException, "Position query misses the trace at " << t);
    for (const GlobalPosition& pos : {GlobalPosition{0.5, 1.0}, GlobalPosition{1.0, 0.5}, GlobalPosition{0.5, 0.5}, GlobalPosition{0.5, 1e-3}})
        if (couplingManager.isCoupledAtPos(pos))
            DUNE_THROW(Dune::InvalidStateException, "Position query finds a trace at " << pos);
    std::cout << "Position query is selective" << std::endl;

    // boundary detection, and the value at a degree of freedom for data per trace cell: the
    // mean over the trace cells containing it, which differ between the legs at the corner
    {
        const auto element = elementAt({0.1, 0.05});
        auto fvGeometry = localView(*gridGeometry);
        fvGeometry.bind(element);
        bool foundCoupled = false, foundUncoupled = false;
        for (const auto& scv : scvs(fvGeometry))
        {
            const bool coupled = Mortar::isOnMortarBoundary(couplingManager, element, scv);
            if (coupled != couplingManager.isCoupled(element, scv))
                DUNE_THROW(Dune::InvalidStateException, "Facade and free function disagree on coupling");
            if (coupled)
            {
                foundCoupled = true;
                const auto value = couplingManager.traceAt(element, scv)[0];
                if (std::abs(value - Mortar::traceVariablesAt(couplingManager, element, scv)[0]) > 1e-12)
                    DUNE_THROW(Dune::MathError, "Facade and free function disagree on the trace value");
                const auto pos = scv.dofPosition();
                const double expected = pos[0] < 1e-10 && pos[1] < 1e-10 ? 150.0
                                      : pos[1] < 1e-10 ? 100.5 : 200.5;
                if (std::abs(value - expected) > 1e-12)
                    DUNE_THROW(Dune::MathError, "Value at the degree of freedom at " << pos << ": " << value << ", expected " << expected);
            }
            else
                foundUncoupled = true;
        }
        if (!foundCoupled || !foundUncoupled)
            DUNE_THROW(Dune::InvalidStateException, "Boundary detection is not selective");
    }
    std::cout << "Mortar-boundary detection is selective" << std::endl;

    // data per trace vertex: the trace function is linear over each trace cell, so a linear
    // field given at the trace vertices is reproduced at any position, at the integration
    // points of the coupled faces and at the coupled degrees of freedom
    {
        couplingManager.setTraceDataOrder(1);
        if (couplingManager.numTraceDofs(0) != traceGridView.size(1))
            DUNE_THROW(Dune::InvalidStateException, "Data per trace vertex has one entry per vertex");
        bool threw = false;
        try { couplingManager.setTraceVariables(0, data); }
        catch (const Dune::InvalidStateException&) { threw = true; }
        if (!threw)
            DUNE_THROW(Dune::InvalidStateException, "Data per cell was accepted for a manager reading data per vertex");

        const auto linear = [] (const GlobalPosition& p) { return 1.0 + 2.0*p[0] + 3.0*p[1]; };
        MortarSolution vertexData(traceGridView.size(1));
        for (const auto& vertex : vertices(traceGridView))
            vertexData[traceGridView.indexSet().index(vertex)] = linear(vertex.geometry().center());
        couplingManager.setTraceVariables(0, vertexData);

        const auto checkLinear = [&] (const GlobalPosition& pos, double value, const std::string& what) {
            if (std::abs(value - linear(pos)) > 1e-12)
                DUNE_THROW(Dune::MathError, what << " at " << pos << ": " << value << ", expected " << linear(pos));
        };
        for (const double x : {0.05, 0.2, 0.45, 0.62, 0.9})
            checkLinear({x, 0.0}, Mortar::traceVariablesAt(couplingManager, elementAt({x, 0.0}), GlobalPosition{x, 0.0})[0], "Pointwise lookup");
        for (const double y : {0.05, 0.2, 0.7})
            checkLinear({0.0, y}, Mortar::traceVariablesAt(couplingManager, elementAt({0.0, y}), GlobalPosition{0.0, y})[0], "Pointwise lookup");

        const auto element = elementAt({0.1, 0.05});
        auto fvGeometry = localView(*gridGeometry);
        fvGeometry.bind(element);
        const auto eIdx = gridGeometry->elementMapper().index(element);
        if (couplingManager.faceCouplingsOf(eIdx).empty())
            DUNE_THROW(Dune::InvalidStateException, "The corner element has coupled faces");
        for (const auto& entry : couplingManager.faceCouplingsOf(eIdx))
        {
            const auto& scvf = fvGeometry.scvf(entry.scvfIndex);
            checkLinear(scvf.ipGlobal(), couplingManager.traceAt(element, scvf)[0], "Face lookup");
        }
        for (const auto& scv : scvs(fvGeometry))
            if (couplingManager.isCoupled(element, scv))
                checkLinear(scv.dofPosition(), couplingManager.traceAt(element, scv)[0], "Lookup at the degree of freedom");
        couplingManager.setTraceDataOrder(0);
        couplingManager.setTraceVariables(0, data);
    }
    std::cout << "Per-vertex trace data is interpolated linearly" << std::endl;

    // area-weighted trace read of a constant field reproduces the constant
    {
        const auto values = Mortar::traceValues<Mortar::TraceEntity::subControlVolumeFace>(couplingManager, 0, ConstantField{});
        for (std::size_t i = 0; i < values.size(); ++i)
            if (std::abs(values[i][0] - 7.0) > 1e-12)
                DUNE_THROW(Dune::MathError, "Trace read of a constant field is " << values[i][0]);
    }
    std::cout << "Trace read reproduces a constant field" << std::endl;

    if (couplingManager.isFloating())
        DUNE_THROW(Dune::InvalidStateException, "Subdomain wrongly detected as floating");

    // the homogeneous flag is coupling state the problem reads, defaulting to off
    if (couplingManager.isHomogeneous())
        DUNE_THROW(Dune::InvalidStateException, "Coupling manager starts out homogeneous");
    couplingManager.setHomogeneous(true);
    if (!couplingManager.isHomogeneous())
        DUNE_THROW(Dune::InvalidStateException, "Homogeneous flag was not taken");
    couplingManager.setHomogeneous(false);
    std::cout << "Homogeneous flag is carried by the coupling manager" << std::endl;

    std::cout << "All coupling-manager checks passed" << std::endl;
    return 0;
}
