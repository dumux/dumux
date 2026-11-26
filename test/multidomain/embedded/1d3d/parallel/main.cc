// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup EmbeddedTests
 * \brief The coupling terms of the embedded coupling on distributed grids equal the sequential ones
 *
 * Every rank sets up the coupled problem on the whole tissue and network grids with a sequential
 * communicator as the reference, and the same problem on the distributed grids. The network is distributed
 * either by the tissue partition or by a partition unrelated to it. For two solution states, the
 * coupling part of the residual (the residual minus the residual without coupling terms) of every owned
 * degree of freedom has to equal the reference at the same position, the coupling terms of the two
 * domains have to cancel over all processes, and the number of point sources over all processes has to
 * equal the sequential number.
 */
#include <config.h>

#include <algorithm>
#include <array>
#include <bitset>
#include <cmath>
#include <iostream>
#include <map>
#include <memory>
#include <span>
#include <string>
#include <vector>

#include <dune/common/exceptions.hh>
#include <dune/common/parallel/mpihelper.hh>
#include <dune/grid/common/gridfactory.hh>
#include <dune/foamgrid/parallel/partition.hh>

#include <dumux/common/initialize.hh>
#include <dumux/common/parameters.hh>
#include <dumux/common/properties.hh>
#include <dumux/discretization/elementscontainingpoint.hh>
#include <dumux/discretization/method.hh>
#include <dumux/geometry/intersectingentities.hh>
#include <dumux/multidomain/fvassembler.hh>
#include <dumux/multidomain/traits.hh>
#include <dumux/multidomain/embedded/circlepoints.hh>
#include <dumux/multidomain/embedded/lowdimpartition.hh>
#include <dumux/parallel/vectorcommdatahandle.hh>

#include "properties.hh"

namespace Dumux {

using BulkTypeTag = Properties::TTag::BULKTYPETAG;
using LowDimTypeTag = Properties::TTag::LOWDIMTYPETAG;

//! a position rounded so that it identifies a degree of freedom in the sequential and the distributed grid
using DofKey = std::array<long long, 3>;

template<class Position>
DofKey dofKey(const Position& x)
{ return {std::llround(x[0]*1e9), std::llround(x[1]*1e9), std::llround(x[2]*1e9)}; }

//! The coupling part of the residual of the owned dofs per domain and its sums over all processes
struct CouplingResidual
{
    std::array<std::map<DofKey, double>, 2> values;
    std::array<double, 2> total = {0.0, 0.0};
    std::array<double, 2> magnitude = {0.0, 0.0};
};

struct CouplingResult
{
    std::vector<CouplingResidual> states;
    std::size_t numBulkPointSources = 0;
    std::size_t numLowDimPointSources = 0;
    //! the bulk point sources the search also finds in elements that do not contain them
    std::size_t numPointSourcesNearFaces = 0;
};

template<class Position>
double bulkPressure(const Position& x, int state)
{ return 1.0 + 0.5*std::sin(2.0*x[0] + state)*std::cos(3.0*x[1]) + 0.25*x[2]*(1.0 + state); }

template<class Position>
double lowDimPressure(const Position& x, int state)
{ return 2.0 + 0.3*std::cos(5.0*x[0] + x[1] - state) - 0.2*x[2]; }

/*!
 * \brief A network of a helix, a branch leaving the unit cube, a straight vessel at generic positions,
 *        a straight vessel on faces and one on edges of the tissue cells
 *
 * Optionally, vessels parallel to the y-axis have quadrature points or points of the circle around them
 * just outside of faces of the tissue cells, within the tolerance of the search for the elements of a point.
 */
template<class Grid>
std::unique_ptr<Grid> makeNetwork()
{
    using Position = Dune::FieldVector<double, 3>;
    Dune::GridFactory<Grid> factory;
    unsigned int numVertices = 0;
    const auto addPolyline = [&](const std::vector<Position>& points, int startVertex)
    {
        int previous = startVertex;
        for (std::size_t i = 0; i < points.size(); ++i)
        {
            factory.insertVertex(points[i]);
            const int current = numVertices++;
            if (previous >= 0)
                factory.insertElement(Dune::GeometryTypes::line, {static_cast<unsigned int>(previous), static_cast<unsigned int>(current)});
            previous = current;
        }
    };

    const int numHelixSegments = 40;
    std::vector<Position> helix;
    for (int i = 0; i <= numHelixSegments; ++i)
    {
        const double t = double(i)/numHelixSegments;
        const double phi = 0.37 + 2.0*M_PI*1.3*t;
        helix.push_back({0.5 + 0.23*std::cos(phi), 0.5 + 0.23*std::sin(phi), 0.06 + 0.88*t});
    }
    addPolyline(helix, -1);

    const auto lerp = [](const Position& a, const Position& b, double t){ auto x = a; x.axpy(t, b - a); return x; };
    std::vector<Position> branch;
    const Position branchEnd{0.93, 0.11, 0.37}, outside{1.07, 0.04, 0.29};
    for (int i = 1; i <= 8; ++i)
        branch.push_back(lerp(helix[20], branchEnd, i/8.0));
    for (int i = 1; i <= 2; ++i)
        branch.push_back(lerp(branchEnd, outside, i/2.0));
    addPolyline(branch, 20);

    std::vector<Position> straight;
    const Position a{0.13, 0.81, 0.07}, b{0.71, 0.27, 0.88};
    for (int i = 0; i <= 15; ++i)
        straight.push_back(lerp(a, b, i/15.0));
    addPolyline(straight, -1);

    // a vessel in the plane x = 0.5, which contains faces of the tissue cells, so its pieces lie in two cells
    std::vector<Position> onFaces;
    const Position c{0.5, 0.17, 0.23}, d{0.5, 0.83, 0.71};
    for (int i = 0; i <= 11; ++i)
        onFaces.push_back(lerp(c, d, i/11.0));
    addPolyline(onFaces, -1);

    // a vessel on the edge x = y = 0.5 of the tissue cells, so its pieces lie in four cells
    std::vector<Position> onEdges;
    const Position e{0.5, 0.5, 0.13}, f{0.5, 0.5, 0.87};
    for (int i = 0; i <= 9; ++i)
        onEdges.push_back(lerp(e, f, i/9.0));
    addPolyline(onEdges, -1);

    if (getParam<bool>("Network.NearFaceVessels", false))
    {
        const auto cells = getParamFromGroup<std::array<int, 3>>("Tissue", "Grid.Cells");
        const auto numCirclePoints = getParam<int>("MixedDimension.NumCircleSegments");
        const auto circle = EmbeddedCoupling::circlePoints(Position(0.0), Position{0.0, 1.0, 0.0}, 1.0, numCirclePoints);
        using SpatialParams = GetPropType<Properties::TTag::LOWDIMTYPETAG, Properties::SpatialParams>;
        for (const int normal : {0, 2})
        {
            const int other = 2 - normal;
            const auto circleOffset = std::max_element(circle.begin(), circle.end(), [&](const auto& x, const auto& y)
            { return std::abs(x[normal]) < std::abs(y[normal]); })->operator[](normal);

            int count = 0;
            for (const double sign : {-1.0, 1.0})
            {
                for (const double fraction : {0.3, 0.6, 0.9})
                {
                    const auto facePosition = 0.5 + sign*fraction*1e-7/cells[normal];
                    for (const bool circleNearFace : {false, true})
                    {
                        Position p;
                        p[other] = (circleNearFace ? 0.633 : 0.133) + 0.0317*count;
                        p[normal] = facePosition;
                        if (circleNearFace)
                            for (int i = 0; i < 100; ++i)
                                p[normal] = facePosition - SpatialParams::radiusAtPos(p)*circleOffset;

                        std::vector<Position> vessel;
                        for (int i = 0; i <= 6; ++i)
                        {
                            p[1] = 0.17 + 0.11*i;
                            vessel.push_back(p);
                        }
                        addPolyline(vessel, -1);
                    }
                    ++count;
                }
            }
        }
    }

    return factory.createGrid();
}

template<class Grid>
inline constexpr bool isTensorProductYaspGrid = false;

template<int dim, class ctype>
inline constexpr bool isTensorProductYaspGrid<Dune::YaspGrid<dim, Dune::TensorProductCoordinates<ctype, dim>>> = true;

/*!
 * \brief The tissue grid, sequential for a self communicator and distributed otherwise
 *
 * The cells are equidistant, or graded symmetrically towards the boundary with a power law if the grid
 * has tensor-product coordinates. Both grids have faces in the plane x = 0.5 for an even number of cells.
 */
template<class Grid>
std::unique_ptr<Grid> makeBulkGrid(Dune::MPIHelper::MPICommunicator mpiComm)
{
    using Position = Dune::FieldVector<double, 3>;
    const auto lowerLeft = getParamFromGroup<Position>("Tissue", "Grid.LowerLeft");
    const auto upperRight = getParamFromGroup<Position>("Tissue", "Grid.UpperRight");
    const auto cells = getParamFromGroup<std::array<int, 3>>("Tissue", "Grid.Cells");
    const auto overlap = getParamFromGroup<int>("Tissue", "Grid.Overlap");
    const typename Grid::Communication comm(mpiComm);

    if constexpr (isTensorProductYaspGrid<Grid>)
    {
        const std::array<double, 3> exponents{1.3, 1.5, 0.8};
        std::array<std::vector<double>, 3> coordinates;
        for (int d = 0; d < 3; ++d)
        {
            const int half = cells[d]/2;
            coordinates[d].resize(cells[d] + 1);
            for (int i = 0; i <= half; ++i)
            {
                const auto t = 0.5*std::pow(double(i)/half, exponents[d]);
                coordinates[d][i] = lowerLeft[d] + t*(upperRight[d] - lowerLeft[d]);
                coordinates[d][cells[d] - i] = upperRight[d] - t*(upperRight[d] - lowerLeft[d]);
            }
        }
        return std::make_unique<Grid>(coordinates, std::bitset<3>(), overlap, comm);
    }
    else
        return std::make_unique<Grid>(lowerLeft, upperRight, cells, std::bitset<3>(), overlap, comm);
}

/*!
 * \brief The coupling part of the residual for two solution states
 * \note Collective over the communication of the bulk grid view.
 */
template<class BulkGridView, class LowDimGridView>
CouplingResult computeCoupling(const BulkGridView& bulkGridView, const LowDimGridView& lowDimGridView)
{
    using Traits = MultiDomainTraits<BulkTypeTag, LowDimTypeTag>;
    constexpr auto bulkIdx = Traits::template SubDomain<0>::Index();
    constexpr auto lowDimIdx = Traits::template SubDomain<1>::Index();

    using BulkGridGeometry = GetPropType<BulkTypeTag, Properties::GridGeometry>;
    using LowDimGridGeometry = GetPropType<LowDimTypeTag, Properties::GridGeometry>;
    auto bulkGridGeometry = std::make_shared<BulkGridGeometry>(bulkGridView);
    auto lowDimGridGeometry = std::make_shared<LowDimGridGeometry>(lowDimGridView);

    using CouplingManager = GetPropType<BulkTypeTag, Properties::CouplingManager>;
    auto couplingManager = std::make_shared<CouplingManager>(bulkGridGeometry, lowDimGridGeometry);
    couplingManager->setLowDimElementData(1, [lowDimGridGeometry](std::size_t eIdx, std::span<double> values)
    { values[0] = 1.0 + lowDimGridGeometry->element(eIdx).geometry().center()[1]; });

    using BulkProblem = GetPropType<BulkTypeTag, Properties::Problem>;
    using LowDimProblem = GetPropType<LowDimTypeTag, Properties::Problem>;
    auto bulkProblem = std::make_shared<BulkProblem>(bulkGridGeometry, couplingManager);
    auto lowDimProblem = std::make_shared<LowDimProblem>(lowDimGridGeometry, couplingManager);

    typename Traits::SolutionVector sol;
    sol[bulkIdx].resize(bulkGridGeometry->numDofs());
    sol[lowDimIdx].resize(lowDimGridGeometry->numDofs());
    const auto setSolution = [&](int state)
    {
        const auto set = [&](const auto& gridGeometry, auto& x, auto&& pressure)
        {
            auto fvGeometry = localView(gridGeometry);
            for (const auto& element : elements(gridGeometry.gridView()))
            {
                fvGeometry.bindElement(element);
                for (const auto& scv : scvs(fvGeometry))
                    x[scv.dofIndex()] = pressure(scv.dofPosition(), state);
            }
        };
        set(*bulkGridGeometry, sol[bulkIdx], [](const auto& x, int s){ return bulkPressure(x, s); });
        set(*lowDimGridGeometry, sol[lowDimIdx], [](const auto& x, int s){ return lowDimPressure(x, s); });
    };
    setSolution(0);

    couplingManager->init(bulkProblem, lowDimProblem, sol);
    bulkProblem->computePointSourceMap();
    lowDimProblem->computePointSourceMap();

    using BulkGridVariables = GetPropType<BulkTypeTag, Properties::GridVariables>;
    using LowDimGridVariables = GetPropType<LowDimTypeTag, Properties::GridVariables>;
    auto bulkGridVariables = std::make_shared<BulkGridVariables>(bulkProblem, bulkGridGeometry);
    auto lowDimGridVariables = std::make_shared<LowDimGridVariables>(lowDimProblem, lowDimGridGeometry);
    bulkGridVariables->init(sol[bulkIdx]);
    lowDimGridVariables->init(sol[lowDimIdx]);

    using Assembler = MultiDomainFVAssembler<Traits, CouplingManager, DiffMethod::numeric>;
    auto assembler = std::make_shared<Assembler>(std::make_tuple(bulkProblem, lowDimProblem),
                                                 std::make_tuple(bulkGridGeometry, lowDimGridGeometry),
                                                 std::make_tuple(bulkGridVariables, lowDimGridVariables),
                                                 couplingManager);

    const auto& comm = bulkGridView.comm();
    CouplingResult result;
    result.numBulkPointSources = comm.sum(couplingManager->bulkPointSources().size());
    result.numLowDimPointSources = comm.sum(couplingManager->lowDimPointSources().size());
    for (const auto& source : couplingManager->bulkPointSources())
        if (intersectingEntities(source.position(), bulkGridGeometry->boundingBoxTree()).size()
            != elementsContainingPoint(*bulkGridGeometry, source.position()).size())
            ++result.numPointSourcesNearFaces;

    for (int state = 0; state < 2; ++state)
    {
        if (state > 0)
        {
            setSolution(state);
            couplingManager->updateSolution(sol);
        }

        typename Assembler::ResidualType residual, uncoupledResidual;
        residual[bulkIdx].resize(bulkGridGeometry->numDofs());
        residual[lowDimIdx].resize(lowDimGridGeometry->numDofs());
        uncoupledResidual = residual;

        bulkProblem->enableCoupling(true);
        lowDimProblem->enableCoupling(true);
        assembler->assembleResidual(residual, sol);
        bulkProblem->enableCoupling(false);
        lowDimProblem->enableCoupling(false);
        assembler->assembleResidual(uncoupledResidual, sol);
        residual -= uncoupledResidual;

        auto& coupling = result.states.emplace_back();
        const auto collect = [&](auto domainIdx, const auto& gridGeometry)
        {
            auto& r = residual[domainIdx];
            const auto gridView = gridGeometry.gridView();
            for (const auto& value : r)
            {
                coupling.total[domainIdx] += value[0];
                coupling.magnitude[domainIdx] += std::abs(value[0]);
            }
            coupling.total[domainIdx] = comm.sum(coupling.total[domainIdx]);
            coupling.magnitude[domainIdx] = comm.sum(coupling.magnitude[domainIdx]);

            using GridGeometry = std::decay_t<decltype(gridGeometry)>;
            if constexpr (GridGeometry::discMethod == DiscretizationMethods::box)
            {
                static constexpr int dim = GridGeometry::GridView::dimension;
                using Vector = std::decay_t<decltype(r)>;
                VectorCommDataHandleSum<typename GridGeometry::VertexMapper, Vector, dim> handle(gridGeometry.vertexMapper(), r);
                if (gridView.comm().size() > 1)
                    gridView.communicate(handle, Dune::InteriorBorder_InteriorBorder_Interface, Dune::ForwardCommunication);

                for (const auto& vertex : vertices(gridView))
                    if (vertex.partitionType() == Dune::InteriorEntity || vertex.partitionType() == Dune::BorderEntity)
                        coupling.values[domainIdx][dofKey(vertex.geometry().center())] = r[gridGeometry.vertexMapper().index(vertex)][0];
            }
            else
            {
                for (const auto& element : elements(gridView, Dune::Partitions::interior))
                    coupling.values[domainIdx][dofKey(element.geometry().center())] = r[gridGeometry.elementMapper().index(element)][0];
            }
        };
        collect(bulkIdx, *bulkGridGeometry);
        collect(lowDimIdx, *lowDimGridGeometry);
    }

    return result;
}

} // end namespace Dumux

int main(int argc, char** argv)
{
    using namespace Dumux;

    initialize(argc, argv);
    const auto& mpiHelper = Dune::MPIHelper::instance();
    Parameters::init(argc, argv);

    using BulkGrid = GetPropType<BulkTypeTag, Properties::Grid>;
    using LowDimGrid = GetPropType<LowDimTypeTag, Properties::Grid>;
    const Dune::Communication<Dune::MPIHelper::MPICommunicator> comm(mpiHelper.getCommunicator());

    // the reference on the whole grids with a sequential communicator
    const auto fullBulkGrid = makeBulkGrid<BulkGrid>(mpiHelper.getLocalCommunicator());
    const auto fullLowDimGrid = makeNetwork<LowDimGrid>();
    const auto reference = computeCoupling(fullBulkGrid->leafGridView(), fullLowDimGrid->leafGridView());

    // the distributed grids
    const auto bulkGrid = makeBulkGrid<BulkGrid>(mpiHelper.getCommunicator());
    const auto bulkGridView = bulkGrid->leafGridView();

    auto lowDimGrid = makeNetwork<LowDimGrid>();
    lowDimGrid->setCommunicator(mpiHelper.getCommunicator());
    const auto partition = getParam<std::string>("Network.Partition");
    const auto owners = partition == "Tissue"
        ? lowDimPartition(bulkGridView, lowDimGrid->leafGridView())
        : Dune::FoamGridParallel::recursiveCoordinateBisection(lowDimGrid->leafGridView(), comm.size());
    lowDimGrid->loadBalance(owners, getParamFromGroup<int>("Network", "Grid.Overlap"));

    const auto distributed = computeCoupling(bulkGridView, lowDimGrid->leafGridView());

    int failures = 0;
    const auto fail = [&](const auto& message){ std::cerr << "rank " << comm.rank() << ": " << message << std::endl; ++failures; };

    if (distributed.numBulkPointSources != reference.numBulkPointSources)
        fail("bulk point sources " + std::to_string(distributed.numBulkPointSources) + " instead of " + std::to_string(reference.numBulkPointSources));
    if (distributed.numLowDimPointSources != reference.numLowDimPointSources)
        fail("lower-dimensional point sources " + std::to_string(distributed.numLowDimPointSources) + " instead of " + std::to_string(reference.numLowDimPointSources));

    using BulkGridGeometry = GetPropType<BulkTypeTag, Properties::GridGeometry>;
    if constexpr (DiscretizationMethods::isCVFE<typename BulkGridGeometry::DiscretizationMethod>)
        if (getParam<bool>("Network.NearFaceVessels", false) && reference.numPointSourcesNearFaces == 0)
            fail("no point source the search also finds in an element that does not contain it");

    static constexpr double tolerance = 1e-12;
    const std::array<std::string, 2> domainNames{"tissue", "network"};
    for (std::size_t state = 0; state < reference.states.size(); ++state)
    {
        const auto& ref = reference.states[state];
        const auto& dist = distributed.states[state];

        for (const auto* r : {&ref, &dist})
        {
            const auto imbalance = std::abs(r->total[0] + r->total[1]);
            const auto scale = r->magnitude[0] + r->magnitude[1];
            if (!(imbalance <= tolerance*scale))
                fail(std::string(r == &ref ? "sequential" : "distributed") + " state " + std::to_string(state)
                     + ": the coupling terms do not cancel, imbalance " + std::to_string(imbalance/scale));
        }

        for (int domain = 0; domain < 2; ++domain)
        {
            double maxRef = 0.0;
            for (const auto& [key, value] : ref.values[domain])
                maxRef = std::max(maxRef, std::abs(value));
            if (maxRef == 0.0)
                fail(domainNames[domain] + ": no coupling terms in the reference");

            double maxDiff = 0.0;
            for (const auto& [key, value] : dist.values[domain])
            {
                const auto it = ref.values[domain].find(key);
                if (it == ref.values[domain].end())
                {
                    fail(domainNames[domain] + ": an owned dof is not in the reference");
                    continue;
                }
                if (!std::isfinite(value))
                    fail(domainNames[domain] + ": a coupling term is not finite");
                maxDiff = std::max(maxDiff, std::abs(value - it->second));
            }

            const auto relDiff = comm.max(maxDiff)/maxRef;
            if (comm.rank() == 0)
                std::cout << "state " << state << ", " << domainNames[domain]
                          << ": max |distributed - sequential| / max |sequential| = " << relDiff << std::endl;
            if (relDiff > tolerance)
                fail(domainNames[domain] + " state " + std::to_string(state) + ": coupling differs from the sequential one by " + std::to_string(relDiff));
        }
    }

    failures = comm.sum(failures);
    if (comm.rank() == 0)
        std::cout << (failures == 0 ? "PASSED" : "FAILED") << " on " << comm.size() << " processes ("
                  << reference.numLowDimPointSources << " point sources)" << std::endl;

    return failures == 0 ? 0 : 1;
}
