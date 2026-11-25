// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \author Timo Koch
 * \ingroup EmbeddedCoupling
 * \brief Coupling manager for low-dimensional domains embedded in the bulk
 *        domain with spatially resolved interface.
 */

#ifndef DUMUX_MULTIDOMAIN_EMBEDDED_COUPLINGMANAGER_1D3D_PROJECTION_HH
#define DUMUX_MULTIDOMAIN_EMBEDDED_COUPLINGMANAGER_1D3D_PROJECTION_HH

#include <vector>
#include <algorithm>
#include <limits>
#include <string>
#include <tuple>
#include <type_traits>
#include <utility>

#include <dune/common/timer.hh>
#include <dune/common/exceptions.hh>
#include <dune/common/fvector.hh>
#include <dune/geometry/quadraturerules.hh>

#include <dumux/common/tag.hh>
#include <dumux/common/properties.hh>
#include <dumux/common/typetraits/griddiscretization.hh>
#include <dumux/common/math.hh>
#include <dumux/common/indextraits.hh>

#include <dumux/geometry/boundingboxtree.hh>
#include <dumux/geometry/distance.hh>
#include <dumux/geometry/refinementquadraturerule.hh>
#include <dumux/geometry/triangulation.hh>
#include <dumux/geometry/distance.hh>
#include <dumux/geometry/grahamconvexhull.hh>

#include <dumux/multidomain/embedded/couplingmanagerbase.hh>
#include <dumux/multidomain/embedded/localrefinementquadrature.hh>

namespace Dumux {

// some implementation details of the projection coupling manager
namespace Detail {

/*!
 * \ingroup EmbeddedCoupling
 * \brief Segment representation of a 1d network grid
 *
 * We use separate segment representation which is fast to iterate over and which
 * can be used to find the unique mapping from points on coupled bulk surface facets
 * to the network centerline (i.e. a 1D dof on the centerline).
 */
template <class GridGeometry>
class SegmentNetwork
{
    using GridView = typename GridGeometry::GridView;
    using GlobalPosition = typename GridView::template Codim<0>::Geometry::GlobalCoordinate;
    using GridIndex = typename IndexTraits<GridView>::GridIndex;

    using BoundingBoxTree = std::decay_t<decltype(std::declval<const GridGeometry&>().boundingBoxTree())>;

    struct Segment
    {
        GlobalPosition a, b;
        double r;
        GridIndex eIdx;
    };

public:
    template <class RadiusFunction>
    SegmentNetwork(const GridGeometry& gg, const RadiusFunction& radiusFunction)
    : tree_(&gg.boundingBoxTree())
    {
        const auto& gridView = gg.gridView();
        segments_.resize(gridView.size(0));
        for (const auto &element : elements(gridView))
        {
            const auto geometry = element.geometry();
            const auto eIdx = gg.elementMapper().index(element);
            const auto radius = radiusFunction(eIdx);

            segments_[eIdx] = Segment{ geometry.corner(0), geometry.corner(1), radius, eIdx };
            maxRadius_ = std::max(maxRadius_, static_cast<double>(radius));
        }

        std::cout << "-- Extracted " << segments_.size() << " segments.\n";
    }

    //! the segment radius
    auto radius(std::size_t eIdx) const
    { return segments_[eIdx].r; }

    //! the segment with index eIdx
    const auto& segment(std::size_t eIdx) const
    { return segments_[eIdx]; }

    //! the segments
    const auto& segments() const
    { return segments_; }

    //! find the segment with the smallest signed distance from its surface to the position (lowest index on ties)
    //! \note the returned distance is the absolute value of that signed distance
    auto findClosestSegmentSurface(const GlobalPosition& globalPos) const
    { return findClosestSegmentSurface(globalPos, closestSegmentSurfaceCandidates(globalPos)); }

    /*!
     * \brief Indices (ascending) of all segments that can be the closest segment surface
     *        for a point within the given distance from the position
     * \note Let c be the segment with the closest centerline (distance d_c, radius r_c). The
     *       segment s that minimizes the signed surface distance has d_s - r_s <= d_c - r_c at the
     *       position. For a point within the margin the distances change by at most the margin.
     *       Hence d_s <= d_c + 2*margin + r_max - r_c at the position.
     */
    std::vector<GridIndex> closestSegmentSurfaceCandidates(const GlobalPosition& globalPos,
                                                           double margin = 0.0) const
    {
        std::vector<GridIndex> candidates;
        if (segments_.empty())
            return candidates;

        using std::sqrt;
        const auto [squaredDist, closestIdx] = closestEntity(globalPos, *tree_);
        auto radius = sqrt(squaredDist) + 2.0*margin + maxRadius_ - segments_[closestIdx].r;

        // slack for the round-off in the distance evaluations
        const auto* rootBox = tree_->getBoundingBoxCoordinates(tree_->numBoundingBoxes() - 1);
        double extent = 0.0;
        for (int i = 0; i < GlobalPosition::dimension; ++i)
            extent = std::max(extent, rootBox[i + GlobalPosition::dimension] - rootBox[i]);
        radius += 1e-9*(radius + extent);

        collectCandidates_(globalPos, tree_->numBoundingBoxes() - 1, radius*radius, candidates);
        std::sort(candidates.begin(), candidates.end());
        return candidates;
    }

    //! same as overload with taking a position but only look at the segments specified by the index vector
    template <class IndexRange>
    auto findClosestSegmentSurface(const GlobalPosition& globalPos,
                                   const IndexRange& segIndices) const
    {
        GridIndex eIdx = 0;
        double minSignedDist = std::numeric_limits<double>::max();

        for (const auto index : segIndices)
        {
            const auto &segment = segments_[index];
            const auto signedDist = computeSignedDistance_(globalPos, segment);
            if (signedDist < minSignedDist)
            {
                minSignedDist = signedDist;
                eIdx = segment.eIdx;
            }
        }

        using std::abs;
        return std::make_tuple(eIdx, abs(minSignedDist));
    }

    // Compute the projection point of p onto the segment connecting a->b
    GlobalPosition projectionPointOnSegment(const GlobalPosition& p, std::size_t idx) const
    {
        const auto& segment = segments_[idx];
        return projectionPointOnSegment_(p, segment.a, segment.b);
    }

private:
    // collect the segments whose bounding box is within the given squared distance of p
    void collectCandidates_(const GlobalPosition& p, std::size_t node, double squaredRadius,
                            std::vector<GridIndex>& candidates) const
    {
        const auto& box = tree_->getBoundingBoxNode(node);
        if (squaredDistancePointBoundingBox(p, tree_->getBoundingBoxCoordinates(node)) > squaredRadius)
            return;

        if (tree_->isLeaf(box, node))
            candidates.push_back(box.child1);
        else
        {
            collectCandidates_(p, box.child0, squaredRadius, candidates);
            collectCandidates_(p, box.child1, squaredRadius, candidates);
        }
    }

    // Compute the projection of p onto the segment
    // returns the closest end point if the orthogonal projection onto the line implied by the segment
    // is outside the segment
    GlobalPosition projectionPointOnSegment_(const GlobalPosition& p, const GlobalPosition& a, const GlobalPosition& b) const
    {
        const auto v = b - a;
        const auto w = p - a;

        const auto proj1 = v*w;
        if (proj1 <= 0.0)
            return a;

        const auto proj2 = v.two_norm2();
        if (proj2 <= proj1)
            return b;

        const auto t = proj1 / proj2;
        auto x = a;
        x.axpy(t, v);
        return x;
    }

    // Compute the signed (!) distance to a cylinder segment surface and the projection point onto the centerline
    auto computeSignedDistance_(const GlobalPosition& p, const Segment& segment) const
    {
        const auto projPoint = projectionPointOnSegment_(p, segment.a, segment.b);
        return (p-projPoint).two_norm()-segment.r;
    }

    const BoundingBoxTree* tree_;
    std::vector<Segment> segments_;
    double maxRadius_ = 0.0;
};

/*!
 * \ingroup EmbeddedCoupling
 * \brief Get the closest segment for a given surface point
 * \note The point does not necessarily be exactly on the virtual surface
 *       implied by the networks radius function
 * \note Only the segments given by the hints are considered
 */
template<class Network, class HintContainer>
class NetworkIndicatorFunction
{
public:
    NetworkIndicatorFunction(const Network& n, const HintContainer& hints)
    : network_(n)
    , hints_(hints)
    {}

    template<class Pos>
    auto operator()(const Pos& pos) const
    {
        const auto& [closestSegmentIdx, _] = network_.findClosestSegmentSurface(pos, hints_);
        return closestSegmentIdx;
    }

private:
    const Network& network_;
    const HintContainer& hints_;
};

/*!
 * \ingroup EmbeddedCoupling
 * \brief The largest distance of a corner of the geometry to the given position
 */
template<class Geometry, class Position>
double maxCornerDistance(const Geometry& geometry, const Position& position)
{
    double maxDistance = 0.0;
    for (int i = 0; i < geometry.corners(); ++i)
        maxDistance = std::max(maxDistance, (geometry.corner(i) - position).two_norm());
    return maxDistance;
}

/*!
 * \ingroup EmbeddedCoupling
 * \brief Simple legacy VTK writer for outputting debug data on the coupling interface
 */
class DebugIntersectionVTKOutput
{
public:
    void push(const std::array<Dune::FieldVector<double, 3>, 3>& corners, double data)
    {
        vtkCells_.emplace_back(std::array<std::size_t, 3>{
            { vtkPoints_.size(), vtkPoints_.size()+1, vtkPoints_.size()+2 }
        });
        for (const auto& c : corners)
            vtkPoints_.push_back(c);
        vtkCellData_.push_back(data);
    }

    void write(const std::string& name)
    {
        std::ofstream intersectionFile(name);
        intersectionFile << "# vtk DataFile Version 2.0\n";
        intersectionFile << "Really cool intersection data\n";
        intersectionFile << "ASCII\n";
        intersectionFile << "DATASET UNSTRUCTURED_GRID\n";
        intersectionFile << "POINTS " << vtkPoints_.size() << " float\n";
        for (const auto& p : vtkPoints_)
            intersectionFile << p[0] << " " << p[1] << " " << p[2] << "\n";
        intersectionFile << "\nCELLS " << vtkCells_.size() << " " << vtkCells_.size()*4 << "\n";
        for (const auto& c : vtkCells_)
            intersectionFile << "3 " << c[0] << " " << c[1] << " " << c[2] << "\n";
        intersectionFile << "\nCELL_TYPES " << vtkCells_.size() << "\n";
        for (int i = 0; i < vtkCells_.size(); ++i)
            intersectionFile << "5\n";
        intersectionFile << "\nCELL_DATA " << vtkCells_.size() << "\nSCALARS index float\nLOOKUP_TABLE default\n";
        for (const auto& c : vtkCellData_)
            intersectionFile << c << "\n";
    }

private:
    std::vector<Dune::FieldVector<double, 3>> vtkPoints_;
    std::vector<std::array<std::size_t, 3>> vtkCells_;
    std::vector<double> vtkCellData_;
};

} // end namespace Detail

namespace Embedded1d3dCouplingMode {
struct Projection : public Utility::Tag<Projection> {
    static std::string name() { return "projection"; }
};

inline constexpr Projection projection{};
} // end namespace Embedded1d3dCouplingMode

// forward declaration
template<class MDTraits, class CouplingMode>
class Embedded1d3dCouplingManager;

/*!
 * \ingroup EmbeddedCoupling
 * \brief Manages the coupling between bulk elements and lower dimensional elements
 *
 * The method is described in detail in the article
 * "Projection-based resolved interface mixed-dimension method for embedded tubular network systems"
 * by Timo Koch (2021) available at https://arxiv.org/abs/2106.06358
 *
 * The network is represented as a line segment network. The bulk domain explicitly resolves the wall
 * of the tube network (e.g. blood vessel wall, outer root wall), so this generally requires an unstructured grid.
 * The coupling term coupling the PDEs on the network and in the bulk domain is integrated over the
 * coupling surface and each integration point couples to quantities evaluated at the closest point on the graph.
 * There is a unique mapping from every point on the virtual tube surface to the tube centerline, therefore
 * this coupling is determined by mapping to the closest point on the virtual surface and then evaluating
 * the 1D quantity on the mapped centerline point. (The inverse mapping may be non-unique.)
 *
 * The original paper also includes an algorithm for computing exact intersection for the surface quadrature
 * in case of simple straight vessels. However, as a comparison in the paper shows that the suggested
 * approximate quadrature is exact enough (configured by parameter MixedDimension.Projection.SimplexIntegrationRefine
 * specifying the number of (adaptive local) virtual refinements of a surface facet for integration), only
 * the approximate general purpose algorithm is included in this implementation. Only one quadrature point
 * will be added for each surface element and coupled 1D dof including all contribution evaluated by
 * local refinement. That means the virtual refinement doesn't increase the number of integration points and
 * additionally is also optimized by using local adaptive refinement only in the places where the index mapping
 * to the 1D degree of freedom is still multivalent.
 *
 * Coupling sources are implement in terms of point sources at quadrature points to make it easy to reuse code
 * written for other 1D-3D coupling manager modes.
 *
 * This algorithm can be configured by several parameters
 * (although this is usually not necessary as there are sensible defaults):
 *
 *   - MixedDimension.Projection.SimplexIntegrationRefine
 *         number of virtual refinement steps to determine coupled surface area
 *   - MixedDimension.Projection.EnableIntersectionOutput
 *         set to true to enable debug VTK output for intersections
 *   - MixedDimension.Projection.EstimateNumberOfPointSources
 *         provide an estimate for the expected number of coupling points for memory allocation
 *   - MixedDimension.Projection.CoupledRadiusFactor
 *         threshold distance in which to search for coupled elements (specified as multiple of radius, default 0.1)
 *   - MixedDimension.Projection.CoupledAngleFactor
 *         threshold for the cosine of the angle between the surface normal vector and the vector
 *         from the face center to its closest point on the centerline, a face is only coupled if the cosine
 *         exceeds this value (default 0.3, i.e. angles up to about 72.5 degrees)
 *   - MixedDimension.Projection.ConsiderFacesWithinBoundingBoxCoupled
 *         determines if all 3D boundary facets within the mesh bounding box should be considered as coupling faces
 *   - MixedDimension.Projection.CoupledBoundingBoxShrinkingFactor
 *         if ConsiderFacesWithinBoundingBoxCoupled=true shrink the bounding box in all directions by this factor
 *
 */
template<class MDTraits>
class Embedded1d3dCouplingManager<MDTraits, Embedded1d3dCouplingMode::Projection>
: public EmbeddedCouplingManagerBase<MDTraits, Embedded1d3dCouplingManager<MDTraits, Embedded1d3dCouplingMode::Projection>>
{
    using ThisType = Embedded1d3dCouplingManager<MDTraits, Embedded1d3dCouplingMode::Projection>;
    using ParentType = EmbeddedCouplingManagerBase<MDTraits, ThisType>;

    using Scalar = typename MDTraits::Scalar;
    using SolutionVector = typename MDTraits::SolutionVector;
    using PointSourceData = typename ParentType::PointSourceTraits::PointSourceData;

    static constexpr auto bulkIdx = typename MDTraits::template SubDomain<0>::Index();
    static constexpr auto lowDimIdx = typename MDTraits::template SubDomain<1>::Index();

    // the sub domain type aliases
    template<std::size_t id> using SubDomainTypeTag = typename MDTraits::template SubDomain<id>::TypeTag;
    template<std::size_t id> using Problem = GetPropType<SubDomainTypeTag<id>, Properties::Problem>;
    template<std::size_t id> using GridGeometry = GetPropType<SubDomainTypeTag<id>, Properties::GridGeometry>;
    template<std::size_t id> using GridView = typename GridGeometry<id>::GridView;
    template<std::size_t id> using Element = typename GridView<id>::template Codim<0>::Entity;
    template<std::size_t id> using GridIndex = typename IndexTraits<GridView<id>>::GridIndex;

    static constexpr int lowDimDim = GridView<lowDimIdx>::dimension;
    static constexpr int bulkDim = GridView<bulkIdx>::dimension;
    static constexpr int dimWorld = GridView<bulkIdx>::dimensionworld;

    template<std::size_t id>
    static constexpr bool isCVFE()
    { return DiscretizationMethods::isCVFE<typename GridGeometry<id>::DiscretizationMethod>; }

    using GlobalPosition = typename Element<bulkIdx>::Geometry::GlobalCoordinate;

    using SegmentNetwork = Detail::SegmentNetwork<GridGeometry<lowDimIdx>>;

public:
    static constexpr Embedded1d3dCouplingMode::Projection couplingMode{};

    using ParentType::ParentType;

    void init(std::shared_ptr<Problem<bulkIdx>> bulkProblem,
              std::shared_ptr<Problem<lowDimIdx>> lowDimProblem,
              const SolutionVector& curSol)
    {
        ParentType::init(bulkProblem, lowDimProblem, curSol);
    }

    /* \brief Compute integration point point sources and associated data
     * This method computes a source term on the explicitly given surface of the network grid
     * by integrating over the elements of the surface, projecting the integration points onto
     * the centerlines to evaluate 1d quantities and set point sources (minimum distance projection)
     *
     * \param order Unused in this algorithm! (passed from default 1D3D coupling manager)
     */
    void computePointSourceData(std::size_t order = 1, bool /*verbose*/ = false)
    {
        if (this->gridView(bulkIdx).comm().size() > 1)
            DUNE_THROW(Dune::NotImplemented, "The projection coupling on distributed grids");

        // initialize the maps
        // do some logging and profiling
        Dune::Timer watch;
        std::cout << "Initializing the coupling manager (projection)" << std::endl;

        // clear all internal members like point source vectors and stencils
        this->clear();

        // debug VTK output
        const bool enableIntersectionOutput = getParam<bool>("MixedDimension.Projection.EnableIntersectionOutput", false);
        std::unique_ptr<Detail::DebugIntersectionVTKOutput> vtkCoupledFaces =
            enableIntersectionOutput ? std::make_unique<Detail::DebugIntersectionVTKOutput>() : nullptr;

        // temporarily represent the network domain as a list of segments
        // a more efficient data structure
        const auto& lowDimProblem = this->problem(lowDimIdx);
        const auto& lowDimFvGridGeometry = gridDiscretization(lowDimProblem);
        const auto& lowDimGridView = lowDimFvGridGeometry.gridView();
        localAreaFactor_.assign(lowDimGridView.size(0), 0.0);
        const auto radiusFunc = [&](const GridIndex<lowDimIdx> eIdx) {
            return lowDimProblem.spatialParams().radius(eIdx);
        };
        const auto network = SegmentNetwork{
            lowDimFvGridGeometry, radiusFunc
        };

        // precompute the dof indices for efficiency (control-volume finite element methods)
        this->precomputeElementDofIndices(bulkIdx);
        this->precomputeElementDofIndices(lowDimIdx);

        const auto& bulkProblem = this->problem(bulkIdx);
        const auto& bulkFvGridGeometry = gridDiscretization(bulkProblem);
        const auto& bulkGridView = bulkFvGridGeometry.gridView();
        bulkElementMarker_.assign(bulkGridView.size(0), 0);
        bulkVertexMarker_.assign(bulkGridView.size(bulkDim), 0);

        // construct a bounding box that may be used to determine coupled faces
        // the shrinking factor shrinks this bounding box
        const auto bBoxShrinking
            = getParam<Scalar>("MixedDimension.Projection.CoupledBoundingBoxShrinkingFactor", 1e-2);
        const GlobalPosition threshold(
            bBoxShrinking*( bulkFvGridGeometry.bBoxMax()-bulkFvGridGeometry.bBoxMin() ).two_norm()
        );
        const auto bBoxMinSmall = bulkFvGridGeometry.bBoxMin() + threshold;
        const auto bBoxMaxSmall = bulkFvGridGeometry.bBoxMax() - threshold;
        auto insideBBox = [bBoxMin=bBoxMinSmall, bBoxMax=bBoxMaxSmall](const GlobalPosition& point) -> bool
        {
            static constexpr Scalar eps_ = 1.0e-7;
            const Scalar eps0 = eps_*(bBoxMax[0] - bBoxMin[0]);
            const Scalar eps1 = eps_*(bBoxMax[1] - bBoxMin[1]);
            const Scalar eps2 = eps_*(bBoxMax[2] - bBoxMin[2]);
            return (bBoxMin[0] - eps0 <= point[0] && point[0] <= bBoxMax[0] + eps0 &&
                    bBoxMin[1] - eps1 <= point[1] && point[1] <= bBoxMax[1] + eps1 &&
                    bBoxMin[2] - eps2 <= point[2] && point[2] <= bBoxMax[2] + eps2);
        };

        // setup simplex "quadrature" rule
        const auto quadSimplexRefineMaxLevel
            = getParam<std::size_t>("MixedDimension.Projection.SimplexIntegrationRefine", 4);

        // estimate number of point sources for memory allocation
        const auto estimateNumberOfPointSources
            = getParam<std::size_t>("MixedDimension.Projection.EstimateNumberOfPointSources", bulkFvGridGeometry.gridView().size(0));

        this->pointSources(bulkIdx).reserve(estimateNumberOfPointSources);
        this->pointSources(lowDimIdx).reserve(estimateNumberOfPointSources);
        this->pointSourceData().reserve(estimateNumberOfPointSources);

        // we go over the bulk surface intersection, check if we are coupled,
        // and if this is the case, we add integration point sources using
        // the local simplex refinement quadrature procedure (see paper Koch 2021)
        for (const auto& element : elements(bulkGridView))
        {
            const auto bulkElementIdx = bulkFvGridGeometry.elementMapper().index(element);
            const auto bulkGeometry = element.geometry();
            const auto refElement = Dune::referenceElement(element);
            for (const auto& intersection : intersections(bulkGridView, element))
            {
                // skip inner intersection
                if (!intersection.boundary())
                    continue;

                // only segments that can be closest for some point on the facet are searched
                const auto interfaceGeometry = intersection.geometry();
                const auto center = interfaceGeometry.center();
                const auto candidates = network.closestSegmentSurfaceCandidates(
                    center, Detail::maxCornerDistance(interfaceGeometry, center)
                );

                // check if we are close enough to the coupling interface
                const auto [isCoupled, closestSegmentIdx, minDist] = isCoupled_(network, candidates, intersection, insideBBox);
                if (!isCoupled)
                    continue;

                // we found an intersection on the coupling interface: mark as coupled element
                bulkElementMarker_[bulkElementIdx] = 1;

                for (int i = 0; i < interfaceGeometry.corners(); ++i)
                {
                    const auto localVIdx = refElement.subEntity(intersection.indexInInside(), 1, i, bulkDim);
                    const auto vIdx = bulkFvGridGeometry.vertexMapper().subIndex(element, localVIdx, bulkDim);
                    bulkVertexMarker_[vIdx] = 1;
                }

                // debug graphical intersection output
                if (enableIntersectionOutput)
                    vtkCoupledFaces->push({
                        { interfaceGeometry.corner(0), interfaceGeometry.corner(1), interfaceGeometry.corner(2)}
                    }, closestSegmentIdx);

                const auto indicator = Detail::NetworkIndicatorFunction { network, candidates };
                const auto quadSimplex = LocalRefinementSimplexQuadrature{ interfaceGeometry, quadSimplexRefineMaxLevel, indicator };

                // iterate over all quadrature points (children simplices)
                for (const auto& qp : quadSimplex)
                {
                    const auto surfacePoint = interfaceGeometry.global(qp.position());
                    const auto closestSegmentIdx = qp.indicator();
                    const auto projPoint = network.projectionPointOnSegment(surfacePoint, closestSegmentIdx);

                    addPointSourceAtIP_(bulkGeometry, interfaceGeometry, qp, bulkElementIdx, closestSegmentIdx, surfacePoint, projPoint);
                    addStencilEntries_(bulkElementIdx, closestSegmentIdx);
                    this->averageDistanceToBulkCell().push_back(averageDistancePointGeometry(surfacePoint, bulkGeometry));
                }
            }
        }

        // make stencils unique
        using namespace Dune::Hybrid;
        forEach(integralRange(Dune::index_constant<2>{}), [&](const auto domainIdx)
        {
            for (auto&& stencil : this->couplingStencils(domainIdx))
            {
                std::sort(stencil.second.begin(), stencil.second.end());
                stencil.second.erase(
                    std::unique(stencil.second.begin(), stencil.second.end()),
                    stencil.second.end()
                );
            }
        });

        // debug VTK output
        if (enableIntersectionOutput)
            vtkCoupledFaces->write("coupledfaces.vtk");

        // compare theoretical cylinder surface to actual discrete surface
        // of the coupled bulk interface facets
        for (const auto& element : elements(lowDimGridView))
        {
            const auto length = element.geometry().volume();
            const auto eIdx = lowDimFvGridGeometry.elementMapper().index(element);
            const auto radius = this->problem(lowDimIdx).spatialParams().radius(eIdx);
            const auto cylinderSurface = 2.0*M_PI*radius*length;
            localAreaFactor_[eIdx] /= cylinderSurface;
            localAreaFactor_[eIdx] -= 1.0;
        }

        this->pointSources(bulkIdx).shrink_to_fit();
        this->pointSources(lowDimIdx).shrink_to_fit();
        this->pointSourceData().shrink_to_fit();

        std::cout << "-- Coupling at " << this->pointSourceData().size() << " integration points." << std::endl;
        std::cout << "-- Finished initialization after " << watch.elapsed() << " seconds." << std::endl;
    }

    /*!
     * \brief Methods to be accessed by the subproblems
     */
    // \{

    //! Returns the radius of the lower dimensional element
    Scalar radius(std::size_t id) const
    {
        const auto& data = this->pointSourceData()[id];
        return this->problem(lowDimIdx).spatialParams().radius(data.lowDimElementIdx());
    }

    // \}

    //! Marker that is non-zero for bulk elements that are coupled
    const std::vector<int>& bulkElementMarker() const
    { return bulkElementMarker_; }

    //! Marker that is non-zero for bulk vertices that belong to coupled surface facets
    const std::vector<int>& bulkVertexMarker() const
    { return bulkVertexMarker_; }

    //! For each 1D element, the ratio of associated discrete bulk surface area
    //! to the theoretical bulk surface area resulting from assuming a cylinder-shaped segment
    //! (excluding cylinder caps).
    const std::vector<Scalar>& localAreaFactor() const
    { return localAreaFactor_; }

private:
    std::vector<int> bulkElementMarker_, bulkVertexMarker_;
    std::vector<double> localAreaFactor_;

    // add a point source containing all data needed to evaluate
    // the coupling source term from both domains
    template<class BulkGeometry, class InterfaceGeometry, class QP>
    void addPointSourceAtIP_(const BulkGeometry& bulkGeometry,
                             const InterfaceGeometry& ifGeometry,
                             const QP& qp,
                             const GridIndex<bulkIdx> bulkElementIdx,
                             const GridIndex<lowDimIdx> lowDimElementIdx,
                             const GlobalPosition& surfacePoint,
                             const GlobalPosition& projPoint)
    {
        // add a new point source id (for the associated data)
        const auto id = this->idCounter_++;

        // add a bulk point source
        const auto qpweight = qp.weight();
        const auto integrationElement = ifGeometry.integrationElement(qp.position());
        this->pointSources(bulkIdx).emplace_back(surfacePoint, id, qpweight, integrationElement, bulkElementIdx);

        // find the corresponding projection point and 1d element
        this->pointSources(lowDimIdx).emplace_back(projPoint, id, qpweight, integrationElement, lowDimElementIdx);
        localAreaFactor_[lowDimElementIdx] += qpweight*integrationElement;

        const auto& bulkFvGridGeometry = gridDiscretization(this->problem(bulkIdx));
        const auto& lowDimFvGridGeometry = gridDiscretization(this->problem(lowDimIdx));

        // pre compute additional data used for the evaluation of
        // the actual solution dependent source term
        PointSourceData psData;

        if constexpr (isCVFE<lowDimIdx>())
        {
            std::vector<Dune::FieldVector<Scalar, 1> > shapeValues;
            const auto lowDimGeometry = gridDiscretization(this->problem(lowDimIdx)).element(lowDimElementIdx).geometry();
            this->getShapeValues(lowDimIdx, lowDimFvGridGeometry, lowDimGeometry, projPoint, shapeValues);
            psData.addLowDimInterpolation(shapeValues, this->elementDofIndices(lowDimIdx, lowDimElementIdx), lowDimElementIdx);
        }
        else
        {
            psData.addLowDimInterpolation(lowDimElementIdx);
        }

        if constexpr (isCVFE<bulkIdx>())
        {
            std::vector<Dune::FieldVector<Scalar, 1> > shapeValues;
            this->getShapeValues(bulkIdx, bulkFvGridGeometry, bulkGeometry, surfacePoint, shapeValues);
            psData.addBulkInterpolation(shapeValues, this->elementDofIndices(bulkIdx, bulkElementIdx), bulkElementIdx);
        }
        else
        {
            psData.addBulkInterpolation(bulkElementIdx);
        }

        // publish point source data in the global vector
        this->pointSourceData().emplace_back(std::move(psData));
    }

    void addStencilEntries_(const GridIndex<bulkIdx> bulkElementIdx, const GridIndex<lowDimIdx> lowDimElementIdx)
    {
        // export the bulk coupling stencil
        if constexpr (isCVFE<lowDimIdx>())
        {
            const auto& dofs = this->elementDofIndices(lowDimIdx, lowDimElementIdx);
            this->couplingStencils(bulkIdx)[bulkElementIdx].insert(
                this->couplingStencils(bulkIdx)[bulkElementIdx].end(),
                dofs.begin(), dofs.end()
            );

        }
        else
        {
            this->couplingStencils(bulkIdx)[bulkElementIdx].push_back(lowDimElementIdx);
        }

        // export the lowdim coupling stencil
        if constexpr (isCVFE<bulkIdx>())
        {
            const auto& dofs = this->elementDofIndices(bulkIdx, bulkElementIdx);
            this->couplingStencils(lowDimIdx)[lowDimElementIdx].insert(
                this->couplingStencils(lowDimIdx)[lowDimElementIdx].end(),
                dofs.begin(), dofs.end()
            );
        }
        else
        {
            this->couplingStencils(lowDimIdx)[lowDimElementIdx].push_back(bulkElementIdx);
        }
    }

    template<class Candidates, class BoundaryIntersection, class BBoxCheck>
    auto isCoupled_(const SegmentNetwork& network, const Candidates& candidates,
                    const BoundaryIntersection& intersection, const BBoxCheck& insideBBox) const
    {
        const auto center = intersection.geometry().center();
        const auto [eIdx, minDist] = network.findClosestSegmentSurface(center, candidates);

        // all elements in the bounding box are coupled if this is enabled
        // this is a simplified check that can be enabled for box-shaped domains
        static const bool considerBBox = getParam<bool>("MixedDimension.Projection.ConsiderFacesWithinBoundingBoxCoupled", false);
        if (considerBBox && insideBBox(center))
            return std::make_tuple(true, eIdx, minDist);

        // general coupling face detection algorithm:
        // if the distance is under a certain threshold we are coupled (and the angle is not very big)
        static const auto rFactor = getParam<Scalar>("MixedDimension.Projection.CoupledRadiusFactor", 0.1);
        static const auto spFactor = getParam<Scalar>("MixedDimension.Projection.CoupledAngleFactor", 0.3);
        const auto projPoint = network.projectionPointOnSegment(center, eIdx);
        const auto distVec = projPoint - center;
        const bool isCoupled = minDist < rFactor * network.radius(eIdx) && distVec * intersection.centerUnitOuterNormal() > distVec.two_norm() * spFactor;
        return std::make_tuple(isCoupled, eIdx, minDist);
    }
};

//! we support multithreaded assembly
template<class MDTraits>
struct CouplingManagerSupportsMultithreadedAssembly<Embedded1d3dCouplingManager<MDTraits, Embedded1d3dCouplingMode::Projection>>
: public std::true_type {};

} // end namespace Dumux

#endif
