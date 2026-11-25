// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup EmbeddedCoupling
 * \brief Coupling manager for low-dimensional domains embedded in the bulk domain.
 */

#ifndef DUMUX_MULTIDOMAIN_EMBEDDED_COUPLINGMANAGER_1D3D_KERNEL_HH
#define DUMUX_MULTIDOMAIN_EMBEDDED_COUPLINGMANAGER_1D3D_KERNEL_HH

#include <array>
#include <type_traits>
#include <vector>

#include <dune/common/timer.hh>
#include <dune/geometry/quadraturerules.hh>
#include <dune/grid/yaspgrid.hh>

#include <dumux/common/tag.hh>
#include <dumux/common/properties.hh>
#include <dumux/common/typetraits/griddiscretization.hh>
#include <dumux/common/indextraits.hh>
#include <dumux/common/typetraits/localdofs_.hh>
#include <dumux/geometry/distance.hh>

#include <dumux/multidomain/embedded/couplingmanager1d3dbase.hh>
#include <dumux/multidomain/embedded/cylinderintegration.hh>
#include <dumux/multidomain/embedded/extendedsourcestencil.hh>

namespace Dumux {

namespace Detail {

//! YaspGrids with equidistant cells admit a direct computation of the cell containing a point
template<class Grid>
inline constexpr bool isEquidistantYaspGrid = false;

template<int dim, class ctype>
inline constexpr bool isEquidistantYaspGrid<Dune::YaspGrid<dim, Dune::EquidistantCoordinates<ctype, dim>>> = true;

template<int dim, class ctype>
inline constexpr bool isEquidistantYaspGrid<Dune::YaspGrid<dim, Dune::EquidistantOffsetCoordinates<ctype, dim>>> = true;

} // end namespace Detail

namespace Embedded1d3dCouplingMode {
struct Kernel : public Utility::Tag<Kernel> {
    static std::string name() { return "kernel"; }
};

inline constexpr Kernel kernel{};
} // end namespace Embedded1d3dCouplingMode

// forward declaration
template<class MDTraits, class CouplingMode>
class Embedded1d3dCouplingManager;

/*!
 * \ingroup EmbeddedCoupling
 * \brief Manages the coupling between bulk elements and lower dimensional elements
 *        Point sources on each integration point are computed by an AABB tree.
 * \note Specialization for coupling method using a distributed kernel source with 3d quantities evaluated on the line
 */
template<class MDTraits>
class Embedded1d3dCouplingManager<MDTraits, Embedded1d3dCouplingMode::Kernel>
: public Embedded1d3dCouplingManagerBase<MDTraits, Embedded1d3dCouplingManager<MDTraits, Embedded1d3dCouplingMode::Kernel>>
{
    using ThisType = Embedded1d3dCouplingManager<MDTraits, Embedded1d3dCouplingMode::Kernel>;
    using ParentType = Embedded1d3dCouplingManagerBase<MDTraits, ThisType>;
    using BaseType = EmbeddedCouplingManagerBase<MDTraits, ParentType>;
    using Scalar = typename MDTraits::Scalar;
    using SolutionVector = typename MDTraits::SolutionVector;
    using PointSourceData = typename BaseType::PointSourceTraits::PointSourceData;

    static constexpr auto bulkIdx = typename MDTraits::template SubDomain<0>::Index();
    static constexpr auto lowDimIdx = typename MDTraits::template SubDomain<1>::Index();

    // the sub domain type aliases
    template<std::size_t id> using SubDomainTypeTag = typename MDTraits::template SubDomain<id>::TypeTag;
    template<std::size_t id> using Problem = GetPropType<SubDomainTypeTag<id>, Properties::Problem>;
    template<std::size_t id> using GridGeometry = GetPropType<SubDomainTypeTag<id>, Properties::GridGeometry>;
    template<std::size_t id> using GridView = typename GridGeometry<id>::GridView;
    template<std::size_t id> using Element = typename GridView<id>::template Codim<0>::Entity;
    template<std::size_t id> using GridIndex = typename IndexTraits<GridView<id>>::GridIndex;

    using GlobalPosition = typename Element<bulkIdx>::Geometry::GlobalCoordinate;

    //! a straight piece of the centerline
    struct Segment
    {
        GlobalPosition a, b;
        const GlobalPosition& corner(int i) const { return i == 0 ? a : b; }
        Scalar volume() const { return (b - a).two_norm(); }
    };

    template<std::size_t id>
    static constexpr bool isCVFE()
    { return DiscretizationMethods::isCVFE<typename GridGeometry<id>::DiscretizationMethod>; }

    //! the maximum number of local dofs of a bulk element
    static constexpr std::size_t maxNumBulkLocalDofs = isCVFE<bulkIdx>()
        ? Detail::LocalDofs::maxNumLocalDofs<typename GridGeometry<bulkIdx>::LocalView>() : 1;

    //! Equidistant YaspGrids allow an O(1) point lookup; otherwise the bounding box tree is used
    static constexpr bool bulkHasCartesianLookup
        = Detail::isEquidistantYaspGrid<typename GridView<bulkIdx>::Grid>
          && std::is_same_v<GridView<bulkIdx>, typename GridView<bulkIdx>::Grid::LeafGridView>;

    enum {
        bulkDim = GridView<bulkIdx>::dimension,
        lowDimDim = GridView<lowDimIdx>::dimension,
        dimWorld = GridView<bulkIdx>::dimensionworld
    };

    // detect if a class (the spatial params) has a kernelWidthFactor() function
    template <typename T, typename ...Ts>
    using VariableKernelWidthDetector = decltype(std::declval<T>().kernelWidthFactor(std::declval<Ts>()...));

    template<class T, typename ...Args>
    static constexpr bool hasKernelWidthFactor()
    { return Dune::Std::is_detected<VariableKernelWidthDetector, T, Args...>::value; }

public:
    static constexpr Embedded1d3dCouplingMode::Kernel couplingMode{};

    using ParentType::ParentType;

    void init(std::shared_ptr<Problem<bulkIdx>> bulkProblem,
              std::shared_ptr<Problem<lowDimIdx>> lowDimProblem,
              const SolutionVector& curSol)
    {
        ParentType::init(bulkProblem, lowDimProblem, curSol);
    }

    /*!
     * \brief extend the jacobian pattern of the diagonal block of domain i
     *        by those entries that are not already in the uncoupled pattern
     * \note Such additional dependencies can arise from the coupling, e.g. if a coupling source
     *       term depends on a non-local average of a quantity of the same domain
     */
    template<std::size_t id, class JacobianPattern>
    void extendJacobianPattern(Dune::index_constant<id> domainI, JacobianPattern& pattern) const
    {
        extendedSourceStencil_.extendJacobianPattern(*this, domainI, pattern);
    }

    /*!
     * \brief evaluate additional derivatives of the element residual of a domain with respect
     *        to dofs in the same domain that are not in the regular stencil (per default this is not the case)
     * \note Such additional dependencies can arise from the coupling, e.g. if a coupling source
     *       term depends on a non-local average of a quantity of the same domain
     * \note This is the same for box and cc
     */
    template<std::size_t i, class LocalAssemblerI, class JacobianMatrixDiagBlock, class GridVariables>
    void evalAdditionalDomainDerivatives(Dune::index_constant<i> domainI,
                                         const LocalAssemblerI& localAssemblerI,
                                         const typename LocalAssemblerI::LocalResidual::ElementResidualVector&,
                                         JacobianMatrixDiagBlock& A,
                                         GridVariables& gridVariables)
    {
        extendedSourceStencil_.evalAdditionalDomainDerivatives(*this, domainI, localAssemblerI, A, gridVariables);
    }

    /* \brief Compute integration point point sources and associated data
     *
     * This method uses grid glue to intersect the given grids. Over each intersection
     * we later need to integrate a source term. This method places point sources
     * at each quadrature point and provides the point source with the necessary
     * information to compute integrals (quadrature weight and integration element).
     * The kernel of a quadrature point is distributed over the cross-section through the point,
     * with the part of the intersection length given by the quadrature weight.
     * The order of the quadrature rule is MixedDimension.KernelIntegrationOrder (default 1).
     * \note The parameters have no effect
     */
    void computePointSourceData(std::size_t /*order*/ = 1, bool /*verbose*/ = false)
    {
        // initialize the maps
        // do some logging and profiling
        Dune::Timer watch;
        std::cout << "[coupling] Initializing the integration point source data structures..." << std::endl;

        // prepare the internal data structures
        prepareDataStructures_();
        std::cout << "[coupling] Resized data structures." << std::endl;

        const auto& bulkGridGeometry = gridDiscretization(this->problem(bulkIdx));
        const auto& lowDimGridGeometry = gridDiscretization(this->problem(lowDimIdx));

        // generate a bunch of random vectors and values for
        // Monte-carlo integration on the cylinder defined by line and radius
        const double characteristicRelativeLength = getParam<double>("MixedDimension.KernelIntegrationCRL", 0.1);
        EmbeddedCoupling::CylinderIntegration<Scalar> cylIntegration(characteristicRelativeLength, 1);
        const int integrationOrder = getParam<int>("MixedDimension.KernelIntegrationOrder", 1);

        const bool writeIntegrationPointsToFile = getParam<bool>("MixedDimension.WriteIntegrationPointsToFile", false);
        if (writeIntegrationPointsToFile)
        {
            std::ofstream ipPointFile("kernel_points.log", std::ios::trunc);
            ipPointFile << "x,y,z\n";
            std::cout << "[coupling] Initialized kernel_points.log." << std::endl;
        }

        for (const auto& is : intersections(this->glue()))
        {
            // all inside elements are identical...
            const auto& inside = is.targetEntity(0);
            // get the intersection geometry for integrating over it
            const auto intersectionGeometry = is.geometry();
            const auto lowDimElementIdx = lowDimGridGeometry.elementMapper().index(inside);

            // for each intersection integrate kernel and add:
            //  * 1d: a new point source
            //  * 3d: a new kernel volume source
            const auto radius = this->problem(lowDimIdx).spatialParams().radius(lowDimElementIdx);
            const auto kernelWidthFactor = kernelWidthFactor_(this->problem(lowDimIdx).spatialParams(), lowDimElementIdx);
            const auto kernelWidth = kernelWidthFactor*radius;
            const auto& quad = Dune::QuadratureRules<Scalar, lowDimDim>::rule(intersectionGeometry.type(), integrationOrder);
            for (const auto& qp : quad)
            {
                // the part of the intersection that belongs to the quadrature point, centered at the point
                const auto center = intersectionGeometry.global(qp.position());
                auto a = intersectionGeometry.corner(0);
                auto b = intersectionGeometry.corner(1);
                if (quad.size() > 1)
                {
                    auto halfSegment = b - a;
                    halfSegment *= 0.5*qp.weight();
                    a = center - halfSegment;
                    b = center + halfSegment;
                }
                const Segment segment{a, b};
                cylIntegration.setGeometry(a, b, kernelWidth);

                // we can have multiple 3d elements per 1d intersection
                for (int outsideIdx = 0; outsideIdx < is.numDomainNeighbors(); ++outsideIdx)
                {
                    // compute source id
                    // for each id we have
                    // (1) a point source for the 1d domain
                    // (2) a bulk source weight for the each element in the 3d domain (the fraction of the total source/sink for the 1d element with the point source)
                    //     TODO: i.e. one integration of the kernel should be enough (for each of the weights it's weight/embeddings)
                    // (3) the flux scaling factor for each outside element, i.e. each id
                    const auto id = this->idCounter_++;

                    // compute the weights for the bulk volume sources
                    const auto& outside = is.domainEntity(outsideIdx);
                    const auto bulkElementIdx = bulkGridGeometry.elementMapper().index(outside);
                    const auto surfaceFactor = computeBulkSource(segment, radius, kernelWidth, id, lowDimElementIdx, bulkElementIdx, cylIntegration, is.numDomainNeighbors());

                    // place a point source at the quadrature point
                    this->pointSources(lowDimIdx).emplace_back(center, id, surfaceFactor, segment.volume(), lowDimElementIdx);
                    this->pointSources(lowDimIdx).back().setEmbeddings(is.numDomainNeighbors());

                    // pre compute additional data used for the evaluation of
                    // the actual solution dependent source term
                    PointSourceData psData;

                    // lowdim interpolation (evaluate at the quadrature point)
                    if constexpr (isCVFE<lowDimIdx>())
                    {
                        using ShapeValues = std::vector<Dune::FieldVector<Scalar, 1> >;
                        const auto lowDimGeometry = gridDiscretization(this->problem(lowDimIdx)).element(lowDimElementIdx).geometry();
                        ShapeValues shapeValues;
                        this->getShapeValues(lowDimIdx, gridDiscretization(this->problem(lowDimIdx)), lowDimGeometry, center, shapeValues);
                        psData.addLowDimInterpolation(shapeValues, this->elementDofIndices(lowDimIdx, lowDimElementIdx), lowDimElementIdx);
                    }
                    else
                    {
                        psData.addLowDimInterpolation(lowDimElementIdx);
                    }

                    // bulk interpolation (evaluate at the quadrature point)
                    if constexpr (isCVFE<bulkIdx>())
                    {
                        using ShapeValues = std::vector<Dune::FieldVector<Scalar, 1> >;
                        const auto bulkGeometry = gridDiscretization(this->problem(bulkIdx)).element(bulkElementIdx).geometry();
                        ShapeValues shapeValues;
                        this->getShapeValues(bulkIdx, gridDiscretization(this->problem(bulkIdx)), bulkGeometry, center, shapeValues);
                        psData.addBulkInterpolation(shapeValues, this->elementDofIndices(bulkIdx, bulkElementIdx), bulkElementIdx);
                    }
                    else
                    {
                        psData.addBulkInterpolation(bulkElementIdx);
                    }

                    // publish point source data in the global vector
                    this->pointSourceData().emplace_back(std::move(psData));

                    // the flux scaling factor assumes the radially symmetric solution around a straight
                    // line, so the distance is measured to the line through the segment. The distance to
                    // the segment itself would grow with the parts of the bulk element beyond the segment
                    // ends, e.g. if several lower-dimensional elements are embedded in one bulk element.
                    const auto avgMinDist = averageDistanceLineGeometry(a, b, outside.geometry());
                    this->averageDistanceToBulkCell().push_back(avgMinDist);

                    // a cell-centered bulk pressure represents the element, an interpolated one the value on the centerline
                    const auto scalingDistance = isCVFE<bulkIdx>() ? Scalar(0.0) : avgMinDist;
                    fluxScalingFactor_.push_back(this->problem(bulkIdx).fluxScalingFactor(scalingDistance, radius, kernelWidth));

                    // export the lowdim coupling stencil
                    // we insert all dofs / elements and make it unique later
                    if constexpr (isCVFE<bulkIdx>())
                    {
                        const auto& dofs = this->elementDofIndices(bulkIdx, bulkElementIdx);
                        this->couplingStencils(lowDimIdx)[lowDimElementIdx].insert(this->couplingStencils(lowDimIdx)[lowDimElementIdx].end(),
                                                                                   dofs.begin(), dofs.end());
                    }
                    else
                    {
                        this->couplingStencils(lowDimIdx)[lowDimElementIdx].push_back(bulkElementIdx);
                    }
                }
            }
        }

        // make the stencils unique
        makeUniqueStencil_();

        if (!this->pointSources(bulkIdx).empty())
            DUNE_THROW(Dune::InvalidStateException, "Kernel method shouldn't have point sources in the bulk domain but only volume sources!");

        std::cout << "[coupling] Finished preparing manager in " << watch.elapsed() << " seconds." << std::endl;
    }

    /*!
     * \brief Methods to be accessed by the subproblems
     */
    // \{


    //! the maximum number of local dofs of a bulk element, i.e. the range of the local dof index
    static constexpr std::size_t maxNumBulkElementDofs()
    { return maxNumBulkLocalDofs; }

    //! return all source ids for the local dof with index localDofIdx of a bulk element
    const std::vector<std::size_t>& bulkSourceIds(GridIndex<bulkIdx> eIdx, int localDofIdx = 0) const
    { return bulkSourceIds_[eIdx][localDofIdx]; }

    //! return the source weights for the local dof with index localDofIdx of a bulk element
    const std::vector<Scalar>& bulkSourceWeights(GridIndex<bulkIdx> eIdx, int localDofIdx = 0) const
    { return bulkSourceWeights_[eIdx][localDofIdx]; }

    //! The flux scaling factor for a source with id
    Scalar fluxScalingFactor(std::size_t id) const
    { return fluxScalingFactor_[id]; }

    // \}

    /*!
     * \brief Extended source stencil (for the bulk domain)
     */
    const typename BaseType::template CouplingStencils<bulkIdx>::mapped_type&
    extendedSourceStencil(std::size_t eIdx) const
    {
        const auto& sourceStencils = extendedSourceStencil_.stencil();
        if (auto stencil = sourceStencils.find(eIdx); stencil != sourceStencils.end())
            return stencil->second;

        return this->emptyStencil(bulkIdx);
    }

private:
    //! compute the kernel distributed sources and add stencils
    template<class Line, class CylIntegration>
    Scalar computeBulkSource(const Line& line, const Scalar radius, const Scalar kernelWidth,
                             std::size_t id, GridIndex<lowDimIdx> lowDimElementIdx, GridIndex<bulkIdx> coupledBulkElementIdx,
                             const CylIntegration& cylIntegration, int embeddings)
    {
        // Monte-carlo integration on the cylinder defined by line and radius
        const auto& min = gridDiscretization(this->problem(bulkIdx)).bBoxMin();
        const auto& max = gridDiscretization(this->problem(bulkIdx)).bBoxMax();

        // the direct cell computation needs the global cell counts, so it is only valid on a single process
        std::array<int, bulkDim> cells{};
        bool useCartesianLookup = false;
        if constexpr (bulkHasCartesianLookup)
        {
            const auto& gridView = this->gridView(bulkIdx);
            useCartesianLookup = gridView.comm().size() == 1;
            if (useCartesianLookup)
                for (int d = 0; d < bulkDim; ++d)
                    cells[d] = gridView.grid().levelSize(gridView.grid().maxLevel(), d);
        }
        const auto cylSamples = cylIntegration.size();
        const auto& a = line.corner(0);
        const auto& b = line.corner(1);

        // optionally write debugging / visual output of the integration points
        static const auto writeIntegrationPointsToFile = getParam<bool>("MixedDimension.WriteIntegrationPointsToFile", false);
        if (writeIntegrationPointsToFile)
        {
            std::ofstream ipPointFile("kernel_points.log", std::ios::app);
            for (int i = 0; i < cylSamples; ++i)
            {
                const auto& point = cylIntegration.integrationPoint(i);
                if (const bool hasIntersection = intersectsPointBoundingBox(point, min, max); hasIntersection)
                    ipPointFile << point[0] << "," << point[1] << "," << point[2] << '\n';
            }
        }

        Scalar integral = 0.0;
        for (int i = 0; i < cylSamples; ++i)
        {
            const auto& point = cylIntegration.integrationPoint(i);
            if (const bool hasIntersection = intersectsPointBoundingBox(point, min, max); hasIntersection)
            {
                const auto bulkElementIdx = bulkElementIndex_(point, min, max, cells, useCartesianLookup);
                if (bulkElementIdx == invalidElement_)
                    continue;
                assert(bulkElementIdx < gridDiscretization(this->problem(bulkIdx)).gridView().size(0));
                const auto localWeight = evalConstKernel_(a, b, point, radius, kernelWidth)*cylIntegration.integrationElement(i)/Scalar(embeddings);
                integral += localWeight;

                if constexpr (isCVFE<bulkIdx>())
                {
                    using ShapeValues = std::vector<Dune::FieldVector<Scalar, 1> >;
                    const auto& gg = gridDiscretization(this->problem(bulkIdx));
                    ShapeValues shapeValues;
                    this->getShapeValues(bulkIdx, gg, gg.element(bulkElementIdx).geometry(), point, shapeValues);

                    bool addedStencil = false;
                    for (int v = 0; v < static_cast<int>(shapeValues.size()); ++v)
                    {
                        const auto w = localWeight*shapeValues[v][0];
                        if (!bulkSourceIds_[bulkElementIdx][v].empty() && id == bulkSourceIds_[bulkElementIdx][v].back())
                            bulkSourceWeights_[bulkElementIdx][v].back() += w;
                        else
                        {
                            bulkSourceIds_[bulkElementIdx][v].emplace_back(id);
                            bulkSourceWeights_[bulkElementIdx][v].emplace_back(w);
                            addedStencil = true;
                        }
                    }
                    if (addedStencil)
                        addBulkSourceStencils_(bulkElementIdx, lowDimElementIdx, coupledBulkElementIdx);
                }
                else
                {
                    if (!bulkSourceIds_[bulkElementIdx][0].empty() && id == bulkSourceIds_[bulkElementIdx][0].back())
                        bulkSourceWeights_[bulkElementIdx][0].back() += localWeight;
                    else
                    {
                        bulkSourceIds_[bulkElementIdx][0].emplace_back(id);
                        bulkSourceWeights_[bulkElementIdx][0].emplace_back(localWeight);
                        addBulkSourceStencils_(bulkElementIdx, lowDimElementIdx, coupledBulkElementIdx);
                    }
                }

                // the sample as a point source for problems with point sources of the element discretization interface
                if constexpr (BaseType::template usesElementVariables<bulkIdx>())
                {
                    this->additionalPointSourceContexts(bulkIdx).emplace_back(point, id, localWeight, bulkElementIdx, 1);
                }
            }
        }

        // the surface factor (which fraction of the source is inside the domain and needs to be considered)
        const auto length = (a-b).two_norm()/Scalar(embeddings);
        return integral/length;
    }

    static constexpr GridIndex<bulkIdx> invalidElement_ = std::numeric_limits<GridIndex<bulkIdx>>::max();

    /*!
     * \brief Locate the bulk element containing a point.
     *
     *        An equidistant YaspGrid on a single process allows an O(1) index computation;
     *        any other grid needs the bounding box tree.
     */
    GridIndex<bulkIdx> bulkElementIndex_(const GlobalPosition& point,
                                        const GlobalPosition& min,
                                        const GlobalPosition& max,
                                        const std::array<int, bulkDim>& cells,
                                        bool useCartesianLookup) const
    {
        if constexpr (bulkHasCartesianLookup)
            if (useCartesianLookup)
                return intersectingEntityCartesianGrid(point, min, max, cells);

        const auto candidates = intersectingEntities(
            point, gridDiscretization(this->problem(bulkIdx)).boundingBoxTree());
        return candidates.empty() ? invalidElement_ : smallestIdElement_(candidates);
    }

    /*!
     * \brief Of the bulk elements containing a point, the one with the smallest global id
     *
     * A point on a face, an edge or a vertex lies in several elements. Choosing by the global id gives the
     * same element on every process that holds the candidates, independent of the order of a tree query.
     */
    GridIndex<bulkIdx> smallestIdElement_(const std::vector<std::size_t>& candidates) const
    {
        const auto& gridGeometry = gridDiscretization(this->problem(bulkIdx));
        const auto& idSet = gridGeometry.gridView().grid().globalIdSet();
        auto best = candidates[0];
        if (candidates.size() > 1)
        {
            auto bestId = idSet.id(gridGeometry.element(best));
            for (std::size_t i = 1; i < candidates.size(); ++i)
            {
                const auto id = idSet.id(gridGeometry.element(candidates[i]));
                if (id < bestId)
                {
                    bestId = id;
                    best = candidates[i];
                }
            }
        }
        return best;
    }

    void prepareDataStructures_()
    {
        // clear all internal members like pointsource vectors and stencils
        // initializes the point source id counter
        this->clear();
        bulkSourceIds_.clear();
        bulkSourceWeights_.clear();
        extendedSourceStencil_.stencil().clear();

        // precompute the dof indices for efficiency for the control-volume finite element methods
        this->precomputeElementDofIndices(bulkIdx);
        this->precomputeElementDofIndices(lowDimIdx);

        bulkSourceIds_.resize(this->gridView(bulkIdx).size(0));
        bulkSourceWeights_.resize(this->gridView(bulkIdx).size(0));

        // intersect the bounding box trees
        this->glueGrids();

        // reserve memory for data
        this->pointSourceData().reserve(this->glue().size());
        this->averageDistanceToBulkCell().reserve(this->glue().size());
        fluxScalingFactor_.reserve(this->glue().size());

        // reserve memory for stencils
        const auto numBulkElements = this->gridView(bulkIdx).size(0);
        for (GridIndex<bulkIdx> bulkElementIdx = 0; bulkElementIdx < numBulkElements; ++bulkElementIdx)
        {
            this->couplingStencils(bulkIdx)[bulkElementIdx].reserve(10);
            extendedSourceStencil_.stencil()[bulkElementIdx].reserve(10);
            for (auto& ids : bulkSourceIds_[bulkElementIdx]) ids.reserve(10);
            for (auto& w : bulkSourceWeights_[bulkElementIdx]) w.reserve(10);
        }
    }

    //! Make the stencils unique
    void makeUniqueStencil_()
    {
        // make extra stencils unique
        for (auto&& stencil : extendedSourceStencil_.stencil())
        {
            std::sort(stencil.second.begin(), stencil.second.end());
            stencil.second.erase(std::unique(stencil.second.begin(), stencil.second.end()), stencil.second.end());

            // remove the dofs of the element (control-volume finite element methods)
            if constexpr (isCVFE<bulkIdx>())
            {
                const auto& indices = this->elementDofIndices(bulkIdx, stencil.first);
                stencil.second.erase(std::remove_if(stencil.second.begin(), stencil.second.end(),
                                                   [&](auto i){ return std::find(indices.begin(), indices.end(), i) != indices.end(); }),
                                     stencil.second.end());
            }
            // remove the own element (cell-centered)
            else
            {
                stencil.second.erase(std::remove_if(stencil.second.begin(), stencil.second.end(),
                                                   [&](auto i){ return i == stencil.first; }),
                                     stencil.second.end());
            }
        }

        // make stencils unique
        using namespace Dune::Hybrid;
        forEach(integralRange(Dune::index_constant<2>{}), [&](const auto domainIdx)
        {
            for (auto&& stencil : this->couplingStencils(domainIdx))
            {
                std::sort(stencil.second.begin(), stencil.second.end());
                stencil.second.erase(std::unique(stencil.second.begin(), stencil.second.end()), stencil.second.end());
            }
        });
    }

    //! add additional stencil entries for the bulk element
    void addBulkSourceStencils_(GridIndex<bulkIdx> bulkElementIdx, GridIndex<lowDimIdx> coupledLowDimElementIdx, GridIndex<bulkIdx> coupledBulkElementIdx)
    {
        // add the lowdim element to the coupling stencil of this bulk element
        if constexpr (isCVFE<lowDimIdx>())
        {
            const auto& dofs = this->elementDofIndices(lowDimIdx, coupledLowDimElementIdx);
            this->couplingStencils(bulkIdx)[bulkElementIdx].insert(this->couplingStencils(bulkIdx)[bulkElementIdx].end(),
                                                                   dofs.begin(), dofs.end());

        }
        else
        {
            auto& s = this->couplingStencils(bulkIdx)[bulkElementIdx];
            s.push_back(coupledLowDimElementIdx);
        }

        // the extended source stencil, every 3d element with a source is coupled to
        // the element/dofs where the 3d quantities are measured
        if constexpr (isCVFE<bulkIdx>())
        {
            const auto& dofs = this->elementDofIndices(bulkIdx, coupledBulkElementIdx);
            extendedSourceStencil_.stencil()[bulkElementIdx].insert(extendedSourceStencil_.stencil()[bulkElementIdx].end(),
                                                                    dofs.begin(), dofs.end());
        }
        else
        {
            auto& s = extendedSourceStencil_.stencil()[bulkElementIdx];
            s.push_back(coupledBulkElementIdx);
        }
    }

    //! a cylindrical kernel around the segment a->b
    Scalar evalConstKernel_(const GlobalPosition& a,
                            const GlobalPosition& b,
                            const GlobalPosition& point,
                            const Scalar R,
                            const Scalar rho) const noexcept
    {
        // projection of point onto line a + t*(b-a)
        const auto ab = b - a;
        const auto t = (point - a)*ab/ab.two_norm2();

        // return 0 if we are outside cylinder
        if (t < 0.0 || t > 1.0)
            return 0.0;

        // compute distance
        auto proj = a; proj.axpy(t, ab);
        const auto radiusSquared = (proj - point).two_norm2();

        if (radiusSquared > rho*rho)
            return 0.0;

        /* The default is the uniform disk. Setting MixedDimension.KernelSecondMoment selects
         * (alpha + beta r^2)(1 - r^2/rho^2)^2 instead, with alpha and beta fixed by unit mass
         * and by a prescribed second moment M_2. Choosing M_2 = R^2 makes the exterior of a
         * transverse source agree with the true inclusion's, which the uniform disk does not:
         * it leaves a spurious doublet of strength (M_2 - R^2)/2 that survives outside the
         * support. The condition forces the profile to change sign, so weights are not positive.
         */
        static const bool useMoment = hasParam("MixedDimension.KernelSecondMoment");
        if (!useMoment)
            return 1.0/(M_PI*rho*rho);

        static const Scalar m2 = getParam<Scalar>("MixedDimension.KernelSecondMoment");
        const auto rho2 = rho*rho;
        const auto alpha = 8.0/(M_PI*rho2) - 20.0*m2/(M_PI*rho2*rho2);
        const auto beta = -20.0/(M_PI*rho2*rho2) + 80.0*m2/(M_PI*rho2*rho2*rho2);
        const auto s = 1.0 - radiusSquared/rho2;
        return (alpha + beta*radiusSquared)*s*s;
    }

    /*!
     * \brief Get the kernel width factor from the spatial params (if possible)
     */
    template<class SpatialParams>
    auto kernelWidthFactor_(const SpatialParams& spatialParams, unsigned int eIdx)
    -> std::enable_if_t<hasKernelWidthFactor<SpatialParams, unsigned int>(), Scalar>
    { return spatialParams.kernelWidthFactor(eIdx); }

    /*!
     * \brief Get the kernel width factor (constant) from the input file (if possible)
     */
    template<class SpatialParams>
    auto kernelWidthFactor_(const SpatialParams& spatialParams, unsigned int eIdx)
    -> std::enable_if_t<!hasKernelWidthFactor<SpatialParams, unsigned int>(), Scalar>
    {
        static const Scalar kernelWidthFactor = getParam<Scalar>("MixedDimension.KernelWidthFactor");
        return kernelWidthFactor;
    }

    //! the extended source stencil object for kernel coupling
    EmbeddedCoupling::ExtendedSourceStencil<ThisType> extendedSourceStencil_;

    //! kernel sources to integrate for each bulk element
    std::vector<std::array<std::vector<std::size_t>, maxNumBulkLocalDofs>> bulkSourceIds_;
    //! the integral of the kernel for each point source / integration point, i.e. weight for the source
    std::vector<std::array<std::vector<Scalar>, maxNumBulkLocalDofs>> bulkSourceWeights_;

    //! the flux scaling factor for the respective source id
    std::vector<Scalar> fluxScalingFactor_;
};

//! we support multithreaded assembly
template<class MDTraits>
struct CouplingManagerSupportsMultithreadedAssembly<Embedded1d3dCouplingManager<MDTraits, Embedded1d3dCouplingMode::Kernel>>
: public std::true_type {};

} // end namespace Dumux

#endif
