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
 * \ingroup EmbeddedCoupling
 * \brief Coupling manager for low-dimensional domains embedded in the bulk domain.
 */

#ifndef DUMUX_MULTIDOMAIN_EMBEDDED_COUPLINGMANAGER_1D3D_KERNEL_HH
#define DUMUX_MULTIDOMAIN_EMBEDDED_COUPLINGMANAGER_1D3D_KERNEL_HH

#include <vector>

#include <dune/common/timer.hh>
#include <dune/geometry/quadraturerules.hh>

#include <dumux/common/tag.hh>
#include <dumux/common/properties.hh>
#include <dumux/common/indextraits.hh>
#include <dumux/geometry/distance.hh>

#include <dumux/multidomain/embedded/couplingmanagerbase.hh>
#include <dumux/multidomain/embedded/cylinderintegration.hh>
#include <dumux/multidomain/embedded/extendedsourcestencil.hh>

namespace Dumux {

namespace Embedded1d3dCouplingMode {
struct Kernel : public Utility::Tag<Kernel> {
    static std::string name() { return "kernel"; }

    [[deprecated("Comparison with enum is deprecated. Removed after 3.3. Use tags.")]]
    friend constexpr bool operator==(Kernel, EmbeddedCouplingMode m) { return m == EmbeddedCouplingMode::kernel; }
    [[deprecated("Comparison with enum is deprecated. Removed after 3.3. Use tags.")]]
    friend constexpr bool operator==(EmbeddedCouplingMode m, Kernel) { return m == EmbeddedCouplingMode::kernel; }
    [[deprecated("Comparison with enum is deprecated. Removed after 3.3. Use tags.")]]
    friend constexpr bool operator!=(Kernel, EmbeddedCouplingMode m) { return m != EmbeddedCouplingMode::kernel; }
    [[deprecated("Comparison with enum is deprecated. Removed after 3.3. Use tags.")]]
    friend constexpr bool operator!=(EmbeddedCouplingMode m, Kernel) { return m != EmbeddedCouplingMode::kernel; }
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
: public EmbeddedCouplingManagerBase<MDTraits, Embedded1d3dCouplingManager<MDTraits, Embedded1d3dCouplingMode::Kernel>>
{
    using ThisType = Embedded1d3dCouplingManager<MDTraits, Embedded1d3dCouplingMode::Kernel>;
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

    using GlobalPosition = typename Element<bulkIdx>::Geometry::GlobalCoordinate;

    template<std::size_t id>
    static constexpr bool isBox()
    { return GridGeometry<id>::discMethod == DiscretizationMethod::box; }

    static_assert(!isBox<bulkIdx>() && !isBox<lowDimIdx>(), "The kernel coupling method is only implemented for the tpfa method");
    static_assert(Dune::Capabilities::isCartesian<typename GridView<bulkIdx>::Grid>::v, "The kernel coupling method is only implemented for structured grids");

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
        computeLowDimVolumeFractions();

        const auto refinement = getParamFromGroup<int>(bulkProblem->paramGroup(), "Grid.Refinement", 0);
        if (refinement > 0)
            DUNE_THROW(Dune::NotImplemented, "The current intersection detection may likely fail for refined grids.");
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
        extendedSourceStencil_.evalAdditionalDomainDerivatives(*this, domainI, localAssemblerI, this->curSol(), A, gridVariables);
    }

    /* \brief Compute integration point point sources and associated data
     *
     * This method uses grid glue to intersect the given grids. Over each intersection
     * we later need to integrate a source term. This method places point sources
     * at each quadrature point and provides the point source with the necessary
     * information to compute integrals (quadrature weight and integration element)
     * \param order The order of the quadrature rule for integration of sources over an intersection
     * \param verbose If the point source computation is verbose
     */
    void computePointSourceData(std::size_t order = 1, bool verbose = false)
    {
        // initialize the maps
        // do some logging and profiling
        Dune::Timer watch;
        std::cout << "[coupling] Initializing the integration point source data structures..." << std::endl;

        // prepare the internal data structures
        prepareDataStructures_();
        std::cout << "[coupling] Resized data structures." << std::endl;

        const auto& bulkGridGeometry = this->problem(bulkIdx).gridGeometry();
        const auto& lowDimGridGeometry = this->problem(lowDimIdx).gridGeometry();

        // generate a bunch of random vectors and values for
        // Monte-carlo integration on the cylinder defined by line and radius
        static const double characteristicRelativeLength = getParam<double>("MixedDimension.KernelIntegrationCRL", 0.1);
        EmbeddedCoupling::CylinderIntegration<Scalar> cylIntegration(characteristicRelativeLength, 1);

        static const bool writeIntegrationPointsToFile = getParam<bool>("MixedDimension.WriteIntegrationPointsToFile", false);
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
            const auto a = intersectionGeometry.corner(0);
            const auto b = intersectionGeometry.corner(1);
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
                const auto surfaceFactor = computeBulkSource(intersectionGeometry, radius, kernelWidth, id, lowDimElementIdx, bulkElementIdx, cylIntegration, is.numDomainNeighbors());

                // place a point source at the intersection center
                const auto center = intersectionGeometry.center();
                this->pointSources(lowDimIdx).emplace_back(center, id, surfaceFactor, intersectionGeometry.volume(), lowDimElementIdx);
                this->pointSources(lowDimIdx).back().setEmbeddings(is.numDomainNeighbors());

                // pre compute additional data used for the evaluation of
                // the actual solution dependent source term
                PointSourceData psData;

                // lowdim interpolation (evaluate at center)
                if constexpr (isBox<lowDimIdx>())
                {
                    using ShapeValues = std::vector<Dune::FieldVector<Scalar, 1> >;
                    const auto lowDimGeometry = this->problem(lowDimIdx).gridGeometry().element(lowDimElementIdx).geometry();
                    ShapeValues shapeValues;
                    this->getShapeValues(lowDimIdx, this->problem(lowDimIdx).gridGeometry(), lowDimGeometry, center, shapeValues);
                    psData.addLowDimInterpolation(shapeValues, this->vertexIndices(lowDimIdx, lowDimElementIdx), lowDimElementIdx);
                }
                else
                {
                    psData.addLowDimInterpolation(lowDimElementIdx);
                }

                // bulk interpolation (evaluate at center)
                if constexpr (isBox<bulkIdx>())
                {
                    using ShapeValues = std::vector<Dune::FieldVector<Scalar, 1> >;
                    const auto bulkGeometry = this->problem(bulkIdx).gridGeometry().element(bulkElementIdx).geometry();
                    ShapeValues shapeValues;
                    this->getShapeValues(bulkIdx, this->problem(bulkIdx).gridGeometry(), bulkGeometry, center, shapeValues);
                    psData.addBulkInterpolation(shapeValues, this->vertexIndices(bulkIdx, bulkElementIdx), bulkElementIdx);
                }
                else
                {
                    psData.addBulkInterpolation(bulkElementIdx);
                }

                // publish point source data in the global vector
                this->pointSourceData().emplace_back(std::move(psData));

                const auto avgMinDist = averageDistanceSegmentGeometry(a, b, outside.geometry());
                this->averageDistanceToBulkCell().push_back(avgMinDist);
                fluxScalingFactor_.push_back(this->problem(bulkIdx).fluxScalingFactor(avgMinDist, radius, kernelWidth));

                // export the lowdim coupling stencil
                // we insert all vertices / elements and make it unique later
                if constexpr (isBox<bulkIdx>())
                {
                    const auto& vertices = this->vertexIndices(bulkIdx, bulkElementIdx);
                    this->couplingStencils(lowDimIdx)[lowDimElementIdx].insert(this->couplingStencils(lowDimIdx)[lowDimElementIdx].end(),
                                                                               vertices.begin(), vertices.end());
                }
                else
                {
                    this->couplingStencils(lowDimIdx)[lowDimElementIdx].push_back(bulkElementIdx);
                }
            }
        }

        // make the stencils unique
        makeUniqueStencil_();

        if (!this->pointSources(bulkIdx).empty())
            DUNE_THROW(Dune::InvalidStateException, "Kernel method shouldn't have point sources in the bulk domain but only volume sources!");

        std::cout << "[coupling] Finished preparing manager in " << watch.elapsed() << " seconds." << std::endl;
    }

    //! Compute the low dim volume fraction in the bulk domain cells
    void computeLowDimVolumeFractions()
    {
        // resize the storage vector
        lowDimVolumeInBulkElement_.resize(this->gridView(bulkIdx).size(0));
        // get references to the grid geometries
        const auto& lowDimGridGeometry = this->problem(lowDimIdx).gridGeometry();
        const auto& bulkGridGeometry = this->problem(bulkIdx).gridGeometry();

        // compute the low dim volume fractions
        for (const auto& is : intersections(this->glue()))
        {
            // all inside elements are identical...
            const auto& inside = is.targetEntity(0);
            const auto intersectionGeometry = is.geometry();
            const auto lowDimElementIdx = lowDimGridGeometry.elementMapper().index(inside);

            // compute the volume the low-dim domain occupies in the bulk domain if it were full-dimensional
            const auto radius = this->problem(lowDimIdx).spatialParams().radius(lowDimElementIdx);
            for (int outsideIdx = 0; outsideIdx < is.numDomainNeighbors(); ++outsideIdx)
            {
                const auto& outside = is.domainEntity(outsideIdx);
                const auto bulkElementIdx = bulkGridGeometry.elementMapper().index(outside);
                lowDimVolumeInBulkElement_[bulkElementIdx] += intersectionGeometry.volume()*M_PI*radius*radius;
            }
        }
    }

    /*!
     * \brief Methods to be accessed by the subproblems
     */
    // \{

    //! Return a reference to the bulk problem
    Scalar radius(std::size_t id) const
    {
        const auto& data = this->pointSourceData()[id];
        return this->problem(lowDimIdx).spatialParams().radius(data.lowDimElementIdx());
    }

    //! The volume the lower dimensional domain occupies in the bulk domain element
    // For one-dimensional low dim domain we assume radial tubes
    Scalar lowDimVolume(const Element<bulkIdx>& element) const
    {
        const auto eIdx = this->problem(bulkIdx).gridGeometry().elementMapper().index(element);
        return lowDimVolumeInBulkElement_[eIdx];
    }

    //! The volume fraction the lower dimensional domain occupies in the bulk domain element
    // For one-dimensional low dim domain we assume radial tubes
    Scalar lowDimVolumeFraction(const Element<bulkIdx>& element) const
    {
        const auto totalVolume = element.geometry().volume();
        return lowDimVolume(element) / totalVolume;
    }

    //! return all source ids for a bulk elements
    const std::vector<std::size_t>& bulkSourceIds(GridIndex<bulkIdx> eIdx, int scvIdx = 0) const
    { return bulkSourceIds_[eIdx][scvIdx]; }

    //! return all source ids for a bulk elements
    const std::vector<Scalar>& bulkSourceWeights(GridIndex<bulkIdx> eIdx, int scvIdx = 0) const
    { return bulkSourceWeights_[eIdx][scvIdx]; }

    //! The flux scaling factor for a source with id
    Scalar fluxScalingFactor(std::size_t id) const
    { return fluxScalingFactor_[id]; }

    // \}

private:
    //! compute the kernel distributed sources and add stencils
    template<class Line, class CylIntegration>
    Scalar computeBulkSource(const Line& line, const Scalar radius, const Scalar kernelWidth,
                             std::size_t id, GridIndex<lowDimIdx> lowDimElementIdx, GridIndex<bulkIdx> coupledBulkElementIdx,
                             const CylIntegration& cylIntegration, int embeddings)
    {
        // Monte-carlo integration on the cylinder defined by line and radius
        static const auto bulkParamGroup = this->problem(bulkIdx).paramGroup();
        static const auto min = getParamFromGroup<GlobalPosition>(bulkParamGroup, "Grid.LowerLeft");
        static const auto max = getParamFromGroup<GlobalPosition>(bulkParamGroup, "Grid.UpperRight");
        static const auto cells = getParamFromGroup<std::array<int, bulkDim>>(bulkParamGroup, "Grid.Cells");
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
            // TODO: below only works for Cartesian grids with ijk numbering (e.g. level 0 YaspGrid (already fails when refined))
            // more general is the bounding box tree solution which always works, however it's much slower
            //const auto bulkIndices = intersectingEntities(point, this->problem(bulkIdx).gridGeometry().boundingBoxTree(), true);
            if (const bool hasIntersection = intersectsPointBoundingBox(point, min, max); hasIntersection)
            {
                const auto bulkElementIdx = intersectingEntityCartesianGrid(point, min, max, cells);
                assert(bulkElementIdx < this->problem(bulkIdx).gridGeometry().gridView().size(0));
                const auto localWeight = evalConstKernel_(a, b, point, radius, kernelWidth)*cylIntegration.integrationElement(i)/Scalar(embeddings);
                integral += localWeight;
                if (!bulkSourceIds_[bulkElementIdx][0].empty() && id == bulkSourceIds_[bulkElementIdx][0].back())
                {
                    bulkSourceWeights_[bulkElementIdx][0].back() += localWeight;
                }
                else
                {
                    bulkSourceIds_[bulkElementIdx][0].emplace_back(id);
                    bulkSourceWeights_[bulkElementIdx][0].emplace_back(localWeight);
                    addBulkSourceStencils_(bulkElementIdx, lowDimElementIdx, coupledBulkElementIdx);
                }
            }
        }

        // the surface factor (which fraction of the source is inside the domain and needs to be considered)
        const auto length = (a-b).two_norm()/Scalar(embeddings);
        return integral/length;
    }

    void prepareDataStructures_()
    {
        // clear all internal members like pointsource vectors and stencils
        // initializes the point source id counter
        this->clear();
        bulkSourceIds_.clear();
        bulkSourceWeights_.clear();
        extendedSourceStencil_.stencil().clear();

        // precompute the vertex indices for efficiency for the box method
        this->precomputeVertexIndices(bulkIdx);
        this->precomputeVertexIndices(lowDimIdx);

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
            bulkSourceIds_[bulkElementIdx][0].reserve(10);
            bulkSourceWeights_[bulkElementIdx][0].reserve(10);
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

            // remove the vertices element (box)
            if constexpr (isBox<bulkIdx>())
            {
                const auto& indices = this->vertexIndices(bulkIdx, stencil.first);
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
        if constexpr (isBox<lowDimIdx>())
        {
            const auto& vertices = this->vertexIndices(lowDimIdx, coupledLowDimElementIdx);
            this->couplingStencils(bulkIdx)[bulkElementIdx].insert(this->couplingStencils(bulkIdx)[bulkElementIdx].end(),
                                                                   vertices.begin(), vertices.end());

        }
        else
        {
            auto& s = this->couplingStencils(bulkIdx)[bulkElementIdx];
            s.push_back(coupledLowDimElementIdx);
        }

        // the extended source stencil, every 3d element with a source is coupled to
        // the element/dofs where the 3d quantities are measured
        if constexpr (isBox<bulkIdx>())
        {
            const auto& vertices = this->vertexIndices(bulkIdx, coupledBulkElementIdx);
            extendedSourceStencil_.stencil()[bulkElementIdx].insert(extendedSourceStencil_.stencil()[bulkElementIdx].end(),
                                                                    vertices.begin(), vertices.end());
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
        const auto r = (proj - point).two_norm();

        if (r > rho)
            return 0.0;

        return 1.0/(M_PI*rho*rho);
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
    //! vector for the volume fraction of the lowdim domain in the bulk domain cells
    std::vector<Scalar> lowDimVolumeInBulkElement_;
    //! kernel sources to integrate for each bulk element
    std::vector<std::array<std::vector<std::size_t>, isBox<bulkIdx>() ? 1<<bulkDim : 1>> bulkSourceIds_;
    //! the integral of the kernel for each point source / integration point, i.e. weight for the source
    std::vector<std::array<std::vector<Scalar>, isBox<bulkIdx>() ? 1<<bulkDim : 1>> bulkSourceWeights_;

    //! the flux scaling factor for the respective source id
    std::vector<Scalar> fluxScalingFactor_;
};

} // end namespace Dumux

#endif
