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
 * \note the kernel per point source is isotropic and it's integral over the domain is one
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

    enum {
        bulkDim = GridView<bulkIdx>::dimension,
        lowDimDim = GridView<lowDimIdx>::dimension,
        dimWorld = GridView<bulkIdx>::dimensionworld
    };
public:
    static constexpr Embedded1d3dCouplingMode::Kernel couplingMode{};

    using ParentType::ParentType;

    void init(std::shared_ptr<Problem<bulkIdx>> bulkProblem,
              std::shared_ptr<Problem<lowDimIdx>> lowDimProblem,
              const SolutionVector& curSol)
    {
        ParentType::init(bulkProblem, lowDimProblem, curSol);
        computeLowDimVolumeFractions();
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
        std::cout << "Initializing the point sources..." << std::endl;

        // clear all internal members like pointsource vectors and stencils
        // initializes the point source id counter
        this->clear();
        bulkSourceIds_.clear();
        bulkSourceWeights_.clear();
        extendedSourceStencil_.clear();

        // precompute the vertex indices for efficiency for the box method
        this->precomputeVertexIndices(bulkIdx);
        this->precomputeVertexIndices(lowDimIdx);

        const auto& bulkGridGeometry = this->problem(bulkIdx).gridGeometry();
        const auto& lowDimGridGeometry = this->problem(lowDimIdx).gridGeometry();

        bulkSourceIds_.resize(this->gridView(bulkIdx).size(0));
        bulkSourceWeights_.resize(this->gridView(bulkIdx).size(0));

        // intersect the bounding box trees
        this->glueGrids();

        this->pointSourceData().reserve(this->glue().size());
        this->averageDistanceToBulkCell().reserve(this->glue().size());
        const Scalar kernelWidth = getParam<Scalar>("MixedDimension.KernelWidth");
        for (const auto& is : intersections(this->glue()))
        {
            // all inside elements are identical...
            const auto& inside = is.targetEntity(0);
            // get the intersection geometry for integrating over it
            const auto intersectionGeometry = is.geometry();

            // get the Gaussian quadrature rule for the local intersection
            const auto& quad = Dune::QuadratureRules<Scalar, lowDimDim>::rule(intersectionGeometry.type(), order);
            const auto lowDimElementIdx = lowDimGridGeometry.elementMapper().index(inside);

            // iterate over all quadrature points and place a source
            // for 1d: make a new point source
            // for 3d: make a new kernel volume source
            for (auto&& qp : quad)
            {
                // compute the coupling stencils
                for (int outsideIdx = 0; outsideIdx < is.numDomainNeighbors(); ++outsideIdx)
                {
                    const auto& outside = is.domainEntity(outsideIdx);
                    const auto bulkElementIdx = bulkGridGeometry.elementMapper().index(outside);

                    // each quadrature point will be a point source for the sub problem
                    const auto globalPos = intersectionGeometry.global(qp.position());
                    const auto id = this->idCounter_++;
                    const auto qpweight = qp.weight();
                    const auto ie = intersectionGeometry.integrationElement(qp.position());
                    this->pointSources(lowDimIdx).emplace_back(globalPos, id, qpweight, ie, std::vector<std::size_t>({lowDimElementIdx}));
                    this->pointSources(lowDimIdx).back().setEmbeddings(is.numDomainNeighbors());
                    computeBulkSource(globalPos, kernelWidth, id, lowDimElementIdx, bulkElementIdx, qpweight*ie/is.numDomainNeighbors());

                    // pre compute additional data used for the evaluation of
                    // the actual solution dependent source term
                    PointSourceData psData;

                    if constexpr (isBox<lowDimIdx>())
                    {
                        using ShapeValues = std::vector<Dune::FieldVector<Scalar, 1> >;
                        const auto lowDimGeometry = lowDimGridGeometry.element(lowDimElementIdx).geometry();
                        ShapeValues shapeValues;
                        this->getShapeValues(lowDimIdx, lowDimGridGeometry, lowDimGeometry, globalPos, shapeValues);
                        psData.addLowDimInterpolation(shapeValues, this->vertexIndices(lowDimIdx, lowDimElementIdx), lowDimElementIdx);
                    }
                    else
                    {
                        psData.addLowDimInterpolation(lowDimElementIdx);
                    }

                    // add data needed to compute integral over the circle
                    if constexpr (isBox<bulkIdx>())
                    {
                        using ShapeValues = std::vector<Dune::FieldVector<Scalar, 1> >;
                        const auto bulkGeometry = bulkGridGeometry.element(bulkElementIdx).geometry();
                        ShapeValues shapeValues;
                        this->getShapeValues(bulkIdx, bulkGridGeometry, bulkGeometry, globalPos, shapeValues);
                        psData.addBulkInterpolation(shapeValues, this->vertexIndices(bulkIdx, bulkElementIdx), bulkElementIdx);
                    }
                    else
                    {
                        psData.addBulkInterpolation(bulkElementIdx);
                    }

                    // publish point source data in the global vector
                    this->pointSourceData().emplace_back(std::move(psData));

                    // compute average distance to bulk cell
                    this->averageDistanceToBulkCell().push_back(averageDistancePointGeometry(globalPos, outside.geometry()));

                    // export the lowdim coupling stencil
                    // we insert all vertices / elements and make it unique later
                    if (isBox<bulkIdx>())
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
        }

        // make extra stencils unique
        for (auto&& stencil : extendedSourceStencil_.stencil())
        {
            std::sort(stencil.second.begin(), stencil.second.end());
            stencil.second.erase(std::unique(stencil.second.begin(), stencil.second.end()), stencil.second.end());

            // remove the vertices element (box)
            if (isBox<bulkIdx>())
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

        if (!this->pointSources(bulkIdx).empty())
            DUNE_THROW(Dune::InvalidStateException, "Kernel method shouldn't have point sources in the bulk domain but only volume sources!");

        std::cout << "took " << watch.elapsed() << " seconds." << std::endl;
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
    const std::vector<std::size_t> bulkSourceIds(GridIndex<bulkIdx> eIdx) const
    { return bulkSourceIds_[eIdx]; }

    //! return all source ids for a bulk elements
    const std::vector<Scalar> bulkSourceWeights(GridIndex<bulkIdx> eIdx) const
    { return bulkSourceWeights_[eIdx]; }

    // \}

private:
    //! TODO: How to optimize this?
    void computeBulkSource(const GlobalPosition& globalPos, const Scalar kernelWidth,
                           std::size_t id, GridIndex<lowDimIdx> lowDimElementIdx, GridIndex<bulkIdx> coupledBulkElementIdx,
                           Scalar pointSourceWeight)
    {
        // make sure it is mass conservative
        // i.e. the point source in the 1d domain needs to have the exact same integral as the distributed
        // kernel source integral in the 3d domain. Correct the integration formula by balancing the error
        // by scaling the kernel with the checksum inverse
        Scalar checkSum = 0.0;
        std::vector<bool> mask(this->gridView(bulkIdx).size(0), false);
        for (const auto& element : elements(this->gridView(bulkIdx)))
        {
            Scalar weight = 0.0;
            const auto geometry = element.geometry();
            const auto& quad = Dune::QuadratureRules<Scalar, bulkDim>::rule(geometry.type(), 3);
            for (auto&& qp : quad)
            {
                const auto qpweight = qp.weight();
                const auto ie = geometry.integrationElement(qp.position());
                weight += evalKernel(globalPos, geometry.global(qp.position()), kernelWidth)*qpweight*ie;
            }

            if (weight > 1e-13)
            {
                const auto& bulkGridGeometry = this->problem(bulkIdx).gridGeometry();
                const auto bulkElementIdx = bulkGridGeometry.elementMapper().index(element);
                bulkSourceIds_[bulkElementIdx].push_back(id);
                bulkSourceWeights_[bulkElementIdx].push_back(weight*pointSourceWeight);
                mask[bulkElementIdx] = true;

                // add lowDim dofs that the source is related to to the bulk stencil
                if (isBox<lowDimIdx>())
                {
                    const auto& vertices = this->vertexIndices(lowDimIdx, lowDimElementIdx);
                    this->couplingStencils(bulkIdx)[bulkElementIdx].insert(this->couplingStencils(bulkIdx)[bulkElementIdx].end(),
                                                                           vertices.begin(), vertices.end());

                }
                else
                {
                    this->couplingStencils(bulkIdx)[bulkElementIdx].push_back(lowDimElementIdx);
                }

                // tpfa
                extendedSourceStencil_.stencil()[bulkElementIdx].push_back(coupledBulkElementIdx);

                // compute check sum -> should sum up to 1.0 to be mass conservative
                checkSum += weight;
            }
        }

        for (GridIndex<bulkIdx> eIdx = 0; eIdx < bulkSourceWeights_.size(); ++eIdx)
            if (mask[eIdx]) bulkSourceWeights_[eIdx].back() /= checkSum;

        // balance error of the quadrature rule -> TODO: what to do at boundaries
        // const auto diff = 1.0 - checkSum;
        // std::cout << "Integrated kernel with integration error of " << diff << std::endl;
    }

    //! an isotropic cubic kernel with derivatives 0 at r=origin and r=width and domain integral 1
    inline Scalar evalKernel(const GlobalPosition& origin,
                             const GlobalPosition& pos,
                             const Scalar width) const noexcept
    {
        const auto r = (pos-origin).two_norm();
        const auto r2 = r*r;
        const auto r3 = r2*r;

        if (r > width)
            return 0.0;

        const Scalar w2 = width*width;
        const Scalar w3 = w2*width;
        const Scalar k = 15.0/(4*M_PI*w3);
        const Scalar a = 2.0/w3;
        const Scalar b = 3.0/w2;

        return k*(a*r3 - b*r2 + 1.0);
    }

    //! the extended source stencil object for kernel coupling
    EmbeddedCoupling::ExtendedSourceStencil<ThisType> extendedSourceStencil_;
    //! vector for the volume fraction of the lowdim domain in the bulk domain cells
    std::vector<Scalar> lowDimVolumeInBulkElement_;
    //! kernel sources to integrate for each bulk element
    std::vector<std::vector<std::size_t>> bulkSourceIds_;
    //! the integral of the kernel for each point source / integration point, i.e. weight for the source
    std::vector<std::vector<Scalar>> bulkSourceWeights_;
};

} // end namespace Dumux

#endif
