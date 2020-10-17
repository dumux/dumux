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

#ifndef DUMUX_MULTIDOMAIN_EMBEDDED_COUPLINGMANAGER_1D3D_AVERAGE_HH
#define DUMUX_MULTIDOMAIN_EMBEDDED_COUPLINGMANAGER_1D3D_AVERAGE_HH

#include <vector>

#include <dune/common/timer.hh>
#include <dune/geometry/quadraturerules.hh>

#include <dumux/common/tag.hh>
#include <dumux/common/properties.hh>
#include <dumux/common/indextraits.hh>

#include <dumux/geometry/distance.hh>

#include <dumux/multidomain/embedded/couplingmanagerbase.hh>
#include <dumux/multidomain/embedded/circlepoints.hh>
#include <dumux/multidomain/embedded/extendedsourcestencil.hh>
#include <dumux/multidomain/embedded/circleaveragepointsourcetraits.hh>

namespace Dumux {

namespace Embedded1d3dCouplingMode {
struct Average : public Utility::Tag<Average> {
    static std::string name() { return "average"; }

    [[deprecated("Comparison with enum is deprecated. Removed after 3.3. Use tags.")]]
    friend constexpr bool operator==(Average, EmbeddedCouplingMode m) { return m == EmbeddedCouplingMode::average; }
    [[deprecated("Comparison with enum is deprecated. Removed after 3.3. Use tags.")]]
    friend constexpr bool operator==(EmbeddedCouplingMode m, Average) { return m == EmbeddedCouplingMode::average; }
    [[deprecated("Comparison with enum is deprecated. Removed after 3.3. Use tags.")]]
    friend constexpr bool operator!=(Average, EmbeddedCouplingMode m) { return m != EmbeddedCouplingMode::average; }
    [[deprecated("Comparison with enum is deprecated. Removed after 3.3. Use tags.")]]
    friend constexpr bool operator!=(EmbeddedCouplingMode m, Average) { return m != EmbeddedCouplingMode::average; }
};

inline constexpr Average average{};
} // end namespace Embedded1d3dCouplingMode

// forward declaration
template<class MDTraits, class CouplingMode>
class Embedded1d3dCouplingManager;

/*!
 * \ingroup EmbeddedCoupling
 * \brief Manages the coupling between bulk elements and lower dimensional elements
 *        Point sources on each integration point are computed by an AABB tree.
 * \note Specialization for coupling method using line sources with 3d quantities averaged on the cylinder surface
 */
template<class MDTraits>
class Embedded1d3dCouplingManager<MDTraits, Embedded1d3dCouplingMode::Average>
: public EmbeddedCouplingManagerBase<MDTraits, Embedded1d3dCouplingManager<MDTraits, Embedded1d3dCouplingMode::Average>,
                                     CircleAveragePointSourceTraits<MDTraits>>
{
    using ThisType = Embedded1d3dCouplingManager<MDTraits, Embedded1d3dCouplingMode::Average>;
    using ParentType = EmbeddedCouplingManagerBase<MDTraits, ThisType, CircleAveragePointSourceTraits<MDTraits>>;
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


public:
    enum {
        bulkDim = GridView<bulkIdx>::dimension,
        lowDimDim = GridView<lowDimIdx>::dimension,
        dimWorld = GridView<bulkIdx>::dimensionworld
    };

    static constexpr Embedded1d3dCouplingMode::Average couplingMode{};

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
        // Initialize the bulk bounding box tree
        const auto& bulkGridGeometry = this->problem(bulkIdx).gridGeometry();
        const auto& lowDimGridGeometry = this->problem(lowDimIdx).gridGeometry();
        const auto& bulkTree = bulkGridGeometry.boundingBoxTree();

        // initialize the maps
        // do some logging and profiling
        Dune::Timer watch;
        std::cout << "Initializing the point sources..." << std::endl;

        // clear all internal members like pointsource vectors and stencils
        // initializes the point source id counter
        this->clear();
        extendedSourceStencil_.clear();

        // precompute the vertex indices for efficiency
        this->precomputeVertexIndices(bulkIdx);
        this->precomputeVertexIndices(lowDimIdx);

        // intersect the bounding box trees
        this->glueGrids();

        // iterate over all intersection and add point sources
        const auto& lowDimProblem = this->problem(lowDimIdx);
        for (const auto& is : intersections(this->glue()))
        {
            // all inside elements are identical...
            const auto& lowDimElement = is.targetEntity(0);
            const auto lowDimElementIdx = lowDimGridGeometry.elementMapper().index(lowDimElement);

            // get the intersection geometry
            const auto intersectionGeometry = is.geometry();
            // get the Gaussian quadrature rule for the local intersection
            const auto& quad = Dune::QuadratureRules<Scalar, lowDimDim>::rule(intersectionGeometry.type(), order);

            // apply the Gaussian quadrature rule and define point sources at each quadrature point
            // note that the approximation is not optimal if
            // (a) the one-dimensional elements are too large,
            // (b) whenever a one-dimensional element is split between two or more elements,
            // (c) when gradients of important quantities in the three-dimensional domain are large.

            // iterate over all quadrature points
            for (auto&& qp : quad)
            {
                // global position of the quadrature point
                const auto globalPos = intersectionGeometry.global(qp.position());

                const auto bulkElementIndices = intersectingEntities(globalPos, bulkTree);

                // do not add a point source if the qp is outside of the 3d grid
                // this is equivalent to having a source of zero for that qp
                if (bulkElementIndices.empty())
                    continue;

                //////////////////////////////////////////////////////////
                // get circle average connectivity and interpolation data
                //////////////////////////////////////////////////////////

                static const auto numIp = getParam<int>("MixedDimension.NumCircleSegments");
                const auto radius = lowDimProblem.spatialParams().radius(lowDimElementIdx);
                const auto normal = intersectionGeometry.corner(1)-intersectionGeometry.corner(0);
                const auto circleAvgWeight = 2*M_PI*radius/numIp;
                const auto integrationElement = intersectionGeometry.integrationElement(qp.position());
                const auto qpweight = qp.weight();

                const auto circlePoints = EmbeddedCoupling::circlePoints(globalPos, normal, radius, numIp);
                std::vector<Scalar> circleIpWeight; circleIpWeight.reserve(circlePoints.size());
                std::vector<GridIndex<bulkIdx>> circleStencil; circleStencil.reserve(circlePoints.size());

                // for box
                std::vector<const std::vector<GridIndex<bulkIdx>>*> circleCornerIndices;
                using ShapeValues = std::vector<Dune::FieldVector<Scalar, 1> >;
                std::vector<ShapeValues> circleShapeValues;

                // go over all points of the average operator
                int insideCirclePoints = 0;
                for (int k = 0; k < circlePoints.size(); ++k)
                {
                    const auto circleBulkElementIndices = intersectingEntities(circlePoints[k], bulkTree);
                    if (circleBulkElementIndices.empty())
                        continue;

                    ++insideCirclePoints;
                    const auto localCircleAvgWeight = circleAvgWeight / circleBulkElementIndices.size();
                    for (const auto bulkElementIdx : circleBulkElementIndices)
                    {
                        circleStencil.push_back(bulkElementIdx);
                        circleIpWeight.push_back(localCircleAvgWeight);

                        // precompute interpolation data for box scheme for each cut bulk element
                        if constexpr (isBox<bulkIdx>())
                        {
                            const auto bulkElement = bulkGridGeometry.element(bulkElementIdx);
                            circleCornerIndices.push_back(&(this->vertexIndices(bulkIdx, bulkElementIdx)));

                            // evaluate shape functions at the integration point
                            const auto bulkGeometry = bulkElement.geometry();
                            ShapeValues shapeValues;
                            this->getShapeValues(bulkIdx, bulkGridGeometry, bulkGeometry, circlePoints[k], shapeValues);
                            circleShapeValues.emplace_back(std::move(shapeValues));
                        }
                    }
                }

                // export low dim circle stencil
                if constexpr (isBox<bulkIdx>())
                {
                    // we insert all vertices and make it unique later
                    for (const auto& vertices : circleCornerIndices)
                    {
                        this->couplingStencils(lowDimIdx)[lowDimElementIdx].insert(this->couplingStencils(lowDimIdx)[lowDimElementIdx].end(),
                                                                                   vertices->begin(), vertices->end());

                    }
                }
                else
                {
                    this->couplingStencils(lowDimIdx)[lowDimElementIdx].insert(this->couplingStencils(lowDimIdx)[lowDimElementIdx].end(),
                                                                               circleStencil.begin(), circleStencil.end());
                }

                // surface fraction that is inside the domain (only different than 1.0 on the boundary)
                const auto surfaceFraction = Scalar(insideCirclePoints)/Scalar(circlePoints.size());

                // loop over the bulk elements at the integration points (usually one except when it is on a face or edge or vertex)
                for (auto bulkElementIdx : bulkElementIndices)
                {
                    const auto id = this->idCounter_++;

                    this->pointSources(bulkIdx).emplace_back(globalPos, id, qpweight, integrationElement*surfaceFraction, bulkElementIdx);
                    this->pointSources(bulkIdx).back().setEmbeddings(bulkElementIndices.size());
                    this->pointSources(lowDimIdx).emplace_back(globalPos, id, qpweight, integrationElement*surfaceFraction, lowDimElementIdx);
                    this->pointSources(lowDimIdx).back().setEmbeddings(bulkElementIndices.size());

                    // pre compute additional data used for the evaluation of
                    // the actual solution dependent source term
                    PointSourceData psData;

                    if constexpr (isBox<lowDimIdx>())
                    {
                        ShapeValues shapeValues;
                        this->getShapeValues(lowDimIdx, lowDimGridGeometry, intersectionGeometry, globalPos, shapeValues);
                        psData.addLowDimInterpolation(shapeValues, this->vertexIndices(lowDimIdx, lowDimElementIdx), lowDimElementIdx);
                    }
                    else
                    {
                        psData.addLowDimInterpolation(lowDimElementIdx);
                    }

                    // add data needed to compute integral over the circle
                    if constexpr (isBox<bulkIdx>())
                    {
                        psData.addCircleInterpolation(circleCornerIndices, circleShapeValues, circleIpWeight, circleStencil);

                        const auto bulkGeometry = bulkGridGeometry.element(bulkElementIdx).geometry();
                        ShapeValues shapeValues;
                        this->getShapeValues(bulkIdx, bulkGridGeometry, bulkGeometry, globalPos, shapeValues);
                        psData.addBulkInterpolation(shapeValues, this->vertexIndices(bulkIdx, bulkElementIdx), bulkElementIdx);
                    }
                    else
                    {
                        psData.addCircleInterpolation(circleIpWeight, circleStencil);
                        psData.addBulkInterpolation(bulkElementIdx);
                    }

                    // publish point source data in the global vector
                    this->pointSourceData().emplace_back(std::move(psData));

                    // mean distance to outside element for source correction schemes
                    const auto outsideGeometry = bulkGridGeometry.element(bulkElementIdx).geometry();
                    this->averageDistanceToBulkCell().push_back(averageDistancePointGeometry(globalPos, outsideGeometry));

                    // export the bulk coupling stencil
                    if constexpr (isBox<lowDimIdx>())
                    {
                        this->couplingStencils(bulkIdx)[bulkElementIdx].insert(this->couplingStencils(bulkIdx)[bulkElementIdx].end(),
                                                                               this->vertexIndices(lowDimIdx, lowDimElementIdx).begin(),
                                                                               this->vertexIndices(lowDimIdx, lowDimElementIdx).end());

                    }
                    else
                    {
                        this->couplingStencils(bulkIdx)[bulkElementIdx].push_back(lowDimElementIdx);
                    }

                    // export bulk circle stencil
                    if constexpr (isBox<bulkIdx>())
                    {
                        // we insert all vertices and make it unique later
                        for (const auto& vertices : circleCornerIndices)
                        {
                            extendedSourceStencil_.stencil()[bulkElementIdx].insert(extendedSourceStencil_.stencil()[bulkElementIdx].end(),
                                                                                    vertices->begin(), vertices->end());

                        }
                    }
                    else
                    {
                        extendedSourceStencil_.stencil()[bulkElementIdx].insert(extendedSourceStencil_.stencil()[bulkElementIdx].end(),
                                                                   circleStencil.begin(), circleStencil.end());
                    }
                }
            }
        }

        // make the circle stencil unique (for source derivatives)
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

    // \}

private:
    //! the extended source stencil object
    EmbeddedCoupling::ExtendedSourceStencil<ThisType> extendedSourceStencil_;

    //! vector for the volume fraction of the lowdim domain in the bulk domain cells
    std::vector<Scalar> lowDimVolumeInBulkElement_;
};

} // end namespace Dumux

#endif
