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

#ifndef DUMUX_MULTIDOMAIN_EMBEDDED_COUPLINGMANAGER_1D3D_SURFACE_HH
#define DUMUX_MULTIDOMAIN_EMBEDDED_COUPLINGMANAGER_1D3D_SURFACE_HH

#include <vector>

#include <dune/common/timer.hh>
#include <dune/geometry/quadraturerules.hh>

#include <dumux/common/tag.hh>
#include <dumux/common/properties.hh>
#include <dumux/common/indextraits.hh>

#include <dumux/multidomain/embedded/couplingmanagerbase.hh>
#include <dumux/multidomain/embedded/circlepoints.hh>

namespace Dumux {

namespace Embedded1d3dCouplingMode {
struct Surface : public Utility::Tag<Surface> {
    static std::string name() { return "surface"; }

    [[deprecated("Comparison with enum is deprecated. Removed after 3.3. Use tags.")]]
    friend constexpr bool operator==(Surface, EmbeddedCouplingMode m) { return m == EmbeddedCouplingMode::cylindersources; }
    [[deprecated("Comparison with enum is deprecated. Removed after 3.3. Use tags.")]]
    friend constexpr bool operator==(EmbeddedCouplingMode m, Surface) { return m == EmbeddedCouplingMode::cylindersources; }
    [[deprecated("Comparison with enum is deprecated. Removed after 3.3. Use tags.")]]
    friend constexpr bool operator!=(Surface, EmbeddedCouplingMode m) { return m != EmbeddedCouplingMode::cylindersources; }
    [[deprecated("Comparison with enum is deprecated. Removed after 3.3. Use tags.")]]
    friend constexpr bool operator!=(EmbeddedCouplingMode m, Surface) { return m != EmbeddedCouplingMode::cylindersources; }
};

inline constexpr Surface surface{};
} // end namespace Embedded1d3dCouplingMode

// forward declaration
template<class MDTraits, class CouplingMode>
class Embedded1d3dCouplingManager;

/*!
 * \ingroup EmbeddedCoupling
 * \brief Manages the coupling between bulk elements and lower dimensional elements
 *        Point sources on each integration point are computed by an AABB tree.
 * \note Specialization for coupling method using cylinder sources with 3d quantities evaluated on the cylinder surface
 * \note the cylinder source is approximated by point sources on the cylinder surface
 */
template<class MDTraits>
class Embedded1d3dCouplingManager<MDTraits, Embedded1d3dCouplingMode::Surface>
: public EmbeddedCouplingManagerBase<MDTraits, Embedded1d3dCouplingManager<MDTraits, Embedded1d3dCouplingMode::Surface>>
{
    using ThisType = Embedded1d3dCouplingManager<MDTraits, Embedded1d3dCouplingMode::Surface>;
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
    static constexpr Embedded1d3dCouplingMode::Surface couplingMode{};

    using ParentType::ParentType;

    void init(std::shared_ptr<Problem<bulkIdx>> bulkProblem,
              std::shared_ptr<Problem<lowDimIdx>> lowDimProblem,
              const SolutionVector& curSol)
    {
        ParentType::init(bulkProblem, lowDimProblem, curSol);
        computeLowDimVolumeFractions();
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

        // precompute the vertex indices for efficiency
        this->precomputeVertexIndices(bulkIdx);
        this->precomputeVertexIndices(lowDimIdx);

        // iterate over all lowdim elements
        const auto& lowDimProblem = this->problem(lowDimIdx);
        for (const auto& lowDimElement : elements(this->gridView(lowDimIdx)))
        {
            // get the Gaussian quadrature rule for the low dim element
            const auto lowDimGeometry = lowDimElement.geometry();
            const auto& quad = Dune::QuadratureRules<Scalar, lowDimDim>::rule(lowDimGeometry.type(), order);

            const auto lowDimElementIdx = lowDimGridGeometry.elementMapper().index(lowDimElement);

            // apply the Gaussian quadrature rule and define point sources at each quadrature point
            // note that the approximation is not optimal if
            // (a) the one-dimensional elements are too large,
            // (b) whenever a one-dimensional element is split between two or more elements,
            // (c) when gradients of important quantities in the three-dimensional domain are large.

            // iterate over all quadrature points
            for (auto&& qp : quad)
            {
                // global position of the quadrature point
                const auto globalPos = lowDimGeometry.global(qp.position());

                const auto bulkElementIndices = intersectingEntities(globalPos, bulkTree);

                // do not add a point source if the qp is outside of the 3d grid
                // this is equivalent to having a source of zero for that qp
                if (bulkElementIndices.empty())
                    continue;

                ////////////////////////////////////////////////////////////////
                // get points on the cylinder surface at the integration point
                ////////////////////////////////////////////////////////////////

                static const auto numIp = getParam<int>("MixedDimension.NumCircleSegments");
                const auto radius = lowDimProblem.spatialParams().radius(lowDimElementIdx);
                const auto normal = lowDimGeometry.corner(1)-lowDimGeometry.corner(0);
                const auto integrationElement = lowDimGeometry.integrationElement(qp.position())*2*M_PI*radius/Scalar(numIp);
                const auto weight = qp.weight()/(2*M_PI*radius);

                const auto circlePoints = EmbeddedCoupling::circlePoints(globalPos, normal, radius, numIp);
                for (int k = 0; k < circlePoints.size(); ++k)
                {
                    const auto& circlePos = circlePoints[k];
                    const auto circleBulkElementIndices = intersectingEntities(circlePos, bulkTree);
                    if (circleBulkElementIndices.empty())
                        continue;

                    // loop over the bulk elements at the integration points (usually one except when it is on a face or edge or vertex)
                    // and add a point source at every point on the circle
                    for (const auto bulkElementIdx : circleBulkElementIndices)
                    {
                        const auto id = this->idCounter_++;
                        this->pointSources(bulkIdx).emplace_back(circlePos, id, weight, integrationElement, std::vector<std::size_t>({bulkElementIdx}));
                        this->pointSources(bulkIdx).back().setEmbeddings(circleBulkElementIndices.size());
                        this->pointSources(lowDimIdx).emplace_back(globalPos, id, weight, integrationElement, std::vector<std::size_t>({lowDimElementIdx}));
                        this->pointSources(lowDimIdx).back().setEmbeddings(circleBulkElementIndices.size());

                        // pre compute additional data used for the evaluation of
                        // the actual solution dependent source term
                        PointSourceData psData;

                        if constexpr (isBox<lowDimIdx>())
                        {
                            using ShapeValues = std::vector<Dune::FieldVector<Scalar, 1> >;
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
                            this->getShapeValues(bulkIdx, bulkGridGeometry, bulkGeometry, circlePos, shapeValues);
                            psData.addBulkInterpolation(shapeValues, this->vertexIndices(bulkIdx, bulkElementIdx), bulkElementIdx);
                        }
                        else
                        {
                            psData.addBulkInterpolation(bulkElementIdx);
                        }

                        // publish point source data in the global vector
                        this->pointSourceData().emplace_back(std::move(psData));

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

                        // export the bulk coupling stencil
                        // we insert all vertices / elements and make it unique later
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

                    }
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

    //! vector for the volume fraction of the lowdim domain in the bulk domain cells
    std::vector<Scalar> lowDimVolumeInBulkElement_;
};

} // end namespace Dumux

#endif
