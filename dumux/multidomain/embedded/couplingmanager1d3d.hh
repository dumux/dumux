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
 * \brief Coupling manager for low-dimensional domains embedded in the bulk
 *        domain. Point sources on each integration point are computed by an AABB tree.
 */

#ifndef DUMUX_MULTIDOMAIN_EMBEDDED_COUPLINGMANAGER_1D3D_HH
#define DUMUX_MULTIDOMAIN_EMBEDDED_COUPLINGMANAGER_1D3D_HH

#include <random>
#include <vector>

#include <dune/common/timer.hh>
#include <dune/common/exceptions.hh>
#include <dune/geometry/quadraturerules.hh>

#include <dumux/common/properties.hh>
#include <dumux/common/math.hh>
#include <dumux/common/integrator.hh>
#include <dumux/multidomain/embedded/pointsourcedata.hh>
#include <dumux/multidomain/embedded/integrationpointsource.hh>
#include <dumux/multidomain/embedded/couplingmanagerbase.hh>
#include <dumux/multidomain/embedded/circlepoints.hh>
#include <dumux/multidomain/embedded/extendedsourcestencil.hh>
#include <dumux/multidomain/embedded/cylinderintegration.hh>

namespace Dumux {

/*!
 * \ingroup EmbeddedCoupling
 * \brief The coupling mode
 */
enum class EmbeddedCouplingMode
{ line, average, cylindersources, kernel };

//! point source traits for the circle average coupling mode
template<class MDTraits>
struct CircleAveragePointSourceTraits
{
private:
    template<std::size_t i> using SubDomainTypeTag = typename MDTraits::template SubDomain<i>::TypeTag;
    template<std::size_t i> using FVGridGeometry = GetPropType<SubDomainTypeTag<i>, Properties::FVGridGeometry>;
    template<std::size_t i> using NumEqVector = GetPropType<SubDomainTypeTag<i>, Properties::NumEqVector>;
public:
    //! export the point source type for domain i
    template<std::size_t i>
    using PointSource = IntegrationPointSource<typename FVGridGeometry<i>::GlobalCoordinate, NumEqVector<i>>;

    //! export the point source helper type  for domain i
    template<std::size_t i>
    using PointSourceHelper = IntegrationPointSourceHelper;

    //! export the point source data type
    using PointSourceData = PointSourceDataCircleAverage<MDTraits>;
};

/*!
 * \ingroup EmbeddedCoupling
 * \brief Manages the coupling between bulk elements and lower dimensional elements
 *        Point sources on each integration point are computed by an AABB tree.
 */
template<class MDTraits, EmbeddedCouplingMode mode>
class EmbeddedCouplingManager1d3d;

/*!
 * \ingroup EmbeddedCoupling
 * \brief Manages the coupling between bulk elements and lower dimensional elements
 *        Point sources on each integration point are computed by an AABB tree.
 * \note Specialization for coupling method using line sources with 3d quantities evaluated on the line
 * \note This is the simplest method but it mathematically not well defined as the 3d quantity is evaluated
 *       where the solution to the continuous problem has a singularity
 */
template<class MDTraits>
class EmbeddedCouplingManager1d3d<MDTraits, EmbeddedCouplingMode::line>
: public EmbeddedCouplingManagerBase<MDTraits, EmbeddedCouplingManager1d3d<MDTraits, EmbeddedCouplingMode::line>>
{
    using ThisType = EmbeddedCouplingManager1d3d<MDTraits, EmbeddedCouplingMode::line>;
    using ParentType = EmbeddedCouplingManagerBase<MDTraits, ThisType>;

    using Scalar = typename MDTraits::Scalar;
    using SolutionVector = typename MDTraits::SolutionVector;

    static constexpr auto bulkIdx = typename MDTraits::template SubDomain<0>::Index();
    static constexpr auto lowDimIdx = typename MDTraits::template SubDomain<1>::Index();

    // the sub domain type aliases
    template<std::size_t id> using SubDomainTypeTag = typename MDTraits::template SubDomain<id>::TypeTag;
    template<std::size_t id> using Problem = GetPropType<SubDomainTypeTag<id>, Properties::Problem>;
    template<std::size_t id> using FVGridGeometry = GetPropType<SubDomainTypeTag<id>, Properties::FVGridGeometry>;
    template<std::size_t id> using GridView = typename FVGridGeometry<id>::GridView;
    template<std::size_t id> using Element = typename GridView<id>::template Codim<0>::Entity;

public:
    static constexpr EmbeddedCouplingMode couplingMode = EmbeddedCouplingMode::line;

    using ParentType::ParentType;

    void init(std::shared_ptr<Problem<bulkIdx>> bulkProblem,
              std::shared_ptr<Problem<lowDimIdx>> lowDimProblem,
              const SolutionVector& curSol)
    {
        ParentType::init(bulkProblem, lowDimProblem, curSol);
        computeLowDimVolumeFractions();
    }

    //! Compute the low dim volume fraction in the bulk domain cells
    void computeLowDimVolumeFractions()
    {
        // resize the storage vector
        lowDimVolumeInBulkElement_.resize(this->gridView(bulkIdx).size(0));

        // compute the low dim volume fractions
        for (const auto& is : intersections(this->glue()))
        {
            // all inside elements are identical...
            const auto& inside = is.inside(0);
            const auto intersectionGeometry = is.geometry();
            const std::size_t lowDimElementIdx = this->problem(lowDimIdx).fvGridGeometry().elementMapper().index(inside);

            // compute the volume the low-dim domain occupies in the bulk domain if it were full-dimensional
            const auto radius = this->problem(lowDimIdx).spatialParams().radius(lowDimElementIdx);
            for (std::size_t outsideIdx = 0; outsideIdx < is.neighbor(0); ++outsideIdx)
            {
                const auto& outside = is.outside(outsideIdx);
                const std::size_t bulkElementIdx = this->problem(bulkIdx).fvGridGeometry().elementMapper().index(outside);
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
        const auto eIdx = this->problem(bulkIdx).fvGridGeometry().elementMapper().index(element);
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

/*!
 * \ingroup EmbeddedCoupling
 * \brief Manages the coupling between bulk elements and lower dimensional elements
 *        Point sources on each integration point are computed by an AABB tree.
 * \note Specialization for coupling method using line sources with 3d quantities averaged on the cylinder surface
 */
template<class MDTraits>
class EmbeddedCouplingManager1d3d<MDTraits, EmbeddedCouplingMode::average>
: public EmbeddedCouplingManagerBase<MDTraits, EmbeddedCouplingManager1d3d<MDTraits, EmbeddedCouplingMode::average>,
                                     CircleAveragePointSourceTraits<MDTraits>>
{
    using ThisType = EmbeddedCouplingManager1d3d<MDTraits, EmbeddedCouplingMode::average>;
    using ParentType = EmbeddedCouplingManagerBase<MDTraits, ThisType, CircleAveragePointSourceTraits<MDTraits>>;
    using Scalar = typename MDTraits::Scalar;
    using SolutionVector = typename MDTraits::SolutionVector;
    using PointSourceData = typename ParentType::PointSourceTraits::PointSourceData;

    static constexpr auto bulkIdx = typename MDTraits::template SubDomain<0>::Index();
    static constexpr auto lowDimIdx = typename MDTraits::template SubDomain<1>::Index();

    // the sub domain type aliases
    template<std::size_t id> using SubDomainTypeTag = typename MDTraits::template SubDomain<id>::TypeTag;
    template<std::size_t id> using Problem = GetPropType<SubDomainTypeTag<id>, Properties::Problem>;
    template<std::size_t id> using FVGridGeometry = GetPropType<SubDomainTypeTag<id>, Properties::FVGridGeometry>;
    template<std::size_t id> using GridView = typename FVGridGeometry<id>::GridView;
    template<std::size_t id> using Element = typename GridView<id>::template Codim<0>::Entity;

    using GlobalPosition = typename Element<bulkIdx>::Geometry::GlobalCoordinate;

    template<std::size_t id>
    static constexpr bool isBox()
    { return FVGridGeometry<id>::discMethod == DiscretizationMethod::box; }


public:
    enum {
        bulkDim = GridView<bulkIdx>::dimension,
        lowDimDim = GridView<lowDimIdx>::dimension,
        dimWorld = GridView<bulkIdx>::dimensionworld
    };

    static constexpr EmbeddedCouplingMode couplingMode = EmbeddedCouplingMode::average;

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
        const auto& bulkTree = this->problem(bulkIdx).fvGridGeometry().boundingBoxTree();

        // initilize the maps
        // do some logging and profiling
        Dune::Timer watch;
        std::cout << "Initializing the point sources..." << std::endl;

        // clear all internal members like pointsource vectors and stencils
        // initializes the point source id counter
        this->clear();
        extendedSourceStencil_.stencil().clear();

        // precompute the vertex indices for efficiency
        this->preComputeVertexIndices(bulkIdx);
        this->preComputeVertexIndices(lowDimIdx);

        // iterate over all lowdim elements
        const auto& lowDimProblem = this->problem(lowDimIdx);
        for (const auto& lowDimElement : elements(this->gridView(lowDimIdx)))
        {
            // get the Gaussian quadrature rule for the low dim element
            const auto lowDimGeometry = lowDimElement.geometry();
            const auto& quad = Dune::QuadratureRules<Scalar, lowDimDim>::rule(lowDimGeometry.type(), order);

            const auto lowDimElementIdx = lowDimProblem.fvGridGeometry().elementMapper().index(lowDimElement);

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

                //////////////////////////////////////////////////////////
                // get circle average connectivity and interpolation data
                //////////////////////////////////////////////////////////

                static const auto numIp = getParam<int>("MixedDimension.NumCircleSegments");
                const auto radius = lowDimProblem.spatialParams().radius(lowDimElementIdx);
                const auto normal = lowDimGeometry.corner(1)-lowDimGeometry.corner(0);
                const auto weight = 2*M_PI*radius/numIp;

                const auto circlePoints = EmbeddedCoupling::circlePoints(globalPos, normal, radius, numIp);
                std::vector<Scalar> circleIpWeight(circlePoints.size());
                std::vector<std::size_t> circleStencil(circlePoints.size());
                // for box
                std::unordered_map<std::size_t, std::vector<std::size_t> > circleCornerIndices;
                using ShapeValues = std::vector<Dune::FieldVector<Scalar, 1> >;
                std::unordered_map<std::size_t, ShapeValues> circleShapeValues;

                for (int k = 0; k < circlePoints.size(); ++k)
                {
                    const auto circleBulkElementIndices = intersectingEntities(circlePoints[k], bulkTree);
                    if (circleBulkElementIndices.empty())
                        continue;

                    const auto bulkElementIdx = circleBulkElementIndices[0];
                    circleStencil[k] = bulkElementIdx;
                    circleIpWeight[k] = weight;

                    if (isBox<bulkIdx>())
                    {
                        if (!static_cast<bool>(circleCornerIndices.count(bulkElementIdx)))
                        {
                            const auto bulkElement = this->problem(bulkIdx).fvGridGeometry().element(bulkElementIdx);
                            circleCornerIndices[bulkElementIdx] = this->vertexIndices(bulkIdx, bulkElementIdx);

                            // evaluate shape functions at the integration point
                            const auto bulkGeometry = bulkElement.geometry();
                            this->getShapeValues(bulkIdx, this->problem(bulkIdx).fvGridGeometry(), bulkGeometry, circlePoints[k], circleShapeValues[bulkElementIdx]);
                        }
                    }
                }

                // export low dim circle stencil
                if (isBox<bulkIdx>())
                {
                    // we insert all vertices and make it unique later
                    for (const auto& vertices : circleCornerIndices)
                    {
                        this->couplingStencils(lowDimIdx)[lowDimElementIdx].insert(this->couplingStencils(lowDimIdx)[lowDimElementIdx].end(),
                                                                                   vertices.second.begin(), vertices.second.end());

                    }
                }
                else
                {
                    this->couplingStencils(lowDimIdx)[lowDimElementIdx].insert(this->couplingStencils(lowDimIdx)[lowDimElementIdx].end(),
                                                                               circleStencil.begin(), circleStencil.end());
                }

                // loop over the bulk elements at the integration points (usually one except when it is on a face or edge or vertex)
                for (auto bulkElementIdx : bulkElementIndices)
                {
                    const auto id = this->idCounter_++;
                    const auto ie = lowDimGeometry.integrationElement(qp.position());
                    const auto qpweight = qp.weight();

                    this->pointSources(bulkIdx).emplace_back(globalPos, id, qpweight, ie, std::vector<std::size_t>({bulkElementIdx}));
                    this->pointSources(bulkIdx).back().setEmbeddings(bulkElementIndices.size());
                    this->pointSources(lowDimIdx).emplace_back(globalPos, id, qpweight, ie, std::vector<std::size_t>({lowDimElementIdx}));
                    this->pointSources(lowDimIdx).back().setEmbeddings(bulkElementIndices.size());

                    // pre compute additional data used for the evaluation of
                    // the actual solution dependent source term
                    PointSourceData psData;

                    if (isBox<lowDimIdx>())
                    {
                        ShapeValues shapeValues;
                        this->getShapeValues(lowDimIdx, this->problem(lowDimIdx).fvGridGeometry(), lowDimGeometry, globalPos, shapeValues);
                        psData.addLowDimInterpolation(shapeValues, this->vertexIndices(lowDimIdx, lowDimElementIdx), lowDimElementIdx);
                    }
                    else
                    {
                        psData.addLowDimInterpolation(lowDimElementIdx);
                    }

                    // add data needed to compute integral over the circle
                    if (isBox<bulkIdx>())
                    {
                        psData.addCircleInterpolation(circleCornerIndices, circleShapeValues, circleIpWeight, circleStencil);

                        const auto bulkGeometry = this->problem(bulkIdx).fvGridGeometry().element(bulkElementIdx).geometry();
                        ShapeValues shapeValues;
                        this->getShapeValues(bulkIdx, this->problem(bulkIdx).fvGridGeometry(), bulkGeometry, globalPos, shapeValues);
                        psData.addBulkInterpolation(shapeValues, this->vertexIndices(bulkIdx, bulkElementIdx), bulkElementIdx);
                    }
                    else
                    {
                        psData.addCircleInterpolation(circleIpWeight, circleStencil);
                        psData.addBulkInterpolation(bulkElementIdx);
                    }

                    // publish point source data in the global vector
                    this->pointSourceData().emplace_back(std::move(psData));

                    // export the bulk coupling stencil
                    if (isBox<lowDimIdx>())
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
                    if (isBox<bulkIdx>())
                    {
                        // we insert all vertices and make it unique later
                        for (const auto& vertices : circleCornerIndices)
                        {
                            extendedSourceStencil_.stencil()[bulkElementIdx].insert(extendedSourceStencil_.stencil()[bulkElementIdx].end(),
                                                                       vertices.second.begin(), vertices.second.end());

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

        std::cout << "took " << watch.elapsed() << " seconds." << std::endl;
    }

    //! Compute the low dim volume fraction in the bulk domain cells
    void computeLowDimVolumeFractions()
    {
        // resize the storage vector
        lowDimVolumeInBulkElement_.resize(this->gridView(bulkIdx).size(0));

        // compute the low dim volume fractions
        for (const auto& is : intersections(this->glue()))
        {
            // all inside elements are identical...
            const auto& inside = is.inside(0);
            const auto intersectionGeometry = is.geometry();
            const std::size_t lowDimElementIdx = this->problem(lowDimIdx).fvGridGeometry().elementMapper().index(inside);

            // compute the volume the low-dim domain occupies in the bulk domain if it were full-dimensional
            const auto radius = this->problem(lowDimIdx).spatialParams().radius(lowDimElementIdx);
            for (std::size_t outsideIdx = 0; outsideIdx < is.neighbor(0); ++outsideIdx)
            {
                const auto& outside = is.outside(outsideIdx);
                const std::size_t bulkElementIdx = this->problem(bulkIdx).fvGridGeometry().elementMapper().index(outside);
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
        const auto eIdx = this->problem(bulkIdx).fvGridGeometry().elementMapper().index(element);
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


/*!
 * \ingroup EmbeddedCoupling
 * \brief Manages the coupling between bulk elements and lower dimensional elements
 *        Point sources on each integration point are computed by an AABB tree.
 * \note Specialization for coupling method using cylinder sources with 3d quantities evaluated on the cylinder surface
 * \note the cylinder source is approximated by point sources on the cylinder surface
 */
template<class MDTraits>
class EmbeddedCouplingManager1d3d<MDTraits, EmbeddedCouplingMode::cylindersources>
: public EmbeddedCouplingManagerBase<MDTraits, EmbeddedCouplingManager1d3d<MDTraits, EmbeddedCouplingMode::cylindersources>>
{
    using ThisType = EmbeddedCouplingManager1d3d<MDTraits, EmbeddedCouplingMode::cylindersources>;
    using ParentType = EmbeddedCouplingManagerBase<MDTraits, ThisType>;
    using Scalar = typename MDTraits::Scalar;
    using SolutionVector = typename MDTraits::SolutionVector;
    using PointSourceData = typename ParentType::PointSourceTraits::PointSourceData;

    static constexpr auto bulkIdx = typename MDTraits::template SubDomain<0>::Index();
    static constexpr auto lowDimIdx = typename MDTraits::template SubDomain<1>::Index();

    // the sub domain type aliases
    template<std::size_t id> using SubDomainTypeTag = typename MDTraits::template SubDomain<id>::TypeTag;
    template<std::size_t id> using Problem = GetPropType<SubDomainTypeTag<id>, Properties::Problem>;
    template<std::size_t id> using FVGridGeometry = GetPropType<SubDomainTypeTag<id>, Properties::FVGridGeometry>;
    template<std::size_t id> using GridView = typename FVGridGeometry<id>::GridView;
    template<std::size_t id> using Element = typename GridView<id>::template Codim<0>::Entity;

    using GlobalPosition = typename Element<bulkIdx>::Geometry::GlobalCoordinate;

    template<std::size_t id>
    static constexpr bool isBox()
    { return FVGridGeometry<id>::discMethod == DiscretizationMethod::box; }

    enum {
        bulkDim = GridView<bulkIdx>::dimension,
        lowDimDim = GridView<lowDimIdx>::dimension,
        dimWorld = GridView<bulkIdx>::dimensionworld
    };
public:
    static constexpr EmbeddedCouplingMode couplingMode = EmbeddedCouplingMode::cylindersources;

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
        const auto& bulkTree = this->problem(bulkIdx).fvGridGeometry().boundingBoxTree();

        // initilize the maps
        // do some logging and profiling
        Dune::Timer watch;
        std::cout << "Initializing the point sources..." << std::endl;

        // clear all internal members like pointsource vectors and stencils
        // initializes the point source id counter
        this->clear();

        // precompute the vertex indices for efficiency
        this->preComputeVertexIndices(bulkIdx);
        this->preComputeVertexIndices(lowDimIdx);

        // iterate over all lowdim elements
        const auto& lowDimProblem = this->problem(lowDimIdx);
        for (const auto& lowDimElement : elements(this->gridView(lowDimIdx)))
        {
            // get the Gaussian quadrature rule for the low dim element
            const auto lowDimGeometry = lowDimElement.geometry();
            const auto& quad = Dune::QuadratureRules<Scalar, lowDimDim>::rule(lowDimGeometry.type(), order);

            const auto lowDimElementIdx = lowDimProblem.fvGridGeometry().elementMapper().index(lowDimElement);

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

                static const auto numIp = getParam<int>("MixedDimension.NumCircleSegments", 25);
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

                        if (isBox<lowDimIdx>())
                        {
                            using ShapeValues = std::vector<Dune::FieldVector<Scalar, 1> >;
                            ShapeValues shapeValues;
                            this->getShapeValues(lowDimIdx, this->problem(lowDimIdx).fvGridGeometry(), lowDimGeometry, globalPos, shapeValues);
                            psData.addLowDimInterpolation(shapeValues, this->vertexIndices(lowDimIdx, lowDimElementIdx), lowDimElementIdx);
                        }
                        else
                        {
                            psData.addLowDimInterpolation(lowDimElementIdx);
                        }

                        // add data needed to compute integral over the circle
                        if (isBox<bulkIdx>())
                        {
                            using ShapeValues = std::vector<Dune::FieldVector<Scalar, 1> >;
                            const auto bulkGeometry = this->problem(bulkIdx).fvGridGeometry().element(bulkElementIdx).geometry();
                            ShapeValues shapeValues;
                            this->getShapeValues(bulkIdx, this->problem(bulkIdx).fvGridGeometry(), bulkGeometry, circlePos, shapeValues);
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

        // compute the low dim volume fractions
        for (const auto& is : intersections(this->glue()))
        {
            // all inside elements are identical...
            const auto& inside = is.inside(0);
            const auto intersectionGeometry = is.geometry();
            const std::size_t lowDimElementIdx = this->problem(lowDimIdx).fvGridGeometry().elementMapper().index(inside);

            // compute the volume the low-dim domain occupies in the bulk domain if it were full-dimensional
            const auto radius = this->problem(lowDimIdx).spatialParams().radius(lowDimElementIdx);
            for (std::size_t outsideIdx = 0; outsideIdx < is.neighbor(0); ++outsideIdx)
            {
                const auto& outside = is.outside(outsideIdx);
                const std::size_t bulkElementIdx = this->problem(bulkIdx).fvGridGeometry().elementMapper().index(outside);
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
        const auto eIdx = this->problem(bulkIdx).fvGridGeometry().elementMapper().index(element);
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


/*!
 * \ingroup EmbeddedCoupling
 * \brief Manages the coupling between bulk elements and lower dimensional elements
 *        Point sources on each integration point are computed by an AABB tree.
 * \note Specialization for coupling method using a distributed kernel source with 3d quantities evaluated on the line
 * \note the kernel per point source is isotropic and it's integral over the domain is one
 */
template<class MDTraits>
class EmbeddedCouplingManager1d3d<MDTraits, EmbeddedCouplingMode::kernel>
: public EmbeddedCouplingManagerBase<MDTraits, EmbeddedCouplingManager1d3d<MDTraits, EmbeddedCouplingMode::kernel>>
{
    using ThisType = EmbeddedCouplingManager1d3d<MDTraits, EmbeddedCouplingMode::kernel>;
    using ParentType = EmbeddedCouplingManagerBase<MDTraits, ThisType>;
    using Scalar = typename MDTraits::Scalar;
    using SolutionVector = typename MDTraits::SolutionVector;
    using PointSourceData = typename ParentType::PointSourceTraits::PointSourceData;

    static constexpr auto bulkIdx = typename MDTraits::template SubDomain<0>::Index();
    static constexpr auto lowDimIdx = typename MDTraits::template SubDomain<1>::Index();

    // the sub domain type aliases
    template<std::size_t id> using SubDomainTypeTag = typename MDTraits::template SubDomain<id>::TypeTag;
    template<std::size_t id> using Problem = GetPropType<SubDomainTypeTag<id>, Properties::Problem>;
    template<std::size_t id> using FVGridGeometry = GetPropType<SubDomainTypeTag<id>, Properties::FVGridGeometry>;
    template<std::size_t id> using GridView = typename FVGridGeometry<id>::GridView;
    template<std::size_t id> using Element = typename GridView<id>::template Codim<0>::Entity;

    using GlobalPosition = typename Element<bulkIdx>::Geometry::GlobalCoordinate;

    template<std::size_t id>
    static constexpr bool isBox()
    { return FVGridGeometry<id>::discMethod == DiscretizationMethod::box; }

    enum {
        bulkDim = GridView<bulkIdx>::dimension,
        lowDimDim = GridView<lowDimIdx>::dimension,
        dimWorld = GridView<bulkIdx>::dimensionworld
    };
public:
    static constexpr EmbeddedCouplingMode couplingMode = EmbeddedCouplingMode::kernel;

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
        // initilize the maps
        // do some logging and profiling
        Dune::Timer watch;
        std::cout << "Initializing the point sources..." << std::endl;

        // clear all internal members like pointsource vectors and stencils
        // initializes the point source id counter
        this->clear();
        bulkSourceIds_.clear();
        bulkSourceWeights_.clear();
        extendedSourceStencil_.stencil().clear();

        // precompute the vertex indices for efficiency for the box method
        this->preComputeVertexIndices(bulkIdx);
        this->preComputeVertexIndices(lowDimIdx);

        const auto& bulkFvGridGeometry = this->problem(bulkIdx).fvGridGeometry();
        const auto& lowDimFvGridGeometry = this->problem(lowDimIdx).fvGridGeometry();

        bulkSourceIds_.resize(this->gridView(bulkIdx).size(0));
        bulkSourceWeights_.resize(this->gridView(bulkIdx).size(0));

        // intersect the bounding box trees
        this->glueGrids();

        // reserve memory for data
        this->pointSourceData().reserve(this->glue().size());
        this->averageDistanceToBulkCell().reserve(this->glue().size());
        fluxScalingFactor_.reserve(this->glue().size());

        // generate a bunch of random vectors and values for
        // Monte-carlo integration on the cylinder defined by line and radius
        static const bool useCylinderIntegration = getParam<bool>("MixedDimension.UseCylinderIntegration");
        const int samples = useCylinderIntegration ? getParam<int>("MixedDimension.NumMCSamples") : 1;
        CylinderIntegration<Scalar, CylinderIntegrationMethod::spaced> cylIntegration(samples);
        static const auto kernelWidthFactor = getParam<Scalar>("MixedDimension.KernelWidthFactor");

        for (const auto& is : intersections(this->glue()))
        {
            // all inside elements are identical...
            const auto& inside = is.inside(0);
            // get the intersection geometry
            const auto intersectionGeometry = is.geometry();
            // get the Gaussian quadrature rule for the local intersection
            const std::size_t lowDimElementIdx = lowDimFvGridGeometry.elementMapper().index(inside);

            // for each intersection integrate kernel and add:
            //  * 1d: a new point source
            //  * 3d: a new kernel volume source
            const auto radius = this->problem(lowDimIdx).spatialParams().radius(lowDimElementIdx);
            const auto kernelWidth = kernelWidthFactor*radius;

            const auto a = intersectionGeometry.corner(0);
            const auto b = intersectionGeometry.corner(1);

            // we can have multiple 3d elements per 1d intersection, for evaluation taking one of them should be fine
            for (std::size_t outsideIdx = 0; outsideIdx < is.neighbor(0); ++outsideIdx)
            {
                // compute source id
                // for each id we have
                // (1) a point source for the 1d domain
                // (2) a bulk source weight for the each element in the 3d domain (the fraction of the total source/sink for the 1d element with the point source)
                //     TODO: i.e. one integration of the kernel should be enough (for each of the weights it's weight/embeddings)
                // (3) the flux scaling factor for each outside element, i.e. each id
                const auto id = this->idCounter_++;

                // place a point source at the intersection center
                const auto center = intersectionGeometry.center();
                this->pointSources(lowDimIdx).emplace_back(center, id, /*weight=*/1.0, intersectionGeometry.volume(), std::vector<std::size_t>({lowDimElementIdx}));
                this->pointSources(lowDimIdx).back().setEmbeddings(is.neighbor(0));

                const auto& outside = is.outside(outsideIdx);
                const std::size_t bulkElementIdx = bulkFvGridGeometry.elementMapper().index(outside);

                // compute the weights for the bulk volume sources
                if (useCylinderIntegration)
                    cylIntegration.setGeometry(a, b, kernelWidth);
                computeBulkSource(intersectionGeometry, radius, kernelWidth, id, lowDimElementIdx, bulkElementIdx, cylIntegration, is.neighbor(0));

                // pre compute additional data used for the evaluation of
                // the actual solution dependent source term
                PointSourceData psData;

                // lowdim interpolation (evaluate at center)
                if (isBox<lowDimIdx>())
                {
                    using ShapeValues = std::vector<Dune::FieldVector<Scalar, 1> >;
                    const auto lowDimGeometry = this->problem(lowDimIdx).fvGridGeometry().element(lowDimElementIdx).geometry();
                    ShapeValues shapeValues;
                    this->getShapeValues(lowDimIdx, this->problem(lowDimIdx).fvGridGeometry(), lowDimGeometry, center, shapeValues);
                    psData.addLowDimInterpolation(shapeValues, this->vertexIndices(lowDimIdx, lowDimElementIdx), lowDimElementIdx);
                }
                else
                {
                    psData.addLowDimInterpolation(lowDimElementIdx);
                }

                // bulk interpolation (evaluate at center)
                if (isBox<bulkIdx>())
                {
                    using ShapeValues = std::vector<Dune::FieldVector<Scalar, 1> >;
                    const auto bulkGeometry = this->problem(bulkIdx).fvGridGeometry().element(bulkElementIdx).geometry();
                    ShapeValues shapeValues;
                    this->getShapeValues(bulkIdx, this->problem(bulkIdx).fvGridGeometry(), bulkGeometry, center, shapeValues);
                    psData.addBulkInterpolation(shapeValues, this->vertexIndices(bulkIdx, bulkElementIdx), bulkElementIdx);
                }
                else
                {
                    psData.addBulkInterpolation(bulkElementIdx);
                }

                // publish point source data in the global vector
                this->pointSourceData().emplace_back(std::move(psData));

                // compute the average flux scaling factor
                const auto outsideGeometry = outside.geometry();
                static const int fluxScalingFactorIntOrder = getParam<int>("MixedDimension.FluxScalingIntegrationOrder", 5);
                const auto& quadBulk = Dune::QuadratureRules<Scalar, bulkDim>::rule(outsideGeometry.type(), fluxScalingFactorIntOrder);
                // iterate over all quadrature points
                Scalar avgMinDist = 0.0;
                for (auto&& qpBulk : quadBulk)
                {
                    const auto globalPosBulk = outsideGeometry.global(qpBulk.position());
                    const auto d = computeMinDistance_(a, b, globalPosBulk);
                    avgMinDist += d*qpBulk.weight();
                }

                // add it to the internal data vector
                this->averageDistanceToBulkCell().push_back(avgMinDist);
                fluxScalingFactor_.push_back(this->problem(bulkIdx).fluxScalingFactor(avgMinDist, radius));

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

        // compute the low dim volume fractions
        for (const auto& is : intersections(this->glue()))
        {
            // all inside elements are identical...
            const auto& inside = is.inside(0);
            const auto intersectionGeometry = is.geometry();
            const std::size_t lowDimElementIdx = this->problem(lowDimIdx).fvGridGeometry().elementMapper().index(inside);

            // compute the volume the low-dim domain occupies in the bulk domain if it were full-dimensional
            const auto radius = this->problem(lowDimIdx).spatialParams().radius(lowDimElementIdx);
            for (std::size_t outsideIdx = 0; outsideIdx < is.neighbor(0); ++outsideIdx)
            {
                const auto& outside = is.outside(outsideIdx);
                const std::size_t bulkElementIdx = this->problem(bulkIdx).fvGridGeometry().elementMapper().index(outside);
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
        const auto eIdx = this->problem(bulkIdx).fvGridGeometry().elementMapper().index(element);
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
    const std::vector<std::size_t> bulkSourceIds(std::size_t eIdx, std::size_t scvIdx = 0) const
    { return bulkSourceIds_[eIdx][scvIdx]; }

    //! return all source ids for a bulk elements
    const std::vector<Scalar> bulkSourceWeights(std::size_t eIdx, std::size_t scvIdx = 0) const
    { return bulkSourceWeights_[eIdx][scvIdx]; }

    //! The flux scaling factor for a source with id
    Scalar fluxScalingFactor(std::size_t id) const
    { return fluxScalingFactor_[id]; }

    // \}

private:
    //! compute the kernel distributed sources and add stencils
    template<class Line, class CylIntegration>
    void computeBulkSource(const Line& line, const Scalar radius, const Scalar kernelWidth,
                           std::size_t id, std::size_t lowDimElementIdx, std::size_t coupledBulkElementIdx,
                           const CylIntegration& cylIntegration, std::size_t embeddings)
    {
        // Monte-carlo integration on the cylinder defined by line and radius
        const auto numBulkElements = this->gridView(bulkIdx).size(0);
        constexpr std::size_t numScv = isBox<bulkIdx>() ? 1<<bulkDim : 1;
        std::vector<std::array<Scalar, numScv>> weights(numBulkElements, std::array<Scalar, numScv>{});
        std::vector<std::bitset<numScv>> visited(numBulkElements, std::bitset<numScv>{});

        Scalar integral = 0.0;
        static const bool useCylinderIntegration = getParam<bool>("MixedDimension.UseCylinderIntegration");
        if (useCylinderIntegration)
        {
            for (int i = 0; i < cylIntegration.size(); ++i)
            {
                const auto point = cylIntegration.getIntegrationPoint(i);
                const auto weight = evalKernel(line.corner(0), line.corner(1), point, radius, kernelWidth)*cylIntegration.integrationElement(i);
                integral += weight;

                static const auto min = getParam<GlobalPosition>("Tissue.Grid.LowerLeft");
                static const auto max = getParam<GlobalPosition>("Tissue.Grid.UpperRight");
                static const auto cells = getParam<GlobalPosition>("Tissue.Grid.Cells");
                const auto is = intersectingEntityCartesian(point, min, max, cells);
                // const auto bulkIndices = intersectingEntities(point, this->problem(bulkIdx).fvGridGeometry().boundingBoxTree(), true);
                if (is.first)
                {
                    const auto bulkElementIdx = is.second;
                    if (isBox<bulkIdx>())
                    {
                        const auto dist = point-min;
                        const auto diag = max-min;
                        std::bitset<3> corner;
                        for (int i = 0; i < 3; ++i)
                            corner[i] = bool(int(std::round(dist[i]/diag[i])));

                        std::size_t scvIdx = corner.to_ulong();
                        visited[bulkElementIdx][scvIdx] = true;
                        weights[bulkElementIdx][scvIdx] += weight;
                    }
                    else
                    {
                        visited[bulkElementIdx][0] = true;
                        weights[bulkElementIdx][0] += weight;
                    }
                }
            }
        }
        // integrate over 3d domain using standard quadrature rules
        else
        {
            static const int maxLevel = getParam<int>("MixedDimension.MaxOctreeLevel");
            CylinderDomain<Scalar> cylinder(line.corner(0), line.corner(1), kernelWidth);
            const auto kernelFunc = [&](const auto& p){ return evalKernel(line.corner(0), line.corner(1), p, radius, kernelWidth); };
            for (const auto& element : elements(this->gridView(bulkIdx)))
            {
                const auto bulkElementIdx = this->problem(bulkIdx).fvGridGeometry().elementMapper().index(element);
                const auto geometry = element.geometry();
                OctreeIntegrator<Scalar> integrator(geometry.corner(0), geometry.corner(7), maxLevel);
                const auto weight = integrator.integrate(kernelFunc, cylinder);
                if (weight > 1e-35)
                {
                    visited[bulkElementIdx][0] = true;
                    weights[bulkElementIdx][0] += weight;
                    integral += weight;
                }
            }
        }

        // const auto length = line.volume();
        // std::cout << "integration error in %: " << std::abs(length-integral)/length*100 << std::endl;
        // const auto correctionFactor = length/integral;
        for (const auto& element : elements(this->problem(bulkIdx).fvGridGeometry().gridView()))
        {
            const auto bulkElementIdx = this->problem(bulkIdx).fvGridGeometry().elementMapper().index(element);
            if (!isBox<bulkIdx>() && !visited[bulkElementIdx][0])
                continue;

            if (isBox<bulkIdx>())
            {
                auto fvGeometry = localView(this->problem(bulkIdx).fvGridGeometry());
                fvGeometry.bindElement(element);
                for (const auto& scv : scvs(fvGeometry))
                {
                    const auto scvIdx = scv.indexInElement();
                    if (!visited[bulkElementIdx][scvIdx])
                        continue;

                    bulkSourceIds_[bulkElementIdx][scvIdx].push_back(id);
                    bulkSourceWeights_[bulkElementIdx][scvIdx].push_back(weights[bulkElementIdx][scvIdx]/embeddings);//*correctionFactor);
                    addBulkSourceStencils_(bulkElementIdx, lowDimElementIdx, coupledBulkElementIdx);
                }

            }
            else
            {
                bulkSourceIds_[bulkElementIdx][0].push_back(id);
                bulkSourceWeights_[bulkElementIdx][0].push_back(weights[bulkElementIdx][0]/embeddings);//*correctionFactor);
                addBulkSourceStencils_(bulkElementIdx, lowDimElementIdx, coupledBulkElementIdx);
            }
        }

    }

    //! add additional stencil entries for the bulk element
    void addBulkSourceStencils_(std::size_t bulkElementIdx, std::size_t coupledLowDimElementIdx, std::size_t coupledBulkElementIdx)
    {
        // add the lowdim element to the coupling stencil of this bulk element
        if (isBox<lowDimIdx>())
        {
            const auto& vertices = this->vertexIndices(lowDimIdx, coupledLowDimElementIdx);
            this->couplingStencils(bulkIdx)[bulkElementIdx].insert(this->couplingStencils(bulkIdx)[bulkElementIdx].end(),
                                                                   vertices.begin(), vertices.end());

        }
        else
        {
            this->couplingStencils(bulkIdx)[bulkElementIdx].push_back(coupledLowDimElementIdx);
        }

        // the extended source stencil, every 3d element with a source is coupled to
        // the element/dofs where the 3d quantities are measured
        if (isBox<bulkIdx>())
        {
            const auto& vertices = this->vertexIndices(bulkIdx, coupledBulkElementIdx);
            extendedSourceStencil_.stencil()[bulkElementIdx].insert(extendedSourceStencil_.stencil()[bulkElementIdx].end(),
                                                                    vertices.begin(), vertices.end());
        }
        else
        {
            extendedSourceStencil_.stencil()[bulkElementIdx].push_back(coupledBulkElementIdx);
        }
    }

    //! a cylindrical kernel around the segment a->b
    inline Scalar evalKernel(const GlobalPosition& a,
                             const GlobalPosition& b,
                             const GlobalPosition& point,
                             const Scalar R,
                             const Scalar rho) const
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

        static const std::string type = getParam<std::string>("MixedDimension.KernelType");
        if (type == "Const")
            return 1.0/(M_PI*rho*rho);
        else if (type == "Cubic")
        {
            const auto r2 = r*r; const auto r3 = r*r*r;
            const auto rho2 = rho*rho; const auto rho3 = rho2*rho;
            return 10.0*(2.0*r3/rho3 - 3.0*r2/rho2 + 1.0)/(3.0*M_PI*rho2);
        }

        DUNE_THROW(Dune::NotImplemented, "Kernel type: " << type);
    }

    inline Scalar computeMinDistance_(const GlobalPosition& a, const GlobalPosition& b, const GlobalPosition& x3d) const
    {
        // r is the minimum distance to the line through a and b
        const auto ab = b - a;
        const auto t = (x3d - a)*ab/ab.two_norm2();
        auto proj = a; proj.axpy(t, ab);
        const auto r = (proj - x3d).two_norm();
        return r;
    }

    //! the extended source stencil object for kernel coupling
    EmbeddedCoupling::ExtendedSourceStencil<ThisType> extendedSourceStencil_;
    //! vector for the volume fraction of the lowdim domain in the bulk domain cells
    std::vector<Scalar> lowDimVolumeInBulkElement_;
    //! kernel sources to integrate for each bulk element
    std::vector<std::array<std::vector<size_t>, isBox<bulkIdx>() ? 1<<bulkDim : 1>> bulkSourceIds_;
    //! the integral of the kernel for each point source / integration point, i.e. weight for the source
    std::vector<std::array<std::vector<Scalar>, isBox<bulkIdx>() ? 1<<bulkDim : 1>> bulkSourceWeights_;
    //! the flux scaling factor for the respective source id
    std::vector<Scalar> fluxScalingFactor_;
};

} // end namespace Dumux

#endif
