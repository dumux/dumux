// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   See the file COPYING for full copying permissions.                      *
 *                                                                           *
 *   This program is free software: you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation, either version 2 of the License, or       *
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

#include <vector>

#include <dune/common/timer.hh>
#include <dune/geometry/quadraturerules.hh>

#include <dumux/common/properties.hh>
#include <dumux/multidomain/embedded/pointsourcedata.hh>
#include <dumux/multidomain/embedded/integrationpointsource.hh>
#include <dumux/multidomain/embedded/couplingmanagerbase.hh>

namespace Dumux {

/*!
 * \ingroup EmbeddedCoupling
 * \brief The coupling mode
 */
enum class EmbeddedCouplingMode
{ simple, average, kernel };

//! point source traits for the circle average coupling mode
template<class MDTraits>
struct CircleAveragePointSourceTraits
{
private:
    template<std::size_t i> using SubDomainTypeTag = typename MDTraits::template SubDomainTypeTag<i>;
    template<std::size_t i> using FVGridGeometry = typename GET_PROP_TYPE(SubDomainTypeTag<i>, FVGridGeometry);
    template<std::size_t i> using NumEqVector = typename GET_PROP_TYPE(SubDomainTypeTag<i>, NumEqVector);
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

    static constexpr auto bulkIdx = typename MDTraits::template DomainIdx<0>();
    static constexpr auto lowDimIdx = typename MDTraits::template DomainIdx<1>();

    // the sub domain type aliases
    template<std::size_t id> using SubDomainTypeTag = typename MDTraits::template SubDomainTypeTag<id>;
    template<std::size_t id> using Problem = typename GET_PROP_TYPE(SubDomainTypeTag<id>, Problem);
    template<std::size_t id> using FVGridGeometry = typename GET_PROP_TYPE(SubDomainTypeTag<id>, FVGridGeometry);
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
    using ParentType::ParentType;

    void init(std::shared_ptr<Problem<bulkIdx>> bulkProblem,
              std::shared_ptr<Problem<lowDimIdx>> lowDimProblem,
              const SolutionVector& curSol)
    {
        ParentType::init(bulkProblem, lowDimProblem, curSol);
        computeLowDimVolumeFractions();
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
        constexpr auto numEq = GET_PROP_TYPE(SubDomainTypeTag<i>, ModelTraits)::numEq();
        const auto& elementI = localAssemblerI.element();

        // only do something if we have an extended stencil
        if (extendedSourceStencil_(domainI, elementI).empty())
            return;

        // compute the undeflected residual (source only!)
        const auto origResidual = localAssemblerI.evalLocalSourceResidual(elementI);

        // compute derivate for all additional dofs in the circle stencil
        for (const auto dofIndex : extendedSourceStencil_(domainI, elementI))
        {
            auto partialDerivs = origResidual;
            const auto origPriVars = this->curSol()[domainI][dofIndex];

            // calculate derivatives w.r.t to the privars at the dof at hand
            for (int pvIdx = 0; pvIdx < numEq; pvIdx++)
            {
                // reset partial derivatives
                partialDerivs = 0.0;

                auto evalResiduals = [&](Scalar priVar)
                {
                    // update the coupling context (solution vector and recompute element residual)
                    auto priVars = origPriVars;
                    priVars[pvIdx] = priVar;
                    this->updateCouplingContext(domainI, localAssemblerI, domainI, dofIndex, priVars, pvIdx);
                    return localAssemblerI.evalLocalSourceResidual(elementI);
                };

                // derive the residuals numerically
                static const int numDiffMethod = getParam<int>("Assembly.NumericDifferenceMethod");
                NumericDifferentiation::partialDerivative(evalResiduals, this->curSol()[domainI][dofIndex][pvIdx],
                                                          partialDerivs, origResidual, numDiffMethod);

                // update the global stiffness matrix with the current partial derivatives
                for (auto&& scvJ : scvs(localAssemblerI.fvGeometry()))
                {
                    for (int eqIdx = 0; eqIdx < numEq; eqIdx++)
                    {
                        // A[i][col][eqIdx][pvIdx] is the rate of change of
                        // the residual of equation 'eqIdx' at dof 'i'
                        // depending on the primary variable 'pvIdx' at dof
                        // 'col'.
                        A[scvJ.dofIndex()][dofIndex][eqIdx][pvIdx] += partialDerivs[scvJ.indexInElement()][eqIdx];
                    }
                }

                // restore the original coupling context
                this->updateCouplingContext(domainI, localAssemblerI, domainI, dofIndex, origPriVars, pvIdx);
            }
        }
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
        // add additional dof dependencies
        for (const auto& element : elements(this->gridView(domainI)))
        {
            const auto& dofs = this->extendedSourceStencil_(domainI, element);

            if (isBox<domainI>())
            {
                for (int i = 0; i < element.subEntities(GridView<domainI>::dimension); ++i)
                    for (const auto globalJ : dofs)
                        pattern.add(this->problem(domainI).fvGridGeometry().vertexMapper().subIndex(element, i, GridView<domainI>::dimension), globalJ);
            }
            else
            {
                const auto globalI = this->problem(domainI).fvGridGeometry().elementMapper().index(element);
                for (const auto globalJ : dofs)
                    pattern.add(globalI, globalJ);
            }
        }
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
        // precompute the bulk vertex indices for efficiency
        this->preComputeBulkVertexIndices();

        // Initialize the bulk bounding box tree
        const auto& bulkTree = this->problem(bulkIdx).fvGridGeometry().boundingBoxTree();

        // initilize the maps
        // do some logging and profiling
        Dune::Timer watch;
        std::cout << "Initializing the point sources..." << std::endl;

        // clear all internal members like pointsource vectors and stencils
        // initializes the point source id counter
        this->clear();
        bulkCircleStencils_.clear();

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
                const auto length = 2*M_PI*radius/numIp;

                const auto circlePoints = getCirclePoints_(globalPos, normal, radius, numIp);
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
                    circleIpWeight[k] = length;

                    if (isBox<bulkIdx>())
                    {
                        if (!static_cast<bool>(circleCornerIndices.count(bulkElementIdx)))
                        {
                            const auto bulkElement = this->problem(bulkIdx).fvGridGeometry().element(bulkElementIdx);
                            circleCornerIndices[bulkElementIdx] = this->bulkVertexIndices(bulkElementIdx);

                            // evaluate shape functions at the integration point
                            const auto bulkGeometry = bulkElement.geometry();
                            this->getShapeValues_(this->problem(bulkIdx).fvGridGeometry(), bulkGeometry, circlePoints[k], circleShapeValues[bulkElementIdx]);
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
                    psData.addLowDimInterpolation(lowDimElementIdx);
                    psData.addBulkInterpolation(bulkElementIdx);
                    // add data needed to compute integral over the circle
                    if (isBox<bulkIdx>())
                        psData.addCircleInterpolation(circleCornerIndices, circleShapeValues, circleIpWeight, circleStencil);
                    else
                        psData.addCircleInterpolation(circleIpWeight, circleStencil);

                    // publish point source data in the global vector
                    this->pointSourceData().push_back(psData);

                    // export the bulk coupling stencil
                    if (isBox<lowDimIdx>())
                    {
                        this->couplingStencils(bulkIdx)[bulkElementIdx].insert(this->couplingStencils(bulkIdx)[bulkElementIdx].end(),
                                                                               this->bulkVertexIndices(bulkElementIdx).begin(),
                                                                               this->bulkVertexIndices(bulkElementIdx).end());

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
                            bulkCircleStencils_[bulkElementIdx].insert(bulkCircleStencils_[bulkElementIdx].end(),
                                                                       vertices.second.begin(), vertices.second.end());

                        }
                    }
                    else
                    {
                        bulkCircleStencils_[bulkElementIdx].insert(bulkCircleStencils_[bulkElementIdx].end(),
                                                                   circleStencil.begin(), circleStencil.end());
                    }
                }
            }
        }

        // make the circle stencil unique (for source derivatives)
        for (auto&& stencil : bulkCircleStencils_)
        {
            std::sort(stencil.second.begin(), stencil.second.end());
            stencil.second.erase(std::unique(stencil.second.begin(), stencil.second.end()), stencil.second.end());

            // remove the vertices element (box)
            if (isBox<bulkIdx>())
            {
                const auto& indices = this->bulkVertexIndices(stencil.first);
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
    //! the extended source stencil for the bulk domain due to the source average
    const std::vector<std::size_t>& extendedSourceStencil_(Dune::index_constant<0> bulkDomain, const Element<0>& bulkElement) const
    {
        const auto bulkElementIdx = this->problem(bulkIdx).fvGridGeometry().elementMapper().index(bulkElement);
        if (bulkCircleStencils_.count(bulkElementIdx))
            return bulkCircleStencils_.at(bulkElementIdx);
        else
            return this->emptyStencil();
    }

    //! the extended source stencil for the low dim domain is empty
    const std::vector<std::size_t>& extendedSourceStencil_(Dune::index_constant<1> bulkDomain, const Element<1>& lowDimElement) const
    { return this->emptyStencil(); }

    //! returns a vector of points on a circle given the circle center
    //! a normal vector for orientaiton of the circle plane and the radius
    std::vector<GlobalPosition> getCirclePoints_(const GlobalPosition& center,
                                                 const GlobalPosition& normal,
                                                 const Scalar radius,
                                                 const std::size_t numIp = 20)
    {
        const Scalar eps = 1.5e-7;
        static_assert(dimWorld == 3, "Only implemented for world dimension 3");

        std::vector<GlobalPosition> points(numIp);

        // make sure n is a unit vector
        auto n = normal;
        n /= n.two_norm();

        // caculate a vector u perpendicular to n
        GlobalPosition u;
        if (std::abs(n[0]) < eps && std::abs(n[1]) < eps)
            if (std::abs(n[2]) < eps)
                DUNE_THROW(Dune::MathError, "The normal vector has to be non-zero!");
            else
                u = {0, 1, 0};
        else
            u = {-n[1], n[0], 0};

        u *= radius/u.two_norm();

        // the circle parameterization is p(t) = r*cos(t)*u + r*sin(t)*(n x u) + c
        auto tangent = crossProduct(u, n);
        tangent *= radius/tangent.two_norm();

        // the parameter with an offset
        Scalar t = 0 + 0.1;
        // insert the vertices
        for (std::size_t i = 0; i < numIp; ++i)
        {
            points[i] = GlobalPosition({u[0]*std::cos(t) + tangent[0]*std::sin(t) + center[0],
                                        u[1]*std::cos(t) + tangent[1]*std::sin(t) + center[1],
                                        u[2]*std::cos(t) + tangent[2]*std::sin(t) + center[2]});

            t += 2*M_PI/numIp;

            // periodic t
            if(t > 2*M_PI)
                t -= 2*M_PI;
        }

        return points;
    }

    //! compute the shape function for a given point and geometry
    template<class FVGG, class Geometry, class ShapeValues, typename std::enable_if_t<FVGG::discMethod == DiscretizationMethod::box, int> = 0>
    void getShapeValues_(const FVGG& fvGridGeometry, const Geometry& geo, const GlobalPosition& globalPos, ShapeValues& shapeValues)
    {
        const auto ipLocal = geo.local(globalPos);
        const auto& localBasis = this->problem(bulkIdx).fvGridGeometry().feCache().get(geo.type()).localBasis();
        localBasis.evaluateFunction(ipLocal, shapeValues);
    }

    //! compute the shape function for a given point and geometry
    template<class FVGG, class Geometry, class ShapeValues, typename std::enable_if_t<FVGG::discMethod != DiscretizationMethod::box, int> = 0>
    void getShapeValues_(const FVGG& fvGridGeometry, const Geometry& geo, const GlobalPosition& globalPos, ShapeValues& shapeValues)
    {
        DUNE_THROW(Dune::InvalidStateException, "Shape values requested for other discretization than box!");
    }

    //! vector for the volume fraction of the lowdim domain in the bulk domain cells
    std::vector<Scalar> lowDimVolumeInBulkElement_;

    //! the additional stencil for the average operator
    typename ParentType::CouplingStencils bulkCircleStencils_;
};

} // end namespace Dumux

#endif
