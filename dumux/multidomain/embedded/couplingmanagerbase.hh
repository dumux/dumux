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

#ifndef DUMUX_MULTIDOMAIN_EMBEDDED_COUPLINGMANAGERBASE_HH
#define DUMUX_MULTIDOMAIN_EMBEDDED_COUPLINGMANAGERBASE_HH

#include <iostream>
#include <fstream>
#include <string>
#include <utility>
#include <unordered_map>

#include <dune/common/timer.hh>
#include <dune/geometry/quadraturerules.hh>

#include <dumux/common/properties.hh>
#include <dumux/geometry/distance.hh>
#include <dumux/geometry/intersectingentities.hh>
#include <dumux/discretization/method.hh>
#include <dumux/multidomain/couplingmanager.hh>
#include <dumux/multidomain/glue.hh>
#include <dumux/multidomain/embedded/pointsourcedata.hh>
#include <dumux/multidomain/embedded/integrationpointsource.hh>

namespace Dumux {

//! the default point source traits
template<class MDTraits>
struct DefaultPointSourceTraits
{
private:
    template<std::size_t i> using SubDomainTypeTag = typename MDTraits::template SubDomain<i>::TypeTag;
    template<std::size_t i> using GridGeometry = GetPropType<SubDomainTypeTag<i>, Properties::GridGeometry>;
    template<std::size_t i> using NumEqVector = GetPropType<SubDomainTypeTag<i>, Properties::NumEqVector>;
public:
    //! export the point source type for domain i
    template<std::size_t i>
    using PointSource = IntegrationPointSource<typename GridGeometry<i>::GlobalCoordinate, NumEqVector<i>>;

    //! export the point source helper type  for domain i
    template<std::size_t i>
    using PointSourceHelper = IntegrationPointSourceHelper;

    //! export the point source data type
    using PointSourceData = Dumux::PointSourceData<MDTraits>;
};

/*!
 * \ingroup EmbeddedCoupling
 * \brief Manages the coupling between bulk elements and lower dimensional elements
 *        Point sources on each integration point are computed by an AABB tree.
 */
template<class MDTraits, class Implementation, class PSTraits = DefaultPointSourceTraits<MDTraits>>
class EmbeddedCouplingManagerBase
: public CouplingManager<MDTraits>
{
    using ParentType = CouplingManager<MDTraits>;
    using Scalar = typename MDTraits::Scalar;
    static constexpr auto bulkIdx = typename MDTraits::template SubDomain<0>::Index();
    static constexpr auto lowDimIdx = typename MDTraits::template SubDomain<1>::Index();
    using SolutionVector = typename MDTraits::SolutionVector;
    using PointSourceData = typename PSTraits::PointSourceData;

    // the sub domain type tags
    template<std::size_t id> using PointSource = typename PSTraits::template PointSource<id>;
    template<std::size_t id> using SubDomainTypeTag = typename MDTraits::template SubDomain<id>::TypeTag;
    template<std::size_t id> using Problem = GetPropType<SubDomainTypeTag<id>, Properties::Problem>;
    template<std::size_t id> using PrimaryVariables = GetPropType<SubDomainTypeTag<id>, Properties::PrimaryVariables>;
    template<std::size_t id> using GridGeometry = GetPropType<SubDomainTypeTag<id>, Properties::GridGeometry>;
    template<std::size_t id> using GridView = typename GridGeometry<id>::GridView;
    template<std::size_t id> using ElementMapper = typename GridGeometry<id>::ElementMapper;
    template<std::size_t id> using Element = typename GridView<id>::template Codim<0>::Entity;
    template<std::size_t id> using GridIndex = typename IndexTraits<GridView<id>>::GridIndex;
    template<std::size_t id> using CouplingStencil = std::vector<GridIndex<id>>;

    static constexpr int bulkDim = GridView<bulkIdx>::dimension;
    static constexpr int lowDimDim = GridView<lowDimIdx>::dimension;
    static constexpr int dimWorld = GridView<bulkIdx>::dimensionworld;

    template<std::size_t id>
    static constexpr bool isBox()
    { return GridGeometry<id>::discMethod == DiscretizationMethod::box; }

    using GlobalPosition = typename Element<bulkIdx>::Geometry::GlobalCoordinate;
    using GlueType = MultiDomainGlue<GridView<bulkIdx>, GridView<lowDimIdx>, ElementMapper<bulkIdx>, ElementMapper<lowDimIdx>>;

public:
    //! export traits
    using MultiDomainTraits = MDTraits;
    //! export the point source traits
    using PointSourceTraits = PSTraits;
    //! export stencil types
    template<std::size_t id> using CouplingStencils = std::unordered_map<GridIndex<id>, CouplingStencil<id>>;

    /*!
    * \brief call this after grid adaption
    */
    void updateAfterGridAdaption(std::shared_ptr<const GridGeometry<bulkIdx>> bulkGridGeometry,
                                 std::shared_ptr<const GridGeometry<lowDimIdx>> lowDimGridGeometry)
    {
        glue_ = std::make_shared<GlueType>();
    }

    /*!
     * \brief Constructor
     */
    EmbeddedCouplingManagerBase(std::shared_ptr<const GridGeometry<bulkIdx>> bulkGridGeometry,
                                std::shared_ptr<const GridGeometry<lowDimIdx>> lowDimGridGeometry)
    {
        updateAfterGridAdaption(bulkGridGeometry, lowDimGridGeometry);
    }

    /*!
     * \brief Methods to be accessed by main
     */
    // \{

    void init(std::shared_ptr<Problem<bulkIdx>> bulkProblem,
              std::shared_ptr<Problem<lowDimIdx>> lowDimProblem,
              const SolutionVector& curSol)
    {
        this->updateSolution(curSol);
        this->setSubProblems(std::make_tuple(bulkProblem, lowDimProblem));

        integrationOrder_ = getParam<int>("MixedDimension.IntegrationOrder", 1);
        asImp_().computePointSourceData(integrationOrder_);
    }

    // \}

    /*!
     * \brief Methods to be accessed by the assembly
     */
    // \{

    /*!
     * \brief returns an iteratable container of all indices of degrees of freedom of domain j
     *        that couple with / influence the element residual of the given element of domain i
     *
     * \param domainI the domain index of domain i
     * \param element the coupled element of domain í
     * \param domainJ the domain index of domain j
     *
     * \note  The element residual definition depends on the discretization scheme of domain i
     *        box: a container of the residuals of all sub control volumes
     *        cc : the residual of the (sub) control volume
     *        fem: the residual of the element
     * \note  This function has to be implemented by all coupling managers for all combinations of i and j
     */
    template<std::size_t i, std::size_t j>
    const CouplingStencil<j>& couplingStencil(Dune::index_constant<i> domainI,
                                              const Element<i>& element,
                                              Dune::index_constant<j> domainJ) const
    {
        static_assert(i != j, "A domain cannot be coupled to itself!");

        const auto eIdx = this->problem(domainI).gridGeometry().elementMapper().index(element);
        if (couplingStencils(domainI).count(eIdx))
            return couplingStencils(domainI).at(eIdx);
        else
            return emptyStencil(domainI);
    }

    /*!
     * \brief evaluates the element residual of a coupled element of domain i which depends on the variables
     *        at the degree of freedom with index dofIdxGlobalJ of domain j
     *
     * \param domainI the domain index of domain i
     * \param localAssemblerI the local assembler assembling the element residual of an element of domain i
     * \param domainJ the domain index of domain j
     * \param dofIdxGlobalJ the index of the degree of freedom of domain j which has an influence on the element residual of domain i
     *
     * \note  we only need to evaluate the source contribution to the residual here as the coupling term is the source
     * \return the element residual
     */
    template<std::size_t i, std::size_t j, class LocalAssemblerI>
    decltype(auto) evalCouplingResidual(Dune::index_constant<i> domainI,
                                        const LocalAssemblerI& localAssemblerI,
                                        Dune::index_constant<j> domainJ,
                                        std::size_t dofIdxGlobalJ)
    {
        static_assert(i != j, "A domain cannot be coupled to itself!");

        typename LocalAssemblerI::LocalResidual::ElementResidualVector residual;

        const auto& element = localAssemblerI.element();
        const auto& fvGeometry = localAssemblerI.fvGeometry();
        const auto& curElemVolVars = localAssemblerI.curElemVolVars();

        residual.resize(fvGeometry.numScv());
        for (const auto& scv : scvs(fvGeometry))
        {
            auto couplingSource = this->problem(domainI).scvPointSources(element, fvGeometry, curElemVolVars, scv);
            couplingSource += this->problem(domainI).source(element, fvGeometry, curElemVolVars, scv);
            couplingSource *= -GridGeometry<i>::Extrusion::volume(scv)*curElemVolVars[scv].extrusionFactor();
            residual[scv.indexInElement()] = couplingSource;
        }
        return residual;
    }

    // \}

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
        clear();

        // precompute the vertex indices for efficiency for the box method
        this->precomputeVertexIndices(bulkIdx);
        this->precomputeVertexIndices(lowDimIdx);

        const auto& bulkGridGeometry = this->problem(bulkIdx).gridGeometry();
        const auto& lowDimGridGeometry = this->problem(lowDimIdx).gridGeometry();

        // intersect the bounding box trees
        glueGrids();

        pointSourceData_.reserve(this->glue().size());
        averageDistanceToBulkCell_.reserve(this->glue().size());
        for (const auto& is : intersections(this->glue()))
        {
            // all inside elements are identical...
            const auto& inside = is.targetEntity(0);
            // get the intersection geometry for integrating over it
            const auto intersectionGeometry = is.geometry();

            // get the Gaussian quadrature rule for the local intersection
            const auto& quad = Dune::QuadratureRules<Scalar, lowDimDim>::rule(intersectionGeometry.type(), order);
            const std::size_t lowDimElementIdx = lowDimGridGeometry.elementMapper().index(inside);

            // iterate over all quadrature points
            for (auto&& qp : quad)
            {
                // compute the coupling stencils
                for (std::size_t outsideIdx = 0; outsideIdx < is.numDomainNeighbors(); ++outsideIdx)
                {
                    const auto& outside = is.domainEntity(outsideIdx);
                    const std::size_t bulkElementIdx = bulkGridGeometry.elementMapper().index(outside);

                    // each quadrature point will be a point source for the sub problem
                    const auto globalPos = intersectionGeometry.global(qp.position());
                    const auto id = idCounter_++;
                    const auto qpweight = qp.weight();
                    const auto ie = intersectionGeometry.integrationElement(qp.position());
                    pointSources(bulkIdx).emplace_back(globalPos, id, qpweight, ie, std::vector<std::size_t>({bulkElementIdx}));
                    pointSources(bulkIdx).back().setEmbeddings(is.numDomainNeighbors());
                    pointSources(lowDimIdx).emplace_back(globalPos, id, qpweight, ie, std::vector<std::size_t>({lowDimElementIdx}));
                    pointSources(lowDimIdx).back().setEmbeddings(is.numDomainNeighbors());

                    // pre compute additional data used for the evaluation of
                    // the actual solution dependent source term
                    PointSourceData psData;

                    if constexpr (isBox<lowDimIdx>())
                    {
                        using ShapeValues = std::vector<Dune::FieldVector<Scalar, 1> >;
                        const auto lowDimGeometry = this->problem(lowDimIdx).gridGeometry().element(lowDimElementIdx).geometry();
                        ShapeValues shapeValues;
                        this->getShapeValues(lowDimIdx, this->problem(lowDimIdx).gridGeometry(), lowDimGeometry, globalPos, shapeValues);
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
                        const auto bulkGeometry = this->problem(bulkIdx).gridGeometry().element(bulkElementIdx).geometry();
                        ShapeValues shapeValues;
                        this->getShapeValues(bulkIdx, this->problem(bulkIdx).gridGeometry(), bulkGeometry, globalPos, shapeValues);
                        psData.addBulkInterpolation(shapeValues, this->vertexIndices(bulkIdx, bulkElementIdx), bulkElementIdx);
                    }
                    else
                    {
                        psData.addBulkInterpolation(bulkElementIdx);
                    }

                    // publish point source data in the global vector
                    this->pointSourceData().emplace_back(std::move(psData));

                    // compute average distance to bulk cell
                    averageDistanceToBulkCell_.push_back(averageDistancePointGeometry(globalPos, outside.geometry()));

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

        // make stencils unique
        using namespace Dune::Hybrid;
        forEach(integralRange(std::integral_constant<std::size_t, MDTraits::numSubDomains>{}), [&](const auto domainIdx)
        {
            for (auto&& stencil : this->couplingStencils(domainIdx))
            {
                std::sort(stencil.second.begin(), stencil.second.end());
                stencil.second.erase(std::unique(stencil.second.begin(), stencil.second.end()), stencil.second.end());
            }
        });

        std::cout << "took " << watch.elapsed() << " seconds." << std::endl;
    }

    /*!
     * \brief Methods to be accessed by the subproblems
     */
    // \{

    //! Return a reference to the pointSource data
    const PointSourceData& pointSourceData(std::size_t id) const
    { return pointSourceData_[id]; }

    //! Return a reference to the bulk problem
    template<std::size_t id>
    const GridView<id>& gridView(Dune::index_constant<id> domainIdx) const
    { return this->problem(domainIdx).gridGeometry().gridView(); }

    //! Return data for a bulk point source with the identifier id
    PrimaryVariables<bulkIdx> bulkPriVars(std::size_t id) const
    { return pointSourceData_[id].interpolateBulk(this->curSol()[bulkIdx]); }

    //! Return data for a low dim point source with the identifier id
    PrimaryVariables<lowDimIdx> lowDimPriVars(std::size_t id) const
    { return pointSourceData_[id].interpolateLowDim(this->curSol()[lowDimIdx]); }

    //! return the average distance to the coupled bulk cell center
    Scalar averageDistance(std::size_t id) const
    { return averageDistanceToBulkCell_[id]; }

    //! Return reference to bulk point sources
    const std::vector<PointSource<bulkIdx>>& bulkPointSources() const
    { return std::get<bulkIdx>(pointSources_); }

    //! Return reference to low dim point sources
    const std::vector<PointSource<lowDimIdx>>& lowDimPointSources() const
    { return std::get<lowDimIdx>(pointSources_); }

    //! Return the point source if domain i
    template<std::size_t i>
    const std::vector<PointSource<i>>& pointSources(Dune::index_constant<i> dom) const
    { return std::get<i>(pointSources_); }

    //! Return reference to bulk coupling stencil member of domain i
    template<std::size_t i>
    const CouplingStencils<i>& couplingStencils(Dune::index_constant<i> dom) const
    { return std::get<i>(couplingStencils_); }

    //! Return reference to point source data vector member
    const std::vector<PointSourceData>& pointSourceData() const
    { return pointSourceData_; }

    //! Return a reference to an empty stencil
    template<std::size_t i>
    const CouplingStencil<i>& emptyStencil(Dune::index_constant<i> dom) const
    { return std::get<i>(emptyStencil_); }

protected:

    //! computes the vertex indices per element for the box method
    template<std::size_t id>
    void precomputeVertexIndices(Dune::index_constant<id> domainIdx)
    {
        // fill helper structure for box discretization
        if constexpr (isBox<domainIdx>())
        {
            this->vertexIndices(domainIdx).resize(gridView(domainIdx).size(0));
            for (const auto& element : elements(gridView(domainIdx)))
            {
                constexpr int dim = GridView<domainIdx>::dimension;
                const auto eIdx = this->problem(domainIdx).gridGeometry().elementMapper().index(element);
                this->vertexIndices(domainIdx, eIdx).resize(element.subEntities(dim));
                for (int i = 0; i < element.subEntities(dim); ++i)
                    this->vertexIndices(domainIdx, eIdx)[i] = this->problem(domainIdx).gridGeometry().vertexMapper().subIndex(element, i, dim);
            }
        }
    }

    //! compute the shape function for a given point and geometry
    template<std::size_t i, class FVGG, class Geometry, class ShapeValues>
    void getShapeValues(Dune::index_constant<i> domainI, const FVGG& gridGeometry, const Geometry& geo, const GlobalPosition& globalPos, ShapeValues& shapeValues)
    {
        if constexpr (FVGG::discMethod == DiscretizationMethod::box)
        {
            const auto ipLocal = geo.local(globalPos);
            const auto& localBasis = this->problem(domainI).gridGeometry().feCache().get(geo.type()).localBasis();
            localBasis.evaluateFunction(ipLocal, shapeValues);
        }
        else
            DUNE_THROW(Dune::InvalidStateException, "Shape values requested for other discretization than box!");
    }

    //! Clear all internal data members
    void clear()
    {
        pointSources(bulkIdx).clear();
        pointSources(lowDimIdx).clear();
        couplingStencils(bulkIdx).clear();
        couplingStencils(lowDimIdx).clear();
        vertexIndices(bulkIdx).clear();
        vertexIndices(lowDimIdx).clear();
        pointSourceData_.clear();
        averageDistanceToBulkCell_.clear();

        idCounter_ = 0;
    }

    //! compute the intersections between the two grids
    void glueGrids()
    {
        const auto& bulkGridGeometry = this->problem(bulkIdx).gridGeometry();
        const auto& lowDimGridGeometry = this->problem(lowDimIdx).gridGeometry();

        // intersect the bounding box trees
        glue_->build(bulkGridGeometry.boundingBoxTree(), lowDimGridGeometry.boundingBoxTree());
    }

    //! Return reference to point source data vector member
    std::vector<PointSourceData>& pointSourceData()
    { return pointSourceData_; }

    //! Return reference to average distances to bulk cell
    std::vector<Scalar>& averageDistanceToBulkCell()
    { return averageDistanceToBulkCell_; }

    //! Return the point source if domain i
    template<std::size_t i>
    std::vector<PointSource<i>>& pointSources(Dune::index_constant<i> dom)
    { return std::get<i>(pointSources_); }

    //! Return reference to bulk coupling stencil member of domain i
    template<std::size_t i>
    CouplingStencils<i>& couplingStencils(Dune::index_constant<i> dom)
    { return std::get<i>(couplingStencils_); }

    //! Return a reference to the vertex indices
    template<std::size_t i>
    std::vector<GridIndex<i>>& vertexIndices(Dune::index_constant<i> dom, GridIndex<i> eIdx)
    { return std::get<i>(vertexIndices_)[eIdx]; }

    //! Return a reference to the vertex indices container
    template<std::size_t i>
    std::vector<std::vector<GridIndex<i>>>& vertexIndices(Dune::index_constant<i> dom)
    { return std::get<i>(vertexIndices_); }

    const GlueType& glue() const
    { return *glue_; }

    //! Returns the implementation of the problem (i.e. static polymorphism)
    Implementation &asImp_()
    { return *static_cast<Implementation *>(this); }

    //! \copydoc asImp_()
    const Implementation &asImp_() const
    { return *static_cast<const Implementation *>(this); }

    //! id generator for point sources
    std::size_t idCounter_ = 0;

private:

    //! the point source in both domains
    std::tuple<std::vector<PointSource<bulkIdx>>, std::vector<PointSource<lowDimIdx>>> pointSources_;
    std::vector<PointSourceData> pointSourceData_;
    std::vector<Scalar> averageDistanceToBulkCell_;

    //! Stencil data
    std::tuple<std::vector<std::vector<GridIndex<bulkIdx>>>,
               std::vector<std::vector<GridIndex<lowDimIdx>>>> vertexIndices_;
    std::tuple<CouplingStencils<bulkIdx>, CouplingStencils<lowDimIdx>> couplingStencils_;
    std::tuple<CouplingStencil<bulkIdx>, CouplingStencil<lowDimIdx>> emptyStencil_;

    //! The glue object
    std::shared_ptr<GlueType> glue_;

    //! integration order for coupling source
    int integrationOrder_ = 1;
};

} // end namespace Dumux

#endif
