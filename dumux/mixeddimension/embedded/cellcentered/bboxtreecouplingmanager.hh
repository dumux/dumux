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
 * \ingroup CCModel
 * \brief Coupling manager for low-dimensional domains embedded in the bulk
 *        domain. Point sources on each integration point are computed by an AABB tree.
 *        Both domain are assumed to be discretized using a cc finite volume scheme.
 */

#ifndef DUMUX_CC_BBOXTREE_EMBEDDEDCOUPLINGMANAGER_HH
#define DUMUX_CC_BBOXTREE_EMBEDDEDCOUPLINGMANAGER_HH

#include <iostream>
#include <fstream>
#include <string>
#include <utility>

#include <dune/common/timer.hh>
#include <dune/geometry/quadraturerules.hh>
#include <dune/geometry/referenceelements.hh>

#include <dumux/common/properties.hh>
#include <dumux/common/geometry/intersectingentities.hh>
#include <dumux/multidomain/couplingmanager.hh>
#include <dumux/mixeddimension/glue/glue.hh>
#include <dumux/mixeddimension/embedded/cellcentered/pointsourcedata.hh>

namespace Dumux {

/*!
 * \ingroup EmbeddedCoupling
 * \ingroup CCModel
 * \brief Manages the coupling between bulk elements and lower dimensional elements
 *        Point sources on each integration point are computed by an AABB tree.
 *        Both domain are assumed to be discretized using a cc finite volume scheme.
 */
template<class MDTraits>
class CCBBoxTreeEmbeddedCouplingManager
: public CouplingManager<MDTraits, CCBBoxTreeEmbeddedCouplingManager<MDTraits>>
{
    using Scalar = typename MDTraits::Scalar;
    static constexpr auto bulkIdx = typename MDTraits::template DomainIdx<0>();
    static constexpr auto lowDimIdx = typename MDTraits::template DomainIdx<1>();
    using SolutionVector = typename MDTraits::SolutionVector;
    using PointSourceData = Dumux::PointSourceDataCircleAverage<MDTraits>;

    // obtain the type tags of the sub problems
    using BulkTypeTag = typename MDTraits::template SubDomainTypeTag<0>;
    using LowDimTypeTag = typename MDTraits::template SubDomainTypeTag<1>;

    using BulkGridView = typename GET_PROP_TYPE(BulkTypeTag, GridView);
    using BulkProblem = typename GET_PROP_TYPE(BulkTypeTag, Problem);
    using BulkPointSource = typename GET_PROP_TYPE(BulkTypeTag, PointSource);
    using BulkPrimaryVariables = typename GET_PROP_TYPE(BulkTypeTag, PrimaryVariables);
    using BulkNumEqVector = typename GET_PROP_TYPE(BulkTypeTag, NumEqVector);
    using BulkElementSolutionVector = typename GET_PROP_TYPE(BulkTypeTag, ElementSolutionVector);
    using BulkVolumeVariables = typename GET_PROP_TYPE(BulkTypeTag, VolumeVariables);
    using BulkElementVolumeVariables = typename GET_PROP_TYPE(BulkTypeTag, ElementVolumeVariables);
    using BulkFVGridGeometry = typename GET_PROP_TYPE(BulkTypeTag, FVGridGeometry);
    using BulkFVElementGeometry = typename BulkFVGridGeometry::LocalView;
    using BulkElementBoundaryTypes = typename GET_PROP_TYPE(BulkTypeTag, ElementBoundaryTypes);
    using BulkElementFluxVariablesCache = typename GET_PROP_TYPE(BulkTypeTag, ElementFluxVariablesCache);
    using BulkElement = typename BulkGridView::template Codim<0>::Entity;

    using LowDimGridView = typename GET_PROP_TYPE(LowDimTypeTag, GridView);
    using LowDimProblem = typename GET_PROP_TYPE(LowDimTypeTag, Problem);
    using LowDimPointSource = typename GET_PROP_TYPE(LowDimTypeTag, PointSource);
    using LowDimPrimaryVariables = typename GET_PROP_TYPE(LowDimTypeTag, PrimaryVariables);
    using LowDimNumEqVector = typename GET_PROP_TYPE(LowDimTypeTag, NumEqVector);
    using LowDimElementSolutionVector = typename GET_PROP_TYPE(LowDimTypeTag, ElementSolutionVector);
    using LowDimVolumeVariables = typename GET_PROP_TYPE(LowDimTypeTag, VolumeVariables);
    using LowDimElementVolumeVariables = typename GET_PROP_TYPE(LowDimTypeTag, ElementVolumeVariables);
    using LowDimFVGridGeometry = typename GET_PROP_TYPE(LowDimTypeTag, FVGridGeometry);
    using LowDimFVElementGeometry = typename LowDimFVGridGeometry::LocalView;
    using LowDimElementBoundaryTypes = typename GET_PROP_TYPE(LowDimTypeTag, ElementBoundaryTypes);
    using LowDimElementFluxVariablesCache = typename GET_PROP_TYPE(LowDimTypeTag, ElementFluxVariablesCache);
    using LowDimElement = typename LowDimGridView::template Codim<0>::Entity;

    enum {
        bulkDim = BulkGridView::dimension,
        lowDimDim = LowDimGridView::dimension,
        dimWorld = BulkGridView::dimensionworld
    };

    using GlobalPosition = Dune::FieldVector<Scalar, dimWorld>;

public:
    /*!
     * \brief Constructor
     */
    CCBBoxTreeEmbeddedCouplingManager(std::shared_ptr<const BulkFVGridGeometry> bulkFvGridGeometry,
                                      std::shared_ptr<const LowDimFVGridGeometry> lowDimFvGridGeometry)
    {

        // Check if we are using the cellcentered method in both domains
        static_assert(lowDimDim == 1, "The bounding box coupling manager only works with one-dimensional low-dim grids");

        using GlueType = CCMixedDimensionGlue<BulkGridView, LowDimGridView>;
        glue_ = std::make_shared<GlueType>();
    }

    /*!
     * \brief Methods to be accessed by main
     */
    // \{

    void init(std::shared_ptr<BulkProblem> bulkProblem,
              std::shared_ptr<LowDimProblem> lowDimProblem,
              const SolutionVector& curSol)
    {
        curSol_ = curSol;
        bulkProblem_ = bulkProblem;
        lowDimProblem_ = lowDimProblem;

        computePointSourceData();
        computeLowDimVolumeFractions();
    }

    /*!
     * \brief Update the solution vector before assembly
     */
    void updateSolution(const SolutionVector& curSol)
    { curSol_ = curSol; }

    // \}

    /*!
     * \brief Methods to be accessed by the assembly
     */
    // \{

    /*!
     * \brief The bulk coupling stencil, i.e. which low-dimensional dofs
     *        the given bulk element dof depends on.
     */
    const std::vector<std::size_t>& couplingStencil(const BulkElement& element, Dune::index_constant<1> lowDimIdx) const
    {
        const std::size_t eIdx = bulkProblem().fvGridGeometry().elementMapper().index(element);
        if (bulkCouplingStencils_.count(eIdx))
            return bulkCouplingStencils_.at(eIdx);
        else
            return emptyStencil_;
    }

    /*!
     * \brief The low dim coupling stencil, i.e. which bulk dofs
     *        the given low dimensional element dof depends on.
     */
    const std::vector<std::size_t>& couplingStencil(const LowDimElement& element, Dune::index_constant<0> bulkIdx) const
    {
        const std::size_t eIdx = lowDimProblem().fvGridGeometry().elementMapper().index(element);
        if (lowDimCouplingStencils_.count(eIdx))
            return lowDimCouplingStencils_.at(eIdx);
        else
            return emptyStencil_;
    }

    //! evaluate coupling residual for the derivative bulk DOF with respect to low dim DOF
    //! we only need to evaluate the part of the residual that will be influence by the low dim DOF
    BulkNumEqVector evalCouplingResidual(const BulkElement& element,
                                         const BulkFVElementGeometry& fvGeometry,
                                         const BulkElementVolumeVariables& curElemVolVars,
                                         const BulkElementBoundaryTypes& elemBcTypes,
                                         const BulkElementFluxVariablesCache& elemFluxVarsCache,
                                         const LowDimElement& lowDimElement)
    {
        const auto bulkElementIdx = bulkProblem().fvGridGeometry().elementMapper().index(element);
        auto&& scv = fvGeometry.scv(bulkElementIdx);
        auto couplingSource = bulkProblem().scvPointSources(element, fvGeometry, curElemVolVars, scv);
        couplingSource *= -scv.volume()*curElemVolVars[scv].extrusionFactor();
        return couplingSource;
    }

    //! evaluate coupling residual for the derivative low dim DOF with respect to bulk DOF
    //! we only need to evaluate the part of the residual that will be influence by the bulk DOF
    LowDimNumEqVector evalCouplingResidual(const LowDimElement& element,
                                           const LowDimFVElementGeometry& fvGeometry,
                                           const LowDimElementVolumeVariables& curElemVolVars,
                                           const LowDimElementBoundaryTypes& elemBcTypes,
                                           const LowDimElementFluxVariablesCache& elemFluxVarsCache,
                                           const BulkElement& bulkElement)
    {
        const auto lowDimElementIdx = lowDimProblem().fvGridGeometry().elementMapper().index(element);
        auto&& scv = fvGeometry.scv(lowDimElementIdx);
        auto couplingSource = lowDimProblem().scvPointSources(element, fvGeometry, curElemVolVars, scv);
        couplingSource *= -scv.volume()*curElemVolVars[scv].extrusionFactor();
        return couplingSource;
    }

    //! evaluate coupling residual for the derivative bulk DOF with respect to low dim DOF
    //! we only need to evaluate the part of the residual that will be influence by the low dim DOF
    template<std::size_t i, class MatrixBlock>
    void addCouplingDerivatives(Dune::index_constant<i> domainId,
                                MatrixBlock& Aij,
                                const BulkElement& element,
                                const BulkFVElementGeometry& fvGeometry,
                                const BulkElementVolumeVariables& curElemVolVars,
                                const LowDimElement& lowDimElement)
    {
        const auto bulkElementIdx = bulkProblem().fvGridGeometry().elementMapper().index(element);

        auto key = std::make_pair(bulkElementIdx, 0);
        if (bulkProblem().pointSourceMap().count(key))
        {
            // call the solDependent function. Herein the user might fill/add values to the point sources
            // we make a copy of the local point sources here
            auto pointSources = bulkProblem().pointSourceMap().at(key);

            // add the point source values to the local residual (negative is sign convention for source term)
            for (const auto& source : pointSources)
                Aij[0][0] -= pointSourceDerivative(source, domainId, lowDimIdx);
        }
    }

    //! evaluate coupling residual for the derivative bulk DOF with respect to low dim DOF
    //! we only need to evaluate the part of the residual that will be influence by the low dim DOF
    template<std::size_t i, class MatrixBlock>
    void addCouplingDerivatives(Dune::index_constant<i> domainId,
                                MatrixBlock& Aij,
                                const LowDimElement& element,
                                const LowDimFVElementGeometry& fvGeometry,
                                const LowDimElementVolumeVariables& curElemVolVars,
                                const BulkElement& bulkElement)
    {
        const auto lowDimElementIdx = lowDimProblem().fvGridGeometry().elementMapper().index(element);
        auto key = std::make_pair(lowDimElementIdx, 0);
        if (lowDimProblem().pointSourceMap().count(key))
        {
            // call the solDependent function. Herein the user might fill/add values to the point sources
            // we make a copy of the local point sources here
            auto pointSources = lowDimProblem().pointSourceMap().at(key);

            // add the point source values to the local residual (negative is sign convention for source term)
            for (const auto& source : pointSources)
                Aij[0][0] -= pointSourceDerivative(source, domainId, bulkIdx);
        }
    }

    //! helper function for the point source derivative (di/dj)
    template<std::size_t i, std::size_t j, class PointSource>
    Scalar pointSourceDerivative(const PointSource& source, Dune::index_constant<i> idI, Dune::index_constant<j> idJ) const
    {
        constexpr Scalar sign = (i == j) ? -1.0 : 1.0;
        // calculate the source derivative
        const Scalar radius = this->radius(source.id());
        const Scalar beta = 2*M_PI/(2*M_PI + std::log(radius));
        return sign*beta*source.quadratureWeight()*source.integrationElement()/source.embeddings();
    }

    /*!
     * \brief Bind the coupling context
     */
    template<class Element, class Assembler>
    void bindCouplingContext(const Element& element, const Assembler& assembler)
    {}

    /*!
     * \brief Update the coupling context
     */
    template<std::size_t i, class Assembler>
    void updateCouplingContext(Dune::index_constant<i> domainId, const BulkElement& element, const BulkPrimaryVariables& priVars, const Assembler& assembler)
    {
        const auto eIdx = bulkProblem().fvGridGeometry().elementMapper().index(element);
        curSol_[bulkIdx][eIdx] = priVars;
    }

    /*!
     * \brief Update the coupling context
     */
    template<std::size_t i, class Assembler>
    void updateCouplingContext(Dune::index_constant<i> domainId, const LowDimElement& element, const LowDimPrimaryVariables& priVars, const Assembler& assembler)
    {
        const auto eIdx = lowDimProblem().fvGridGeometry().elementMapper().index(element);
        curSol_[lowDimIdx][eIdx] = priVars;
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
    void computePointSourceData(std::size_t order = 1,
                                bool verbose = false)
    {
        // Initialize the bulk bounding box tree
        const auto& bulkTree = bulkProblem().fvGridGeometry().boundingBoxTree();

        // initilize the maps
        // do some logging and profiling
        Dune::Timer watch;
        std::cout << "Initializing the point sources..." << std::endl;

        // clear all internal members like pointsource vectors and stencils
        // initializes the point source id counter
        clear();

        // iterate over all lowdim elements
        for (const auto& lowDimElement : elements(lowDimGridView()))
        {
            // get the Gaussian quadrature rule for the low dim element
            const auto lowDimGeometry = lowDimElement.geometry();
            const auto& quad = Dune::QuadratureRules<Scalar, lowDimDim>::rule(lowDimGeometry.type(), order);

            const auto lowDimElementIdx = lowDimProblem().fvGridGeometry().elementMapper().index(lowDimElement);

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
                const auto radius = this->lowDimProblem().spatialParams().radius(lowDimElementIdx);
                const auto normal = lowDimGeometry.corner(1)-lowDimGeometry.corner(0);
                const auto length = 2*M_PI*radius/numIp;

                const auto circlePoints = getCirclePoints_(globalPos, normal, radius, numIp);
                std::vector<Scalar> circleIpWeight;
                std::vector<std::size_t> circleStencil;
                for (const auto& p : circlePoints)
                {
                    const auto circleBulkElementIndices = intersectingEntities(p, bulkTree);
                    if (circleBulkElementIndices.empty())
                        continue;

                    for (auto bulkElementIdx : circleBulkElementIndices)
                    {
                        circleStencil.push_back(bulkElementIdx);
                        circleIpWeight.push_back(length);
                    }
                }

                // export low dim circle stencil
                lowDimCouplingStencils_[lowDimElementIdx].insert(lowDimCouplingStencils_[lowDimElementIdx].end(),
                                                                 circleStencil.begin(),
                                                                 circleStencil.end());

                for (auto bulkElementIdx : bulkElementIndices)
                {
                    const auto id = idCounter_++;
                    const auto ie = lowDimGeometry.integrationElement(qp.position());
                    const auto qpweight = qp.weight();

                    bulkPointSources_.emplace_back(globalPos, id, qpweight, ie, std::vector<std::size_t>({bulkElementIdx}));
                    bulkPointSources_.back().setEmbeddings(bulkElementIndices.size());
                    lowDimPointSources_.emplace_back(globalPos, id, qpweight, ie, std::vector<std::size_t>({lowDimElementIdx}));
                    lowDimPointSources_.back().setEmbeddings(bulkElementIndices.size());

                    // pre compute additional data used for the evaluation of
                    // the actual solution dependent source term
                    PointSourceData psData;
                    psData.addLowDimInterpolation(lowDimElementIdx);
                    psData.addBulkInterpolation(bulkElementIdx);
                    // add data needed to compute integral over the circle
                    psData.addCircleInterpolation(circleIpWeight, circleStencil);

                    // publish point source data in the global vector
                    pointSourceData_.push_back(psData);

                    // compute the coupling stencils
                    bulkCouplingStencils_[bulkElementIdx].push_back(lowDimElementIdx);

                    // export bulk circle stencil
                    bulkCircleStencils_[bulkElementIdx].insert(bulkCircleStencils_[bulkElementIdx].end(),
                                                               circleStencil.begin(),
                                                               circleStencil.end());
                }
            }
        }

        // make the circle stencil unique
        for (auto&& stencil : bulkCircleStencils_)
        {
            std::sort(stencil.second.begin(), stencil.second.end());
            stencil.second.erase(std::unique(stencil.second.begin(), stencil.second.end()), stencil.second.end());
            stencil.second.erase(std::remove_if(stencil.second.begin(), stencil.second.end(),
                                               [&](auto i){ return i == stencil.first; }),
                                 stencil.second.end());
        }

        // low dim coupling stencil
        for (auto&& stencil : lowDimCouplingStencils_)
        {
            std::sort(stencil.second.begin(), stencil.second.end());
            stencil.second.erase(std::unique(stencil.second.begin(), stencil.second.end()), stencil.second.end());
        }

        // bulk coupling stencils
        for (auto&& stencil : bulkCouplingStencils_)
        {
            std::sort(stencil.second.begin(), stencil.second.end());
            stencil.second.erase(std::unique(stencil.second.begin(), stencil.second.end()), stencil.second.end());
        }

        std::cout << "took " << watch.elapsed() << " seconds." << std::endl;
    }

    //! Compute the low dim volume fraction in the bulk domain cells
    void computeLowDimVolumeFractions()
    {
        // resize the storage vector
        lowDimVolumeInBulkElement_.resize(bulkGridView().size(0));

        // compute the low dim volume fractions
        for (const auto& is : intersections(*glue_))
        {
            // all inside elements are identical...
            const auto& inside = is.inside(0);
            const auto intersectionGeometry = is.geometry();
            const std::size_t lowDimElementIdx = lowDimProblem().fvGridGeometry().elementMapper().index(inside);

            // compute the volume the low-dim domain occupies in the bulk domain if it were full-dimensional
            const auto radius = lowDimProblem().spatialParams().radius(lowDimElementIdx);
            for (std::size_t outsideIdx = 0; outsideIdx < is.neighbor(0); ++outsideIdx)
            {
                const auto& outside = is.outside(outsideIdx);
                const std::size_t bulkElementIdx = bulkProblem().fvGridGeometry().elementMapper().index(outside);
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
        const auto& data = pointSourceData_[id];
        return lowDimProblem().spatialParams().radius(data.lowDimElementIdx());
    }

    //! The volume the lower dimensional domain occupies in the bulk domain element
    // For one-dimensional low dim domain we assume radial tubes
    Scalar lowDimVolume(const BulkElement& element) const
    {
        auto eIdx = bulkProblem().elementMapper().index(element);
        return lowDimVolumeInBulkElement_[eIdx];
    }

    //! The volume fraction the lower dimensional domain occupies in the bulk domain element
    // For one-dimensional low dim domain we assume radial tubes
    Scalar lowDimVolumeFraction(const BulkElement& element) const
    {
        auto totalVolume = element.geometry().volume();
        return lowDimVolume(element) / totalVolume;
    }

    //! Return a reference to the pointSource data
    const PointSourceData& pointSourceData(std::size_t id) const
    {
        return pointSourceData_[id];
    }

    //! Return a reference to the bulk problem
    const BulkProblem& bulkProblem() const
    {
        return *bulkProblem_;
    }

    //! Return a reference to the low dimensional problem
    const LowDimProblem& lowDimProblem() const
    {
        return *lowDimProblem_;
    }

    //! Return a reference to the bulk problem
    BulkProblem& bulkProblem()
    {
        return *bulkProblem_;
    }

    //! Return a reference to the low dimensional problem
    LowDimProblem& lowDimProblem()
    {
        return *lowDimProblem_;
    }

    //! Return a reference to the bulk problem
    const BulkGridView& bulkGridView() const
    {
        return bulkProblem().fvGridGeometry().gridView();
    }

    //! Return a reference to the low dimensional problem
    const LowDimGridView& lowDimGridView() const
    {
        return lowDimProblem().fvGridGeometry().gridView();
    }

    //! Return a reference to the point sources
    const std::vector<BulkPointSource>& bulkPointSources() const
    {
        return bulkPointSources_;
    }

    //! Return a reference to the low dimensional problem
    const std::vector<LowDimPointSource>& lowDimPointSources() const
    {
        return lowDimPointSources_;
    }

    //! Return data for a bulk point source with the identifier id
    BulkPrimaryVariables bulkPriVars(std::size_t id) const
    {
        auto& data = pointSourceData_[id];
        return data.interpolateBulk(curSol_[bulkIdx]);
    }

    //! Return data for a low dim point source with the identifier id
    LowDimPrimaryVariables lowDimPriVars(std::size_t id) const
    {
        auto& data = pointSourceData_[id];
        return data.interpolateLowDim(curSol_[lowDimIdx]);
    }

    // //! Compute bulk volume variables for the identifier id
    // BulkVolumeVariables bulkVolVars(std::size_t id) const
    // {
    //     // use circle interpolated data and construct volVar object for the interpolated privars
    //     auto& data = pointSourceData_[id];
    //     auto bulkPriVars = data.interpolateBulk(bulkProblem().model().curSol());
    //
    //     const auto element = bulkProblem().model().fvGridGeometry().element(data.bulkElementIdx());
    //     auto fvGeometry = localView(bulkProblem().model().fvGridGeometry());
    //     fvGeometry.bindElement(element);
    //
    //     BulkVolumeVariables volVars;
    //     volVars.update(BulkElementSolutionVector(bulkPriVars),
    //                    bulkProblem(),
    //                    element,
    //                    fvGeometry.scv(data.bulkElementIdx()));
    //
    //     return volVars;
    // }

    // //! Compute lowDim volume variables for the identifier id
    // LowDimVolumeVariables lowDimVolVars(std::size_t id) const
    // {
    //     const auto& data = pointSourceData_[id];
    //     const auto element = lowDimProblem().model().fvGridGeometry().element(data.lowDimElementIdx());
    //     auto fvGeometry = localView(lowDimProblem().model().fvGridGeometry());
    //     fvGeometry.bindElement(element);
    //
    //     const auto& curSol = lowDimProblem().model().curSol();
    //     LowDimVolumeVariables volVars;
    //     volVars.update(lowDimProblem().model().elementSolution(element, curSol),
    //                    lowDimProblem(),
    //                    element,
    //                    fvGeometry.scv(data.lowDimElementIdx()));
    //
    //     return volVars;
    // }
    // \}

    //! Clear all internal data members
    void clear()
    {
        bulkPointSources_.clear();
        lowDimPointSources_.clear();
        pointSourceData_.clear();
        bulkCouplingStencils_.clear();
        lowDimCouplingStencils_.clear();
        idCounter_ = 0;
    }

    const std::vector<std::size_t>& getAdditionalDofDependencies(Dune::index_constant<0> bulkDomain, std::size_t bulkElementIdx) const
    { return bulkCircleStencils_.count(bulkElementIdx) ? bulkCircleStencils_.at(bulkElementIdx) : emptyStencil_; }


    const std::vector<std::size_t>& getAdditionalDofDependencies(Dune::index_constant<1> lowDimDomain, std::size_t lowDimElementIdx) const
    { return emptyStencil_; }

protected:
    //! Return reference to point source data vector member
    std::vector<PointSourceData>& pointSourceData()
    { return pointSourceData_; }

    //! Return reference to bulk point sources
    std::vector<BulkPointSource>& bulkPointSources()
    { return bulkPointSources_; }

    //! Return reference to low dim point sources
    std::vector<LowDimPointSource>& lowDimPointSources()
    { return lowDimPointSources_; }

    //! Return reference to bulk coupling stencil member
    std::unordered_map<std::size_t, std::vector<std::size_t> >& bulkCouplingStencils()
    { return bulkCouplingStencils_; }

    //! Return reference to low dim coupling stencil member
    std::unordered_map<std::size_t, std::vector<std::size_t> >& lowDimCouplingStencils()
    { return lowDimCouplingStencils_; }

    //! Return a reference to an empty stencil
    std::vector<std::size_t>& emptyStencil()
    { return emptyStencil_; }

private:

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

    std::shared_ptr<BulkProblem> bulkProblem_;
    std::shared_ptr<LowDimProblem> lowDimProblem_;

    std::vector<BulkPointSource> bulkPointSources_;
    std::vector<LowDimPointSource> lowDimPointSources_;

    mutable std::vector<PointSourceData> pointSourceData_;

    std::unordered_map<std::size_t, std::vector<std::size_t> > bulkCouplingStencils_;
    std::unordered_map<std::size_t, std::vector<std::size_t> > lowDimCouplingStencils_;
    std::vector<std::size_t> emptyStencil_;

     // circle stencils
    std::unordered_map<std::size_t, std::vector<std::size_t> > bulkCircleStencils_;

    //! vector for the volume fraction of the lowdim domain in the bulk domain cells
    std::vector<Scalar> lowDimVolumeInBulkElement_;

    //! id generator for point sources
    std::size_t idCounter_ = 0;

    //! The glue object
    std::shared_ptr<CCMixedDimensionGlue<BulkGridView, LowDimGridView>> glue_;

    ////////////////////////////////////////////////////////////////////////////
    //! The coupling context
    ////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////
    //! TODO: this is the simplest context -> just the solutionvector
    ////////////////////////////////////////////////////////////////////////////
    SolutionVector curSol_;
};

} // end namespace Dumux

#endif
