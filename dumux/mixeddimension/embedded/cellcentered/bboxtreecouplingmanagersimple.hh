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

#ifndef DUMUX_CC_BBOXTREE_EMBEDDEDCOUPLINGMANAGER_SIMPLE_HH
#define DUMUX_CC_BBOXTREE_EMBEDDEDCOUPLINGMANAGER_SIMPLE_HH

#include <iostream>
#include <fstream>
#include <string>
#include <utility>

#include <dune/common/timer.hh>
#include <dune/geometry/quadraturerules.hh>
#include <dune/geometry/referenceelements.hh>

#include <dumux/mixeddimension/glue/glue.hh>
#include <dumux/mixeddimension/embedded/cellcentered/pointsourcedata.hh>

namespace Dumux
{

namespace Properties
{
// Property forward declarations
NEW_PROP_TAG(BulkProblemTypeTag);
NEW_PROP_TAG(LowDimProblemTypeTag);
NEW_PROP_TAG(GridView);
} // namespace Properties

/*!
 * \ingroup EmbeddedCoupling
 * \ingroup CCModel
 * \brief Manages the coupling between bulk elements and lower dimensional elements
 *        Point sources on each integration point are computed by an AABB tree.
 *        Both domain are assumed to be discretized using a cc finite volume scheme.
 */
template<class TypeTag>
class CCBBoxTreeEmbeddedCouplingManagerSimple
{
    using Implementation = typename GET_PROP_TYPE(TypeTag, CouplingManager);
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using Problem = typename GET_PROP_TYPE(TypeTag, Problem);
    using PointSourceData = Dumux::PointSourceData<TypeTag>;

    // obtain the type tags of the sub problems
    using BulkProblemTypeTag = typename GET_PROP_TYPE(TypeTag, BulkProblemTypeTag);
    using LowDimProblemTypeTag = typename GET_PROP_TYPE(TypeTag, LowDimProblemTypeTag);

    using BulkGridView = typename GET_PROP_TYPE(BulkProblemTypeTag, GridView);
    using BulkProblem = typename GET_PROP_TYPE(BulkProblemTypeTag, Problem);
    using BulkPointSource = typename GET_PROP_TYPE(BulkProblemTypeTag, PointSource);
    using BulkPrimaryVariables = typename GET_PROP_TYPE(BulkProblemTypeTag, PrimaryVariables);
    using BulkElementSolutionVector = typename GET_PROP_TYPE(BulkProblemTypeTag, ElementSolutionVector);
    using BulkVolumeVariables = typename GET_PROP_TYPE(BulkProblemTypeTag, VolumeVariables);
    using BulkElementVolumeVariables = typename GET_PROP_TYPE(BulkProblemTypeTag, ElementVolumeVariables);
    using BulkFVElementGeometry = typename GET_PROP_TYPE(BulkProblemTypeTag, FVElementGeometry);
    using BulkElementBoundaryTypes = typename GET_PROP_TYPE(BulkProblemTypeTag, ElementBoundaryTypes);
    using BulkElementFluxVariablesCache = typename GET_PROP_TYPE(BulkProblemTypeTag, ElementFluxVariablesCache);
    using BulkElement = typename BulkGridView::template Codim<0>::Entity;

    using LowDimGridView = typename GET_PROP_TYPE(LowDimProblemTypeTag, GridView);
    using LowDimProblem = typename GET_PROP_TYPE(LowDimProblemTypeTag, Problem);
    using LowDimPointSource = typename GET_PROP_TYPE(LowDimProblemTypeTag, PointSource);
    using LowDimPrimaryVariables = typename GET_PROP_TYPE(LowDimProblemTypeTag, PrimaryVariables);
    using LowDimVolumeVariables = typename GET_PROP_TYPE(BulkProblemTypeTag, VolumeVariables);
    using LowDimElementVolumeVariables = typename GET_PROP_TYPE(LowDimProblemTypeTag, ElementVolumeVariables);
    using LowDimFVElementGeometry = typename GET_PROP_TYPE(LowDimProblemTypeTag, FVElementGeometry);
    using LowDimElementBoundaryTypes = typename GET_PROP_TYPE(LowDimProblemTypeTag, ElementBoundaryTypes);
    using LowDimElementFluxVariablesCache = typename GET_PROP_TYPE(LowDimProblemTypeTag, ElementFluxVariablesCache);
    using LowDimElement = typename LowDimGridView::template Codim<0>::Entity;

    enum {
        bulkDim = BulkGridView::dimension,
        lowDimDim = LowDimGridView::dimension,
        dimWorld = BulkGridView::dimensionworld
    };

public:

    /*!
     * \brief Constructor
     */
    CCBBoxTreeEmbeddedCouplingManagerSimple(BulkProblem& bulkProblem, LowDimProblem& lowDimProblem)
    : bulkProblem_(bulkProblem), lowDimProblem_(lowDimProblem)
    {
        // Check if we are using the cellcentered method in both domains
        static_assert(lowDimDim == 1, "The bounding box coupling manager only works with one-dimensional low-dim grids");
    }

    /*!
     * \brief Methods to be accessed by the coupled problem
     */
    // \{

    /*!
     * \brief Called by the coupled problem before
     *        initializing the subproblems / models
     */
    void preInit()
    {
        asImp_().computePointSourceData();
    }

    /*!
     * \brief Called by the coupled problem after
     *        initializing the subproblems / models
     */
    void postInit()
    {}

    /*!
     * \brief The bulk coupling stencil, i.e. which low-dimensional dofs
     *        the given bulk element dof depends on.
     */
    const std::vector<unsigned int>& couplingStencil(const BulkElement& element) const
    {
        const unsigned int eIdx = bulkProblem().elementMapper().index(element);
        if (bulkCouplingStencils_.count(eIdx))
            return bulkCouplingStencils_.at(eIdx);
        else
            return emptyStencil_;
    }

    /*!
     * \brief The low dim coupling stencil, i.e. which bulk dofs
     *        the given low dimensional element dof depends on.
     */
    const std::vector<unsigned int>& couplingStencil(const LowDimElement& element) const
    {
        const unsigned int eIdx = lowDimProblem().elementMapper().index(element);
        if (lowDimCouplingStencils_.count(eIdx))
            return lowDimCouplingStencils_.at(eIdx);
        else
            return emptyStencil_;
    }

    //! evaluate coupling residual for the derivative bulk DOF with respect to low dim DOF
    //! we only need to evaluate the part of the residual that will be influence by the low dim DOF
    BulkPrimaryVariables evalCouplingResidual(const BulkElement& element,
                                              const BulkFVElementGeometry& fvGeometry,
                                              const BulkElementVolumeVariables& curElemVolVars,
                                              const BulkElementBoundaryTypes& elemBcTypes,
                                              const BulkElementFluxVariablesCache& elemFluxVarsCache,
                                              const LowDimElement& lowDimElement)
    {
        const auto bulkElementIdx = bulkProblem().elementMapper().index(element);
        auto&& scv = fvGeometry.scv(bulkElementIdx);
        auto couplingSource = bulkProblem().scvPointSources(element, fvGeometry, curElemVolVars, scv);
        couplingSource *= -scv.volume()*curElemVolVars[scv].extrusionFactor();
        return couplingSource;
    }

    //! evaluate coupling residual for the derivative low dim DOF with respect to bulk DOF
    //! we only need to evaluate the part of the residual that will be influence by the bulk DOF
    LowDimPrimaryVariables evalCouplingResidual(const LowDimElement& element,
                                                const LowDimFVElementGeometry& fvGeometry,
                                                const LowDimElementVolumeVariables& curElemVolVars,
                                                const LowDimElementBoundaryTypes& elemBcTypes,
                                                const LowDimElementFluxVariablesCache& elemFluxVarsCache,
                                                const BulkElement& bulkElement)
    {
        const auto lowDimElementIdx = lowDimProblem().elementMapper().index(element);
        auto&& scv = fvGeometry.scv(lowDimElementIdx);
        auto couplingSource = lowDimProblem().scvPointSources(element, fvGeometry, curElemVolVars, scv);
        couplingSource *= -scv.volume()*curElemVolVars[scv].extrusionFactor();
        return couplingSource;
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
    void computePointSourceData(unsigned int order = 1, bool verbose = false)
    {
        // initilize the maps
        // do some logging and profiling
        Dune::Timer watch;
        std::cout << "Initializing the point sources..." << std::endl;

        // clear all internal members like pointsource vectors and stencils
        // initializes the point source id counter
        clear();

        // a vector for the volumes the low-dim domain occupies in the bulk domain if it were full-dimensional
        lowDimVolumeInBulkElement_.resize(this->bulkGridView().size(0));

        // intersect the bounding box trees
        Dumux::CCMultiDimensionGlue<TypeTag> glue(bulkProblem(), lowDimProblem());
        glue.build();

        for (const auto& is : intersections(glue))
        {
            // all inside elements are identical...
            const auto& inside = is.inside(0);
            // get the intersection geometry for integrating over it
            const auto intersectionGeometry = is.geometry();

            // get the Gaussian quadrature rule for the local intersection
            const auto& quad = Dune::QuadratureRules<Scalar, lowDimDim>::rule(intersectionGeometry.type(), order);
            const unsigned int lowDimElementIdx = lowDimProblem().elementMapper().index(inside);

            // compute the volume the low-dim domain occupies in the bulk domain if it were full-dimensional
            const auto radius = lowDimProblem().spatialParams().radius(lowDimElementIdx);
            for (unsigned int outsideIdx = 0; outsideIdx < is.neighbor(0); ++outsideIdx)
            {
                const auto& outside = is.outside(outsideIdx);
                const unsigned int bulkElementIdx = bulkProblem().elementMapper().index(outside);
                lowDimVolumeInBulkElement_[bulkElementIdx] += intersectionGeometry.volume()*M_PI*radius*radius;
            }

            // iterate over all quadrature points
            for (auto&& qp : quad)
            {
                // compute the coupling stencils
                for (unsigned int outsideIdx = 0; outsideIdx < is.neighbor(0); ++outsideIdx)
                {
                    const auto& outside = is.outside(outsideIdx);
                    const unsigned int bulkElementIdx = bulkProblem().elementMapper().index(outside);

                    // each quadrature point will be a point source for the sub problem
                    const auto globalPos = intersectionGeometry.global(qp.position());
                    const auto id = idCounter_++;
                    const auto qpweight = qp.weight();
                    const auto ie = intersectionGeometry.integrationElement(qp.position());
                    bulkPointSources_.emplace_back(globalPos, id, qpweight, ie, std::vector<unsigned int>({bulkElementIdx}));
                    bulkPointSources_.back().setEmbeddings(is.neighbor(0));
                    lowDimPointSources_.emplace_back(globalPos, id, qpweight, ie, std::vector<unsigned int>({lowDimElementIdx}));
                    lowDimPointSources_.back().setEmbeddings(is.neighbor(0));

                    // pre compute additional data used for the evaluation of
                    // the actual solution dependent source term
                    PointSourceData psData;
                    psData.addLowDimInterpolation(lowDimElementIdx);
                    psData.addBulkInterpolation(bulkElementIdx);

                    // publish point source data in the global vector
                    pointSourceData_.push_back(psData);

                    // compute the coupling stencils
                    bulkCouplingStencils_[bulkElementIdx].push_back(lowDimElementIdx);
                    // add this bulk element to the low dim coupling stencil
                    lowDimCouplingStencils_[lowDimElementIdx].push_back(bulkElementIdx);
                }
            }
        }

        // coupling stencils
        for (auto&& stencil : bulkCouplingStencils_)
        {
            std::sort(stencil.second.begin(), stencil.second.end());
            stencil.second.erase(std::unique(stencil.second.begin(), stencil.second.end()), stencil.second.end());
        }
        // the low dim coupling stencil
        for (auto&& stencil : lowDimCouplingStencils_)
        {
            std::sort(stencil.second.begin(), stencil.second.end());
            stencil.second.erase(std::unique(stencil.second.begin(), stencil.second.end()), stencil.second.end());
        }

        std::cout << "took " << watch.elapsed() << " seconds." << std::endl;
    }

    /*!
     * \brief Methods to be accessed by the subproblems
     */
    // \{

    //! Return a reference to the bulk problem
    Scalar radius(unsigned int id) const
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
        return asImp_().lowDimVolume(element) / totalVolume;
    }

    //! Return a reference to the pointSource data
    const PointSourceData& pointSourceData(unsigned int id) const
    {
        return pointSourceData_[id];
    }

    //! Return a reference to the bulk problem
    const BulkProblem& bulkProblem() const
    {
        return bulkProblem_;
    }

    //! Return a reference to the low dimensional problem
    const LowDimProblem& lowDimProblem() const
    {
        return lowDimProblem_;
    }

    //! Return a reference to the bulk problem
    BulkProblem& bulkProblem()
    {
        return bulkProblem_;
    }

    //! Return a reference to the low dimensional problem
    LowDimProblem& lowDimProblem()
    {
        return lowDimProblem_;
    }

    //! Return a reference to the bulk problem
    const BulkGridView& bulkGridView() const
    {
        return bulkProblem().gridView();
    }

    //! Return a reference to the low dimensional problem
    const LowDimGridView& lowDimGridView() const
    {
        return lowDimProblem().gridView();
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
    BulkPrimaryVariables bulkPriVars(unsigned int id) const
    {
        auto& data = pointSourceData_[id];
        return data.interpolateBulk(bulkProblem().model().curSol());
    }

    //! Return data for a low dim point source with the identifier id
    LowDimPrimaryVariables lowDimPriVars(unsigned int id) const
    {
        auto& data = pointSourceData_[id];
        return data.interpolateLowDim(lowDimProblem().model().curSol());
    }

    //! Compute bulk volume variables for the identifier id
    BulkVolumeVariables bulkVolVars(unsigned int id) const
    {
        // use circle interpolated data and construct volVar object for the interpolated privars
        auto& data = pointSourceData_[id];
        auto bulkPriVars = data.interpolateBulk(bulkProblem().model().curSol());

        const auto element = bulkProblem().model().globalFvGeometry().element(data.bulkElementIdx());
        auto fvGeometry = localView(bulkProblem().model().globalFvGeometry());
        fvGeometry.bindElement(element);

        BulkVolumeVariables volVars;
        volVars.update(BulkElementSolutionVector({bulkPriVars}),
                       bulkProblem(),
                       element,
                       fvGeometry.scv(data.bulkElementIdx()));

        return volVars;
    }

    //! Compute lowDim volume variables for the identifier id
    LowDimVolumeVariables lowDimVolVars(unsigned int id) const
    {
        const auto& data = pointSourceData_[id];
        const auto element = lowDimProblem().model().globalFvGeometry().element(data.lowDimElementIdx());
        auto fvGeometry = localView(lowDimProblem().model().globalFvGeometry());
        fvGeometry.bindElement(element);

        const auto& curSol = lowDimProblem().model().curSol();
        LowDimVolumeVariables volVars;
        volVars.update(lowDimProblem().model().elementSolution(element, curSol),
                       lowDimProblem(),
                       element,
                       fvGeometry.scv(data.lowDimElementIdx()));

        return volVars;
    }
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

    const std::vector<unsigned int>& getAdditionalDofDependencies(unsigned int bulkElementIdx) const
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
    std::unordered_map<unsigned int, std::vector<unsigned int> >& bulkCouplingStencils()
    { return bulkCouplingStencils_; }

    //! Return reference to low dim coupling stencil member
    std::unordered_map<unsigned int, std::vector<unsigned int> >& lowDimCouplingStencils()
    { return lowDimCouplingStencils_; }

    //! Return a reference to an empty stencil
    std::vector<unsigned int>& emptyStencil()
    { return emptyStencil_; }

private:
    //! Returns the implementation of the problem (i.e. static polymorphism)
    Implementation &asImp_()
    { return *static_cast<Implementation *>(this); }

    //! Returns the implementation of the problem (i.e. static polymorphism)
    const Implementation &asImp_() const
    { return *static_cast<const Implementation *>(this); }

    BulkProblem& bulkProblem_;
    LowDimProblem& lowDimProblem_;

    std::vector<BulkPointSource> bulkPointSources_;
    std::vector<LowDimPointSource> lowDimPointSources_;

    mutable std::vector<PointSourceData> pointSourceData_;

    std::unordered_map<unsigned int, std::vector<unsigned int> > bulkCouplingStencils_;
    std::unordered_map<unsigned int, std::vector<unsigned int> > lowDimCouplingStencils_;
    std::vector<unsigned int> emptyStencil_;

    //! vector for the volume fraction of the lowdim domain in the bulk domain cells
    std::vector<Scalar> lowDimVolumeInBulkElement_;

    //! id generator for point sources
    unsigned int idCounter_;
};

} // end namespace Dumux

#endif
