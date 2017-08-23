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
 *   along with this program. If not, see <http://www.gnu.org/licenses/>.    *
 *****************************************************************************/
/*!
 * \file
 * \ingroup FacetCoupling
 * \brief Manages the coupling between bulk elements and lower dimensional elements on the element facets
 */

#ifndef DUMUX_MIXEDDIMENSION_CCMPFA_FACET_COUPLINGMANAGER_HH
#define DUMUX_MIXEDDIMENSION_CCMPFA_FACET_COUPLINGMANAGER_HH

#include <dumux/implicit/properties.hh>
#include <dumux/mixeddimension/facet/mpfa/couplingmapper.hh>

namespace Dumux
{

/*!
 * \ingroup FacetCoupling
 * \brief Manages the interaction between the bulk and lower dimensional (facet) elements
 */
template<class TypeTag>
class CCMpfaFacetCouplingManager
{
    using Implementation = typename GET_PROP_TYPE(TypeTag, CouplingManager);
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using GlobalProblem = typename GET_PROP_TYPE(TypeTag, Problem);
    using GridCreator = typename GET_PROP_TYPE(TypeTag, GridCreator);

    using BulkProblemTypeTag = typename GET_PROP_TYPE(TypeTag, BulkProblemTypeTag);
    using LowDimProblemTypeTag = typename GET_PROP_TYPE(TypeTag, LowDimProblemTypeTag);

    using BulkGridView = typename GET_PROP_TYPE(BulkProblemTypeTag, GridView);
    using BulkElement = typename BulkGridView::template Codim<0>::Entity;
    using BulkIntersection = typename BulkGridView::Intersection;
    using BulkIndexType = typename BulkGridView::IndexSet::IndexType;
    using BulkProblem = typename GET_PROP_TYPE(BulkProblemTypeTag, Problem);
    using BulkMpfaHelper = typename GET_PROP_TYPE(BulkProblemTypeTag, MpfaHelper);
    using BulkLocalResidual = typename GET_PROP_TYPE(BulkProblemTypeTag, LocalResidual);
    using BulkPrimaryVariables = typename GET_PROP_TYPE(BulkProblemTypeTag, PrimaryVariables);
    using BulkFVElementGeometry = typename GET_PROP_TYPE(BulkProblemTypeTag, FVElementGeometry);
    using BulkSubControlVolumeFace = typename GET_PROP_TYPE(BulkProblemTypeTag, SubControlVolumeFace);
    using BulkElementBoundaryTypes = typename GET_PROP_TYPE(BulkProblemTypeTag, ElementBoundaryTypes);
    using BulkElementVolumeVariables = typename GET_PROP_TYPE(BulkProblemTypeTag, ElementVolumeVariables);
    using BulkElementFluxVariablesCache = typename GET_PROP_TYPE(BulkProblemTypeTag, ElementFluxVariablesCache);

    using LowDimGridView = typename GET_PROP_TYPE(LowDimProblemTypeTag, GridView);
    using LowDimElement = typename LowDimGridView::template Codim<0>::Entity;
    using LowDimIndexType = typename LowDimGridView::IndexSet::IndexType;
    using LowDimProblem = typename GET_PROP_TYPE(LowDimProblemTypeTag, Problem);
    using LowDimLocalResidual = typename GET_PROP_TYPE(LowDimProblemTypeTag, LocalResidual);
    using LowDimPrimaryVariables = typename GET_PROP_TYPE(LowDimProblemTypeTag, PrimaryVariables);
    using LowDimSubControlVolume = typename GET_PROP_TYPE(LowDimProblemTypeTag, SubControlVolume);
    using LowDimFVElementGeometry = typename GET_PROP_TYPE(LowDimProblemTypeTag, FVElementGeometry);
    using LowDimElementBoundaryTypes = typename GET_PROP_TYPE(LowDimProblemTypeTag, ElementBoundaryTypes);
    using LowDimVolumeVariables = typename GET_PROP_TYPE(LowDimProblemTypeTag, VolumeVariables);
    using LowDimElementVolumeVariables = typename GET_PROP_TYPE(LowDimProblemTypeTag, ElementVolumeVariables);
    using LowDimElementFluxVariablesCache = typename GET_PROP_TYPE(LowDimProblemTypeTag, ElementFluxVariablesCache);

    static constexpr int dimWorld = BulkGridView::dimensionworld;
    static constexpr int lowDimDimWorld = LowDimGridView::dimensionworld;
    static_assert(dimWorld == lowDimDimWorld, "dimWorld cannot be different for the two sub-domains!");

    using GlobalPosition = Dune::FieldVector<Scalar, dimWorld>;

public:

    /*!
     * \brief Constructor
     */
    CCMpfaFacetCouplingManager(BulkProblem& bulkProblem, LowDimProblem& lowDimProblem)
    : bulkProblemPtr_(&bulkProblem),
      lowDimProblemPtr_(&lowDimProblem),
      couplingMapper_(bulkProblem, lowDimProblem)
    {
        // initialize the local residuals
        bulkLocalResidual_.init(bulkProblem);
        lowDimLocalResidual_.init(lowDimProblem);
    }

    void preInit() {}

    void postInit()
    {
        // abfrage ob no flow tip??
        couplingMapper_.init();
    }

    //! evaluates if an intersection is on an interior boundary
    bool isInteriorBoundary(const BulkElement& element, const BulkIntersection& is) const
    {
        // TODO: Facet elements that are on the bulk domain boundary
        if (is.boundary())
            return false;

        const auto insideIdx = bulkProblem().elementMapper().index(element);
        const auto outsideIdx = bulkProblem().elementMapper().index(is.outside());
        return isInteriorBoundary_(insideIdx, outsideIdx);
    }

    //! evaluates if an scvf is on an interior boundary
    bool isInteriorBoundary(const BulkElement& element, const BulkSubControlVolumeFace& scvf) const
    {
        // TODO: Facet elements that are on the bulk domain boundary
        if (scvf.boundary())
            return false;

        if (couplingMapper_.isInitialized())
        {
            const auto& couplingData = couplingMapper_.getBulkCouplingData(element);

            // if the element is coupled, check if the scvf is so as well
            if (!couplingData.isCoupled)
                return false;
            return couplingData.getScvfCouplingData(scvf).first;
        }

        // If the mapper has not been initialized yet, use the private function
        return isInteriorBoundary_(scvf.insideScvIdx(), scvf.outsideScvIdx());
    }

    const std::vector<LowDimIndexType>& couplingStencil(const BulkElement& element)
    { return couplingMapper_.getBulkCouplingData(element).couplingStencil; }

    const std::vector<BulkIndexType>& couplingStencil(const LowDimElement& element)
    { return couplingMapper_.getLowDimCouplingData(element).couplingStencil; }

    //! evaluate coupling residual for the derivative bulk DOF with respect to low dim DOF
    //! we only need to evaluate the part of the residual that will be influence by the low dim DOF
    BulkPrimaryVariables evalCouplingResidual(const BulkElement& element,
                                              const BulkFVElementGeometry& fvGeometry,
                                              const BulkElementVolumeVariables& curElemVolVars,
                                              const BulkElementBoundaryTypes& elemBcTypes,
                                              const BulkElementFluxVariablesCache& elemFluxVarsCache,
                                              const LowDimElement& lowDimElement)
    {
        const auto& couplingData = couplingMapper_.getBulkCouplingData(element);

        if (!couplingData.isCoupled)
            return BulkPrimaryVariables(0.0);

        // calculate the local residual of the bulk element
        auto prevElemVolVars = localView(bulkProblem().model().prevGlobalVolVars());
        prevElemVolVars.bindElement(element, fvGeometry, bulkProblem().model().prevSol());
        bulkLocalResidual_.eval(element, fvGeometry, prevElemVolVars, curElemVolVars, elemBcTypes, elemFluxVarsCache);
        return bulkLocalResidual_.residual(0);
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
        const auto& couplingData = couplingMapper_.getLowDimCouplingData(element);

        if (!couplingData.isCoupled)
            return LowDimPrimaryVariables(0.0);

        // calculate the source term of the low dim element
        const auto eIdx = lowDimProblem().elementMapper().index(element);
        const auto& scv = fvGeometry.scv(eIdx);

        auto source = lowDimLocalResidual_.computeSource(element, fvGeometry, curElemVolVars, scv);
        source *= -scv.volume()*curElemVolVars[scv].extrusionFactor();
        return source;
    }

    //! evaluates the sources in the low dim domain coming from the bulk domain
    //! i.e. the fluxes from the bulk domain into the given low dim element. We return bulk priVars here,
    //! the low dim problem itself needs to transform that into suitable information
    BulkPrimaryVariables evalSourcesFromBulk(const LowDimElement& element,
                                             const LowDimFVElementGeometry& fvGeometry,
                                             const LowDimElementVolumeVariables& curElemVolVars,
                                             const LowDimSubControlVolume& scv) const
    {
        const auto& couplingData = couplingMapper_.getLowDimCouplingData(element);

        if (!couplingData.isCoupled)
            return BulkPrimaryVariables(0.0);

        // evaluate the fluxes over each scvf in each bulk neighbor. Prepare the local views
        // of the first element, the necessary quantities for all the scvfs will be in there,
        // as the scvfs in the other elements are the "flipped" scvfs of the ones in the first element
        const auto& firstPair = couplingData.elementScvfList[0];
        const auto firstBulkElement = bulkProblem().model().fvGridGeometry().element(firstPair.first);

        auto bulkFvGeometry = localView(bulkProblem().model().fvGridGeometry());
        bulkFvGeometry.bind(firstBulkElement);

        auto bulkElemVolVars = localView(bulkProblem().model().curGlobalVolVars());
        bulkElemVolVars.bind(firstBulkElement, bulkFvGeometry, bulkProblem().model().curSol());

        auto bulkElemFluxVarsCache = localView(bulkProblem().model().globalFluxVarsCache());
        bulkElemFluxVarsCache.bind(firstBulkElement, bulkFvGeometry, bulkElemVolVars);

        BulkPrimaryVariables flux(0.0);
        for (const auto& pair : couplingData.elementScvfList)
        {
            const auto bulkElement = bulkProblem().model().fvGridGeometry().element(pair.first);
            for (auto scvfIdx : pair.second)
                flux += bulkLocalResidual_.computeFlux(bulkElement,
                                                       bulkFvGeometry,
                                                       bulkElemVolVars,
                                                       bulkFvGeometry.scvf(scvfIdx),
                                                       bulkElemFluxVarsCache);
        }

        // The source has to be formulated per volume
        flux /= scv.volume()*curElemVolVars[scv].extrusionFactor();
        return flux;
    }

    //! Compute lowDim volume variables for given bulk element and scvf
    LowDimVolumeVariables lowDimVolVars(const BulkElement& element,
                                        const BulkFVElementGeometry& fvGeometry,
                                        const BulkSubControlVolumeFace& scvf) const
    {
        const auto& scvfCouplingData = couplingMapper_.getBulkCouplingData(element).getScvfCouplingData(scvf);

        assert(scvfCouplingData.first && "no coupled facet element found for given scvf!");
        const auto lowDimElement = lowDimProblem().model().fvGridGeometry().element(scvfCouplingData.second);
        auto lowDimFvGeometry = localView(lowDimProblem().model().fvGridGeometry());
        lowDimFvGeometry.bindElement(lowDimElement);

        const auto& lowDimCurSol = lowDimProblem().model().curSol();
        LowDimVolumeVariables lowDimVolVars;
        lowDimVolVars.update(lowDimProblem().model().elementSolution(lowDimElement, lowDimCurSol),
                             lowDimProblem(),
                             lowDimElement,
                             lowDimFvGeometry.scv(scvfCouplingData.second));
        return lowDimVolVars;
    }

    const BulkProblem& bulkProblem() const
    { return *bulkProblemPtr_; }

    const LowDimProblem& lowDimProblem() const
    { return *lowDimProblemPtr_; }

    const CCMpfaFacetCouplingMapper<TypeTag>& couplingMapper() const
    { return couplingMapper_; }

private:

    //! evaluates if there is an interior boundary between two bulk elements
    bool isInteriorBoundary_(BulkIndexType eIdxInside, BulkIndexType eIdxOutside) const
    {
        for (const auto& lowDimElement : elements(lowDimProblem().gridView()))
        {
            // The indices of bulk elements in which the actual low dim element is a facet
            const auto& bulkElementIndices = GridCreator::getCoupledBulkElementIndices(lowDimElement);

            // if inside & outside element are contained in them, this is a "coupling" intersection
            if (BulkMpfaHelper::contains(bulkElementIndices, eIdxInside) && BulkMpfaHelper::contains(bulkElementIndices, eIdxOutside))
                return true;
        }

        return false;
    }

    BulkProblem& bulkProblem()
    { return *bulkProblemPtr_; }

    LowDimProblem& lowDimProblem()
    { return *lowDimProblemPtr_; }

    //! Returns the implementation of the problem (i.e. static polymorphism)
    Implementation &asImp_()
    { return *static_cast<Implementation *>(this); }

    //! \copydoc asImp_()
    const Implementation &asImp_() const
    { return *static_cast<const Implementation *>(this); }

    BulkProblem *bulkProblemPtr_;
    LowDimProblem *lowDimProblemPtr_;

    mutable BulkLocalResidual bulkLocalResidual_;
    mutable LowDimLocalResidual lowDimLocalResidual_;

    CCMpfaFacetCouplingMapper<TypeTag> couplingMapper_;
};

} // end namespace

#endif // DUMUX_FACETCOUPLINGMANAGER_HH
