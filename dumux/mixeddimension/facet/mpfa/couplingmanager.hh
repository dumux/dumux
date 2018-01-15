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

    using GlobalPosition = Dune::FieldVector<Scalar, dimWorld>;

    //! The bulk coupling context stores all data required
    //! for the assembly of a bulk element
    struct BulkCouplingContext
    {
        BulkIndexType bulkElementIndex;
        std::vector<LowDimFVElementGeometry> lowDimFvGeometries;
        std::vector<LowDimVolumeVariables> lowDimVolVars;

        BulkCouplingContext(const LowDimProblem& lowDimProblem)
        : lowDimFvGeometries(1, lowDimProblem.model().globalFvGeometry()) {}

        void clear()
        {
            lowDimFvGeometries.clear();
            lowDimVolVars.clear();
        }
    };

    //! The low dim coupling context stores all data required
    //! for the assembly of a bulk element
    struct LowDimCouplingContext
    {
        LowDimIndexType lowDimElementIndex;
        BulkElement bulkElement;
        BulkFVElementGeometry bulkFVGeometry;
        BulkElementVolumeVariables bulkElemVolVars;
        BulkElementFluxVariablesCache bulkElemFluxVarsCache;

        //! The constructor. Sets the global pointers in the local views
        LowDimCouplingContext(const BulkProblem& bulkProblem)
        : bulkFVGeometry(bulkProblem.model().globalFvGeometry()),
          bulkElemVolVars(bulkProblem.model().curGlobalVolVars()),
          bulkElemFluxVarsCache(bulkProblem.model().globalFluxVarsCache()) {}
    };

public:

    /*!
     * \brief Constructor
     */
    CCMpfaFacetCouplingManager(BulkProblem& bulkProblem, LowDimProblem& lowDimProblem)
    : bulkProblemPtr_(&bulkProblem),
      lowDimProblemPtr_(&lowDimProblem),
      couplingMapper_(bulkProblem, lowDimProblem),
      bulkLocalResidual_(),
      lowDimLocalResidual_(),
      bulkCouplingContext_(lowDimProblem),
      lowDimCouplingContext_(bulkProblem)
    {
        // initialize the local residuals
        bulkLocalResidual_.init(bulkProblem);
        lowDimLocalResidual_.init(lowDimProblem);
    }

    //! Called before the initialization of the sub problems
    void preInit()
    {
        couplingMapper_.setInitializationStatus(false);
    }

    //! convenience function to do pre- and postInit simultaneously
    void init()
    {
        preInit();
        postInit();
    }

    //! Called after the sub problems have been initialized
    void postInit()
    {
        // in case of not using no flow tips, we should prepare this here??
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
    { return bulkProblem().model().globalFvGeometry().isOnInteriorBoundary(scvf); }

    //! returns the low dim DOFs a given bulk element is coupled with
    const std::vector<LowDimIndexType>& couplingStencil(const BulkElement& element)
    { return couplingMapper_.getBulkCouplingData(element).couplingStencil; }

    //! returns the bulk DOFs a given low dim element is coupled with
    const std::vector<BulkIndexType>& couplingStencil(const LowDimElement& element)
    { return couplingMapper_.getLowDimCouplingData(element).couplingStencil; }

    //! evaluates the coupling residual of the bulk element with respect to low dim DOF
    //! we only need to evaluate the part of the residual that will be influenced by the low dim DOF
    BulkPrimaryVariables evalCouplingResidual(const BulkElement& element,
                                              const BulkFVElementGeometry& fvGeometry,
                                              const BulkElementVolumeVariables& elemVolVars,
                                              const BulkElementBoundaryTypes& elemBcTypes,
                                              BulkElementFluxVariablesCache& elemFluxVarsCache,
                                              const LowDimElement& lowDimElement)
    {
        // ensure this is only called for the element the context is set for
        assert(bulkCouplingContext_.bulkElementIndex == bulkProblem().elementMapper().index(element) && "Context has not been set to the given element!");

        const auto& couplingData = couplingMapper_.getBulkCouplingData(bulkCouplingContext_.bulkElementIndex);

        // This should never be called for non-coupled elements. But if so, derivative is zero
        if (!couplingData.isCoupled)
            DUNE_THROW(Dune::InvalidStateException, "Bulk element is not coupled");

        // obtain the corresponding volvars in the coupling context (store copy)
        const auto lowDimElementIndex = lowDimProblem().elementMapper().index(lowDimElement);
        const auto idxInContext = findIndexInVector(couplingData.couplingStencil, lowDimElementIndex);
        // const auto origVolVars = bulkCouplingContext_.lowDimVolVars[idxInContext];
        const auto& lowDimCurSol = lowDimProblem().model().curSol();

        // update volvars of the coupling context
        bulkCouplingContext_.lowDimVolVars[idxInContext].update(lowDimProblem().model().elementSolution(lowDimElement, lowDimCurSol),
                                                                lowDimProblem(),
                                                                lowDimElement,
                                                                bulkCouplingContext_.lowDimFvGeometries[idxInContext].scv(lowDimElementIndex));

        // after updating the vol vars, we have to update the flux vars cache
        elemFluxVarsCache.update(element, fvGeometry, elemVolVars);

        // calculate the sum of the fluxes of the scvfs that touch an interior boundary
        static const bool bulkUseTpfaBoundary = GET_PROP_VALUE(BulkProblemTypeTag, UseTpfaBoundary);
        BulkPrimaryVariables flux(0.0);
        for (const auto& scvf : scvfs(fvGeometry))
        {
            // if the scvf does not touch an interior boundary, skip the rest
            if (!bulkProblem().model().globalFvGeometry().touchesInteriorBoundary(scvf))
                continue;

            if (bulkUseTpfaBoundary)
            {
                // do not calculate fluxes on neumann boundary scvfs using tpfa
                if (scvf.boundary())
                {
                    const auto bcTypes = bulkProblem().boundaryTypes(element, scvf);
                    if (!bcTypes.hasOnlyNeumann())
                        flux += bulkLocalResidual_.computeFlux(element, fvGeometry, elemVolVars, scvf, elemFluxVarsCache);
                }
                else
                    flux += bulkLocalResidual_.computeFlux(element, fvGeometry, elemVolVars, scvf, elemFluxVarsCache);
            }
            else
                flux += bulkLocalResidual_.computeFlux(element, fvGeometry, elemVolVars, scvf, elemFluxVarsCache);
        }

        // restore the vol vars
        // bulkCouplingContext_.lowDimVolVars[idxInContext] = origVolVars;

        // return the sum of the fluxes
        return flux;
    }

    //! evaluate coupling residual of the low dim element with respect to bulk DOF
    //! we only need to evaluate the part of the residual that will be influence by the bulk DOF
    LowDimPrimaryVariables evalCouplingResidual(const LowDimElement& element,
                                                const LowDimFVElementGeometry& fvGeometry,
                                                const LowDimElementVolumeVariables& curElemVolVars,
                                                const LowDimElementBoundaryTypes& elemBcTypes,
                                                const LowDimElementFluxVariablesCache& elemFluxVarsCache,
                                                const BulkElement& bulkElement)
    {
        // ensure this is only called for the element the context is set for
        assert(lowDimCouplingContext_.lowDimElementIndex == lowDimProblem().elementMapper().index(element) && "Context has not been set to the given element!");
        const auto& couplingData = couplingMapper_.getLowDimCouplingData(lowDimCouplingContext_.lowDimElementIndex);

        // This should never be called for non-coupled elements. But if so, derivative is zero
        if (!couplingData.isCoupled)
            DUNE_THROW(Dune::InvalidStateException, "low dim element is not coupled");

        // update the corresponding vol vars in the coupling context
        const auto bulkElementIndex = bulkProblem().elementMapper().index(bulkElement);
        const auto& bulkCurSol = bulkProblem().model().curSol();

        auto& bulkVolVars = lowDimCouplingContext_.bulkElemVolVars[bulkElementIndex];
        // const auto origVolVars = bulkVolVars;
        bulkVolVars.update(bulkProblem().model().elementSolution(bulkElement, bulkCurSol),
                           bulkProblem(),
                           bulkElement,
                           lowDimCouplingContext_.bulkFVGeometry.scv(bulkElementIndex));

        // update the elem flux vars cache
        lowDimCouplingContext_.bulkElemFluxVarsCache.update(lowDimCouplingContext_.bulkElement,
                                                            lowDimCouplingContext_.bulkFVGeometry,
                                                            lowDimCouplingContext_.bulkElemVolVars);

        // update the volvars of the bulk coupling context (store original ones for reset)
        const auto& couplingStencil = couplingMapper_.getBulkCouplingData(bulkCouplingContext_.bulkElementIndex).couplingStencil;
        const auto idxInContext = findIndexInVector(couplingStencil, lowDimCouplingContext_.lowDimElementIndex);

        // const auto origLowDimVolVars = bulkCouplingContext_.lowDimVolVars[idxInContext];
        // const auto& lowDimCurSol = lowDimProblem().model().curSol();
        // bulkCouplingContext_.lowDimVolVars[idxInContext].update(lowDimProblem().model().elementSolution(element, lowDimCurSol),
        //                                                         lowDimProblem(),
        //                                                         element,
        //                                                         bulkCouplingContext_.lowDimFvGeometries[idxInContext].scv(lowDimCouplingContext_.lowDimElementIndex));

        // Calculate the sources stemming from the bulk domain.
        // We call the private routine directly, as we know we do not have to update
        // the bulk coupling context, because the low dim solution has not been deflected.
        auto sources = evalSourcesFromBulk_(element, fvGeometry, curElemVolVars, fvGeometry.scv(lowDimCouplingContext_.lowDimElementIndex));
        sources *= -1.0;

        // reset the corresponding volume variables in the coupling context
        // bulkVolVars = origVolVars;

        // restore the low dim vol vars
        // bulkCouplingContext_.lowDimVolVars[idxInContext] = origLowDimVolVars;

        // return sources
        return sources;
    }

    //! evaluates the sources in the low dim domain coming from the bulk domain
    //! updates the respective vol vars in the bulk coupling context in case
    //! the solution has been deflected and the forwards to private routine
    BulkPrimaryVariables evalSourcesFromBulk(const LowDimElement& element,
                                             const LowDimFVElementGeometry& fvGeometry,
                                             const LowDimElementVolumeVariables& curElemVolVars,
                                             const LowDimSubControlVolume& scv) const
    {
        // ensure this is only called for the element the context is set for
        assert(lowDimCouplingContext_.lowDimElementIndex == lowDimProblem().elementMapper().index(element) && "Context has not been set to the given element!");
        const auto& couplingData = couplingMapper_.getLowDimCouplingData(lowDimCouplingContext_.lowDimElementIndex);
        if (!couplingData.isCoupled) return BulkPrimaryVariables(0.0);

        // update the volvars of the coupling context (store original ones for reset)
        const auto& couplingStencil = couplingMapper_.getBulkCouplingData(bulkCouplingContext_.bulkElementIndex).couplingStencil;
        const auto idxInContext = findIndexInVector(couplingStencil, lowDimCouplingContext_.lowDimElementIndex);

        const auto origVolVars = bulkCouplingContext_.lowDimVolVars[idxInContext];
        const auto& lowDimCurSol = lowDimProblem().model().curSol();
        bulkCouplingContext_.lowDimVolVars[idxInContext].update(lowDimProblem().model().elementSolution(element, lowDimCurSol),
                                                                lowDimProblem(),
                                                                element,
                                                                scv);

        // after updating the low dim vol vars, update the bulk flux variables cache
        auto origFluxVarCache = lowDimCouplingContext_.bulkElemFluxVarsCache;
        lowDimCouplingContext_.bulkElemFluxVarsCache.update(lowDimCouplingContext_.bulkElement,
                                                            lowDimCouplingContext_.bulkFVGeometry,
                                                            lowDimCouplingContext_.bulkElemVolVars);

        // let the private method evaluate the actual sources
        const auto flux = evalSourcesFromBulk_(element, fvGeometry, curElemVolVars, scv);

        // restore the vol vars
        bulkCouplingContext_.lowDimVolVars[idxInContext] = origVolVars;
        lowDimCouplingContext_.bulkElemFluxVarsCache = origFluxVarCache;

        return flux;
    }

    //! return the lowDim volume variables (from the coupling context) for given bulk element and scvf
    const LowDimVolumeVariables& lowDimVolVars(const BulkElement& element,
                                               const BulkFVElementGeometry& fvGeometry,
                                               const BulkSubControlVolumeFace& scvf) const
    {
        const auto lowDimElementIndex = couplingMapper_.getBulkCouplingData(element).getCouplingLowDimElementIndex(scvf);
        const auto& contextCouplingStencil = couplingMapper_.getBulkCouplingData(bulkCouplingContext_.bulkElementIndex).couplingStencil;
        return bulkCouplingContext_.lowDimVolVars[findIndexInVector(contextCouplingStencil, lowDimElementIndex)];
    }

    //! return the lowDim FVElementGeometry (from the coupling context) for given bulk element and scvf
    const LowDimFVElementGeometry& lowDimFvGeometry(const BulkElement& element,
                                                    const BulkFVElementGeometry& fvGeometry,
                                                    const BulkSubControlVolumeFace& scvf) const
    {
        const auto lowDimElementIndex = couplingMapper_.getBulkCouplingData(element).getCouplingLowDimElementIndex(scvf);
        const auto& contextCouplingStencil = couplingMapper_.getBulkCouplingData(bulkCouplingContext_.bulkElementIndex).couplingStencil;
        return bulkCouplingContext_.lowDimFvGeometries[findIndexInVector(contextCouplingStencil, lowDimElementIndex)];
    }

    //! Returns the lowDim element's extrusion factor (obtained from volvars in context) for given bulk element and scvf
    Scalar lowDimExtrusionFactor(const BulkElement& element,
                                 const BulkSubControlVolumeFace& scvf) const
    {
        const auto lowDimElementIndex = couplingMapper_.getBulkCouplingData(element).getCouplingLowDimElementIndex(scvf);
        const auto& contextCouplingStencil = couplingMapper_.getBulkCouplingData(bulkCouplingContext_.bulkElementIndex).couplingStencil;
        return bulkCouplingContext_.lowDimVolVars[findIndexInVector(contextCouplingStencil, lowDimElementIndex)].extrusionFactor();
    }

    void setCouplingContext(const LowDimElement& element)
    {
        // first set the index of the low dim element
        lowDimCouplingContext_.lowDimElementIndex = lowDimProblem().elementMapper().index(element);

        // Prepare the local views. We only need to prepare the local view of the first meighboring
        // element as the necessary quantities for all the scvfs will be in there. This is because the
        // scvfs in the other elements are the "flipped" scvfs of the ones in the first element
        const auto& couplingData = couplingMapper_.getLowDimCouplingData(lowDimCouplingContext_.lowDimElementIndex);

        // if low dim element is not coupled, do nothing
        if (couplingData.elementList.size() == 0) return;

        // set the bulk element in the context
        lowDimCouplingContext_.bulkElement = bulkProblem().model().globalFvGeometry().element(couplingData.elementList[0]);

        // set the bulk coupling context to prepare the low dim data required for the neighboring bulk element
        // this data will be required in updating the elem flux var cache below
        setCouplingContext(lowDimCouplingContext_.bulkElement);

        // set the remaining variables in the low dim context
        lowDimCouplingContext_.bulkFVGeometry = localView(bulkProblem().model().globalFvGeometry());
        lowDimCouplingContext_.bulkElemVolVars = localView(bulkProblem().model().curGlobalVolVars());
        lowDimCouplingContext_.bulkElemFluxVarsCache = localView(bulkProblem().model().globalFluxVarsCache());

        lowDimCouplingContext_.bulkFVGeometry.bind(lowDimCouplingContext_.bulkElement);
        lowDimCouplingContext_.bulkElemVolVars.bind(lowDimCouplingContext_.bulkElement,
                                                    lowDimCouplingContext_.bulkFVGeometry,
                                                    bulkProblem().model().curSol());
        lowDimCouplingContext_.bulkElemFluxVarsCache.bind(lowDimCouplingContext_.bulkElement,
                                                          lowDimCouplingContext_.bulkFVGeometry,
                                                          lowDimCouplingContext_.bulkElemVolVars);
    }

    void setCouplingContext(const BulkElement& element)
    {
        bulkCouplingContext_.bulkElementIndex = bulkProblem().elementMapper().index(element);

        // prepare all the vol vars of the low dim elements in the stencil
        const auto& couplingStencil = couplingMapper_.getBulkCouplingData(bulkCouplingContext_.bulkElementIndex).couplingStencil;
        const auto numLowDimElements = couplingStencil.size();

        bulkCouplingContext_.clear();
        bulkCouplingContext_.lowDimFvGeometries.reserve(numLowDimElements);
        bulkCouplingContext_.lowDimVolVars.resize(numLowDimElements);

        for (unsigned int i = 0; i < numLowDimElements; ++i)
        {
            const auto lowDimElementIndex = couplingStencil[i];
            const auto lowDimElement = lowDimProblem().model().globalFvGeometry().element(lowDimElementIndex);
            bulkCouplingContext_.lowDimFvGeometries.push_back(localView(lowDimProblem().model().globalFvGeometry()));
            bulkCouplingContext_.lowDimFvGeometries.back().bindElement(lowDimElement);

            const auto& lowDimCurSol = lowDimProblem().model().curSol();
            bulkCouplingContext_.lowDimVolVars[i].update(lowDimProblem().model().elementSolution(lowDimElement, lowDimCurSol),
                                                         lowDimProblem(),
                                                         lowDimElement,
                                                         bulkCouplingContext_.lowDimFvGeometries[i].scv(lowDimElementIndex));
        }
    }

    const BulkProblem& bulkProblem() const
    { return *bulkProblemPtr_; }

    const LowDimProblem& lowDimProblem() const
    { return *lowDimProblemPtr_; }

    const CCMpfaFacetCouplingMapper<TypeTag>& couplingMapper() const
    { return couplingMapper_; }

private:

    //! evaluates the sources in the low dim domain coming from the bulk domain
    //! i.e. the fluxes from the bulk domain into the given low dim element. We return bulk priVars here,
    //! the low dim problem itself needs to transform that into suitable information
    BulkPrimaryVariables evalSourcesFromBulk_(const LowDimElement& element,
                                              const LowDimFVElementGeometry& fvGeometry,
                                              const LowDimElementVolumeVariables& curElemVolVars,
                                              const LowDimSubControlVolume& scv) const
    {
        const auto& couplingData = couplingMapper_.getLowDimCouplingData(lowDimCouplingContext_.lowDimElementIndex);

        BulkPrimaryVariables flux(0.0);
        for (unsigned int i = 0; i < couplingData.elementList.size(); ++i)
        {
            const auto bulkElement = bulkProblem().model().globalFvGeometry().element(couplingData.elementList[i]);
            for (const auto scvfIdx : couplingData.scvfLists[i])
                flux += bulkLocalResidual_.computeFlux(bulkElement,
                                                       lowDimCouplingContext_.bulkFVGeometry,
                                                       lowDimCouplingContext_.bulkElemVolVars,
                                                       lowDimCouplingContext_.bulkFVGeometry.scvf(scvfIdx),
                                                       lowDimCouplingContext_.bulkElemFluxVarsCache);
        }

        return flux;
    }

    //! returns the local index in a vector for a given global index
    template<typename IdxType1, typename IdxType2>
    unsigned int findIndexInVector(const std::vector<IdxType1>& vector, const IdxType2 globalIdx) const
    {
        auto it = std::find(vector.begin(), vector.end(), globalIdx);
        assert(it != vector.end() && "could not find local index in the vector for the given global index!");
        return std::distance(vector.begin(), it);
    }

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

    CCMpfaFacetCouplingMapper<TypeTag> couplingMapper_;

    BulkLocalResidual bulkLocalResidual_;
    LowDimLocalResidual lowDimLocalResidual_;

    mutable BulkCouplingContext bulkCouplingContext_;
    mutable LowDimCouplingContext lowDimCouplingContext_;
};

} // end namespace

#endif // DUMUX_FACETCOUPLINGMANAGER_HH
