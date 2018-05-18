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
 * \ingroup MixedDimension
 * \ingroup MixedDimensionFacet
 * \copydoc Dumux::CCTpfaFacetCouplingManager
 */

#ifndef DUMUX_CCTPFA_FACETCOUPLING_MANAGER_HH
#define DUMUX_CCTPFA_FACETCOUPLING_MANAGER_HH

#include <algorithm>

#include <dumux/common/properties.hh>
#include <dumux/discretization/elementsolution.hh>
#include <dumux/multidomain/couplingmanager.hh>

#include "couplingmapper.hh"

namespace Dumux {

//! Forward declaration of the manager coupling three domains
template<class MDTraits, class CouplingMapper>
class CCTpfaFacetCouplingThreeDomainManager;

/*!
 * \ingroup MixedDimension
 * \ingroup MixedDimensionFacet
 * \brief Manages the coupling between bulk elements and lower dimensional elements
 *        where the coupling occurs across the facets of the bulk grid. This implementation
 *        is to be used in conjunction with models using the cell-centered tpfa scheme.
 *
 * \tparam MDTraits The multidomain traits containing the types on all sub-domains
 * \tparam CouplingMapper Class containing maps on the coupling between dofs of different grids
 * \tparam idOffset The offset to be used for the domain ids. This is used to specify
 *                  which of the present grids on the hierarchy are to be managed by this
 *                  class. For instance, if a bulk, a facet and an edge grid is present and
 *                  you want to use this coupling manager for the coupling between the facet
 *                  and the edge grid, you should provide an offet of 1.
 */
template<class MDTraits, class CouplingMapper, std::size_t idOffset = 0>
class CCTpfaFacetCouplingManager : public CouplingManager< MDTraits >
{
    using Scalar = typename MDTraits::Scalar;
    using BaseCouplingManager = CouplingManager< MDTraits >;
    using SolutionVector = typename MDTraits::SolutionVector;

    using BulkIdType = typename MDTraits::template DomainIdx<idOffset+0>;
    using LowDimIdType = typename MDTraits::template DomainIdx<idOffset+1>;

    static constexpr auto bulkId = BulkIdType();
    static constexpr auto lowDimId = LowDimIdType();

    // the sub-domain type tags
    template<std::size_t id> using SubDomainTypeTag = typename MDTraits::template SubDomainTypeTag<id>;

    // further types specific to the sub-problems
    template<std::size_t id> using ElementBoundaryTypes = typename GET_PROP_TYPE(SubDomainTypeTag<id>, ElementBoundaryTypes);
    template<std::size_t id> using PrimaryVariables = typename GET_PROP_TYPE(SubDomainTypeTag<id>, PrimaryVariables);
    template<std::size_t id> using Problem = typename GET_PROP_TYPE(SubDomainTypeTag<id>, Problem);
    template<std::size_t id> using NumEqVector = typename GET_PROP_TYPE(SubDomainTypeTag<id>, NumEqVector);
    template<std::size_t id> using LocalResidual = typename GET_PROP_TYPE(SubDomainTypeTag<id>, LocalResidual);

    template<std::size_t id> using FVGridGeometry = typename GET_PROP_TYPE(SubDomainTypeTag<id>, FVGridGeometry);
    template<std::size_t id> using FVElementGeometry = typename FVGridGeometry<id>::LocalView;
    template<std::size_t id> using SubControlVolume = typename FVGridGeometry<id>::SubControlVolume;
    template<std::size_t id> using SubControlVolumeFace = typename FVGridGeometry<id>::SubControlVolumeFace;
    template<std::size_t id> using GridView = typename FVGridGeometry<id>::GridView;
    template<std::size_t id> using Element = typename GridView<id>::template Codim<0>::Entity;
    template<std::size_t id> using IndexType = typename GridView<id>::IndexSet::IndexType;

    template<std::size_t id> using GridVariables = typename GET_PROP_TYPE(SubDomainTypeTag<id>, GridVariables);
    template<std::size_t id> using GridVolumeVariables = typename GridVariables<id>::GridVolumeVariables;
    template<std::size_t id> using ElementVolumeVariables = typename GridVolumeVariables<id>::LocalView;
    template<std::size_t id> using VolumeVariables = typename ElementVolumeVariables<id>::VolumeVariables;
    template<std::size_t id> using GridFluxVariablesCache = typename GridVariables<id>::GridFluxVariablesCache;
    template<std::size_t id> using ElementFluxVariablesCache = typename GridFluxVariablesCache<id>::LocalView;

    //! The coupling context of the bulk domain
    struct BulkCouplingContext
    {
        std::vector< FVElementGeometry<lowDimId> > lowDimFvGeometries;
        std::vector< VolumeVariables<lowDimId> > lowDimVolVars;
        IndexType< bulkId > elementIdx;
        bool isSet;

        void clear()
        {
            lowDimFvGeometries.clear();
            lowDimVolVars.clear();
            isSet = false;
        }
    };

    //! The coupling context of the lowdim domain
    struct LowDimCouplingContext
    {
        //! The local views of one of the neighboring bulk elements
        //! will be enough to calculate all fluxes. We need unique
        //! ptrs because the local views have no default constructor
        std::unique_ptr< FVElementGeometry<bulkId> > bulkFvGeometry;
        std::unique_ptr< ElementVolumeVariables<bulkId> > bulkElemVolVars;
        std::unique_ptr< ElementFluxVariablesCache<bulkId> > bulkElemFluxVarsCache;
        std::unique_ptr< LocalResidual<bulkId> > bulkLocalResidual;
        IndexType< lowDimId > elementIdx;
        bool isSet;

        void clear()
        {
            bulkFvGeometry.reset(nullptr);
            bulkElemVolVars.reset(nullptr);
            bulkElemFluxVarsCache.reset(nullptr);
            isSet = false;
        }
    };

public:

    //! types used for coupling stencils
    //! the poro-mechanical domain elements only couple to themselves on the cc grid
    template<std::size_t id>
    using CouplingStencilType = typename std::conditional< id == bulkId,
                                                           typename CouplingMapper::template Stencil<bulkId>,
                                                           typename CouplingMapper::template Stencil<lowDimId> >::type;

    /*!
     * \brief Methods to be accessed by main
     */
    // \{

    /*!
     * \brief Initialize the coupling manager.
     *
     * \param bulkProblem The problem to be solved on the bulk domain
     * \param lowDimProblem The problem to be solved on the lower-dimensional domain
     * \param curSol The current solution
     * \param couplingMapper The mapper object containing the connectivity between the domains
     */
    void init(std::shared_ptr< Problem<bulkId> > bulkProblem,
              std::shared_ptr< Problem<lowDimId> > lowDimProblem,
              std::shared_ptr< CouplingMapper > couplingMapper,
              const SolutionVector& curSol)
    {
        problemTuple_ = std::make_tuple(bulkProblem, lowDimProblem);
        couplingMapperPtr_ = couplingMapper;
        curSol_ = curSol;

        bulkElemIsCoupled_.resize(bulkProblem->fvGridGeometry().gridView().size(0), false);
        bulkScvfIsCoupled_.resize(bulkProblem->fvGridGeometry().numScvf(), false);

        const auto& bulkMap = couplingMapperPtr_->couplingMap(bulkId, lowDimId);
        for (const auto& entry : bulkMap)
        {
            bulkElemIsCoupled_[entry.first] = true;
            for (const auto& scvfs : entry.second.couplingScvfs)
                bulkScvfIsCoupled_[scvfs[0]] = true;
        }
    }

    /*!
     * \brief Update after the grid has changed.
     */
    void update()
    { DUNE_THROW(Dune::NotImplemented, "Grid updates for models using facet coupling"); }

    // \}

    /*!
     * \brief Methods to be accessed by the assembly
     */
    // \{

    /*!
     * \brief The coupling stencil of the bulk domain with the lowdim domain.
     */
    const typename CouplingMapper::template Stencil<lowDimId>&
    couplingStencil(const Element<bulkId>& element, BulkIdType, LowDimIdType) const
    {
        const auto eIdx = problem<bulkId>().fvGridGeometry().elementMapper().index(element);

        if (bulkElemIsCoupled_[eIdx])
        {
            const auto& map = couplingMapperPtr_->couplingMap(bulkId, lowDimId);
            auto it = map.find(eIdx);
            assert(it != map.end());
            return it->second.couplingStencil;
        }

        return getEmptyStencil(lowDimId);
    }

    /*!
     * \brief The coupling stencil of the lower-dimensional domain with the bulk domain.
     */
    const typename CouplingMapper::template Stencil<bulkId>&
    couplingStencil(const Element<lowDimId>& element, LowDimIdType, BulkIdType) const
    {
        const auto eIdx = problem<lowDimId>().fvGridGeometry().elementMapper().index(element);

        const auto& map = couplingMapperPtr_->couplingMap(lowDimId, bulkId);
        auto it = map.find(eIdx);
        if (it != map.end()) return it->second.couplingStencil;
        else return getEmptyStencil(bulkId);
    }

    /*!
     * \brief updates the current solution.
     */
    void updateSolution(const SolutionVector& sol)
    {
        curSol_ = sol;
    }

    /*!
     * \brief returns true if a bulk scvf coincides with a facet element.
     */
    const bool isCoupled(const Element<bulkId>& element, const SubControlVolumeFace<bulkId>& scvf) const
    { return bulkScvfIsCoupled_[scvf.index()]; }

    /*!
     * \brief returns the low-dim vol vars coinciding with a bulk scvf.
     */
    const VolumeVariables<lowDimId>& getLowDimVolVars(const Element<bulkId>& element,
                                                      const SubControlVolumeFace<bulkId>& scvf) const
    {
        assert(bulkContext_.isSet);
        assert(bulkScvfIsCoupled_[scvf.index()]);
        assert(scvf.insideScvIdx() == problem<bulkId>().fvGridGeometry().elementMapper().index(element));

        const auto& map = couplingMapperPtr_->couplingMap(bulkId, lowDimId);
        const auto& couplingData = map.find(scvf.insideScvIdx())->second;
        const auto lowDimElemIdx = couplingData.getCoupledFacetElementIdx(scvf.index());

        const auto& s = map.find(bulkContext_.elementIdx)->second.couplingStencil;
        const auto& idxInContext = std::distance( s.begin(), std::find(s.begin(), s.end(), lowDimElemIdx) );
        assert(std::find(s.begin(), s.end(), lowDimElemIdx) != s.end());
        return bulkContext_.lowDimVolVars[idxInContext];
    }

    //! evaluate coupling residual for the derivative residual i with respect to privars of dof j
    //! we only need to evaluate the part of the residual that will be influenced by the the privars of dof j
    //! i.e. the source term for the lower-dimensional and the fluxes for the bulk domain.
    //! This implementation here is the overload for the bulk domain.
    typename LocalResidual<bulkId>::ElementResidualVector
    evalCouplingResidual(Dune::index_constant<bulkId> domainI,
                         const Element<bulkId>& elementI,
                         const FVElementGeometry<bulkId>& fvGeometry,
                         const ElementVolumeVariables<bulkId>& elemVolVars,
                         const ElementBoundaryTypes<bulkId>& elemBcTypes,
                         const ElementFluxVariablesCache<bulkId>& elemFluxVarsCache,
                         LocalResidual<bulkId> localResidual,
                         Dune::index_constant<lowDimId> domainJ,
                         IndexType<lowDimId> globalJ)
    {
        typename LocalResidual<bulkId>::ElementResidualVector res(1);
        res = 0.0;
        res[0] = evalFluxToFacetElement_(elementI, fvGeometry, elemVolVars, elemFluxVarsCache, localResidual, globalJ);
        return res;
    }

    //! Overload for the lower-dimensional domain.
    typename LocalResidual<lowDimId>::ElementResidualVector
    evalCouplingResidual(Dune::index_constant<lowDimId> domainI,
                         const Element<lowDimId>& elementI,
                         const FVElementGeometry<lowDimId>& fvGeometry,
                         const ElementVolumeVariables<lowDimId>& elemVolVars,
                         const ElementBoundaryTypes<lowDimId>& elemBcTypes,
                         const ElementFluxVariablesCache<lowDimId>& elemFluxVarsCache,
                         LocalResidual<lowDimId> localResidual,
                         Dune::index_constant<bulkId> domainJ,
                         IndexType<bulkId> globalJ)
    {
        assert(problem<bulkId>().fvGridGeometry().elementMapper().index(elementI) == lowDimContext_.elementIdx);
        const auto elementJ = problem<bulkId>().fvGridGeometry().element(globalJ);
        auto source = evalFluxToFacetElement_(elementJ,
                                              *lowDimContext_.bulkFvGeometry,
                                              *lowDimContext_.bulkElemVolVars,
                                              *lowDimContext_.bulkElemFluxVarsCache,
                                              *lowDimContext_.bulkLocalResidual,
                                              lowDimContext_.elementIdx);
        source *= -1.0;

        // pass to element residual vector and return
        typename LocalResidual<lowDimId>::ElementResidualVector res(1);
        res = 0.0;
        res[0] = source;
        return res;
    }

    //! Computes the sources in a lower-dimensional element stemming from the bulk domain
    NumEqVector<lowDimId> evalSourcesFromBulk(const Element<lowDimId>& element,
                                              const FVElementGeometry<lowDimId>& fvGeometry,
                                              const ElementVolumeVariables<lowDimId>& elemVolVars,
                                              const SubControlVolume<lowDimId>& scv)
    {
        assert(lowDimContext_.isSet);
        NumEqVector<lowDimId> sources(0.0);
        const auto eIdx = problem<lowDimId>().fvGridGeometry().elementMapper().index(element);
        const auto& map = couplingMapperPtr_->couplingMap(lowDimId, bulkId);

        auto it = map.find(eIdx);
        if (it == map.end())
            return sources;

        for (const auto& embedment : it->second.embedments)
            sources += evalFluxToFacetElement_(problem<bulkId>().fvGridGeometry().element(embedment.first),
                                               *lowDimContext_.bulkFvGeometry,
                                               *lowDimContext_.bulkElemVolVars,
                                               *lowDimContext_.bulkElemFluxVarsCache,
                                               *lowDimContext_.bulkLocalResidual,
                                               eIdx);

        return sources;
    }

    /*!
     * \brief Bind the coupling context of the bulk domain
     */
    template< class Assembler >
    void bindCouplingContext(BulkIdType, const Element<bulkId>& element, const Assembler& assembler)
    {
        // clear context
        bulkContext_.clear();

        const auto bulkElemIdx = problem<bulkId>().fvGridGeometry().elementMapper().index(element);

        if (bulkElemIsCoupled_[bulkElemIdx])
        {
            const auto& map = couplingMapperPtr_->couplingMap(bulkId, lowDimId);

            auto it = map.find(bulkElemIdx); assert(it != map.end());
            const auto stencilSize = it->second.couplingStencil.size();
            bulkContext_.lowDimFvGeometries.reserve(stencilSize);
            bulkContext_.lowDimVolVars.reserve(stencilSize);

            for (const auto lowDimIdx : it->second.couplingStencil)
            {
                const auto& ldGridGeometry = problem<lowDimId>().fvGridGeometry();

                const auto elemJ = ldGridGeometry.element(lowDimIdx);
                auto fvGeom = localView(ldGridGeometry);
                fvGeom.bindElement(elemJ);

                const auto elemSol = elementSolution(elemJ, curSol_[lowDimId], ldGridGeometry);
                VolumeVariables<lowDimId> volVars;
                volVars.update(elemSol, problem<lowDimId>(), elemJ, fvGeom.scv(lowDimIdx));

                bulkContext_.isSet = true;
                bulkContext_.elementIdx = bulkElemIdx;
                bulkContext_.lowDimFvGeometries.emplace_back( std::move(fvGeom) );
                bulkContext_.lowDimVolVars.emplace_back( std::move(volVars) );
            }
        }
    }

    /*!
     * \brief Bind the coupling context of the lower-dimensional domain
     * \note Since the low-dim coupling residua are fluxes stemming from
     *       the bulk domain, we have to prepare the bulk coupling context
     *       for the neighboring element (where fluxes are calculated) as well.
     */
    template< class Assembler >
    void bindCouplingContext(LowDimIdType, const Element<lowDimId>& element, const Assembler& assembler)
    {
        // clear contexts
        bulkContext_.clear();
        lowDimContext_.clear();

        const auto lowDimElemIdx = problem<lowDimId>().fvGridGeometry().elementMapper().index(element);
        const auto& map = couplingMapperPtr_->couplingMap(lowDimId, bulkId);
        auto it = map.find(lowDimElemIdx);

        if (it != map.end())
        {
            const auto& bulkGridGeom = problem<bulkId>().fvGridGeometry();
            const auto bulkElem = bulkGridGeom.element(it->second.couplingStencil[0]);

            // first bind the low dim context for the first neighboring bulk element
            bindCouplingContext(bulkId, bulkElem, assembler);

            // then simply bind the local views of that first neighbor
            auto bulkFvGeom = localView(bulkGridGeom);
            auto bulkElemVolVars = localView(assembler.gridVariables(bulkId).curGridVolVars());
            auto bulkElemFluxVarsCache = localView(assembler.gridVariables(bulkId).gridFluxVarsCache());

            bulkFvGeom.bind(bulkElem);
            bulkElemVolVars.bind(bulkElem, bulkFvGeom, curSol_[bulkId]);
            bulkElemFluxVarsCache.bind(bulkElem, bulkFvGeom, bulkElemVolVars);

            lowDimContext_.isSet = true;
            lowDimContext_.elementIdx = lowDimElemIdx;
            lowDimContext_.bulkFvGeometry = std::make_unique< FVElementGeometry<bulkId> >(bulkFvGeom);
            lowDimContext_.bulkElemVolVars = std::make_unique< ElementVolumeVariables<bulkId> >(bulkElemVolVars);
            lowDimContext_.bulkElemFluxVarsCache = std::make_unique< ElementFluxVariablesCache<bulkId> >(bulkElemFluxVarsCache);
            lowDimContext_.bulkLocalResidual = std::make_unique< LocalResidual<bulkId> >(assembler.localResidual(bulkId));
        }
    }

    /*!
     * \brief Update the coupling context for a derivative bulk -> lowDim
     */
    template< class Assembler >
    void updateCouplingContext(BulkIdType,
                               LowDimIdType,
                               IndexType<lowDimId> globalJ,
                               const PrimaryVariables<lowDimId>& priVarsJ,
                               unsigned int pvIdxJ,
                               const Assembler& assembler)
    {
        // do nothing if context is empty
        if (bulkContext_.isSet)
        {
            const auto& lowDimGridGeom = problem<lowDimId>().fvGridGeometry();
            const auto elementJ = lowDimGridGeom.element(globalJ);

            // update the solution (cc->elemsol has only one entry)
            curSol_[lowDimId][globalJ][pvIdxJ] = priVarsJ[pvIdxJ];

            // update vol vars in context
            const auto& map = couplingMapperPtr_->couplingMap(bulkId, lowDimId);
            const auto& couplingStencil = map.find(bulkContext_.elementIdx)->second.couplingStencil;
            auto it = std::find(couplingStencil.begin(), couplingStencil.end(), globalJ);

            assert(it != map.end());
            const auto idxInContext = std::distance(couplingStencil.begin(), it);
            const auto& lowDimScv = bulkContext_.lowDimFvGeometries[idxInContext].scv(globalJ);
            const auto elemSol = elementSolution(elementJ, curSol_[lowDimId], lowDimGridGeom);
            bulkContext_.lowDimVolVars[idxInContext].update(elemSol, problem<lowDimId>(), elementJ, lowDimScv);
        }
    }

    /*!
     * \brief Update the coupling context for a derivative bulk -> bulk.
     */
    template< class Assembler >
    void updateCouplingContext(BulkIdType,
                               BulkIdType,
                               IndexType<bulkId> globalI,
                               const PrimaryVariables<bulkId>& priVarsJ,
                               unsigned int pvIdxJ,
                               const Assembler& assembler)
    { /* do nothing here */ }

    /*!
     * \brief Update the coupling context for a derivative lowDim -> bulk
     */
    template< class Assembler >
    void updateCouplingContext(LowDimIdType,
                               BulkIdType,
                               IndexType<bulkId> globalJ,
                               const PrimaryVariables<bulkId>& priVarsJ,
                               unsigned int pvIdxJ,
                               const Assembler& assembler)
    {
        // do nothing if context is empty
        if (lowDimContext_.isSet)
        {
            const auto& bulkGridGeom = problem<bulkId>().fvGridGeometry();
            const auto elementJ = bulkGridGeom.element(globalJ);

            // update the solution (cc->elemsol has only one entry)
            curSol_[bulkId][globalJ][pvIdxJ] = priVarsJ[pvIdxJ];

            // update corresponding vol vars in context
            const auto& scv = lowDimContext_.bulkFvGeometry->scv(globalJ);
            const auto elemSol = elementSolution(elementJ, curSol_[bulkId], bulkGridGeom);
            (*lowDimContext_.bulkElemVolVars)[globalJ].update(elemSol, problem<bulkId>(), elementJ, scv);

            // update the element flux variables cache (tij depend on low dim values)
            if (globalJ == bulkContext_.elementIdx)
                lowDimContext_.bulkElemFluxVarsCache->update( elementJ, *lowDimContext_.bulkFvGeometry, *lowDimContext_.bulkElemVolVars);
            else
                lowDimContext_.bulkElemFluxVarsCache->update( problem<bulkId>().fvGridGeometry().element(bulkContext_.elementIdx),
                                                              *lowDimContext_.bulkFvGeometry,
                                                              *lowDimContext_.bulkElemVolVars );
        }
    }

    /*!
     * \brief Update the coupling context for a derivative lowDim -> lowDim
     */
    template< class Assembler >
    void updateCouplingContext(LowDimIdType,
                               LowDimIdType,
                               IndexType<lowDimId> globalI,
                               const PrimaryVariables<lowDimId>& priVarsI,
                               unsigned int pvIdxI,
                               const Assembler& assembler)
    {
        // do nothing if context is empty
        if (lowDimContext_.isSet)
        {
            // update the solution (cc->elemsol has only one entry)
            assert(lowDimContext_.elementIdx == globalI);
            curSol_[lowDimId][lowDimContext_.elementIdx][pvIdxI] = priVarsI[pvIdxI];

            // update the corresponding vol vars in the bulk context
            const auto& lowDimGridGeom = problem<lowDimId>().fvGridGeometry();
            const auto& bulkMap = couplingMapperPtr_->couplingMap(bulkId, lowDimId);
            const auto& couplingStencil = bulkMap.find(bulkContext_.elementIdx)->second.couplingStencil;
            auto it = std::find(couplingStencil.begin(), couplingStencil.end(), lowDimContext_.elementIdx);

            assert(it != couplingStencil.end());
            const auto idxInContext = std::distance(couplingStencil.begin(), it);
            const auto& lowDimScv = bulkContext_.lowDimFvGeometries[idxInContext].scv(lowDimContext_.elementIdx);
            const auto lowDimElement = lowDimGridGeom.element(lowDimContext_.elementIdx);
            const auto elemSol = elementSolution(lowDimElement, curSol_[lowDimId], lowDimGridGeom);
            bulkContext_.lowDimVolVars[idxInContext].update(elemSol, problem<lowDimId>(), lowDimElement, lowDimScv);

            // update the element flux variables cache (tij depend on low dim values)
            const auto contextElem = problem<bulkId>().fvGridGeometry().element(bulkContext_.elementIdx);
            lowDimContext_.bulkElemFluxVarsCache->update(contextElem, *lowDimContext_.bulkFvGeometry, *lowDimContext_.bulkElemVolVars);
        }
    }

    /*!
     * \brief Update the local views of the bulk domain after the coupling context changed
     * \note Specialization of the function for deactivated grid-wide caching
     */
    template< class FluxVarsCacheContainer >
    void updateSelf(BulkIdType,
                    const Element<bulkId>& element,
                    const FVElementGeometry<bulkId>& fvGeometry,
                    ElementVolumeVariables<bulkId>& elemVolVars,
                    FluxVarsCacheContainer& fluxVarsCache)
    {
        // update transmissibilities after low dim context has changed
        fluxVarsCache.update(element, fvGeometry, elemVolVars);
    }

    /*!
     * \brief Update the local views of the bulk domain after the coupling context changed
     * \note Specialization of the function for deactivated grid-wide caching
     */
    template< class FluxVarsCacheContainer >
    void updateSelf(BulkIdType,
                    const Element<bulkId>& element,
                    const FVElementGeometry<bulkId>& fvGeometry,
                    GridVolumeVariables<bulkId>& gridVolVars,
                    FluxVarsCacheContainer& fluxVarsCache)
    {
        // update transmissibilities after low dim context has changed
        auto elemVolVars = localView(gridVolVars);
        elemVolVars.bind(element, fvGeometry, curSol_[bulkId]);
        fluxVarsCache.update(element, fvGeometry, elemVolVars);
    }

    /*!
     * \brief In the low dim domain the local views do not depend on coupling data
     */
    template< class VolVarsContainer, class FluxVarsCacheContainer >
    void updateSelf(LowDimIdType,
                    const Element<lowDimId>& element,
                    const FVElementGeometry<lowDimId>& fvGeometry,
                    VolVarsContainer& gridVolVars,
                    FluxVarsCacheContainer& fluxVarsCache)
    {}

    // \}

    /*!
     * \brief Methods to be accessed by the subproblems
     */
    // \{

    //! Return a const reference to one of the problems
    template<std::size_t id, std::enable_if_t<(id == bulkId || id == lowDimId), int> = 0>
    const Problem<id>& problem() const { return *std::get<id-idOffset>(problemTuple_); }

    //! Return a reference to one of the problems
    template<std::size_t id, std::enable_if_t<(id == bulkId || id == lowDimId), int> = 0>
    Problem<id>& problem() { return *std::get<id-idOffset>(problemTuple_); }

    //! We do not have additional dof dependencies here
    template<std::size_t id, std::enable_if_t<(id == bulkId || id == lowDimId), int> = 0>
    std::vector< typename CouplingMapper::template Stencil<id> >
    additionalDofDependencies(Dune::index_constant<id>) const
    { return std::vector< typename CouplingMapper::template Stencil<id> >(); }

    //! We do not have additional dof dependencies here
    template<std::size_t id, std::enable_if_t<(id == bulkId || id == lowDimId), int> = 0>
    const typename CouplingMapper::template Stencil<id>&
    getAdditionalDofDependencies(Dune::index_constant<id>, std::size_t) const
    { return std::get<id-idOffset>(emptyStencilTuple_); }

    //! We do not have additional dof dependencies here
    template<std::size_t id, std::enable_if_t<(id == bulkId || id == lowDimId), int> = 0>
    const typename CouplingMapper::template Stencil<id>&
    getAdditionalDofDependenciesInverse(Dune::index_constant<id>, std::size_t) const
    { return std::get<id-idOffset>(emptyStencilTuple_); }

    //! Empty stencil to be returned for elements that aren't coupled
    template<std::size_t id, std::enable_if_t<(id == bulkId || id == lowDimId), int> = 0>
    const typename CouplingMapper::template Stencil<id>&
    getEmptyStencil(Dune::index_constant<id>) const
    { return std::get<id-idOffset>(emptyStencilTuple_); }

private:
    //! evaluates the bulk-facet exchange fluxes for a given facet element
    NumEqVector<bulkId> evalFluxToFacetElement_(const Element<bulkId>& elementI,
                                                const FVElementGeometry<bulkId>& fvGeometry,
                                                const ElementVolumeVariables<bulkId>& elemVolVars,
                                                const ElementFluxVariablesCache<bulkId>& elemFluxVarsCache,
                                                LocalResidual<bulkId>& localResidual,
                                                IndexType<lowDimId> globalJ)
    {
        const auto bulkElemIdx = problem<bulkId>().fvGridGeometry().elementMapper().index(elementI);

        assert(bulkContext_.isSet);
        assert(bulkElemIsCoupled_[bulkElemIdx]);

        NumEqVector<bulkId> coupledFluxes(0.0);
        const auto& map = couplingMapperPtr_->couplingMap(bulkId, lowDimId);
        const auto& couplingScvfs = map.find(bulkElemIdx)->second.getCoupledScvfs(globalJ);

        for (const auto& scvfIdx : couplingScvfs)
            coupledFluxes += localResidual.evalFlux(problem<bulkId>(),
                                                    elementI,
                                                    fvGeometry,
                                                    elemVolVars,
                                                    elemFluxVarsCache,
                                                    fvGeometry.scvf(scvfIdx));

        return coupledFluxes;
    }

    using BulkProblemPtr = std::shared_ptr< Problem<bulkId> >;
    using LowDimProblemPtr = std::shared_ptr< Problem<lowDimId> >;
    std::tuple<BulkProblemPtr, LowDimProblemPtr> problemTuple_;
    std::shared_ptr<CouplingMapper> couplingMapperPtr_;

    //! store bools for all bulk elements/scvfs that indicate if they
    //! are coupled, so that we don't have to search in the map every time
    std::vector<bool> bulkElemIsCoupled_;
    std::vector<bool> bulkScvfIsCoupled_;

    //! empty stencil to return for non-coupled elements
    using BulkStencil = typename CouplingMapper::template Stencil<bulkId>;
    using LowDimStencil = typename CouplingMapper::template Stencil<lowDimId>;
    std::tuple<BulkStencil, LowDimStencil> emptyStencilTuple_;

    //! Store copy of the current solution
    SolutionVector curSol_;

    //! The coupling contexts of the two domains
    BulkCouplingContext bulkContext_;
    LowDimCouplingContext lowDimContext_;
};

/*!
 * \ingroup MixedDimension
 * \ingroup MixedDimensionFacet
 * \brief Manages the coupling between bulk elements and lower dimensional elements
 *        where the coupling occurs across the facets of the bulk grid. This implementation is
 *        to be used in conjunction with models using the cell-centered tpfa scheme and in problems
 *        where grids and coupling along the hierarchy from 3d cells to 1d edges is considered.
 *
 * \tparam MDTraits The multidomain traits containing the types on all sub-domains
 * \tparam CouplingMapper Class containing maps on the coupling between dofs of different grids
 */
template<class MDTraits, class CouplingMapper>
class CCTpfaFacetCouplingThreeDomainManager
      : public CCTpfaFacetCouplingManager<MDTraits, CouplingMapper, 0>,
        public CCTpfaFacetCouplingManager<MDTraits, CouplingMapper, 1>
{
    using BulkFacetManager = CCTpfaFacetCouplingManager<MDTraits, CouplingMapper, 0>;
    using FacetEdgeManager = CCTpfaFacetCouplingManager<MDTraits, CouplingMapper, 1>;

    using BulkIdType = typename MDTraits::template DomainIdx<0>;
    using FacetIdType = typename MDTraits::template DomainIdx<1>;
    using EdgeIdType = typename MDTraits::template DomainIdx<2>;
    static constexpr auto bulkId = BulkIdType();
    static constexpr auto facetId = FacetIdType();
    static constexpr auto edgeId = EdgeIdType();

    using SolutionVector = typename MDTraits::SolutionVector;

    // the sub-domain type tags
    template<std::size_t id> using SubDomainTypeTag = typename MDTraits::template SubDomainTypeTag<id>;

    // further types specific to the sub-problems
    template<std::size_t id> using ElementBoundaryTypes = typename GET_PROP_TYPE(SubDomainTypeTag<id>, ElementBoundaryTypes);
    template<std::size_t id> using LocalResidual = typename GET_PROP_TYPE(SubDomainTypeTag<id>, LocalResidual);
    template<std::size_t id> using PrimaryVariables = typename GET_PROP_TYPE(SubDomainTypeTag<id>, PrimaryVariables);
    template<std::size_t id> using Problem = typename GET_PROP_TYPE(SubDomainTypeTag<id>, Problem);

    template<std::size_t id> using FVGridGeometry = typename GET_PROP_TYPE(SubDomainTypeTag<id>, FVGridGeometry);
    template<std::size_t id> using FVElementGeometry = typename FVGridGeometry<id>::LocalView;
    template<std::size_t id> using GridView = typename FVGridGeometry<id>::GridView;
    template<std::size_t id> using IndexType = typename GridView<id>::IndexSet::IndexType;
    template<std::size_t id> using Element = typename GridView<id>::template Codim<0>::Entity;

    template<std::size_t id> using GridVariables = typename GET_PROP_TYPE(SubDomainTypeTag<id>, GridVariables);
    template<std::size_t id> using ElementVolumeVariables = typename GridVariables<id>::GridVolumeVariables::LocalView;
    template<std::size_t id> using ElementFluxVariablesCache = typename GridVariables<id>::GridFluxVariablesCache::LocalView;

public:
    /*!
     * \brief Initialize the coupling manager.
     *
     * \param bulkProblem The problem to be solved on the (3d) bulk domain
     * \param facetProblem The problem to be solved on the (2d) facet domain
     * \param edgeProblem The problem to be solved on the (1d) edge domain
     * \tparam curSol The current solution
     * \param couplingMapper The mapper object containing the connectivity between the domains
     */
    void init(std::shared_ptr< Problem<bulkId> > bulkProblem,
              std::shared_ptr< Problem<facetId> > facetProblem,
              std::shared_ptr< Problem<edgeId> > edgeProblem,
              std::shared_ptr< CouplingMapper > couplingMapper,
              const SolutionVector& curSol)
    {
        problemTuple_ = std::make_tuple(bulkProblem, facetProblem, edgeProblem);
        BulkFacetManager::init(bulkProblem, facetProblem, couplingMapper, curSol);
        FacetEdgeManager::init(facetProblem, edgeProblem, couplingMapper, curSol);
    }

    using BulkFacetManager::couplingStencil;
    using FacetEdgeManager::couplingStencil;
    /*!
     * \brief The coupling stencil of the bulk with the edge domain (empty stencil).
     */
    const typename CouplingMapper::template Stencil<edgeId>&
    couplingStencil(const Element<bulkId>& element, BulkIdType, EdgeIdType) const
    { return FacetEdgeManager::getEmptyStencil(edgeId); }

    /*!
     * \brief The coupling stencil of the edge with the bulk domain (empty stencil).
     */
    const typename CouplingMapper::template Stencil<bulkId>&
    couplingStencil(const Element<edgeId>& element, EdgeIdType, BulkIdType) const
    { return BulkFacetManager::getEmptyStencil(bulkId); }

    /*!
     * \brief updates the current solution.
     */
    void updateSolution(const SolutionVector& sol)
    {
        BulkFacetManager::updateSolution(sol);
        FacetEdgeManager::updateSolution(sol);
    }

    using BulkFacetManager::isCoupled;
    using FacetEdgeManager::isCoupled;

    using BulkFacetManager::getLowDimVolVars;
    using FacetEdgeManager::getLowDimVolVars;

    using BulkFacetManager::evalSourcesFromBulk;
    using FacetEdgeManager::evalSourcesFromBulk;

    using BulkFacetManager::evalCouplingResidual;
    using FacetEdgeManager::evalCouplingResidual;

    /*!
     * \brief Interface for evaluating the coupling residual between
     *        the bulk and the edge domain (is always zero).
     */
    template<std::size_t i,
             std::size_t j,
             std::enable_if_t<((i==bulkId && j==edgeId) || ((i==edgeId && j==bulkId))), int> = 0>
    typename LocalResidual<i>::ElementResidualVector
    evalCouplingResidual(Dune::index_constant<i> domainI,
                         const Element<i>& elementI,
                         const FVElementGeometry<i>& fvGeometry,
                         const ElementVolumeVariables<i>& elemVolVars,
                         const ElementBoundaryTypes<i>& elemBcTypes,
                         const ElementFluxVariablesCache<i>& elemFluxVarsCache,
                         LocalResidual<i> localResidual,
                         Dune::index_constant<j> domainJ,
                         IndexType<j> globalJ)
    {
        typename LocalResidual<i>::ElementResidualVector res(1);
        res = 0.0;
        return res;
    }

    using BulkFacetManager::bindCouplingContext;
    using FacetEdgeManager::bindCouplingContext;

    /*!
     * \brief Interface for binding the coupling context for the facet domain. In this case
     *        we have to bind both the facet -> bulk and the facet -> edge coupling context.
     */
    template<class Assembler>
    void bindCouplingContext(FacetIdType, const Element<facetId>& element, const Assembler& assembler)
    {
        BulkFacetManager::bindCouplingContext(facetId, element, assembler);
        FacetEdgeManager::bindCouplingContext(facetId, element, assembler);
    }

    using BulkFacetManager::updateCouplingContext;
    using FacetEdgeManager::updateCouplingContext;

    /*!
     * \brief Interface for updating the coupling context of the facet domain. In this case
     *        we have to update both the facet -> bulk and the facet -> edge coupling context.
     */
    template< class Assembler >
    void updateCouplingContext(FacetIdType,
                               FacetIdType,
                               IndexType<facetId> globalI,
                               const PrimaryVariables<facetId>& priVarsI,
                               unsigned int pvIdxI,
                               const Assembler& assembler)
    {
        BulkFacetManager::updateCouplingContext(facetId, facetId, globalI, priVarsI, pvIdxI, assembler);
        FacetEdgeManager::updateCouplingContext(facetId, facetId, globalI, priVarsI, pvIdxI, assembler);
    }

    /*!
     * \brief Interface for updating the coupling context between
     *        the bulk and the edge domain (do nothing).
     */
    template<std::size_t i,
             std::size_t j,
             class Assembler,
             std::enable_if_t<((i==bulkId && j==edgeId) || (i==edgeId && j==bulkId)), int> = 0>
    void updateCouplingContext(Dune::index_constant<i> domainI,
                               Dune::index_constant<j> domainJ,
                               IndexType<j> globalJ,
                               const PrimaryVariables<j>& priVarsJ,
                               unsigned int pvIdxJ,
                               const Assembler& assembler)
    { /*do nothing here*/ }

    using BulkFacetManager::updateSelf;
    using FacetEdgeManager::updateSelf;

    /*!
     * \brief Interface for updating the local views of the facet domain after updateCouplingContext
     *        the coupling context. In this case we have to forward the both managers as the facet
     *        domain is a part in both.
     */
    template<class VolVarsContainer, class FluxVarsCacheContainer>
    void updateSelf(FacetIdType,
                    const IndexType<facetId>& element,
                    const FVElementGeometry<facetId>& fvGeometry,
                    VolVarsContainer& volVarsContainer,
                    FluxVarsCacheContainer& fluxVarsCacheContainer)
    {
        BulkFacetManager::updateSelf(facetId, element, fvGeometry, volVarsContainer, fluxVarsCacheContainer);
        FacetEdgeManager::updateSelf(facetId, element, fvGeometry, volVarsContainer, fluxVarsCacheContainer);
    }

    //! Return a reference to one of the problems
    template<std::size_t id> const Problem<id>& problem() const { return *std::get<id>(problemTuple_); }
    template<std::size_t id> Problem<id>& problem() { return *std::get<id>(problemTuple_); }

    //! Returns grid-wide stencils with additional dof dependencies
    template<std::size_t id, std::enable_if_t<(id == bulkId || id == facetId), int> = 0>
    std::vector< typename CouplingMapper::template Stencil<id> >
    additionalDofDependencies(Dune::index_constant<id> domainI) const
    { return BulkFacetManager::additionalDofDependencies(domainI); }

    //! Returns grid-wide stencils with additional dof dependencies
    template<std::size_t id, std::enable_if_t<(id == edgeId), int> = 0>
    std::vector< typename CouplingMapper::template Stencil<id> >
    additionalDofDependencies(Dune::index_constant<id> domainI) const
    { return FacetEdgeManager::additionalDofDependencies(domainI); }

    //! Returns stencils with additional dof dependencies
    template<std::size_t id, std::enable_if_t<(id == bulkId || id == facetId), int> = 0>
    const typename CouplingMapper::template Stencil<id>&
    getAdditionalDofDependencies(Dune::index_constant<id> domainI, std::size_t eIdx) const
    { return BulkFacetManager::getAdditionalDofDependencies(domainI, eIdx); }

    //! Returns stencils with additional dof dependencies
    template<std::size_t id, std::enable_if_t<id == edgeId, int> = 0>
    const typename CouplingMapper::template Stencil<id>&
    getAdditionalDofDependencies(Dune::index_constant<id> domainI, std::size_t eIdx) const
    { return FacetEdgeManager::getAdditionalDofDependencies(domainI, eIdx); }

    //! Returns stencils with additional dof dependencies
    template<std::size_t id, std::enable_if_t<(id == bulkId || id == facetId), int> = 0>
    const typename CouplingMapper::template Stencil<id>&
    getAdditionalDofDependenciesInverse(Dune::index_constant<id> domainI, std::size_t eIdx) const
    { return BulkFacetManager::getAdditionalDofDependenciesInverse(domainI, eIdx); }

    //! Returns stencils with additional dof dependencies
    template<std::size_t id, std::enable_if_t<id == edgeId, int> = 0>
    const typename CouplingMapper::template Stencil<id>&
    getAdditionalDofDependenciesInverse(Dune::index_constant<id> domainI, std::size_t eIdx) const
    { return FacetEdgeManager::getAdditionalDofDependenciesInverse(domainI, eIdx); }

private:
    using BulkProblemPtr = std::shared_ptr< Problem<bulkId> >;
    using FacetProblemPtr = std::shared_ptr< Problem<facetId> >;
    using EdgeProblemPtr = std::shared_ptr< Problem<edgeId> >;
    std::tuple<BulkProblemPtr, FacetProblemPtr, EdgeProblemPtr> problemTuple_;
};

} // end namespace Dumux

#endif
