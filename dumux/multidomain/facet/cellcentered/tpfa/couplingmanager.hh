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
 * \ingroup FacetCoupling
 * \copydoc Dumux::FacetCouplingManager
 */
#ifndef DUMUX_CCTPFA_FACETCOUPLING_MANAGER_HH
#define DUMUX_CCTPFA_FACETCOUPLING_MANAGER_HH

#include <algorithm>
#include <cassert>

#include <dumux/common/properties.hh>
#include <dumux/common/indextraits.hh>
#include <dumux/discretization/method.hh>
#include <dumux/discretization/elementsolution.hh>
#include <dumux/multidomain/couplingmanager.hh>
#include <dumux/multidomain/facet/couplingmanager.hh>

namespace Dumux {

/*!
 * \ingroup FacetCoupling
 * \brief Manages the coupling between bulk elements and lower dimensional elements
 *        where the coupling occurs across the facets of the bulk grid. This implementation
 *        is to be used in conjunction with models using the cell-centered tpfa scheme.
 *
 * \tparam MDTraits The multidomain traits containing the types on all sub-domains
 * \tparam CouplingMapper Class containing maps on the coupling between dofs of different grids
 * \tparam bulkDomainId The domain id of the bulk problem
 * \tparam lowDimDomainId The domain id of the lower-dimensional problem
 */
template<class MDTraits, class CouplingMapper, std::size_t bulkDomainId, std::size_t lowDimDomainId>
class FacetCouplingManager<MDTraits, CouplingMapper, bulkDomainId, lowDimDomainId, DiscretizationMethod::cctpfa>
: public virtual CouplingManager< MDTraits >
{
    using ParentType = CouplingManager< MDTraits >;

    // convenience aliases and instances of the two domain ids
    using BulkIdType = typename MDTraits::template SubDomain<bulkDomainId>::Index;
    using LowDimIdType = typename MDTraits::template SubDomain<lowDimDomainId>::Index;
    static constexpr auto bulkId = BulkIdType();
    static constexpr auto lowDimId = LowDimIdType();

    // the sub-domain type tags
    template<std::size_t id> using SubDomainTypeTag = typename MDTraits::template SubDomain<id>::TypeTag;

    // further types specific to the sub-problems
    template<std::size_t id> using PrimaryVariables = GetPropType<SubDomainTypeTag<id>, Properties::PrimaryVariables>;
    template<std::size_t id> using Problem = GetPropType<SubDomainTypeTag<id>, Properties::Problem>;
    template<std::size_t id> using NumEqVector = GetPropType<SubDomainTypeTag<id>, Properties::NumEqVector>;
    template<std::size_t id> using LocalResidual = GetPropType<SubDomainTypeTag<id>, Properties::LocalResidual>;

    template<std::size_t id> using GridGeometry = GetPropType<SubDomainTypeTag<id>, Properties::GridGeometry>;
    template<std::size_t id> using FVElementGeometry = typename GridGeometry<id>::LocalView;
    template<std::size_t id> using SubControlVolume = typename GridGeometry<id>::SubControlVolume;
    template<std::size_t id> using SubControlVolumeFace = typename GridGeometry<id>::SubControlVolumeFace;
    template<std::size_t id> using GridView = typename GridGeometry<id>::GridView;
    template<std::size_t id> using Element = typename GridView<id>::template Codim<0>::Entity;
    template<std::size_t id> using GridIndexType = typename IndexTraits< GridView<id> >::GridIndex;

    template<std::size_t id> using GridVariables = GetPropType<SubDomainTypeTag<id>, Properties::GridVariables>;
    template<std::size_t id> using GridVolumeVariables = typename GridVariables<id>::GridVolumeVariables;
    template<std::size_t id> using ElementVolumeVariables = typename GridVolumeVariables<id>::LocalView;
    template<std::size_t id> using VolumeVariables = typename ElementVolumeVariables<id>::VolumeVariables;
    template<std::size_t id> using GridFluxVariablesCache = typename GridVariables<id>::GridFluxVariablesCache;
    template<std::size_t id> using ElementFluxVariablesCache = typename GridFluxVariablesCache<id>::LocalView;

    // this currently does not work for some grid-wide caches being active
    static_assert(!getPropValue<SubDomainTypeTag<bulkId>, Properties::EnableGridFluxVariablesCache>(),
                  "Grid flux variables caching currently not supported in the bulk domain of cc-facet coupling models");
    static_assert(!getPropValue<SubDomainTypeTag<lowDimId>, Properties::EnableGridVolumeVariablesCache>(),
                  "Grid volume variables caching currently not supported in the lower-dimensional domain of cc-facet coupling models");
    static_assert(!getPropValue<SubDomainTypeTag<bulkId>, Properties::EnableGridVolumeVariablesCache>(),
                  "Grid volume variables caching currently not supported in the bulk domain of cc-facet coupling models");

    // extract corresponding grid ids from the mapper
    static constexpr int bulkDim = GridView<bulkDomainId>::dimension;
    static constexpr int lowDimDim = GridView<lowDimDomainId>::dimension;
    static constexpr auto bulkGridId = CouplingMapper::template gridId<bulkDim>();
    static constexpr auto lowDimGridId = CouplingMapper::template gridId<lowDimDim>();

    static constexpr bool lowDimUsesBox = GridGeometry<lowDimId>::discMethod == DiscretizationMethod::box;

    /*!
     * \brief The coupling context of the bulk domain. Contains all data of the lower-
     *        dimensional domain which is required for the computation of a bulk element
     *        residual. This boils down to the geometries and volume variables of all
     *        lower-dimensional elements connected to a given bulk element.
     */
    struct BulkCouplingContext
    {
        bool isSet;
        GridIndexType< bulkId > elementIdx;
        std::vector< FVElementGeometry<lowDimId> > lowDimFvGeometries;
        std::vector< VolumeVariables<lowDimId> > lowDimVolVars;

        void reset()
        {
            lowDimFvGeometries.clear();
            lowDimVolVars.clear();
            isSet = false;
        }
    };

    /*!
     * \brief The coupling context of the lower-dimensional (codim 1) domain. Contains
     *        all data of the bulk domain which is required for computation of element
     *        residuals in the lower-dimensional domain. This is essentially everything
     *        that is necessary to compute the fluxes in bulk domain entering a given
     *        lower-dimensional element. Thus, we store and bind the local views of one
     *        of the neighboring elements, which will be enough to compute the fluxes
     *        stemming from all neighboring elements.
     *
     * \note We need unique ptrs here because the local views have no default constructor.
     */
    struct LowDimCouplingContext
    {
        bool isSet;
        GridIndexType< lowDimId > elementIdx;
        std::unique_ptr< FVElementGeometry<bulkId> > bulkFvGeometry;
        std::unique_ptr< ElementVolumeVariables<bulkId> > bulkElemVolVars;
        std::unique_ptr< ElementFluxVariablesCache<bulkId> > bulkElemFluxVarsCache;
        std::unique_ptr< LocalResidual<bulkId> > bulkLocalResidual;

        void reset()
        {
            bulkFvGeometry.reset(nullptr);
            bulkElemVolVars.reset(nullptr);
            bulkElemFluxVarsCache.reset(nullptr);
            isSet = false;
        }
    };

public:

    //! types used for coupling stencils
    template<std::size_t i, std::size_t j = (i == bulkId) ? lowDimId : bulkId>
    using CouplingStencilType = typename CouplingMapper::template Stencil< CouplingMapper::template gridId<GridView<j>::dimension>() >;

    //! the type of the solution vector
    using SolutionVector = typename MDTraits::SolutionVector;

    /*!
     * \brief Initialize the coupling manager.
     *
     * \param bulkProblem The problem to be solved on the bulk domain
     * \param lowDimProblem The problem to be solved on the lower-dimensional domain
     * \param couplingMapper The mapper object containing the connectivity between the domains
     * \param curSol The current solution
     */
    void init(std::shared_ptr< Problem<bulkId> > bulkProblem,
              std::shared_ptr< Problem<lowDimId> > lowDimProblem,
              std::shared_ptr< CouplingMapper > couplingMapper,
              const SolutionVector& curSol)
    {
        couplingMapperPtr_ = couplingMapper;

        // set the sub problems
        this->setSubProblem(bulkProblem, bulkId);
        this->setSubProblem(lowDimProblem, lowDimId);

        // copy the solution vector
        ParentType::updateSolution(curSol);

        // determine all bulk elements/scvfs that couple to low dim elements
        bulkElemIsCoupled_.assign(bulkProblem->gridGeometry().gridView().size(0), false);
        bulkScvfIsCoupled_.assign(bulkProblem->gridGeometry().numScvf(), false);

        const auto& bulkMap = couplingMapperPtr_->couplingMap(bulkGridId, lowDimGridId);
        for (const auto& entry : bulkMap)
        {
            bulkElemIsCoupled_[entry.first] = true;
            for (const auto& couplingEntry : entry.second.dofToCouplingScvfMap)
                for (const auto& scvfIdx : couplingEntry.second)
                    bulkScvfIsCoupled_[scvfIdx] = true;
        }
    }

    /*!
     * \brief The coupling stencil of a given bulk domain element.
     */
    const CouplingStencilType<bulkId>& couplingStencil(BulkIdType domainI,
                                                       const Element<bulkId>& element,
                                                       LowDimIdType domainJ) const
    {
        const auto eIdx = this->problem(bulkId).gridGeometry().elementMapper().index(element);

        if (bulkElemIsCoupled_[eIdx])
        {
            const auto& map = couplingMapperPtr_->couplingMap(bulkGridId, lowDimGridId);
            auto it = map.find(eIdx);
            assert(it != map.end());
            return it->second.couplingStencil;
        }

        return getEmptyStencil(lowDimId);
    }

    /*!
     * \brief The coupling stencil of the lower-dimensional domain with the bulk domain.
     */
    const CouplingStencilType<lowDimId>& couplingStencil(LowDimIdType domainI,
                                                         const Element<lowDimId>& element,
                                                         BulkIdType domainJ) const
    {
        const auto eIdx = this->problem(lowDimId).gridGeometry().elementMapper().index(element);

        const auto& map = couplingMapperPtr_->couplingMap(lowDimGridId, bulkGridId);
        auto it = map.find(eIdx);
        if (it != map.end()) return it->second.couplingStencil;
        else return getEmptyStencil(bulkId);
    }

    /*!
     * \brief returns true if a bulk scvf flux depends on data in the facet domain.
     */
    bool isCoupled(const Element<bulkId>& element,
                   const SubControlVolumeFace<bulkId>& scvf) const
    { return bulkScvfIsCoupled_[scvf.index()]; }

    /*!
     * \brief returns true if a bulk scvf coincides with a facet element.
     * \note for tpfa, this is always true for coupled scvfs
     */
    bool isOnInteriorBoundary(const Element<bulkId>& element,
                              const SubControlVolumeFace<bulkId>& scvf) const
    { return isCoupled(element, scvf); }

    /*!
     * \brief returns the vol vars of a lower-dimensional element coinciding with a bulk scvf.
     */
    const VolumeVariables<lowDimId>& getLowDimVolVars(const Element<bulkId>& element,
                                                      const SubControlVolumeFace<bulkId>& scvf) const
    {
        assert(bulkContext_.isSet);

        const auto lowDimElemIdx = getLowDimElementIndex(element, scvf);
        const auto& map = couplingMapperPtr_->couplingMap(bulkGridId, lowDimGridId);
        const auto& s = map.find(bulkContext_.elementIdx)->second.couplingElementStencil;
        const auto& idxInContext = std::distance( s.begin(), std::find(s.begin(), s.end(), lowDimElemIdx) );

        assert(std::find(s.begin(), s.end(), lowDimElemIdx) != s.end());
        return bulkContext_.lowDimVolVars[idxInContext];
    }

    /*!
     * \brief returns the lower-dimensional element coinciding with a bulk scvf.
     */
    const Element<lowDimId> getLowDimElement(const Element<bulkId>& element,
                                             const SubControlVolumeFace<bulkId>& scvf) const
    {
        const auto lowDimElemIdx = getLowDimElementIndex(element, scvf);
        return this->problem(lowDimId).gridGeometry().element(lowDimElemIdx);
    }

    /*!
     * \brief returns the index of the lower-dimensional element coinciding with a bulk scvf.
     */
    const GridIndexType<lowDimId> getLowDimElementIndex(const Element<bulkId>& element,
                                                        const SubControlVolumeFace<bulkId>& scvf) const
    {
        assert(bulkScvfIsCoupled_[scvf.index()]);

        const auto& map = couplingMapperPtr_->couplingMap(bulkGridId, lowDimGridId);
        const auto& couplingData = map.at(scvf.insideScvIdx());

        // search the low dim element idx this scvf is embedded in
        auto it = std::find_if( couplingData.elementToScvfMap.begin(),
                                couplingData.elementToScvfMap.end(),
                                [&scvf] (auto& dataPair)
                                {
                                    const auto& scvfs = dataPair.second;
                                    return std::find(scvfs.begin(), scvfs.end(), scvf.index()) != scvfs.end();
                                } );

        assert(it != couplingData.elementToScvfMap.end());
        return it->first;
    }

    /*!
     * \brief Evaluates the coupling element residual of a bulk domain element with respect
     *        to a dof in the lower-dimensional domain (dofIdxGlobalJ). This is essentially
     *        the fluxes across the bulk element facets that depend on dofIdxGlobalJ.
     */
    template< class BulkLocalAssembler >
    typename LocalResidual<bulkId>::ElementResidualVector
    evalCouplingResidual(BulkIdType,
                         const BulkLocalAssembler& bulkLocalAssembler,
                         LowDimIdType,
                         GridIndexType<lowDimId> dofIdxGlobalJ)
    {
        const auto& map = couplingMapperPtr_->couplingMap(bulkGridId, lowDimGridId);

        assert(bulkContext_.isSet);
        assert(bulkElemIsCoupled_[bulkContext_.elementIdx]);
        assert(map.find(bulkContext_.elementIdx) != map.end());
        assert(bulkContext_.elementIdx == this->problem(bulkId).gridGeometry().elementMapper().index(bulkLocalAssembler.element()));

        typename LocalResidual<bulkId>::ElementResidualVector res(1);
        res = 0.0;
        res[0] = evalBulkFluxes(bulkLocalAssembler.element(),
                                bulkLocalAssembler.fvGeometry(),
                                bulkLocalAssembler.curElemVolVars(),
                                bulkLocalAssembler.elemFluxVarsCache(),
                                bulkLocalAssembler.localResidual(),
                                map.find(bulkContext_.elementIdx)->second.dofToCouplingScvfMap.at(dofIdxGlobalJ));
        return res;
    }

    /*!
     * \brief Evaluates the coupling element residual of a lower-dimensional domain element
     *        with respect to a dof in the bulk domain (dofIdxGlobalJ). This is essentially
     *        the fluxes across the facets of the neighboring bulk elements that coincide with
     *        the given element.
     *
     * \note The coupling residual in this case is always the entire transfer flux from bulk
     *       to the lowDim domain. It is therefore independent of the given dof index in the
     *       bulk domain, which is why we directly forward to the index-independent function.
     */
    template< class LowDimLocalAssembler >
    typename LocalResidual<lowDimId>::ElementResidualVector
    evalCouplingResidual(LowDimIdType,
                         const LowDimLocalAssembler& lowDimLocalAssembler,
                         BulkIdType,
                         GridIndexType<bulkId> dofIdxGlobalJ)
    { return evalCouplingResidual(lowDimId, lowDimLocalAssembler, bulkId); }

    /*!
     * \brief Evaluates the coupling element residual of a lower-dimensional domain element
     *        with respect to a dof in the bulk domain (dofIdxGlobalJ). This is essentially
     *        the fluxes across the facets of the neighboring bulk elements.
     */
    template< class LowDimLocalAssembler >
    typename LocalResidual<lowDimId>::ElementResidualVector
    evalCouplingResidual(LowDimIdType, const LowDimLocalAssembler& lowDimLocalAssembler, BulkIdType)
    {
        // make sure this is called for the element for which the context was set
        assert(lowDimContext_.isSet);
        assert(this->problem(lowDimId).gridGeometry().elementMapper().index(lowDimLocalAssembler.element()) == lowDimContext_.elementIdx);

        // evaluate sources for the first scv
        // the sources are element-wise & scv-independent since we use tpfa in bulk domain
        const auto source = evalSourcesFromBulk(lowDimLocalAssembler.element(),
                                                lowDimLocalAssembler.fvGeometry(),
                                                lowDimLocalAssembler.curElemVolVars(),
                                                *scvs(lowDimLocalAssembler.fvGeometry()).begin());

        // fill element residual vector with the sources
        typename LocalResidual<lowDimId>::ElementResidualVector res(lowDimLocalAssembler.fvGeometry().numScv());
        res = 0.0;
        for (const auto& scv : scvs(lowDimLocalAssembler.fvGeometry()))
             res[scv.localDofIndex()] -= source;

        return res;
    }

    /*!
     * \brief Computes the sources in a lower-dimensional element stemming from the bulk domain.
     */
    NumEqVector<lowDimId> evalSourcesFromBulk(const Element<lowDimId>& element,
                                              const FVElementGeometry<lowDimId>& fvGeometry,
                                              const ElementVolumeVariables<lowDimId>& elemVolVars,
                                              const SubControlVolume<lowDimId>& scv)
    {
        // make sure the this is called for the element of the context
        assert(this->problem(lowDimId).gridGeometry().elementMapper().index(element) == lowDimContext_.elementIdx);

        NumEqVector<lowDimId> sources(0.0);

        const auto& map = couplingMapperPtr_->couplingMap(lowDimGridId, bulkGridId);
        auto it = map.find(lowDimContext_.elementIdx);
        if (it == map.end())
            return sources;

        assert(lowDimContext_.isSet);
        for (const auto& embedment : it->second.embedments)
            sources += evalBulkFluxes(this->problem(bulkId).gridGeometry().element(embedment.first),
                                      *lowDimContext_.bulkFvGeometry,
                                      *lowDimContext_.bulkElemVolVars,
                                      *lowDimContext_.bulkElemFluxVarsCache,
                                      *lowDimContext_.bulkLocalResidual,
                                      embedment.second);

        // if lowdim domain uses box, we distribute the sources equally among the scvs
        if (lowDimUsesBox)
            sources /= fvGeometry.numScv();

        return sources;
    }

    /*!
     * \brief For the assembly of the element residual of a bulk domain element
     *        we need to prepare all variables of lower-dimensional domain elements
     *        that are coupled to the given bulk element
     */
    template< class Assembler >
    void bindCouplingContext(BulkIdType, const Element<bulkId>& element, const Assembler& assembler)
    {
        // clear context
        bulkContext_.reset();

        // set index in context in any case
        const auto bulkElemIdx = this->problem(bulkId).gridGeometry().elementMapper().index(element);
        bulkContext_.elementIdx = bulkElemIdx;

        // if element is coupled, actually set the context
        if (bulkElemIsCoupled_[bulkElemIdx])
        {
            const auto& map = couplingMapperPtr_->couplingMap(bulkGridId, lowDimGridId);

            auto it = map.find(bulkElemIdx); assert(it != map.end());
            const auto& elementStencil = it->second.couplingElementStencil;
            bulkContext_.lowDimFvGeometries.reserve(elementStencil.size());
            bulkContext_.lowDimVolVars.reserve(elementStencil.size());

            for (const auto lowDimElemIdx : elementStencil)
            {
                const auto& ldSol = Assembler::isImplicit() ? this->curSol()[lowDimId] : assembler.prevSol()[lowDimId];
                const auto& ldProblem = this->problem(lowDimId);
                const auto& ldGridGeometry = this->problem(lowDimId).gridGeometry();

                const auto elemJ = ldGridGeometry.element(lowDimElemIdx);
                auto fvGeom = localView(ldGridGeometry);
                fvGeom.bindElement(elemJ);

                VolumeVariables<lowDimId> volVars;

                // if low dim domain uses the box scheme, we have to create interpolated vol vars
                if (lowDimUsesBox)
                {
                    const auto elemGeom = elemJ.geometry();
                    FacetCoupling::makeInterpolatedVolVars(volVars, ldProblem, ldSol, fvGeom, elemJ, elemGeom, elemGeom.center());
                }
                // if low dim domain uses a cc scheme we can directly update the vol vars
                else
                    volVars.update( elementSolution(elemJ, ldSol, ldGridGeometry),
                                    ldProblem,
                                    elemJ,
                                    fvGeom.scv(lowDimElemIdx) );

                bulkContext_.isSet = true;
                bulkContext_.lowDimFvGeometries.emplace_back( std::move(fvGeom) );
                bulkContext_.lowDimVolVars.emplace_back( std::move(volVars) );
            }
        }
    }

    /*!
     * \brief For the assembly of the element residual of a bulk domain element
     *        we need to prepare the local views of one of the neighboring bulk
     *        domain elements. These are used later to compute the fluxes across
     *        the faces over which the coupling occurs
     * \note Since the low-dim coupling residua are fluxes stemming from
     *       the bulk domain, we have to prepare the bulk coupling context
     *       for the neighboring element (where fluxes are calculated) as well.
     */
    template< class Assembler >
    void bindCouplingContext(LowDimIdType, const Element<lowDimId>& element, const Assembler& assembler)
    {
        // reset contexts
        bulkContext_.reset();
        lowDimContext_.reset();

        // set index in context in any case
        const auto lowDimElemIdx = this->problem(lowDimId).gridGeometry().elementMapper().index(element);
        lowDimContext_.elementIdx = lowDimElemIdx;

        const auto& map = couplingMapperPtr_->couplingMap(lowDimGridId, bulkGridId);
        auto it = map.find(lowDimElemIdx);

        // if element is coupled, actually set the context
        if (it != map.end())
        {
            // first bind the low dim context for the first neighboring bulk element
            const auto& bulkGridGeom = this->problem(bulkId).gridGeometry();
            const auto bulkElem = bulkGridGeom.element(it->second.embedments[0].first);
            bindCouplingContext(bulkId, bulkElem, assembler);

            // then simply bind the local views of that first neighbor
            auto bulkFvGeom = localView(bulkGridGeom);
            auto bulkElemVolVars = Assembler::isImplicit() ? localView(assembler.gridVariables(bulkId).curGridVolVars())
                                                           : localView(assembler.gridVariables(bulkId).prevGridVolVars());
            auto bulkElemFluxVarsCache = localView(assembler.gridVariables(bulkId).gridFluxVarsCache());

            // evaluate variables on old/new time level depending on time disc scheme
            const auto& bulkSol = Assembler::isImplicit() ? this->curSol()[bulkId] : assembler.prevSol()[bulkId];
            bulkFvGeom.bind(bulkElem);
            bulkElemVolVars.bind(bulkElem, bulkFvGeom, bulkSol);
            bulkElemFluxVarsCache.bind(bulkElem, bulkFvGeom, bulkElemVolVars);

            lowDimContext_.isSet = true;
            lowDimContext_.bulkFvGeometry = std::make_unique< FVElementGeometry<bulkId> >( std::move(bulkFvGeom) );
            lowDimContext_.bulkElemVolVars = std::make_unique< ElementVolumeVariables<bulkId> >( std::move(bulkElemVolVars) );
            lowDimContext_.bulkElemFluxVarsCache = std::make_unique< ElementFluxVariablesCache<bulkId> >( std::move(bulkElemFluxVarsCache) );
            lowDimContext_.bulkLocalResidual = std::make_unique< LocalResidual<bulkId> >(assembler.localResidual(bulkId));
        }
    }

    /*!
     * \brief After deflecting the solution of the lower-dimensional domain,
     *        we have to update the element volume variables object if the context.
     */
    template< class BulkLocalAssembler >
    void updateCouplingContext(BulkIdType domainI,
                               const BulkLocalAssembler& bulkLocalAssembler,
                               LowDimIdType domainJ,
                               GridIndexType<lowDimId> dofIdxGlobalJ,
                               const PrimaryVariables<lowDimId>& priVarsJ,
                               unsigned int pvIdxJ)
    {
        // communicate deflected solution
        ParentType::updateCouplingContext(domainI, bulkLocalAssembler, domainJ, dofIdxGlobalJ, priVarsJ, pvIdxJ);

        // Since coupling only occurs via the fluxes, the context does not
        // have to be updated in explicit time discretization schemes, where
        // they are strictly evaluated on the old time level
        if (!BulkLocalAssembler::isImplicit())
            return;

        // skip the rest if context is empty
        if (bulkContext_.isSet)
        {
            const auto& map = couplingMapperPtr_->couplingMap(bulkGridId, lowDimGridId);
            const auto& couplingElemStencil = map.find(bulkContext_.elementIdx)->second.couplingElementStencil;
            const auto& ldSol = this->curSol()[lowDimId];
            const auto& ldProblem = this->problem(lowDimId);
            const auto& ldGridGeometry = this->problem(lowDimId).gridGeometry();

            // find the low-dim elements in coupling stencil, where this dof is contained in
            const auto couplingElements = [&] ()
            {
                if (lowDimUsesBox)
                {
                    std::vector< Element<lowDimId> > lowDimElems;
                    std::for_each( couplingElemStencil.begin(), couplingElemStencil.end(),
                                   [&] (auto lowDimElemIdx)
                                   {
                                       auto element = ldGridGeometry.element(lowDimElemIdx);
                                       for (int i = 0; i < element.geometry().corners(); ++i)
                                       {
                                           const auto dofIdx = ldGridGeometry.vertexMapper().subIndex(element, i, lowDimDim);
                                           if (dofIdxGlobalJ == dofIdx) { lowDimElems.emplace_back( std::move(element) ); break; }
                                       }
                                   } );
                    return lowDimElems;
                }
                // dof index = element index for cc schemes
                else
                    return std::vector<Element<lowDimId>>( {ldGridGeometry.element(dofIdxGlobalJ)} );
            } ();

            // update all necessary vol vars in context
            for (const auto& element : couplingElements)
            {
                // find index in coupling context
                const auto eIdxGlobal = ldGridGeometry.elementMapper().index(element);
                auto it = std::find(couplingElemStencil.begin(), couplingElemStencil.end(), eIdxGlobal);
                const auto idxInContext = std::distance(couplingElemStencil.begin(), it);
                assert(it != couplingElemStencil.end());

                auto& volVars = bulkContext_.lowDimVolVars[idxInContext];
                const auto& fvGeom = bulkContext_.lowDimFvGeometries[idxInContext];
                // if low dim domain uses the box scheme, we have to create interpolated vol vars
                if (lowDimUsesBox)
                {
                    const auto elemGeom = element.geometry();
                    FacetCoupling::makeInterpolatedVolVars(volVars, ldProblem, ldSol, fvGeom, element, elemGeom, elemGeom.center());
                }
                // if low dim domain uses a cc scheme we can directly update the vol vars
                else
                    volVars.update( elementSolution(element, ldSol, ldGridGeometry),
                                    ldProblem,
                                    element,
                                    fvGeom.scv(eIdxGlobal) );
            }
        }
    }

    /*!
     * \brief Update the coupling context for a derivative bulk -> bulk.
     *        Here, we simply have to update the solution.
     */
    template< class BulkLocalAssembler >
    void updateCouplingContext(BulkIdType domainI,
                               const BulkLocalAssembler& bulkLocalAssembler,
                               BulkIdType domainJ,
                               GridIndexType<bulkId> dofIdxGlobalJ,
                               const PrimaryVariables<bulkId>& priVarsJ,
                               unsigned int pvIdxJ)
    {
        // communicate deflected solution
        ParentType::updateCouplingContext(domainI, bulkLocalAssembler, domainJ, dofIdxGlobalJ, priVarsJ, pvIdxJ);
    }

    /*!
     * \brief After deflecting the solution of the bulk domain, we have to update
     *        the element volume variables and transmissibilities of the neighboring
     *        bulk element stored in the context.
     */
    template< class LowDimLocalAssembler >
    void updateCouplingContext(LowDimIdType domainI,
                               const LowDimLocalAssembler& lowDimLocalAssembler,
                               BulkIdType domainJ,
                               GridIndexType<bulkId> dofIdxGlobalJ,
                               const PrimaryVariables<bulkId>& priVarsJ,
                               unsigned int pvIdxJ)
    {
        // communicate deflected solution
        ParentType::updateCouplingContext(domainI, lowDimLocalAssembler, domainJ, dofIdxGlobalJ, priVarsJ, pvIdxJ);

        // Since coupling only occurs via the fluxes, the context does not
        // have to be updated in explicit time discretization schemes, where
        // they are strictly evaluated on the old time level
        if (!LowDimLocalAssembler::isImplicit())
            return;

        // skip the rest if context is empty
        if (lowDimContext_.isSet)
        {
            assert(lowDimContext_.elementIdx == this->problem(lowDimId).gridGeometry().elementMapper().index(lowDimLocalAssembler.element()));

            // since we use cc scheme in bulk domain: dof index = element index
            const auto& bulkGridGeom = this->problem(bulkId).gridGeometry();
            const auto elementJ = bulkGridGeom.element(dofIdxGlobalJ);

            // update corresponding vol vars in context
            const auto& scv = lowDimContext_.bulkFvGeometry->scv(dofIdxGlobalJ);
            const auto elemSol = elementSolution(elementJ, this->curSol()[bulkId], bulkGridGeom);
            (*lowDimContext_.bulkElemVolVars)[dofIdxGlobalJ].update(elemSol, this->problem(bulkId), elementJ, scv);

            // update the element flux variables cache (tij might be solution-dependent)
            if (dofIdxGlobalJ == bulkContext_.elementIdx)
                lowDimContext_.bulkElemFluxVarsCache->update( elementJ, *lowDimContext_.bulkFvGeometry, *lowDimContext_.bulkElemVolVars);
            else
                lowDimContext_.bulkElemFluxVarsCache->update( this->problem(bulkId).gridGeometry().element(bulkContext_.elementIdx),
                                                              *lowDimContext_.bulkFvGeometry,
                                                              *lowDimContext_.bulkElemVolVars );
        }
    }

    /*!
     * \brief After deflecting the solution of the lower-dimensional domain has been deflected
     *        during the assembly of the element residual of a lower-dimensional element, we
     *        have to communicate this to the volume variables stored in the context as well
     *        as the transmissibilities.
     */
    template< class LowDimLocalAssembler >
    void updateCouplingContext(LowDimIdType domainI,
                               const LowDimLocalAssembler& lowDimLocalAssembler,
                               LowDimIdType domainJ,
                               GridIndexType<lowDimId> dofIdxGlobalJ,
                               const PrimaryVariables<lowDimId>& priVarsJ,
                               unsigned int pvIdxJ)
    {
        // communicate deflected solution
        ParentType::updateCouplingContext(domainI, lowDimLocalAssembler, domainJ, dofIdxGlobalJ, priVarsJ, pvIdxJ);

        // Since coupling only occurs via the fluxes, the context does not
        // have to be updated in explicit time discretization schemes, where
        // they are strictly evaluated on the old time level
        if (!LowDimLocalAssembler::isImplicit())
            return;

        // skip the rest if context is empty
        if (lowDimContext_.isSet)
        {
            const auto& ldSol = this->curSol()[lowDimId];
            const auto& ldProblem = this->problem(lowDimId);
            const auto& ldGridGeometry = this->problem(lowDimId).gridGeometry();

            assert(bulkContext_.isSet);
            assert(lowDimContext_.elementIdx == ldGridGeometry.elementMapper().index(lowDimLocalAssembler.element()));

            // update the corresponding vol vars in the bulk context
            const auto& bulkMap = couplingMapperPtr_->couplingMap(bulkGridId, lowDimGridId);
            const auto& couplingElementStencil = bulkMap.find(bulkContext_.elementIdx)->second.couplingElementStencil;
            auto it = std::find(couplingElementStencil.begin(), couplingElementStencil.end(), lowDimContext_.elementIdx);
            assert(it != couplingElementStencil.end());
            const auto idxInContext = std::distance(couplingElementStencil.begin(), it);

            auto& volVars = bulkContext_.lowDimVolVars[idxInContext];
            const auto& fvGeom = bulkContext_.lowDimFvGeometries[idxInContext];
            const auto& element = lowDimLocalAssembler.element();
            // if low dim domain uses the box scheme, we have to create interpolated vol vars
            if (lowDimUsesBox)
            {
                const auto elemGeom = element.geometry();
                FacetCoupling::makeInterpolatedVolVars(volVars, ldProblem, ldSol, fvGeom, element, elemGeom, elemGeom.center());
            }
            // if low dim domain uses a cc scheme we can directly update the vol vars
            else
                volVars.update( elementSolution(element, ldSol, ldGridGeometry),
                                ldProblem,
                                element,
                                fvGeom.scv(lowDimContext_.elementIdx) );

            // update the element flux variables cache (tij depend on low dim values in context)
            const auto contextElem = this->problem(bulkId).gridGeometry().element(bulkContext_.elementIdx);
            lowDimContext_.bulkElemFluxVarsCache->update(contextElem, *lowDimContext_.bulkFvGeometry, *lowDimContext_.bulkElemVolVars);
        }
    }

    //! pull up empty update function to be used for low-dim domain
    using ParentType::updateCoupledVariables;

    /*!
     * \brief Update the transmissibilities in the bulk domain after the coupling context changed
     * \note Specialization of the function for deactivated grid-wide volume variables caching
     */
    template< class BulkLocalAssembler, class UpdatableFluxVarCache >
    void updateCoupledVariables(BulkIdType domainI,
                                const BulkLocalAssembler& bulkLocalAssembler,
                                ElementVolumeVariables<bulkId>& elemVolVars,
                                UpdatableFluxVarCache& fluxVarsCache)
    {
        // update transmissibilities after low dim context has changed (implicit only)
        if (BulkLocalAssembler::isImplicit())
            fluxVarsCache.update(bulkLocalAssembler.element(),
                                 bulkLocalAssembler.fvGeometry(),
                                 bulkLocalAssembler.curElemVolVars());
    }

    /*!
     * \brief Update the transmissibilities in the bulk domain after the coupling context changed
     * \note Specialization of the function for enabled grid-wide volume variables caching
     */
    template< class BulkLocalAssembler, class UpdatableFluxVarCache >
    void updateCoupledVariables(BulkIdType domainI,
                                const BulkLocalAssembler& bulkLocalAssembler,
                                GridVolumeVariables<bulkId>& gridVolVars,
                                UpdatableFluxVarCache& fluxVarsCache)
    {
        // update transmissibilities after low dim context has changed (implicit only)
        if (BulkLocalAssembler::isImplicit())
        {
            auto elemVolVars = localView(gridVolVars);
            elemVolVars.bind(bulkLocalAssembler.element(), bulkLocalAssembler.fvGeometry(), this->curSol()[bulkId]);
            fluxVarsCache.update(bulkLocalAssembler.element(), bulkLocalAssembler.fvGeometry(), elemVolVars);
        }
    }

    //! Empty stencil to be returned for elements that aren't coupled
    template<std::size_t id, std::enable_if_t<(id == bulkId || id == lowDimId), int> = 0>
    const typename CouplingMapper::template Stencil<id>&
    getEmptyStencil(Dune::index_constant<id>) const
    { return std::get<(id == bulkId ? 0 : 1)>(emptyStencilTuple_); }

protected:
    //! Return const references to the bulk coupling contexts
    const BulkCouplingContext& bulkCouplingContext() const { return bulkContext_; }
    const LowDimCouplingContext& lowDimCouplingContext() const { return lowDimContext_; }

    //! Return references to the bulk coupling contexts
    BulkCouplingContext& bulkCouplingContext() { return bulkContext_; }
    LowDimCouplingContext& lowDimCouplingContext() { return lowDimContext_; }

    //! evaluates the bulk-facet exchange fluxes for a given facet element
    template<class BulkScvfIndices>
    NumEqVector<bulkId> evalBulkFluxes(const Element<bulkId>& elementI,
                                       const FVElementGeometry<bulkId>& fvGeometry,
                                       const ElementVolumeVariables<bulkId>& elemVolVars,
                                       const ElementFluxVariablesCache<bulkId>& elemFluxVarsCache,
                                       const LocalResidual<bulkId>& localResidual,
                                       const BulkScvfIndices& scvfIndices) const
    {

        NumEqVector<bulkId> coupledFluxes(0.0);
        for (const auto& scvfIdx : scvfIndices)
            coupledFluxes += localResidual.evalFlux(this->problem(bulkId),
                                                    elementI,
                                                    fvGeometry,
                                                    elemVolVars,
                                                    elemFluxVarsCache,
                                                    fvGeometry.scvf(scvfIdx));
        return coupledFluxes;
    }

private:
    std::shared_ptr<CouplingMapper> couplingMapperPtr_;

    //! store bools for all bulk elements/scvfs that indicate if they
    //! are coupled, so that we don't have to search in the map every time
    std::vector<bool> bulkElemIsCoupled_;
    std::vector<bool> bulkScvfIsCoupled_;

    //! empty stencil to return for non-coupled elements
    using BulkStencil = typename CouplingMapper::template Stencil<bulkId>;
    using LowDimStencil = typename CouplingMapper::template Stencil<lowDimId>;
    std::tuple<BulkStencil, LowDimStencil> emptyStencilTuple_;

    //! The coupling contexts of the two domains
    BulkCouplingContext bulkContext_;
    LowDimCouplingContext lowDimContext_;
};

} // end namespace Dumux

#endif
