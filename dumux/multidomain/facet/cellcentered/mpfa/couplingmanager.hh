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
 * \ingroup FacetCoupling
 * \copydoc Dumux::FacetCouplingManager
 */
#ifndef DUMUX_CCMPFA_FACETCOUPLING_MANAGER_HH
#define DUMUX_CCMPFA_FACETCOUPLING_MANAGER_HH

#include <algorithm>
#include <cassert>

#include <dumux/common/properties.hh>
#include <dumux/common/indextraits.hh>
#include <dumux/common/numericdifferentiation.hh>
#include <dumux/assembly/numericepsilon.hh>

#include <dumux/discretization/method.hh>
#include <dumux/discretization/elementsolution.hh>
#include <dumux/multidomain/couplingmanager.hh>
#include <dumux/multidomain/facet/cellcentered/tpfa/couplingmanager.hh>

namespace Dumux {

/*!
 * \ingroup FacetCoupling
 * \brief Manages the coupling between bulk elements and lower dimensional elements
 *        where the coupling occurs across the facets of the bulk grid. This implementation
 *        is to be used in conjunction with models using the cell-centered mpfa scheme.
 *
 * \tparam MDTraits The multidomain traits containing the types on all sub-domains
 * \tparam CouplingMapper Class containing maps on the coupling between dofs of different grids
 * \tparam bulkDomainId The domain id of the bulk problem
 * \tparam lowDimDomainId The domain id of the lower-dimensional problem
 */
template<class MDTraits, class CouplingMapper, std::size_t bulkDomainId, std::size_t lowDimDomainId>
class FacetCouplingManager<MDTraits, CouplingMapper, bulkDomainId, lowDimDomainId, DiscretizationMethod::ccmpfa>
: public FacetCouplingManager<MDTraits, CouplingMapper, bulkDomainId, lowDimDomainId, DiscretizationMethod::cctpfa>
{
    using ParentType = FacetCouplingManager<MDTraits, CouplingMapper, bulkDomainId, lowDimDomainId, DiscretizationMethod::cctpfa>;

    // domain id instances
    using BulkIdType = typename MDTraits::template SubDomain<bulkDomainId>::Index;
    using LowDimIdType = typename MDTraits::template SubDomain<lowDimDomainId>::Index;
    static constexpr auto bulkId = BulkIdType();
    static constexpr auto lowDimId = LowDimIdType();

    // the sub-domain type tags
    template<std::size_t id> using SubDomainTypeTag = typename MDTraits::template SubDomain<id>::TypeTag;

    // further types specific to the sub-problems
    template<std::size_t id> using Scalar = GetPropType<SubDomainTypeTag<id>, Properties::Scalar>;
    template<std::size_t id> using Problem = GetPropType<SubDomainTypeTag<id>, Properties::Problem>;
    template<std::size_t id> using NumEqVector = GetPropType<SubDomainTypeTag<id>, Properties::NumEqVector>;
    template<std::size_t id> using GridGeometry = GetPropType<SubDomainTypeTag<id>, Properties::GridGeometry>;
    template<std::size_t id> using FVElementGeometry = typename GridGeometry<id>::LocalView;
    template<std::size_t id> using SubControlVolume = typename GridGeometry<id>::SubControlVolume;
    template<std::size_t id> using SubControlVolumeFace = typename GridGeometry<id>::SubControlVolumeFace;
    template<std::size_t id> using GridView = typename GridGeometry<id>::GridView;
    template<std::size_t id> using GridIndexType = typename IndexTraits< GridView<id> >::GridIndex;
    template<std::size_t id> using Element = typename GridView<id>::template Codim<0>::Entity;
    template<std::size_t id> using LocalResidual = GetPropType<SubDomainTypeTag<id>, Properties::LocalResidual>;

    template<std::size_t id> using GridVariables = GetPropType<SubDomainTypeTag<id>, Properties::GridVariables>;
    template<std::size_t id> using GridVolumeVariables = typename GridVariables<id>::GridVolumeVariables;
    template<std::size_t id> using ElementVolumeVariables = typename GridVolumeVariables<id>::LocalView;

    // grid ids
    static constexpr int bulkDim = GridView<bulkDomainId>::dimension;
    static constexpr int lowDimDim = GridView<lowDimDomainId>::dimension;
    static constexpr auto bulkGridId = CouplingMapper::template gridId<bulkDim>();
    static constexpr auto lowDimGridId = CouplingMapper::template gridId<lowDimDim>();

    static constexpr bool lowDimUsesBox = GridGeometry<lowDimId>::discMethod == DiscretizationMethod::box;

public:

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
        // Initialize the parent class
        ParentType::init(bulkProblem, lowDimProblem, couplingMapper, curSol);

        // determine all bulk scvfs that coincide with low dim elements
        bulkScvfIsOnFacetElement_.assign(bulkProblem->gridGeometry().numScvf(), false);
        const auto& bulkMap = couplingMapper->couplingMap(bulkGridId, lowDimGridId);
        for (const auto& entry : bulkMap)
            for (const auto& couplingEntry : entry.second.elementToScvfMap)
                for (const auto& scvfIdx : couplingEntry.second)
                    bulkScvfIsOnFacetElement_[scvfIdx] = true;

        // store pointer to mapper
        couplingMapperPtr_ = couplingMapper;
    }

    /*!
     * \brief returns true if a bulk scvf coincides with a facet element.
     */
    bool isOnInteriorBoundary(const Element<bulkId>& element,
                              const SubControlVolumeFace<bulkId>& scvf) const
    { return bulkScvfIsOnFacetElement_[scvf.index()]; }

    using ParentType::evalCouplingResidual;
    /*!
     * \brief Evaluates the coupling element residual of a lower-dimensional domain element
     *        with respect to a dof in the bulk domain (dofIdxGlobalJ). This is essentially
     *        the fluxes across the facets of the neighboring bulk elements.
     * \note  The coupling residual is independent of w.r.t. which bulk dof it is computed
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
        assert(this->lowDimCouplingContext().isSet);
        assert(this->problem(lowDimId).gridGeometry().elementMapper().index(lowDimLocalAssembler.element()) == this->lowDimCouplingContext().elementIdx);

        // fill element residual vector with the sources
        typename LowDimLocalAssembler::LocalResidual::ElementResidualVector res(lowDimLocalAssembler.fvGeometry().numScv());
        res = 0.0;
        for (const auto& scv : scvs(lowDimLocalAssembler.fvGeometry()))
             res[scv.localDofIndex()] -= evalSourcesFromBulk(lowDimLocalAssembler.element(),
                                                             lowDimLocalAssembler.fvGeometry(),
                                                             lowDimLocalAssembler.curElemVolVars(),
                                                             scv);
        return res;
    }

    /*!
     * \brief Computes the sources in a lower-dimensional sub-control volume stemming from the bulk domain.
     */
    NumEqVector<lowDimId> evalSourcesFromBulk(const Element<lowDimId>& element,
                                              const FVElementGeometry<lowDimId>& fvGeometry,
                                              const ElementVolumeVariables<lowDimId>& elemVolVars,
                                              const SubControlVolume<lowDimId>& scv)
    {
        // make sure the this is called for the element of the context
        assert(this->problem(lowDimId).gridGeometry().elementMapper().index(element) == this->lowDimCouplingContext().elementIdx);

        NumEqVector<lowDimId> sources(0.0);

        const auto& map = couplingMapperPtr_->couplingMap(lowDimGridId, bulkGridId);
        auto it = map.find(this->lowDimCouplingContext().elementIdx);
        if (it == map.end())
            return sources;

        assert(this->lowDimCouplingContext().isSet);
        for (const auto& embedment : it->second.embedments)
        {
            // list of scvfs in the bulk domain whose fluxes enter this scv
            // if low dim domain uses a cc scheme, this is all scvfs lying on this element
            // if it uses box, it is the one scvf coinciding with the given scv
            const auto& coincidingScvfs = embedment.second;
            const auto& scvfList = lowDimUsesBox ? std::vector<GridIndexType<lowDimId>>{ coincidingScvfs[scv.localDofIndex()] }
                                                 : coincidingScvfs;

            sources += this->evalBulkFluxes(this->problem(bulkId).gridGeometry().element(embedment.first),
                                            *this->lowDimCouplingContext().bulkFvGeometry,
                                            *this->lowDimCouplingContext().bulkElemVolVars,
                                            *this->lowDimCouplingContext().bulkElemFluxVarsCache,
                                            *this->lowDimCouplingContext().bulkLocalResidual,
                                            scvfList);
        }

        return sources;
    }

    /*!
     * \brief Extend the jacobian pattern of the diagonal block of the lowdim domain
     *        by the elements that are in the coupling stencil of the neighboring bulk elements
     */
    template<class JacobianPattern>
    void extendJacobianPattern(LowDimIdType, JacobianPattern& pattern) const
    {
        const auto& lowDimFVGridGeometry = this->problem(lowDimId).gridGeometry();
        for (const auto& element : elements(lowDimFVGridGeometry.gridView()))
        {

            const auto eIdx = lowDimFVGridGeometry.elementMapper().index(element);
            const auto& map = couplingMapperPtr_->couplingMap(lowDimGridId, bulkGridId);
            auto it = map.find(eIdx);

            // if element is coupled, take one of the neighbors and add coupling stencil to pattern
            if (it != map.end())
            {
                // coupling stencil of the first neighbor
                const auto bulkElemIdx = it->second.embedments[0].first;
                const auto& bulkMapEntry = couplingMapperPtr_->couplingMap(bulkGridId, lowDimGridId).at(bulkElemIdx);
                const auto& couplingStencil = bulkMapEntry.couplingStencil;

                for (auto globalJ : couplingStencil)
                {
                    if (lowDimUsesBox)
                    {
                        for (int i = 0; i < element.subEntities(lowDimDim); ++i)
                            pattern.add(lowDimFVGridGeometry.vertexMapper().subIndex(element, i, lowDimDim), globalJ);
                    }
                    else
                        pattern.add(eIdx, globalJ);
                }
            }
        }
    }

    //! The bulk domain has no extended jacobian pattern
    template<class JacobianPattern>
    void extendJacobianPattern(BulkIdType, JacobianPattern& pattern) const
    {}

    /*!
     * \brief evaluate additional derivatives of the element residual of the low-dim domain with respect
     *        to dofs in the same domain that are not in the regular stencil (see extendJacobianPattern)
     * \note Here, this is the change of the source term with respect to changes in the variables of the
     *       other elements in the coupling stencil of the neighboring bulk elements.
     */
    template<class LowDimLocalAssembler, class JacobianMatrixDiagBlock, class GridVariables>
    void evalAdditionalDomainDerivatives(LowDimIdType,
                                         const LowDimLocalAssembler& lowDimLocalAssembler,
                                         const typename LowDimLocalAssembler::LocalResidual::ElementResidualVector&,
                                         JacobianMatrixDiagBlock& A,
                                         GridVariables& gridVariables)
    {
        // Since coupling only occurs via the fluxes, there are no
        // additional derivatives for explicit time integration schemes
        if (!LowDimLocalAssembler::isImplicit())
            return;

        // lambda to update the coupling context for a given lowDim element/dofIdx
        auto updateContext = [&] (auto elemIdx, auto dofIdx, auto priVars, auto pvIdx)
        {
            // deflect the solution
            auto& ldSol = this->curSol()[lowDimId];
            ldSol[dofIdx][pvIdx] = priVars[pvIdx];

            // update the corresponding vol vars in the bulk context
            assert(this->bulkCouplingContext().isSet);
            const auto& bulkMap = couplingMapperPtr_->couplingMap(bulkGridId, lowDimGridId);
            const auto& couplingElementStencil = bulkMap.find(this->bulkCouplingContext().elementIdx)->second.couplingElementStencil;

            auto it = std::find(couplingElementStencil.begin(), couplingElementStencil.end(), elemIdx);
            assert(it != couplingElementStencil.end());
            const auto idxInContext = std::distance(couplingElementStencil.begin(), it);

            auto& volVars = this->bulkCouplingContext().lowDimVolVars[idxInContext];
            const auto& fvGeom = this->bulkCouplingContext().lowDimFvGeometries[idxInContext];
            const auto& element = this->problem(lowDimId).gridGeometry().element(elemIdx);

            // if low dim domain uses the box scheme, we have to create interpolated vol vars
            if (lowDimUsesBox)
            {
                const auto elemGeom = element.geometry();
                FacetCoupling::makeInterpolatedVolVars(volVars, this->problem(lowDimId), ldSol, fvGeom, element, elemGeom, elemGeom.center());
            }
            // if low dim domain uses a cc scheme we can directly update the vol vars
            else
                volVars.update( elementSolution(element, ldSol, this->problem(lowDimId).gridGeometry()),
                                this->problem(lowDimId),
                                element,
                                fvGeom.scv(elemIdx) );

            // update the element flux variables cache (tij depend on low dim values in context)
            const auto contextElem = this->problem(bulkId).gridGeometry().element(this->bulkCouplingContext().elementIdx);
            this->lowDimCouplingContext().bulkElemFluxVarsCache->update(contextElem,
                                                                        *this->lowDimCouplingContext().bulkFvGeometry,
                                                                        *this->lowDimCouplingContext().bulkElemVolVars);
        };

        const auto eIdx = this->problem(lowDimId).gridGeometry().elementMapper().index(lowDimLocalAssembler.element());

        // bug tracking
        assert(this->lowDimCouplingContext().isSet);
        assert(this->lowDimCouplingContext().elementIdx == eIdx);

        // if the element is coupled, evaluate additional source derivatives
        const auto& map = couplingMapperPtr_->couplingMap(lowDimGridId, bulkGridId);
        auto it = map.find(eIdx);
        if (it != map.end())
            evalLowDimSourceDerivatives_(updateContext, lowDimLocalAssembler, A);
    }

    //! The bulk domain has no additional derivatives
    template<class LocalAssemblerI, class JacobianMatrixDiagBlock, class GridVariables>
    void evalAdditionalDomainDerivatives(BulkIdType,
                                         const LocalAssemblerI& localAssemblerI,
                                         const typename LocalAssemblerI::LocalResidual::ElementResidualVector& origResiduals,
                                         JacobianMatrixDiagBlock& A,
                                         GridVariables& gridVariables)
    {}

private:
    //! evaluates the additional source derivatives w.r.t. to neighboring elements
    template<class UpdateContext, class LowDimLocalAssembler, class JacobianMatrixDiagBlock>
    void evalLowDimSourceDerivatives_(const UpdateContext& updateContext,
                                      const LowDimLocalAssembler& lowDimLocalAssembler,
                                      JacobianMatrixDiagBlock& A)
    {
        const auto& lowDimFVGridGeometry = this->problem(lowDimId).gridGeometry();
        const auto eIdx = lowDimFVGridGeometry.elementMapper().index(lowDimLocalAssembler.element());

        // coupling stencil of the first neighbor
        const auto bulkElemIdx = this->bulkCouplingContext().elementIdx;
        const auto& bulkMapEntry = couplingMapperPtr_->couplingMap(bulkGridId, lowDimGridId).at(bulkElemIdx);
        const auto& couplingStencil = bulkMapEntry.couplingStencil;
        const auto& couplingElementStencil = bulkMapEntry.couplingElementStencil;

        // compute the undeflected residual (reuse coupling residual function)
        const auto origResidual = evalCouplingResidual(lowDimId, lowDimLocalAssembler, bulkId);

        // container of dofs within this element
        std::vector< std::decay_t<decltype(couplingStencil[0])> > elemDofs;
        elemDofs.reserve(lowDimLocalAssembler.fvGeometry().numScv());
        for (const auto& scv : scvs(lowDimLocalAssembler.fvGeometry()))
            elemDofs.push_back(scv.dofIndex());

        // compute derivate for all additional dofs in the stencil
        for (const auto couplingElemIdx : couplingElementStencil)
        {
            // skip the same element
            if (couplingElemIdx == eIdx)
                continue;

            // container of dofs within the other element
            std::vector< std::decay_t<decltype(couplingStencil[0])> > elemDofsJ;
            if (lowDimUsesBox)
            {
                const auto& elemJ = lowDimFVGridGeometry.element(couplingElemIdx);
                for (int i = 0; i < elemJ.subEntities(lowDimDim); ++i)
                    elemDofsJ.push_back(lowDimFVGridGeometry.vertexMapper().subIndex(elemJ, i, lowDimDim));
            }
            else
                elemDofsJ.push_back(couplingElemIdx);

            for (auto dofIndex : elemDofsJ)
            {
                auto partialDerivs = origResidual;
                const auto origPriVars = this->curSol()[lowDimId][dofIndex];

                // calculate derivatives w.r.t to the privars at the dof at hand
                static constexpr auto numEq = std::decay_t<decltype(origPriVars)>::dimension;
                for (int pvIdx = 0; pvIdx < numEq; pvIdx++)
                {
                    // reset partial derivatives
                    partialDerivs = 0.0;

                    auto evalResiduals = [&](Scalar<lowDimId> priVar)
                    {
                        auto priVars = origPriVars;
                        priVars[pvIdx] = priVar;

                        // Update context to deflected solution and reevaluate residual
                        updateContext(couplingElemIdx, dofIndex, priVars, pvIdx);
                        return this->evalCouplingResidual(lowDimId, lowDimLocalAssembler, bulkId);
                    };

                    static const int numDiffMethod = getParamFromGroup<int>(this->problem(lowDimId).paramGroup(), "Assembly.NumericDifferenceMethod");
                    static const NumericEpsilon< Scalar<lowDimId>, numEq > eps{this->problem(lowDimId).paramGroup()};
                    NumericDifferentiation::partialDerivative(evalResiduals, origPriVars[pvIdx], partialDerivs,
                                                              origResidual, eps(origPriVars[pvIdx], pvIdx), numDiffMethod);

                    // update the global stiffness matrix with the current partial derivatives
                    // A[i][col][eqIdx][pvIdx] is the rate of change of the residual of equation
                    // 'eqIdx' at dof 'i' depending on the primary variable 'pvIdx' at dof 'col'.
                    for (const auto& scv : scvs(lowDimLocalAssembler.fvGeometry()))
                        for (int eqIdx = 0; eqIdx < numEq; eqIdx++)
                            A[scv.dofIndex()][dofIndex][eqIdx][pvIdx] += partialDerivs[scv.indexInElement()][eqIdx];

                    // restore the original coupling context
                    updateContext(couplingElemIdx, dofIndex, origPriVars, pvIdx);
                }
            }
        }
    }

    //! store which scvfs coincide with facet element
    std::vector<bool> bulkScvfIsOnFacetElement_;

    //! store shared_ptr to coupling mapper
    std::shared_ptr< CouplingMapper > couplingMapperPtr_;
};

} // end namespace Dumux

#endif
