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
 * \ingroup Assembly
 * \ingroup CCDiscretization
 * \brief An assembler for Jacobian and residual contribution per element (cell-centered methods)
 * \tparam TypeTag the TypeTag
 * \tparam DM the differentiation method to residual compute derivatives
 * \tparam implicit if to use an implicit or explicit time discretization
 */
#ifndef DUMUX_MULTIDOMAIN_CC_LOCAL_ASSEMBLER_HH
#define DUMUX_MULTIDOMAIN_CC_LOCAL_ASSEMBLER_HH

#include <dune/common/reservedvector.hh>
#include <dune/grid/common/gridenums.hh> // for GhostEntity
#include <dune/istl/matrixindexset.hh>

#include <dumux/common/reservedblockvector.hh>
#include <dumux/common/properties.hh>
#include <dumux/common/parameters.hh>
#include <dumux/common/numericdifferentiation.hh>
#include <dumux/assembly/diffmethod.hh>

namespace Dumux {

/*!
 * \ingroup Assembly
 * \ingroup CCDiscretization
 * \brief A base class for all local assemblers
 * \tparam id the id of the sub domain
 * \tparam TypeTag the TypeTag
 * \tparam Assembler the assembler type
 */
template<std::size_t id, class TypeTag, class Assembler, class Implementation>
class SubDomainCCLocalAssemblerBase
{
    using Problem = typename GET_PROP_TYPE(TypeTag, Problem);
    using LocalResidualValues = typename GET_PROP_TYPE(TypeTag, NumEqVector);
    using LocalResidual = typename GET_PROP_TYPE(TypeTag, LocalResidual);
    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using JacobianMatrix = typename GET_PROP_TYPE(TypeTag, JacobianMatrix);
    using GridVariables = typename GET_PROP_TYPE(TypeTag, GridVariables);
    using SubSolutionVector = typename GET_PROP_TYPE(TypeTag, SolutionVector);
    using SolutionVector = typename Assembler::SolutionVector;
    using ElementBoundaryTypes = typename GET_PROP_TYPE(TypeTag, ElementBoundaryTypes);
    using FVGridGeometry = typename GET_PROP_TYPE(TypeTag, FVGridGeometry);
    using FVElementGeometry = typename FVGridGeometry::LocalView;
    using SubControlVolume = typename FVGridGeometry::SubControlVolume;
    using SubControlVolumeFace = typename FVGridGeometry::SubControlVolumeFace;
    using GridVolumeVariables = typename GET_PROP_TYPE(TypeTag, GridVolumeVariables);
    using ElementVolumeVariables = typename GET_PROP_TYPE(TypeTag, ElementVolumeVariables);
    using VolumeVariables = typename GET_PROP_TYPE(TypeTag, VolumeVariables);
    using ElementFluxVariablesCache = typename GET_PROP_TYPE(TypeTag, ElementFluxVariablesCache);
    using Element = typename GridView::template Codim<0>::Entity;
    using CouplingManager = typename Assembler::CouplingManager;

public:
    static constexpr auto domainId = typename Dune::index_constant<id>();

    explicit SubDomainCCLocalAssemblerBase(const Assembler& assembler,
                                           const Element& element,
                                           const SolutionVector& curSol,
                                           CouplingManager& couplingManager)
    : assembler_(assembler)
    , element_(element)
    , curSol_(curSol)
    , couplingManager_(couplingManager)
    , fvGeometry_(localView(assembler.fvGridGeometry(domainId)))
    , curElemVolVars_(localView(assembler.gridVariables(domainId).curGridVolVars()))
    , prevElemVolVars_(localView(assembler.gridVariables(domainId).prevGridVolVars()))
    , elemFluxVarsCache_(localView(assembler.gridVariables(domainId).gridFluxVarsCache()))
    , localResidual_(assembler.localResidual(domainId))
    , elementIsGhost_((element.partitionType() == Dune::GhostEntity))
    {}

    /*!
     * \brief Computes the derivatives with respect to the given element and adds them
     *        to the global matrix. The element residual is written into the right hand side.
     */
    template<class JacobianMatrixRow, class GridVariablesTuple>
    void assembleJacobianAndResidual(JacobianMatrixRow& jacRow, SubSolutionVector& res, GridVariablesTuple& gridVariables)
    {
        asImp_().bindLocalViews();

        // for the diagonal jacobian block
        const auto globalI = this->fvGeometry().fvGridGeometry().elementMapper().index(this->element());
        // forward to the internal implementation
        res[globalI] = asImp_().assembleJacobianAndResidualImpl(jacRow[domainId], *std::get<domainId>(gridVariables));
        // for the coupling blocks
        using namespace Dune::Hybrid;
        forEach(integralRange(Dune::Hybrid::size(jacRow)), [&](auto&& i)
        {
            if (i != domainId)
                assembleJacobianCoupling(i, jacRow, res[globalI], gridVariables);
        });
    }

    template<std::size_t otherId, class JacRow, class GridVariables,
             typename std::enable_if_t<(otherId == id), int> = 0>
    void assembleJacobianCoupling(Dune::index_constant<otherId> domainIdJ, JacRow& jacRow,
                                  const LocalResidualValues& res, GridVariables& gridVariables)
    {}

    template<std::size_t otherId, class JacRow, class GridVariables,
             typename std::enable_if_t<(otherId != id), int> = 0>
    void assembleJacobianCoupling(Dune::index_constant<otherId> domainIdJ, JacRow& jacRow,
                                  const LocalResidualValues& res, GridVariables& gridVariables)
    {
        asImp_().assembleJacobianCoupling(domainIdJ, jacRow[domainIdJ], res, *std::get<otherId>(gridVariables));
    }

    /*!
     * \brief Assemble the residual only
     */
    void assembleResidual(SubSolutionVector& res)
    {
        asImp_().bindLocalViews();
        const auto globalI = this->fvGeometry().fvGridGeometry().elementMapper().index(this->element());
        res[globalI] = asImp_().assembleResidualImpl(); // forward to the internal implementation
    }

    LocalResidualValues evalLocalResidual(const ElementVolumeVariables& elemVolVars) const
    {
        if (!assembler().isStationaryProblem())
        {
            auto residual = evalLocalFluxAndSourceResidual(elemVolVars);
            residual += evalLocalStorageResidual();
            return residual;
        }
        else
            return evalLocalFluxAndSourceResidual(elemVolVars);
    }

    LocalResidualValues evalLocalFluxAndSourceResidual(const ElementVolumeVariables& elemVolVars) const
    {
        return localResidual_.evalFluxAndSource(element_, fvGeometry_, elemVolVars, elemFluxVarsCache_, elemBcTypes_)[0];
    }

    LocalResidualValues evalLocalStorageResidual() const
    {
        return localResidual_.evalStorage(element_, fvGeometry_, prevElemVolVars_, curElemVolVars_)[0];
    }

    LocalResidualValues evalFluxResidual(const Element& neigbor,
                                         const SubControlVolumeFace& scvf) const
    {
        return localResidual_.evalFlux(problem(), neigbor, fvGeometry_, curElemVolVars_, elemFluxVarsCache_, scvf);
    }

    const Problem& problem() const
    { return assembler_.problem(domainId); }

    const Assembler& assembler() const
    { return assembler_; }

    const Element& element() const
    { return element_; }

    bool elementIsGhost() const
    { return elementIsGhost_; }

    const SolutionVector& curSol() const
    { return curSol_; }

    FVElementGeometry& fvGeometry()
    { return fvGeometry_; }

    ElementVolumeVariables& curElemVolVars()
    { return curElemVolVars_; }

    ElementVolumeVariables& prevElemVolVars()
    { return prevElemVolVars_; }

    ElementFluxVariablesCache& elemFluxVarsCache()
    { return elemFluxVarsCache_; }

    LocalResidual& localResidual()
    { return localResidual_; }

    ElementBoundaryTypes& elemBcTypes()
    { return elemBcTypes_; }

    const FVElementGeometry& fvGeometry() const
    { return fvGeometry_; }

    const ElementVolumeVariables& curElemVolVars() const
    { return curElemVolVars_; }

    const ElementVolumeVariables& prevElemVolVars() const
    { return prevElemVolVars_; }

    const ElementFluxVariablesCache& elemFluxVarsCache() const
    { return elemFluxVarsCache_; }

    const ElementBoundaryTypes& elemBcTypes() const
    { return elemBcTypes_; }

    const LocalResidual& localResidual() const
    { return localResidual_; }

    CouplingManager& couplingManager()
    { return couplingManager_; }

    template<class T = TypeTag, typename std::enable_if_t<!GET_PROP_VALUE(T, EnableGridVolumeVariablesCache), int> = 0>
    VolumeVariables& getVolVarAccess(GridVolumeVariables& gridVolVars, ElementVolumeVariables& elemVolVars, const SubControlVolume& scv)
    { return elemVolVars[scv]; }

    template<class T = TypeTag, typename std::enable_if_t<GET_PROP_VALUE(T, EnableGridVolumeVariablesCache), int> = 0>
    VolumeVariables& getVolVarAccess(GridVolumeVariables& gridVolVars, ElementVolumeVariables& elemVolVars, const SubControlVolume& scv)
    { return gridVolVars.volVars(scv); }

private:
    Implementation &asImp_()
    { return *static_cast<Implementation*>(this); }

    const Implementation &asImp_() const
    { return *static_cast<const Implementation*>(this); }

    const Assembler& assembler_; //!< access pointer to assembler instance
    const Element& element_; //!< the element whose residual is assembled
    const SolutionVector& curSol_; //!< the current solution
    CouplingManager& couplingManager_; //!< the coupling manager

    FVElementGeometry fvGeometry_;
    ElementVolumeVariables curElemVolVars_;
    ElementVolumeVariables prevElemVolVars_;
    ElementFluxVariablesCache elemFluxVarsCache_;
    ElementBoundaryTypes elemBcTypes_;

    LocalResidual localResidual_; //!< the local residual evaluating the equations per element
    bool elementIsGhost_; //!< whether the element's partitionType is ghost
};

/*!
 * \ingroup Assembly
 * \ingroup CCDiscretization
 * \brief A base class for all implicit local assemblers
 * \tparam TypeTag the TypeTag
 * \tparam Assembler the assembler type
 */
template<std::size_t id, class TypeTag, class Assembler, class Implementation>
class SubDomainCCLocalAssemblerImplicitBase : public SubDomainCCLocalAssemblerBase<id, TypeTag, Assembler, Implementation>
{
    using ParentType = SubDomainCCLocalAssemblerBase<id, TypeTag, Assembler, Implementation>;
    using LocalResidualValues = typename GET_PROP_TYPE(TypeTag, NumEqVector);
    static constexpr auto domainId = Dune::index_constant<id>();
public:
    using ParentType::ParentType;

    void bindLocalViews()
    {
        // get some references for convenience
        auto& couplingManager = this->couplingManager();
        const auto& element = this->element();
        const auto& curSol = this->curSol()[domainId];
        auto&& fvGeometry = this->fvGeometry();
        auto&& curElemVolVars = this->curElemVolVars();
        auto&& elemFluxVarsCache = this->elemFluxVarsCache();

        // bind the caches
        couplingManager.bindCouplingContext(element, this->assembler());
        fvGeometry.bind(element);
        curElemVolVars.bind(element, fvGeometry, curSol);
        elemFluxVarsCache.bind(element, fvGeometry, curElemVolVars);
        if (!this->assembler().isStationaryProblem())
            this->prevElemVolVars().bindElement(element, fvGeometry, this->assembler().prevSol()[domainId]);
    }

    using ParentType::evalLocalResidual;
    LocalResidualValues evalLocalResidual() const
    { return this->evalLocalResidual(this->curElemVolVars()); }

    using ParentType::evalLocalFluxAndSourceResidual;
    LocalResidualValues evalLocalFluxAndSourceResidual() const
    { return this->evalLocalFluxAndSourceResidual(this->curElemVolVars()); }

    /*!
     * \brief Computes the residual
     * \return The element residual at the current solution.
     */
    LocalResidualValues assembleResidualImpl()
    { return this->elementIsGhost() ? LocalResidualValues(0.0) : evalLocalResidual(); }
};

/*!
 * \ingroup Assembly
 * \ingroup CCDiscretization
 * \brief An assembler for Jacobian and residual contribution per element (cell-centered methods)
 * \tparam TypeTag the TypeTag
 * \tparam DM the differentiation method to residual compute derivatives
 * \tparam implicit if to use an implicit or explicit time discretization
 */
template<std::size_t id, class TypeTag, class Assembler, DiffMethod DM = DiffMethod::numeric, bool implicit = true>
class SubDomainCCLocalAssembler;

/*!
 * \ingroup Assembly
 * \ingroup CCDiscretization
 * \brief Cell-centered scheme local assembler using numeric differentiation and implicit time discretization
 */
template<std::size_t id, class TypeTag, class Assembler>
class SubDomainCCLocalAssembler<id, TypeTag, Assembler, DiffMethod::numeric, /*implicit=*/true>
: public SubDomainCCLocalAssemblerImplicitBase<id, TypeTag, Assembler,
            SubDomainCCLocalAssembler<id, TypeTag, Assembler, DiffMethod::numeric, true> >
{
    using ThisType = SubDomainCCLocalAssembler<id, TypeTag, Assembler, DiffMethod::numeric, /*implicit=*/true>;
    using ParentType = SubDomainCCLocalAssemblerImplicitBase<id, TypeTag, Assembler, ThisType>;
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using LocalResidualValues = typename GET_PROP_TYPE(TypeTag, NumEqVector);
    using Element = typename GET_PROP_TYPE(TypeTag, GridView)::template Codim<0>::Entity;

    enum { numEq = GET_PROP_VALUE(TypeTag, NumEq) };
    enum { dim = GET_PROP_TYPE(TypeTag, GridView)::dimension };

    static constexpr bool enableGridFluxVarsCache = GET_PROP_VALUE(TypeTag, EnableGridFluxVariablesCache);
    static constexpr int maxNeighbors = 4*(2*dim);
    static constexpr auto domainIdI = Dune::index_constant<id>();

public:
    using ParentType::ParentType;

    /*!
     * \brief Computes the derivatives with respect to the given element and adds them
     *        to the global matrix.
     *
     * \return The element residual at the current solution.
     */
    template<class JacobianMatrixDiagBlock, class GridVariables>
    LocalResidualValues assembleJacobianAndResidualImpl(JacobianMatrixDiagBlock& A, GridVariables& gridVariables)
    {
        //////////////////////////////////////////////////////////////////////////////////////////////////
        // Calculate derivatives of all dofs in stencil with respect to the dofs in the element. In the //
        // neighboring elements we do so by computing the derivatives of the fluxes which depend on the //
        // actual element. In the actual element we evaluate the derivative of the entire residual.     //
        //////////////////////////////////////////////////////////////////////////////////////////////////

        // get some aliases for convenience
        const auto& element = this->element();
        const auto& fvGeometry = this->fvGeometry();
        const auto& fvGridGeometry = fvGeometry.fvGridGeometry();
        auto&& curElemVolVars = this->curElemVolVars();
        auto&& elemFluxVarsCache = this->elemFluxVarsCache();

        // get stencil informations
        const auto globalI = fvGridGeometry.elementMapper().index(element);
        const auto& connectivityMap = fvGridGeometry.connectivityMap();
        const auto numNeighbors = connectivityMap[globalI].size();

        // container to store the neighboring elements
        Dune::ReservedVector<Element, maxNeighbors+1> neighborElements;
        neighborElements.resize(numNeighbors);

        // assemble the undeflected residual
        using Residuals = ReservedBlockVector<LocalResidualValues, maxNeighbors+1>;
        Residuals origResiduals(numNeighbors + 1); origResiduals = 0.0;
        origResiduals[0] = this->assembleResidualImpl();

        // get the elements in which we need to evaluate the fluxes
        // and calculate these in the undeflected state
        unsigned int j = 1;
        for (const auto& dataJ : connectivityMap[globalI])
        {
            neighborElements[j-1] = fvGridGeometry.element(dataJ.globalJ);
            for (const auto scvfIdx : dataJ.scvfsJ)
                origResiduals[j] += this->evalFluxResidual(neighborElements[j-1], fvGeometry.scvf(scvfIdx));

            ++j;
        }

        // reference to the element's scv (needed later) and corresponding vol vars
        const auto& scv = fvGeometry.scv(globalI);
        auto& curVolVars = ParentType::getVolVarAccess(gridVariables.curGridVolVars(), curElemVolVars, scv);

        // save a copy of the original privars and vol vars in order
        // to restore the original solution after deflection
        const auto& curSol = this->curSol()[domainIdI];
        const auto origPriVars = curSol[globalI];
        const auto origVolVars = curVolVars;

        // element solution container to be deflected
        using ElementSolutionVector = typename GET_PROP_TYPE(TypeTag, ElementSolutionVector);
        ElementSolutionVector elemSol(origPriVars);

        // derivatives in the neighbors with repect to the current elements
        // in index 0 we save the derivative of the element residual with respect to it's own dofs
        Residuals partialDerivs(numNeighbors + 1);

        for (int pvIdx = 0; pvIdx < numEq; ++pvIdx)
        {

            // for ghost elements we assemble a 1.0 where the primary variable and zero everywhere else
            // as we always solve for a delta of the solution with repect to the initial solution this
            // results in a delta of zero for ghosts, we still need to do the neighbor derivatives though
            // so we are not done yet here.
            partialDerivs = 0.0;
            if (this->elementIsGhost()) partialDerivs[0][pvIdx] = 1.0;

            auto evalResiduals = [&](Scalar priVar)
            {
                Residuals partialDerivsTmp(numNeighbors + 1);
                partialDerivsTmp = 0.0;
                // update the volume variables and the flux var cache
                elemSol[0][pvIdx] = priVar;
                this->couplingManager().updateCouplingContext(domainIdI, element, elemSol[0], this->assembler());
                curVolVars.update(elemSol, this->problem(), element, scv);
                if (enableGridFluxVarsCache)
                    gridVariables.gridFluxVarsCache().updateElement(element, fvGeometry, curElemVolVars);
                else
                    elemFluxVarsCache.update(element, fvGeometry, curElemVolVars);

                // calculate the residual with the deflected primary variables
                if (!this->elementIsGhost()) partialDerivsTmp[0] = this->evalLocalResidual();

                // calculate the fluxes in the neighbors with the deflected primary variables
                for (std::size_t k = 0; k < numNeighbors; ++k)
                    for (auto scvfIdx : connectivityMap[globalI][k].scvfsJ)
                        partialDerivsTmp[k+1] += this->evalFluxResidual(neighborElements[k], fvGeometry.scvf(scvfIdx));

                return partialDerivsTmp;
            };

            // derive the residuals numerically
            NumericDifferentiation::partialDerivative(evalResiduals, elemSol[0][pvIdx], partialDerivs, origResiduals);

            // add the current partial derivatives to the global jacobian matrix
            for (int eqIdx = 0; eqIdx < numEq; eqIdx++)
            {
                // the diagonal entries
                A[globalI][globalI][eqIdx][pvIdx] += partialDerivs[0][eqIdx];

                // off-diagonal entries
                j = 1;
                for (const auto& dataJ : connectivityMap[globalI])
                    A[dataJ.globalJ][globalI][eqIdx][pvIdx] += partialDerivs[j++][eqIdx];
            }

            // restore the original state of the scv's volume variables
            curVolVars = origVolVars;

            // restore the current element solution
            elemSol[0][pvIdx] = origPriVars[pvIdx];

            // restore the undeflected state of the coupling context
            this->couplingManager().updateCouplingContext(domainIdI, element, elemSol[0], this->assembler());
        }

        // restore original state of the flux vars cache in case of global caching.
        // This has to be done in order to guarantee that everything is in an undeflected
        // state before the assembly of another element is called. In the case of local caching
        // this is obsolete because the elemFluxVarsCache used here goes out of scope after this.
        // We only have to do this for the last primary variable, for all others the flux var cache
        // is updated with the correct element volume variables before residual evaluations
        if (enableGridFluxVarsCache)
            gridVariables.gridFluxVarsCache().updateElement(element, fvGeometry, curElemVolVars);

        // evaluate additional derivatives that can be caused by coupling conditions
        const auto additionalDofDeps = this->couplingManager().getAdditionalDofDependencies(domainIdI, globalI);
        if (!additionalDofDeps.empty())
            evalAdditionalDerivatives(additionalDofDeps, A, gridVariables);

        // return the original residual
        return origResiduals[0];
    }

    /*!
     * \brief Computes the derivatives with respect to the given element and adds them
     *        to the global matrix.
     *
     * \return The element residual at the current solution.
     */
    template<std::size_t otherId, class JacobianBlock, class GridVariables>
    void assembleJacobianCoupling(Dune::index_constant<otherId> domainIdJ, JacobianBlock& A,
                                  const LocalResidualValues& res, GridVariables& gridVariables)
    {
        ////////////////////////////////////////////////////////////////////////////////////////////////////////
        // Calculate derivatives of all dofs in the element with respect to all dofs in the coupling stencil. //
        ////////////////////////////////////////////////////////////////////////////////////////////////////////

        // get some aliases for convenience
        const auto& element = this->element();
        const auto& fvGeometry = this->fvGeometry();
        const auto& fvGridGeometry = fvGeometry.fvGridGeometry();
        auto&& curElemVolVars = this->curElemVolVars();
        auto&& elemFluxVarsCache = this->elemFluxVarsCache();

        // get stencil informations
        const auto globalI = fvGridGeometry.elementMapper().index(element);
        const auto& stencil = this->couplingManager().couplingStencil(element, domainIdJ);

        for (const auto globalJ : stencil)
        {
            const auto& elementJ = this->assembler().fvGridGeometry(domainIdJ).element(globalJ);

            const auto& curSol = this->curSol()[domainIdJ];
            const auto origPriVarsJ = curSol[globalJ];

            // element solution container to be deflected
            using CoupledDomainTypeTag = typename Assembler::Traits::template SubDomainTypeTag<domainIdJ>;
            using ElementSolutionVector = typename GET_PROP_TYPE(CoupledDomainTypeTag, ElementSolutionVector);
            ElementSolutionVector elemSolJ(origPriVarsJ);

            const auto origResidual = this->couplingManager().evalCouplingResidual(element, fvGeometry, curElemVolVars, this->elemBcTypes(), elemFluxVarsCache, elementJ);

            for (int pvIdx = 0; pvIdx < JacobianBlock::block_type::cols; ++pvIdx)
            {
                auto evalCouplingResidual = [&](Scalar priVar)
                {
                    LocalResidualValues partialDerivTmp;
                    partialDerivTmp = 0.0;

                    // update the volume variables and the flux var cache
                    elemSolJ[0][pvIdx] = priVar;
                    this->couplingManager().updateCouplingContext(domainIdI, elementJ, elemSolJ[0], this->assembler());

                    if (enableGridFluxVarsCache)
                        this->couplingManager().updateSelf(domainIdI, element, fvGeometry, curElemVolVars, gridVariables.gridFluxVarsCache());
                    else
                        this->couplingManager().updateSelf(domainIdI, element, fvGeometry, curElemVolVars, elemFluxVarsCache);

                    // calculate the residual with the deflected coupling neighbor primary variables
                    partialDerivTmp = this->couplingManager().evalCouplingResidual(element, fvGeometry, curElemVolVars, this->elemBcTypes(), elemFluxVarsCache, elementJ);

                    return partialDerivTmp;
                };

                // derive the residuals numerically
                LocalResidualValues partialDeriv = 0.0;
                NumericDifferentiation::partialDerivative(evalCouplingResidual, elemSolJ[0][pvIdx], partialDeriv, origResidual);

                // add the current partial derivatives to the global jacobian matrix
                for (int eqIdx = 0; eqIdx < numEq; eqIdx++)
                    A[globalI][globalJ][eqIdx][pvIdx] += partialDeriv[eqIdx];

                // restore the current element solution
                elemSolJ[0][pvIdx] = origPriVarsJ[pvIdx];

                // restore the undeflected state of the coupling context
                this->couplingManager().updateCouplingContext(domainIdI, elementJ, elemSolJ[0], this->assembler());
            }
        }

        // restore original state of the flux vars cache in case of global caching.
        // This has to be done in order to guarantee that everything is in an undeflected
        // state before the assembly of another element is called. In the case of local caching
        // this is obsolete because the elemFluxVarsCache used here goes out of scope after this.
        // We only have to do this for the last primary variable, for all others the flux var cache
        // is updated with the correct element volume variables before residual evaluations
        if (enableGridFluxVarsCache)
            this->couplingManager().updateSelf(domainIdI, element, fvGeometry, curElemVolVars, gridVariables.gridFluxVarsCache());
        else
            this->couplingManager().updateSelf(domainIdI, element, fvGeometry, curElemVolVars, elemFluxVarsCache);
    }

    template<class JacobianMatrixDiagBlock, class GridVariables>
    void evalAdditionalDerivatives(const std::vector<std::size_t>& additionalDofDependencies,
                                   JacobianMatrixDiagBlock& A, GridVariables& gridVariables)
    {
        const auto& fvGeometry = this->fvGeometry();
        const auto& fvGridGeometry = fvGeometry.fvGridGeometry();
        auto&& curElemVolVars = this->curElemVolVars();
        const auto& element = this->element();
        const auto globalI = fvGridGeometry.elementMapper().index(element);
        const auto& scv = fvGeometry.scv(globalI);

        const auto& curVolVarsI = curElemVolVars[scv];
        auto source = this->localResidual().computeSource(this->problem(), element, fvGeometry, curElemVolVars, scv);
        source *= -scv.volume()*curVolVarsI.extrusionFactor();

        for (const auto globalJ : additionalDofDependencies)
        {
            const auto& scvJ = fvGeometry.scv(globalJ);
            auto& curVolVarsJ = ParentType::getVolVarAccess(gridVariables.curGridVolVars(), curElemVolVars, scvJ);
            const auto& elementJ = fvGridGeometry.element(globalJ);

            // save a copy of the original privars and vol vars in order
            // to restore the original solution after deflection
            const auto& curSol = this->curSol()[domainIdI];
            const auto origPriVarsJ = curSol[globalJ];
            const auto origVolVarsJ = curVolVarsJ;

            // element solution container to be deflected
            using ElementSolutionVector = typename GET_PROP_TYPE(TypeTag, ElementSolutionVector);
            ElementSolutionVector elemSolJ(origPriVarsJ);

            // derivatives with repect to the additional DOF we depend on
            for (int pvIdx = 0; pvIdx < numEq; pvIdx++)
            {
                auto evalResiduals = [&](Scalar priVar)
                {
                    LocalResidualValues partialDerivTmp(0.0);
                    // update the volume variables and the flux var cache
                    elemSolJ[0][pvIdx] = priVar;
                    this->couplingManager().updateCouplingContext(domainIdI, elementJ, elemSolJ[0], this->assembler());
                    curVolVarsJ.update(elemSolJ, this->problem(), elementJ, scvJ);

                    // calculate the residual with the deflected primary variables
                    if (!this->elementIsGhost())
                    {
                        partialDerivTmp = this->localResidual().computeSource(this->problem(), element, fvGeometry, curElemVolVars, scv);
                        partialDerivTmp *= -scv.volume()*curVolVarsI.extrusionFactor();
                    }

                    return partialDerivTmp;
                };

                // derive the residuals numerically
                LocalResidualValues partialDeriv(0.0);
                NumericDifferentiation::partialDerivative(evalResiduals, elemSolJ[0][pvIdx], partialDeriv, source);

                // add the current partial derivatives to the global jacobian matrix
                for (int eqIdx = 0; eqIdx < numEq; eqIdx++)
                    A[globalI][globalJ][eqIdx][pvIdx] += partialDeriv[eqIdx];

                // restore the original state of the dofs privars and the volume variables
                elemSolJ[0][pvIdx] = origPriVarsJ[pvIdx];
                curVolVarsJ = origVolVarsJ;

                // restore the undeflected state of the coupling context
                this->couplingManager().updateCouplingContext(domainIdI, elementJ, elemSolJ[0], this->assembler());
            }
        }
    }
};


/*!
 * \ingroup Assembly
 * \ingroup CCDiscretization
 * \brief Cell-centered scheme local assembler using numeric differentiation and implicit time discretization
 */
template<std::size_t id, class TypeTag, class Assembler>
class SubDomainCCLocalAssembler<id, TypeTag, Assembler, DiffMethod::analytic, /*implicit=*/true>
: public SubDomainCCLocalAssemblerImplicitBase<id, TypeTag, Assembler,
            SubDomainCCLocalAssembler<id, TypeTag, Assembler, DiffMethod::analytic, true> >
{
    using ThisType = SubDomainCCLocalAssembler<id, TypeTag, Assembler, DiffMethod::analytic, /*implicit=*/true>;
    using ParentType = SubDomainCCLocalAssemblerImplicitBase<id, TypeTag, Assembler, ThisType>;
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using LocalResidualValues = typename GET_PROP_TYPE(TypeTag, NumEqVector);
    using Element = typename GET_PROP_TYPE(TypeTag, GridView)::template Codim<0>::Entity;

    enum { numEq = GET_PROP_VALUE(TypeTag, NumEq) };
    enum { dim = GET_PROP_TYPE(TypeTag, GridView)::dimension };

    static constexpr bool enableGridFluxVarsCache = GET_PROP_VALUE(TypeTag, EnableGridFluxVariablesCache);
    static constexpr int maxNeighbors = 4*(2*dim);
    static constexpr auto domainIdI = Dune::index_constant<id>();

public:
    using ParentType::ParentType;

    /*!
     * \brief Computes the derivatives with respect to the given element and adds them
     *        to the global matrix.
     *
     * \return The element residual at the current solution.
     */
    template<class JacobianMatrixDiagBlock, class GridVariables>
    LocalResidualValues assembleJacobianAndResidualImpl(JacobianMatrixDiagBlock& A, GridVariables& gridVariables)
    {
        // get some aliases for convenience
        const auto& problem = this->problem();
        const auto& element = this->element();
        const auto& fvGeometry = this->fvGeometry();
        const auto& curElemVolVars = this->curElemVolVars();
        const auto& elemFluxVarsCache = this->elemFluxVarsCache();

        // get reference to the element's current vol vars
        const auto globalI = this->assembler().fvGridGeometry(domainIdI).elementMapper().index(element);
        const auto& scv = fvGeometry.scv(globalI);
        const auto& volVars = curElemVolVars[scv];

        // if the problem is instationary, add derivative of storage term
        if (!this->assembler().isStationaryProblem())
            this->localResidual().addStorageDerivatives(A[globalI][globalI], problem, element, fvGeometry, volVars, scv);

        // add source term derivatives
        this->localResidual().addSourceDerivatives(A[globalI][globalI], problem, element, fvGeometry, volVars, scv);

        // add flux derivatives for each scvf
        for (const auto& scvf : scvfs(fvGeometry))
        {
            // inner faces
            if (!scvf.boundary())
                this->localResidual().addFluxDerivatives(A[globalI], problem, element, fvGeometry, curElemVolVars, elemFluxVarsCache, scvf);

            // boundary faces
            else
            {
                const auto& bcTypes = problem.boundaryTypes(element, scvf);

                // add Dirichlet boundary flux derivatives
                if (bcTypes.hasDirichlet() && !bcTypes.hasNeumann())
                    this->localResidual().addCCDirichletFluxDerivatives(A[globalI], problem, element, fvGeometry, curElemVolVars, elemFluxVarsCache, scvf);

                // add Robin ("solution dependent Neumann") boundary flux derivatives
                else if (bcTypes.hasNeumann() && !bcTypes.hasDirichlet())
                    this->localResidual().addRobinFluxDerivatives(A[globalI], problem, element, fvGeometry, curElemVolVars, elemFluxVarsCache, scvf);

                else
                    DUNE_THROW(Dune::NotImplemented, "Mixed boundary conditions. Use pure boundary conditions by converting Dirichlet BCs to Robin BCs");
            }
        }

        // return element residual
        return this->assembleResidualImpl();
    }

    /*!
     * \brief Computes the derivatives with respect to the given element and adds them
     *        to the global matrix.
     *
     * \return The element residual at the current solution.
     */
    template<std::size_t otherId, class JacobianBlock, class GridVariables>
    void assembleJacobianCoupling(Dune::index_constant<otherId> domainIdJ, JacobianBlock& A,
                                  const LocalResidualValues& res, GridVariables& gridVariables)
    {
        ////////////////////////////////////////////////////////////////////////////////////////////////////////
        // Calculate derivatives of all dofs in the element with respect to all dofs in the coupling stencil. //
        ////////////////////////////////////////////////////////////////////////////////////////////////////////

        // get some aliases for convenience
        const auto& element = this->element();
        const auto& fvGeometry = this->fvGeometry();
        const auto& fvGridGeometry = fvGeometry.fvGridGeometry();
        auto&& curElemVolVars = this->curElemVolVars();
        // auto&& elemFluxVarsCache = this->elemFluxVarsCache();

        // get stencil informations
        const auto globalI = fvGridGeometry.elementMapper().index(element);
        const auto& stencil = this->couplingManager().couplingStencil(element, domainIdJ);

        for (const auto globalJ : stencil)
        {
            const auto& elementJ = this->assembler().fvGridGeometry(domainIdJ).element(globalJ);
            this->couplingManager().addCouplingDerivatives(domainIdI, A[globalI][globalJ], element, fvGeometry, curElemVolVars, elementJ);
        }
    }
};

} // end namespace Dumux

#endif
