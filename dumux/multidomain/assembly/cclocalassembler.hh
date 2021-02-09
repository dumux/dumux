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
 * \ingroup Assembly
 * \brief An assembler for Jacobian and residual contribution per element
 *        for cell-centered schemes in the context of multidomain models.
 */
#ifndef DUMUX_MULTIDOMAIN_CC_LOCAL_ASSEMBLER_HH
#define DUMUX_MULTIDOMAIN_CC_LOCAL_ASSEMBLER_HH

#include <cassert>
#include <memory>

#include <dune/grid/common/gridenums.hh>
#include <dune/common/reservedvector.hh>

#include <dumux/common/reservedblockvector.hh>
#include <dumux/assembly/diffmethod.hh>
#include <dumux/assembly/numericepsilon.hh>

#include <dumux/timestepping/multistagetimestepper.hh>

namespace Dumux {

/*!
 * \ingroup Assembly
 * \brief An assembler for Jacobian and residual contribution
 *        per element for cell-centered schemes.
 * \tparam Assembler The grid-wide assembler type
 */
template<std::size_t id, class Assembler, DiffMethod diffMethod>
class SubDomainCCLocalAssembler
{
    using CouplingManager = typename Assembler::CouplingManager;
    using GridVariables = typename Assembler::template SubDomainGridVariables<id>;
    using GridGeometry = typename Assembler::template SubDomainGridGeometry<id>;
    using LocalOperator = typename Assembler::template SubDomainLocalOperator<id>;

    using Operators = typename LocalOperator::Operators;
    using NumEqVector = typename Operators::NumEqVector;

    using FVElementGeometry = typename GridGeometry::LocalView;
    using ElementVariables = typename GridVariables::LocalView;
    using PrimaryVariables = typename GridVariables::PrimaryVariables;
    using Scalar = typename GridVariables::Scalar;

    using GridView = typename GridGeometry::GridView;
    using Element = typename GridView::template Codim<0>::Entity;
    using Context = typename CouplingManager::template CouplingContext<id>;

    using JacobianMatrix = typename Assembler::JacobianMatrix;
    using ResidualVector = typename Assembler::SolutionVector;

    // TODO: remove this assert? Get numEq from somewhere else?
    static constexpr int numEq = PrimaryVariables::size();
    static constexpr int maxElementStencilSize = GridGeometry::maxElementStencilSize;
    static_assert(diffMethod == DiffMethod::numeric, "Analytical assembly not implemented");
    // static_assert(numEq == JacobianMatrix::block_type::cols, "Matrix block size doesn't match privars size");

    //! the parameters of a stage in time integration
    using StageParams = MultiStageParams<Scalar>;

public:
    static constexpr auto domainId = typename Dune::index_constant<id>();
    using ElementResidualVector = typename LocalOperator::ElementResidualVector;

    /*!
     * \brief Constructor for stationary problems.
     */
    explicit SubDomainCCLocalAssembler(const Element& element,
                                       const FVElementGeometry& fvGeometry,
                                       std::vector<std::shared_ptr<Context>>& contexts,
                                       std::vector<ElementVariables>& elemVars,
                                       std::shared_ptr<CouplingManager> cm)
    : element_(element)
    , fvGeometry_(fvGeometry)
    , contexts_(contexts)
    , elementVariables_(elemVars)
    , elementIsGhost_(element.partitionType() == Dune::GhostEntity)
    , stageParams_(nullptr)
    , cm_(cm)
    {
        assert(elemVars.size() == 1);
        assert(fvGeometry_.numScv() == 1);
    }

    /*!
     * \brief Constructor for instationary problems.
     * \note Using this constructor, we assemble one stage within
     *       a time integration step using multi-stage methods.
     */
    explicit SubDomainCCLocalAssembler(const Element& element,
                                       const FVElementGeometry& fvGeometry,
                                       std::vector<std::shared_ptr<Context>>& contexts,
                                       std::vector<ElementVariables>& elemVars,
                                       std::shared_ptr<const StageParams> stageParams,
                                       std::shared_ptr<CouplingManager> cm)
    : element_(element)
    , fvGeometry_(fvGeometry)
    , contexts_(contexts)
    , elementVariables_(elemVars)
    , elementIsGhost_(element.partitionType() == Dune::GhostEntity)
    , stageParams_(stageParams)
    , cm_(cm)
    { assert(fvGeometry_.numScv() == 1); }

    /*!
     * \brief Computes the derivatives with respect to the given element and adds
     *        them to the global matrix. The element residual is written into the
     *        right hand side.
     */
    template<class JacRow, class SubRes>
    void assembleJacobianAndResidual(JacRow& jacRow, SubRes& res)
    {
        // TODO: treat ghost elements for parallelism

        // assemble the diagonal block
        const auto residual = assembleJacobianAndResidual_(jacRow[domainId]);
        res[scvs(fvGeometry()).begin()->dofIndex()] += residual;

        // assemble the coupling blocks
        using namespace Dune::Hybrid;
        forEach(integralRange(Dune::Hybrid::size(jacRow)), [&](auto&& j)
        {
            if constexpr (j != id)
            {
                static constexpr auto domainJ = Dune::index_constant<j>();
                assembleJacobianCoupling_(j, jacRow[domainJ]);
            }
        });
    }

    /*!
     * \brief Computes the derivatives with respect to the given element and adds them
     *        to the global matrix.
     */
    template<class JacobianBlock>
    void assembleJacobian(JacobianBlock& jac)
    {
        assembleJacobianAndResidual_(jac);
    }

    /*!
     * \brief Assemble the residual only
     */
    template<class SubRes>
    void assembleResidual(SubRes& res)
    {
        const auto residual = evalLocalResidual();
        for (const auto& scv : scvs(fvGeometry()))
            res[scv.dofIndex()] += residual[scv.localDofIndex()];
    }

    /*!
     * \brief Evaluate the complete local residual for the current element.
     */
    ElementResidualVector evalLocalResidual() const
    {
        if (isStationary())
        {
            LocalOperator localOperator(element(), fvGeometry(), elemVariables());
            return elementIsGhost_ ? localOperator.getEmptyResidual()
                                   : localOperator.evalFluxesAndSources();
        }
        else
        {
            ElementResidualVector residual(fvGeometry().numScv());
            residual = 0.0;

            for (std::size_t k = 0; k < stageParams_->size(); ++k)
            {
                cm_->setCouplingContext(contexts_[k]);
                LocalOperator localOperator(element(), fvGeometry(), elementVariables_[k]);

                if (!stageParams_->skipTemporal(k))
                    residual.axpy(stageParams_->temporalWeight(k), localOperator.evalStorage());
                if (!stageParams_->skipSpatial(k))
                    residual.axpy(stageParams_->spatialWeight(k), localOperator.evalFluxesAndSources());
            }

            return residual;
        }
    }

    //! Return references to the local views
    const Element& element() const { return element_; }
    const FVElementGeometry& fvGeometry() const { return fvGeometry_; }
    const ElementVariables& elemVariables() const { return elementVariables_.back(); }

    //! Still needed by some coupling managers
    //! \todo TODO: Double check if still necessary after complete
    //!             integration of new time integration schemes.
    const auto& curElemVolVars() const { return elemVariables().elemVolVars(); }
    const auto& elemFluxVarsCache() const { return elemVariables().elemFluxVarsCache(); }

    //! Returns if a stationary problem is assembled
    bool isStationary() const { return !stageParams_; }

    // TODO: This is currently required by some coupling managers.
    //       We should get rid of this coupling probabaly.
    static constexpr bool isImplicit() { return /*this is buggy*/ true; }

    //! Return a reference to the underlying problem
    //! TODO: Should grid vars return problem directly!?
    const auto& problem() const
    { return elemVariables().gridVariables().gridVolVars().problem(); }

private:
    ElementVariables& elemVariables_()
    { return elementVariables_.back(); }

    /*!
     * \brief Computes the derivatives with respect to the dofs of the given
     *        element and adds them to the global matrix.
     * \return The element residual at the current solution.
     */
    template<class JacobianBlock>
    NumEqVector assembleJacobianAndResidual_(JacobianBlock& A)
    {
        if constexpr (diffMethod == DiffMethod::numeric)
            return assembleJacobianAndResidualNumeric_(A);
        else
            DUNE_THROW(Dune::NotImplemented, "Analytic assembler for box");
    }

    /*!
     * \brief Computes the derivatives by means of numeric differentiation
     *        and adds them to the global matrix.
     * \return The element residual at the current solution.
     * \note This specialization is for the box scheme with numeric differentiation
     */
    template<class JacobianBlock>
    NumEqVector assembleJacobianAndResidualNumeric_(JacobianBlock& A)
    {
        using Problem = std::decay_t<decltype(problem())>;

        // get the variables of the current stage
        auto& curVariables = elemVariables_();
        auto& curElemVolVars = curVariables.elemVolVars();
        auto& curElemFluxVarsCache = curVariables.elemFluxVarsCache();
        const auto& x = curVariables.gridVariables().dofs();

        // get stencil informations
        const auto& gridGeometry = fvGeometry().gridGeometry();
        const auto& connectivityMap = gridGeometry.connectivityMap();
        const auto globalI = gridGeometry.elementMapper().index(element_);
        const auto numNeighbors = connectivityMap[globalI].size();

        // container to store the neighboring elements
        Dune::ReservedVector<Element, maxElementStencilSize> neighborElements;
        neighborElements.resize(numNeighbors);

        // assemble the undeflected residual
        using Residuals = ReservedBlockVector<NumEqVector, maxElementStencilSize>;
        Residuals origResiduals(numNeighbors + 1); origResiduals = 0.0;
        origResiduals[0] = evalLocalResidual()[0];

        // lambda for convenient evaluation of the fluxes across scvfs in the neighbors
        // if the neighbor is a ghost we don't want to add anything to their residual
        // so we return 0 and omit computing the flux
        auto evalNeighborFlux = [&] (const auto& neighbor, const auto& scvf)
        {
            if (neighbor.partitionType() == Dune::GhostEntity)
                return NumEqVector(0.0);
            else
                return Operators::flux(problem(), neighbor, fvGeometry(),
                                       curElemVolVars, curElemFluxVarsCache, scvf);
        };

        // get the elements in which we need to evaluate the fluxes
        // and calculate these in the undeflected state
        unsigned int j = 1;
        for (const auto& dataJ : connectivityMap[globalI])
        {
            neighborElements[j-1] = gridGeometry.element(dataJ.globalJ);
            for (const auto scvfIdx : dataJ.scvfsJ)
                origResiduals[j] += evalNeighborFlux(neighborElements[j-1], fvGeometry().scvf(scvfIdx));
            ++j;
        }

        // reference to the element's scv (needed later) and corresponding vol vars
        // TODO: Support the case of caching?!
        const auto& scv = fvGeometry().scv(globalI);
        auto& curVolVars = curElemVolVars[scv];

        // save a copy of the original privars and vol vars in order
        // to restore the original solution after deflection
        const auto origPriVars = x[globalI];
        const auto origVolVars = curVolVars;

        // element solution container to be deflected
        auto elemSol = elementSolution(element_, x, gridGeometry);

        // derivatives in the neighbors with repect to the current elements
        // in index 0 we save the derivative of the element residual with respect to it's own dofs
        Residuals partialDerivs(numNeighbors + 1);
        for (int pvIdx = 0; pvIdx < numEq; ++pvIdx)
        {
            partialDerivs = 0.0;

            auto evalResiduals = [&](Scalar priVar)
            {
                Residuals partialDerivsTmp(numNeighbors + 1);
                partialDerivsTmp = 0.0;

                // update the volume variables and the flux var cache
                elemSol[0][pvIdx] = priVar;
                curVolVars.update(elemSol, problem(), element_, scv);
                curElemFluxVarsCache.update(element_, fvGeometry(), curElemVolVars);
                cm_->updateCouplingContext(domainId, *this, domainId, globalI, elemSol[0], pvIdx);

                // TODO: UPDATE GLOBAL FLUX VARS CACHE HERE?

                // calculate the residual with the deflected primary variables
                partialDerivsTmp[0] = evalLocalResidual()[0];

                // calculate the fluxes in the neighbors with the deflected primary variables
                for (std::size_t k = 0; k < numNeighbors; ++k)
                    for (auto scvfIdx : connectivityMap[globalI][k].scvfsJ)
                        partialDerivsTmp[k+1] += evalNeighborFlux(neighborElements[k], fvGeometry().scvf(scvfIdx));

                return partialDerivsTmp;
            };

            // derive the residuals numerically
            static const NumericEpsilon<Scalar, numEq> eps_{problem().paramGroup()};
            static const int numDiffMethod = getParamFromGroup<int>(problem().paramGroup(), "Assembly.NumericDifferenceMethod");
            NumericDifferentiation::partialDerivative(evalResiduals, elemSol[0][pvIdx], partialDerivs, origResiduals,
                                                      eps_(elemSol[0][pvIdx], pvIdx), numDiffMethod);

            // Correct derivative for ghost elements, i.e. set a 1 for the derivative w.r.t. the
            // current primary variable and a 0 elsewhere. As we always solve for a delta of the
            // solution with repect to the initial one, this results in a delta of 0 for ghosts.
            if (elementIsGhost_)
            {
                partialDerivs[0] = 0.0;
                partialDerivs[0][pvIdx] = 1.0;
            }

            // add the current partial derivatives to the global jacobian matrix
            // no special treatment is needed if globalJ is a ghost because then derivatives have been assembled to 0 above
            if constexpr (Problem::enableInternalDirichletConstraints())
            {
                // check if own element has internal Dirichlet constraint
                const auto internalDirichletConstraintsOwnElement = this->problem().hasInternalDirichletConstraint(element_, scv);
                const auto dirichletValues = this->problem().internalDirichlet(element_, scv);

                for (int eqIdx = 0; eqIdx < numEq; ++eqIdx)
                {
                    if (internalDirichletConstraintsOwnElement[eqIdx])
                    {
                        origResiduals[0][eqIdx] = origVolVars.priVars()[eqIdx] - dirichletValues[eqIdx];
                        A[globalI][globalI][eqIdx][pvIdx] = (eqIdx == pvIdx) ? 1.0 : 0.0;
                    }
                    else
                        A[globalI][globalI][eqIdx][pvIdx] += partialDerivs[0][eqIdx];
                }

                // off-diagonal entries
                j = 1;
                for (const auto& dataJ : connectivityMap[globalI])
                {
                    const auto& neighborElement = neighborElements[j-1];
                    const auto& neighborScv = fvGeometry().scv(dataJ.globalJ);
                    const auto internalDirichletConstraintsNeighbor = problem().hasInternalDirichletConstraint(neighborElement, neighborScv);

                    for (int eqIdx = 0; eqIdx < numEq; ++eqIdx)
                    {
                        if (internalDirichletConstraintsNeighbor[eqIdx])
                            A[dataJ.globalJ][globalI][eqIdx][pvIdx] = 0.0;
                        else
                            A[dataJ.globalJ][globalI][eqIdx][pvIdx] += partialDerivs[j][eqIdx];
                    }

                    ++j;
                }
            }
            else // no internal Dirichlet constraints specified
            {
                for (int eqIdx = 0; eqIdx < numEq; eqIdx++)
                {
                    // the diagonal entries
                    A[globalI][globalI][eqIdx][pvIdx] += partialDerivs[0][eqIdx];

                    // off-diagonal entries
                    j = 1;
                    for (const auto& dataJ : connectivityMap[globalI])
                        A[dataJ.globalJ][globalI][eqIdx][pvIdx] += partialDerivs[j++][eqIdx];
                }
            }

            // restore the original state of the scv's volume variables
            curVolVars = origVolVars;

            // restore the current element solution
            elemSol[0][pvIdx] = origPriVars[pvIdx];

            // restore the undeflected state of the coupling context
            cm_->updateCouplingContext(domainId, *this, domainId, globalI, elemSol[0], pvIdx);
        }

        // restore original state of the flux vars cache in case of global caching.
        // This has to be done in order to guarantee that everything is in an undeflected
        // state before the assembly of another element is called. In the case of local caching
        // this is obsolete because the elemFluxVarsCache used here goes out of scope after this.
        // We only have to do this for the last primary variable, for all others the flux var cache
        // is updated with the correct element volume variables before residual evaluations
        // TODO: RESET GLOBAL FLUX VARS CACHE HERE

        // return the original residual
        return origResiduals[0];
    }

    /*!
     * \brief Computes the derivatives with respect to the given element and adds them
     *        to the global matrix.
     */
    template<std::size_t otherId, class JacobianBlock>
    void assembleJacobianCoupling_(Dune::index_constant<otherId> domainJ, JacobianBlock& A)
    {
        if constexpr (diffMethod == DiffMethod::numeric)
            return assembleJacobianCouplingNumeric_(domainJ, A);
        else
            DUNE_THROW(Dune::NotImplemented, "Analytic assembler for box");
    }

    /*!
     * \brief Computes the derivatives with respect to the given element by
     *        means of numeric differentiation and adds them to the global matrix.
     */
    template<std::size_t otherId, class JacobianBlock>
    void assembleJacobianCouplingNumeric_(Dune::index_constant<otherId> domainJ, JacobianBlock& A)
    {
        using Problem = std::decay_t<decltype(problem())>;

        const auto& stencil = cm_->couplingStencil(domainId, element(), domainJ);
        const auto& dofs = cm_->dofs(domainJ);

        auto& elemVolVars = elemVariables_().elemVolVars();
        auto& elemFluxVarsCache = elemVariables_().elemFluxVarsCache();

        // TODO: How to handle the case of caching?
        auto updateCoupledVariables = [&] ()
        {
            // Update ourself after the context has been modified. Depending on the
            // type of caching, other objects might have to be updated. All ifs can be optimized away.
            cm_->updateCoupledVariables(domainId, *this, elemVolVars, elemFluxVarsCache);
        };

        for (const auto globalJ : stencil)
        {
            // undeflected privars and privars to be deflected
            const auto origPriVarsJ = dofs[globalJ];
            auto priVarsJ = origPriVarsJ;

            // the undeflected coupling residual
            const auto origResidual = cm_->evalCouplingResidual(domainId, *this, domainJ, globalJ);

            for (int pvIdx = 0; pvIdx < JacobianBlock::block_type::cols; ++pvIdx)
            {
                auto evalCouplingResidual = [&](Scalar priVar)
                {
                    priVarsJ[pvIdx] = priVar;
                    cm_->updateCouplingContext(domainId, *this, domainJ, globalJ, priVarsJ, pvIdx);
                    updateCoupledVariables();
                    return cm_->evalCouplingResidual(domainId, *this, domainJ, globalJ);
                };

                // derive the residuals numerically
                decltype(evalCouplingResidual(0)) partialDerivs(element().subEntities(GridView::dimension));

                const auto& paramGroup = cm_->problem(domainJ).paramGroup();
                static const int numDiffMethod = getParamFromGroup<int>(paramGroup, "Assembly.NumericDifferenceMethod");
                static const auto epsCoupl = cm_->numericEpsilon(domainJ, paramGroup);

                NumericDifferentiation::partialDerivative(evalCouplingResidual, origPriVarsJ[pvIdx], partialDerivs, origResidual,
                                                          epsCoupl(origPriVarsJ[pvIdx], pvIdx), numDiffMethod);

                if constexpr (Problem::enableInternalDirichletConstraints())
                    DUNE_THROW(Dune::NotImplemented, "Internal Dirichlet");

                // update the global stiffness matrix with the current partial derivatives
                for (auto&& scv : scvs(fvGeometry()))
                {
                    for (int eqIdx = 0; eqIdx < numEq; eqIdx++)
                    {

                        // A[i][col][eqIdx][pvIdx] is the rate of change of
                        // the residual of equation 'eqIdx' at dof 'i'
                        // depending on the primary variable 'pvIdx' at dof
                        // 'col'.
                        A[scv.dofIndex()][globalJ][eqIdx][pvIdx] += partialDerivs[scv.localDofIndex()][eqIdx];
                    }
                }

                // restore the current element solution
                priVarsJ[pvIdx] = origPriVarsJ[pvIdx];

                // restore the undeflected state of the coupling context
                cm_->updateCouplingContext(domainId, *this, domainJ, globalJ, priVarsJ, pvIdx);
            }

            // Restore original state of the flux vars cache and/or vol vars.
            // This has to be done in case they depend on variables of domainJ before
            // we continue with the numeric derivative w.r.t the next globalJ. Otherwise,
            // the next "origResidual" will be incorrect.
            updateCoupledVariables();
        }
    }

private:
    const Element& element_;
    const FVElementGeometry& fvGeometry_;
    std::vector<std::shared_ptr<Context>>& contexts_;
    std::vector<ElementVariables>& elementVariables_;

    bool elementIsGhost_;
    std::shared_ptr<const StageParams> stageParams_;
    std::shared_ptr<CouplingManager> cm_;
};

} // end namespace Dumux

#endif
