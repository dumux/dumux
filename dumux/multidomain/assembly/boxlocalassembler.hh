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
 *        for the box scheme in the context of multidomain models.
 */
#ifndef DUMUX_MULTIDOMAIN_BOX_LOCAL_ASSEMBLER_HH
#define DUMUX_MULTIDOMAIN_BOX_LOCAL_ASSEMBLER_HH

#include <cassert>
#include <memory>

#include <dune/grid/common/gridenums.hh>

#include <dumux/assembly/diffmethod.hh>
#include <dumux/assembly/numericepsilon.hh>

#include <dumux/timestepping/multistagetimestepper.hh>

namespace Dumux {

/*!
 * \ingroup Assembly
 * \brief An assembler for Jacobian and residual contribution
 *        per element for the box scheme.
 * \tparam Assembler The grid-wide assembler type
 */
template<std::size_t id, class Assembler, DiffMethod diffMethod>
class SubDomainBoxLocalAssembler
{
    using CouplingManager = typename Assembler::CouplingManager;
    using GridVariables = typename Assembler::template SubDomainGridVariables<id>;
    using GridGeometry = typename Assembler::template SubDomainGridGeometry<id>;
    using LocalOperator = typename Assembler::template SubDomainLocalOperator<id>;

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
    explicit SubDomainBoxLocalAssembler(const Element& element,
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
    }

    /*!
     * \brief Constructor for instationary problems.
     * \note Using this constructor, we assemble one stage within
     *       a time integration step using multi-stage methods.
     */
    explicit SubDomainBoxLocalAssembler(const Element& element,
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
    {}

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
        for (const auto& scv : scvs(fvGeometry()))
            res[scv.dofIndex()] += residual[scv.localDofIndex()];

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

        auto applyDirichlet = [&] (const auto& scvI,
                                   const auto& dirichletValues,
                                   const auto eqIdx,
                                   const auto pvIdx)
        {
            res[scvI.dofIndex()][eqIdx] = elemVariables().elemVolVars()[scvI].priVars()[pvIdx] - dirichletValues[pvIdx];

            // TODO: Skip this for explicit schemes?
            auto& jac = jacRow[domainId];
            auto& row = jac[scvI.dofIndex()];
            for (auto col = row.begin(); col != row.end(); ++col)
                row[col.index()][eqIdx] = 0.0;

            jac[scvI.dofIndex()][scvI.dofIndex()][eqIdx][pvIdx] = 1.0;

            // TODO: Periodicity
            // if (fvGeometry().gridGeometry().dofOnPeriodicBoundary(scvI.dofIndex())) {}
        };

        enforceDirichletConstraints(applyDirichlet);
    }

    /*!
     * \brief Computes the derivatives with respect to the given element and adds them
     *        to the global matrix.
     */
    template<class JacobianBlock>
    void assembleJacobian(JacobianBlock& jac)
    {
        assembleJacobianAndResidual_(jac);

        auto applyDirichlet = [&] (const auto& scvI,
                                   const auto& dirichletValues,
                                   const auto eqIdx,
                                   const auto pvIdx)
        {
            // TODO: Skip this for explicit schemes?
            auto& row = jac[scvI.dofIndex()];
            for (auto col = row.begin(); col != row.end(); ++col)
                row[col.index()][eqIdx] = 0.0;

            jac[scvI.dofIndex()][scvI.dofIndex()][eqIdx][pvIdx] = 1.0;
        };

        enforceDirichletConstraints(applyDirichlet);
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

        auto applyDirichlet = [&] (const auto& scvI,
                                   const auto& dirichletValues,
                                   const auto eqIdx,
                                   const auto pvIdx)
        {
            res[scvI.dofIndex()][eqIdx] = elemVariables().elemVolVars()[scvI].priVars()[pvIdx] - dirichletValues[pvIdx];
        };

        enforceDirichletConstraints(applyDirichlet);
    }

    //! Enforce Dirichlet constraints
    template<typename ApplyFunction>
    void enforceDirichletConstraints(const ApplyFunction& applyDirichlet)
    {
        // enforce Dirichlet boundary conditions
        evalDirichletBoundaries(applyDirichlet);
        // TODO: take care of internal Dirichlet constraints (if enabled)
        // enforceInternalDirichletConstraints(applyDirichlet);
    }

    /*!
     * \brief Evaluates Dirichlet boundaries
     */
    template< typename ApplyDirichletFunctionType >
    void evalDirichletBoundaries(ApplyDirichletFunctionType applyDirichlet)
    {
        for (const auto& scvI : scvs(fvGeometry()))
        {
            if (!fvGeometry().gridGeometry().dofOnBoundary(scvI.dofIndex()))
                continue;

            const auto bcTypes = problem().boundaryTypes(element(), scvI);
            if (bcTypes.hasDirichlet())
            {
                const auto dirichletValues = problem().dirichlet(element(), scvI);

                // set the Dirichlet conditions in residual and jacobian
                for (int eqIdx = 0; eqIdx < numEq; ++eqIdx)
                {
                    if (bcTypes.isDirichlet(eqIdx))
                    {
                        const auto pvIdx = bcTypes.eqToDirichletIndex(eqIdx);
                        assert(0 <= pvIdx && pvIdx < numEq);
                        applyDirichlet(scvI, dirichletValues, eqIdx, pvIdx);
                    }
                }
            }
        }
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
    ElementResidualVector assembleJacobianAndResidual_(JacobianBlock& A)
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
    ElementResidualVector assembleJacobianAndResidualNumeric_(JacobianBlock& A)
    {
        // get the variables of the current stage
        auto& curVariables = elemVariables_();
        auto& curElemVolVars = curVariables.elemVolVars();
        const auto& x = curVariables.gridVariables().dofs();

        const auto origResiduals = evalLocalResidual();
        const auto origElemSol = elementSolution(element(), x, fvGeometry().gridGeometry());
        auto elemSol = origElemSol;

        //////////////////////////////////////////////////////////////////////////////////////////////
        // Calculate derivatives of the residual of all dofs in element with respect to themselves. //
        //////////////////////////////////////////////////////////////////////////////////////////////

        ElementResidualVector partialDerivs(fvGeometry().numScv());
        for (const auto& scvI : scvs(fvGeometry()))
        {
            // dof index and corresponding actual pri vars
            const auto globalI = scvI.dofIndex();
            const auto localI = scvI.localDofIndex();

            const auto origCurVolVars = curElemVolVars[scvI];
            auto& curVolVars = curElemVolVars[scvI];

            // calculate derivatives w.r.t to the privars at the dof at hand
            for (int pvIdx = 0; pvIdx < numEq; pvIdx++)
            {
                partialDerivs = 0.0;
                auto evalResiduals = [&](Scalar priVar)
                {
                    // update the volume variables and compute element residual
                    elemSol[localI][pvIdx] = priVar;
                    cm_->updateCouplingContext(domainId, *this, domainId, scvI.dofIndex(), elemSol[localI], pvIdx);
                    curVolVars.update(elemSol, problem(), element(), scvI);
                    return evalLocalResidual();
                };

                // derive the residuals numerically
                static const NumericEpsilon<Scalar, numEq> eps_{problem().paramGroup()};
                static const int numDiffMethod = getParamFromGroup<int>(problem().paramGroup(), "Assembly.NumericDifferenceMethod");
                NumericDifferentiation::partialDerivative(evalResiduals, elemSol[localI][pvIdx], partialDerivs,
                                                          origResiduals, eps_(elemSol[localI][pvIdx], pvIdx),
                                                          numDiffMethod);

                // update the global stiffness matrix with the current partial derivatives
                for (const auto& scvJ : scvs(fvGeometry()))
                {
                    const auto globalJ = scvJ.dofIndex();
                    const auto localJ = scvJ.localDofIndex();

                    for (int eqIdx = 0; eqIdx < numEq; eqIdx++)
                    {
                        // A[i][col][eqIdx][pvIdx] is the rate of change of the
                        // the residual of equation 'eqIdx' at dof 'i'
                        // depending on the primary variable 'pvIdx' at dof 'col'
                        A[globalJ][globalI][eqIdx][pvIdx] += partialDerivs[localJ][eqIdx];
                    }
                }

                // restore the original element solution & volume variables
                elemSol[localI][pvIdx] = origElemSol[localI][pvIdx];
                curVolVars = origCurVolVars;
                cm_->updateCouplingContext(domainId, *this, domainId, globalI, elemSol[localI], pvIdx);

                // TODO additional dof dependencies
            }
        }

        return origResiduals;
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

                // update the global stiffness matrix with the current partial derivatives
                for (auto&& scv : scvs(fvGeometry()))
                {
                    for (int eqIdx = 0; eqIdx < numEq; eqIdx++)
                    {
                        // A[i][col][eqIdx][pvIdx] is the rate of change of
                        // the residual of equation 'eqIdx' at dof 'i'
                        // depending on the primary variable 'pvIdx' at dof
                        // 'col'.

                        // If the dof is coupled by a Dirichlet condition,
                        // set the derived value only once (i.e. overwrite existing values).
                        // For other dofs, add the contribution of the partial derivative.
                        if (fvGeometry().gridGeometry().dofOnBoundary(scv.dofIndex()))
                        {
                            const auto bcTypes = problem().boundaryTypes(element(), scv);
                            if (bcTypes.isCouplingDirichlet(eqIdx))
                                A[scv.dofIndex()][globalJ][eqIdx][pvIdx] = partialDerivs[scv.localDofIndex()][eqIdx];
                            else if (bcTypes.isDirichlet(eqIdx))
                                A[scv.dofIndex()][globalJ][eqIdx][pvIdx] = 0.0;
                        }
                        else
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
