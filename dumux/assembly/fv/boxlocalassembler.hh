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
 * \brief An assembler for Jacobian and residual contribution
 *        per element for the box scheme.
 */
#ifndef DUMUX_BOX_LOCAL_ASSEMBLER_HH
#define DUMUX_BOX_LOCAL_ASSEMBLER_HH

#include <cassert>
#include <memory>

#include <dune/grid/common/gridenums.hh>

#include <dumux/assembly/diffmethod.hh>
#include <dumux/assembly/numericepsilon.hh>

#include <dumux/discretization/method.hh>
#include <dumux/timestepping/multistagetimestepper.hh>


namespace Dumux {

/*!
 * \ingroup Assembly
 * \brief An assembler for Jacobian and residual contribution
 *        per element for the box scheme.
 * \tparam Assembler The grid-wide assembler type
 */
template<class Assembler, DiffMethod diffMethod>
class BoxLocalAssembler
{
    using GridVariables = typename Assembler::GridVariables;
    using GridGeometry = typename GridVariables::GridGeometry;

    using FVElementGeometry = typename GridGeometry::LocalView;
    using ElementVariables = typename GridVariables::LocalView;
    using PrimaryVariables = typename GridVariables::PrimaryVariables;
    using Scalar = typename GridVariables::Scalar;

    using GridView = typename GridGeometry::GridView;
    using Element = typename GridView::template Codim<0>::Entity;

    using LocalOperator = typename Assembler::LocalOperator;
    using JacobianMatrix = typename Assembler::JacobianMatrix;
    using ResidualVector = typename Assembler::ResidualVector;

    static constexpr int numEq = PrimaryVariables::size();
    static_assert(diffMethod == DiffMethod::numeric, "Analytical assembly not implemented");
    static_assert(numEq == JacobianMatrix::block_type::cols, "Matrix block size doesn't match privars size");

    //! the parameters of a stage in time integration
    using StageParams = MultiStageParams<Scalar>;

public:
    using ElementResidualVector = typename LocalOperator::ElementResidualVector;

    /*!
     * \brief Constructor for stationary problems.
     */
    explicit BoxLocalAssembler(const Element& element,
                               const FVElementGeometry& fvGeometry,
                               std::vector<ElementVariables>& elemVars)
    : element_(element)
    , fvGeometry_(fvGeometry)
    , elemVariables_(elemVars)
    , elementIsGhost_(element.partitionType() == Dune::GhostEntity)
    , stageParams_(nullptr)
    {
        assert(elemVars.size() == 1);
    }

    /*!
     * \brief Constructor for instationary problems.
     * \note Using this constructor, we assemble one stage within
     *       a time integration step using multi-stage methods.
     */
    explicit BoxLocalAssembler(const Element& element,
                               const FVElementGeometry& fvGeometry,
                               std::vector<ElementVariables>& elemVars,
                               std::shared_ptr<const StageParams> stageParams)
    : element_(element)
    , fvGeometry_(fvGeometry)
    , elemVariables_(elemVars)
    , elementIsGhost_(element.partitionType() == Dune::GhostEntity)
    , stageParams_(stageParams)
    {}

    /*!
     * \brief Computes the derivatives with respect to the given element and adds
     *        them to the global matrix. The element residual is written into the
     *        right hand side.
     */
    template <class PartialReassembler = DefaultPartialReassembler>
    void assembleJacobianAndResidual(JacobianMatrix& jac,
                                     ResidualVector& res,
                                     const PartialReassembler* partialReassembler = nullptr)
    {
        const auto eIdxGlobal = fvGeometry().gridGeometry().elementMapper().index(element());

        if (partialReassembler && partialReassembler->elementColor(eIdxGlobal) == EntityColor::green)
        {
            const auto residual = evalLocalResidual();
            for (const auto& scv : scvs(fvGeometry()))
                res[scv.dofIndex()] += residual[scv.localDofIndex()];
        }
        else if (!elementIsGhost_)
        {
            const auto residual = assembleJacobianAndResidual_(jac, partialReassembler);
            for (const auto& scv : scvs(fvGeometry()))
                res[scv.dofIndex()] += residual[scv.localDofIndex()];
        }
        else
        {
            static constexpr int dim = GridView::dimension;

            for (const auto& scv : scvs(fvGeometry()))
            {
                const auto vIdxLocal = scv.indexInElement();
                const auto& v = element().template subEntity<dim>(vIdxLocal);

                // do not change the non-ghost vertices
                if (v.partitionType() == Dune::InteriorEntity ||
                    v.partitionType() == Dune::BorderEntity)
                    continue;

                const auto dofIdx = scv.dofIndex();
                res[dofIdx] = 0;
                for (int pvIdx = 0; pvIdx < numEq; ++pvIdx)
                    jac[dofIdx][dofIdx][pvIdx][pvIdx] = 1.0;
            }
        }

        auto applyDirichlet = [&] (const auto& scvI,
                                   const auto& dirichletValues,
                                   const auto eqIdx,
                                   const auto pvIdx)
        {
            res[scvI.dofIndex()][eqIdx] = elemVariables().elemVolVars()[scvI].priVars()[pvIdx] - dirichletValues[pvIdx];

            auto& row = jac[scvI.dofIndex()];
            for (auto col = row.begin(); col != row.end(); ++col)
                row[col.index()][eqIdx] = 0.0;

            jac[scvI.dofIndex()][scvI.dofIndex()][eqIdx][pvIdx] = 1.0;

            // if a periodic dof has Dirichlet values also apply the same Dirichlet values to the other dof
            if (fvGeometry().gridGeometry().dofOnPeriodicBoundary(scvI.dofIndex()))
            {
                const auto periodicDof = fvGeometry().gridGeometry().periodicallyMappedDof(scvI.dofIndex());
                res[periodicDof][eqIdx] = elemVariables().elemVolVars()[scvI].priVars()[pvIdx] - dirichletValues[pvIdx];
                const auto end = jac[periodicDof].end();
                for (auto it = jac[periodicDof].begin(); it != end; ++it)
                    (*it) = periodicDof != it.index() ? 0.0 : 1.0;
            }
        };

        enforceDirichletConstraints(applyDirichlet);
    }

    /*!
     * \brief Computes the derivatives with respect to the given element and adds them
     *        to the global matrix.
     */
    void assembleJacobian(JacobianMatrix& jac)
    {
        assembleJacobianAndResidual_(jac);

        auto applyDirichlet = [&] (const auto& scvI,
                                   const auto& dirichletValues,
                                   const auto eqIdx,
                                   const auto pvIdx)
        {
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
    void assembleResidual(ResidualVector& res)
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

protected:

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
                LocalOperator localOperator(element(), fvGeometry(), elemVariables_[k]);

                if (!stageParams_->skipTemporal(k))
                    residual.axpy(stageParams_->temporalWeight(k), localOperator.evalStorage());
                if (!stageParams_->skipSpatial(k))
                    residual.axpy(stageParams_->spatialWeight(k), localOperator.evalFluxesAndSources());
            }

            return residual;
        }
    }

    /*!
     * \brief Computes the derivatives with respect to the dofs of the given
     *        element and adds them to the global matrix.
     * \return The element residual at the current solution.
     */
    template< class PartialReassembler = DefaultPartialReassembler >
    ElementResidualVector assembleJacobianAndResidual_(JacobianMatrix& A,
                                                       const PartialReassembler* partialReassembler = nullptr)
    {
        if constexpr (diffMethod == DiffMethod::numeric)
            return assembleJacobianAndResidualNumeric_(A, partialReassembler);
        else
            DUNE_THROW(Dune::NotImplemented, "Analytic assembler for box");
    }

    /*!
     * \brief Computes the derivatives by means of numeric differentiation
     *        and adds them to the global matrix.
     * \return The element residual at the current solution.
     * \note This specialization is for the box scheme with numeric differentiation
     */
    template< class PartialReassembler = DefaultPartialReassembler >
    ElementResidualVector assembleJacobianAndResidualNumeric_(JacobianMatrix& A,
                                                              const PartialReassembler* partialReassembler = nullptr)
    {
        // get the variables of the current stage
        auto& curVariables = elemVariables();
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
                    elemSol[scvI.localDofIndex()][pvIdx] = priVar;
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

                    // don't add derivatives for green entities
                    if (!partialReassembler || partialReassembler->dofColor(globalJ) != EntityColor::green)
                    {
                        for (int eqIdx = 0; eqIdx < numEq; eqIdx++)
                        {
                            // A[i][col][eqIdx][pvIdx] is the rate of change of the
                            // the residual of equation 'eqIdx' at dof 'i'
                            // depending on the primary variable 'pvIdx' at dof 'col'
                            A[globalJ][globalI][eqIdx][pvIdx] += partialDerivs[localJ][eqIdx];
                        }
                    }
                }

                // restore the original element solution & volume variables
                elemSol[localI][pvIdx] = origElemSol[localI][pvIdx];
                curVolVars = origCurVolVars;

                // TODO additional dof dependencies
            }
        }

        return origResiduals;
    }

    //! Return references to the local views
    const Element& element() const { return element_; }
    const FVElementGeometry& fvGeometry() const { return fvGeometry_; }
    const ElementVariables& elemVariables() const { return elemVariables_.back(); }
    ElementVariables& elemVariables() { return elemVariables_.back(); }

    //! Returns if a stationary problem is assembled
    bool isStationary() const { return !stageParams_; }

    //! Return a reference to the underlying problem
    //! TODO: Should grid vars return problem directly!?
    const auto& problem() const
    { return elemVariables().gridVariables().gridVolVars().problem(); }

private:
    const Element& element_;
    const FVElementGeometry& fvGeometry_;
    std::vector<ElementVariables>& elemVariables_;

    bool elementIsGhost_;
    std::shared_ptr<const StageParams> stageParams_;
};

} // end namespace Dumux

#endif
