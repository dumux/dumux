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
 *        per element for cell-centered schemes.
 */
#ifndef DUMUX_CC_LOCAL_ASSEMBLER_HH
#define DUMUX_CC_LOCAL_ASSEMBLER_HH

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
template<class Assembler, DiffMethod diffMethod>
class CCLocalAssembler
{
    using GridVariables = typename Assembler::GridVariables;
    using GridGeometry = typename GridVariables::GridGeometry;

    using FVElementGeometry = typename GridGeometry::LocalView;
    using ElementVariables = typename GridVariables::LocalView;
    using PrimaryVariables = typename GridVariables::PrimaryVariables;
    using Scalar = typename GridVariables::Scalar;

    using GridView = typename GridGeometry::GridView;
    using Element = typename GridView::template Codim<0>::Entity;

    using JacobianMatrix = typename Assembler::JacobianMatrix;
    using ResidualVector = typename Assembler::ResidualVector;
    using LocalOperator = typename Assembler::LocalOperator;
    using Operators  = typename LocalOperator::Operators;
    using NumEqVector = typename Operators::NumEqVector;

    static constexpr int numEq = PrimaryVariables::size();
    static constexpr int maxElementStencilSize = GridGeometry::maxElementStencilSize;
    static_assert(diffMethod == DiffMethod::numeric, "Analytical assembly not implemented");
    static_assert(numEq == JacobianMatrix::block_type::cols, "Matrix block size doesn't match privars size");

    //! the parameters of a stage in time integration
    using StageParams = MultiStageParams<Scalar>;

public:
    using ElementResidualVector = typename LocalOperator::ElementResidualVector;

    /*!
     * \brief Constructor for stationary problems.
     */
    explicit CCLocalAssembler(const Element& element,
                              const FVElementGeometry& fvGeometry,
                              std::vector<ElementVariables>& elemVars)
    : element_(element)
    , fvGeometry_(fvGeometry)
    , elementVariables_(elemVars)
    , elementIsGhost_(element.partitionType() == Dune::GhostEntity)
    , stageParams_(nullptr)
    {
        assert(elemVars.size() == 1);
        assert(fvGeometry_.numScv() == 1);
    }

    /*!
     * \brief Constructor for instationary problems.
     * \note Using this constructor, we assemble one stage within
     *       a time integration step using multi-stage methods.
     */
    explicit CCLocalAssembler(const Element& element,
                              const FVElementGeometry& fvGeometry,
                              std::vector<ElementVariables>& elemVars,
                              std::shared_ptr<const StageParams> stageParams)
    : element_(element)
    , fvGeometry_(fvGeometry)
    , elementVariables_(elemVars)
    , elementIsGhost_(element.partitionType() == Dune::GhostEntity)
    , stageParams_(stageParams)
    { assert(fvGeometry_.numScv() == 1); }

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
        const auto globalI = fvGeometry().gridGeometry().elementMapper().index(element());
        if (partialReassembler
            && partialReassembler->elementColor(globalI) == EntityColor::green)
        {
            res[globalI] = evalLocalResidual()[0]; // forward to the internal implementation
        }
        else
        {
            res[globalI] = assembleJacobianAndResidual_(jac); // forward to the internal implementation
        }
    }

    /*!
     * \brief Computes the derivatives with respect to the given element and adds them
     *        to the global matrix.
     */
    void assembleJacobian(JacobianMatrix& jac)
    {
        assembleJacobianAndResidual_(jac);
    }

    /*!
     * \brief Assemble the residual only
     */
    void assembleResidual(ResidualVector& res)
    {
        const auto residual = evalLocalResidual();
        for (const auto& scv : scvs(fvGeometry()))
            res[scv.dofIndex()] += residual[scv.localDofIndex()];
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
                LocalOperator localOperator(element(), fvGeometry(), elementVariables_[k]);

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
    NumEqVector assembleJacobianAndResidual_(JacobianMatrix& A)
    {
        if constexpr (diffMethod == DiffMethod::numeric)
            return assembleJacobianAndResidualNumeric_(A);
        else
            DUNE_THROW(Dune::NotImplemented, "Analytic assembler for cc schemes");
    }

    /*!
     * \brief Computes the derivatives by means of numeric differentiation
     *        and adds them to the global matrix.
     * \return The element residual at the current solution.
     * \note This specialization is for numeric differentiation
     */
    NumEqVector assembleJacobianAndResidualNumeric_(JacobianMatrix& A)
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

            // For instationary simulations, scale the coupling
            // fluxes of the current stage correctly
            if (stageParams_)
            {
                for (std::size_t k = 0; k < numNeighbors; ++k)
                    partialDerivs[k+1] *= stageParams_->spatialWeight(stageParams_->size()-1);
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

    //! Return references to the local views
    const Element& element() const { return element_; }
    const FVElementGeometry& fvGeometry() const { return fvGeometry_; }
    const ElementVariables& elemVariables() const { return elementVariables_.back(); }
    ElementVariables& elemVariables_() { return elementVariables_.back(); }

    //! Returns if a stationary problem is assembled
    bool isStationary() const { return !stageParams_; }

    //! Return a reference to the underlying problem
    //! TODO: Should grid vars return problem directly!?
    const auto& problem() const
    { return elemVariables().gridVariables().gridVolVars().problem(); }

private:
    const Element& element_;
    const FVElementGeometry& fvGeometry_;
    std::vector<ElementVariables>& elementVariables_;

    bool elementIsGhost_;
    std::shared_ptr<const StageParams> stageParams_;
};

} // end namespace Dumux

#endif
