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
 *        per element for cell-centered finite volume schemes.
 */
#ifndef DUMUX_CC_LOCAL_ASSEMBLER_HH
#define DUMUX_CC_LOCAL_ASSEMBLER_HH

#include <cassert>
#include <memory>

#include <dune/grid/common/gridenums.hh>
#include <dune/common/reservedvector.hh>

#include <dumux/assembly/diffmethod.hh>
#include <dumux/assembly/numericepsilon.hh>

#include <dumux/common/numeqvector.hh>
#include <dumux/common/reservedblockvector.hh>

#include <dumux/discretization/method.hh>
#include <dumux/discretization/elementsolution.hh>
#include <dumux/discretization/solutionstateview.hh>

#include <dumux/timestepping/multistagetimestepper.hh>


namespace Dumux::Experimental {

/*!
 * \ingroup Assembly
 * \brief An assembler for Jacobian and residual contribution
 *        per element for cell-centered finite volume schemes.
 * \tparam LO The element-local operator
 */
template<class LO>
class CCLocalAssembler
{
    using Operators = typename LO::Operators;
    using LocalContext = typename LO::LocalContext;
    using ElementVariables = typename LocalContext::ElementVariables;
    using ElementGridGeometry = typename LocalContext::ElementGridGeometry;

    using GG = typename ElementGridGeometry::GridGeometry;
    using GV = typename ElementVariables::GridVariables;

    using FVElementGeometry = typename GG::LocalView;
    using SubControlVolume = typename GG::SubControlVolume;
    using Element = typename GG::GridView::template Codim<0>::Entity;

    using PrimaryVariables = typename GV::PrimaryVariables;
    using NumEqVector = Dumux::NumEqVector<PrimaryVariables>;
    using Scalar = typename GV::Scalar;

    static constexpr int numEq = NumEqVectorTraits<PrimaryVariables>::numEq;
    static constexpr int maxElementStencilSize = GG::maxElementStencilSize;
    static constexpr bool enableGridFluxVarsCache = GV::GridFluxVariablesCache::cachingEnabled;

public:

    //! the parameters of a stage in time integration
    using StageParams = MultiStageParams<Scalar>;

    //! export the grid variables type on which to operate
    using GridVariables = GV;

    //! export the underlying local operator
    using LocalOperator = LO;

    //! export the vector storing the residuals of all dofs of the element
    using ElementResidualVector = typename LO::ElementResidualVector;

    /*!
     * \brief Constructor for stationary problems.
     */
    explicit CCLocalAssembler(const Element& element,
                              GridVariables& gridVariables,
                              DiffMethod dm = DiffMethod::numeric)
    : diffMethod_(dm)
    , gridVariables_(gridVariables)
    , fvGeometry_(localView(gridVariables.gridGeometry()))
    , elemVariables_(localView(gridVariables))
    , prevElemVariables_()
    , elementIsGhost_(element.partitionType() == Dune::GhostEntity)
    , stageParams_(nullptr)
    {
        if (diffMethod_ != DiffMethod::numeric)
            DUNE_THROW(Dune::NotImplemented, "Provided differentiation method");

        fvGeometry_.bind(element);
        elemVariables_.bind(element, fvGeometry_);
    }

    /*!
     * \brief Constructor for instationary problems.
     * \note Using this constructor, we assemble one stage within
     *       a time integration step using multi-stage methods.
     */
    explicit CCLocalAssembler(const Element& element,
                              const std::vector<GridVariables>& prevGridVars,
                              GridVariables& gridVariables,
                              std::shared_ptr<const StageParams> stageParams,
                              DiffMethod dm = DiffMethod::numeric)
    : diffMethod_(dm)
    , gridVariables_(gridVariables)
    , fvGeometry_(localView(gridVariables.gridGeometry()))
    , elemVariables_(localView(gridVariables))
    , prevElemVariables_()
    , elementIsGhost_(element.partitionType() == Dune::GhostEntity)
    , stageParams_(stageParams)
    {
        if (diffMethod_ != DiffMethod::numeric)
            DUNE_THROW(Dune::NotImplemented, "Provided differentiation method");

        if (prevGridVars.size() != stageParams->size() - 1)
            DUNE_THROW(Dune::InvalidStateException, "Grid Variables for all stages needed");

        fvGeometry_.bind(element);
        for (const auto& gridVars : prevGridVars)
        {
            prevElemVariables_.emplace_back(localView(gridVars));
            prevElemVariables_.back().bind(element, fvGeometry_);
        }

        elemVariables_.bind(element, fvGeometry_);
    }

    /*!
     * \brief Computes the derivatives with respect to the given element and adds
     *        them to the global matrix. The element residual is written into the
     *        right hand side.
     */
    template<class JacobianMatrix,
             class ResidualVector,
             class PartialReassembler = DefaultPartialReassembler>
    void assembleJacobianAndResidual(JacobianMatrix& jac,
                                     ResidualVector& res,
                                     const PartialReassembler* partialReassembler = nullptr)
    {
        const auto& element = fvGeometry_.element();
        const auto globalI = fvGeometry_.gridGeometry().elementMapper().index(element);
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
     * \brief Computes the derivatives with respect to the dofs of the given
     *        element and adds them to the given jacobian matrix.
     */
    template<class JacobianMatrix>
    void assembleJacobian(JacobianMatrix& jac)
    {
        assembleJacobianAndResidual_(jac);
    }

    //! Assemble the residual only
    template<class ResidualVector>
    void assembleResidual(ResidualVector& res)
    {
        const auto residual = evalLocalResidual();
        for (const auto& scv : scvs(fvGeometry_))
            res[scv.dofIndex()] += residual[scv.localDofIndex()];
    }

    //! Evaluate the complete local residual for the current element.
    ElementResidualVector evalLocalResidual() const
    {
        if (isStationary())
        {
            const auto context = makeLocalContext(fvGeometry_, elemVariables_);
            LocalOperator localOperator(context);
            return elementIsGhost_ ? localOperator.getEmptyResidual()
                                   : localOperator.evalFluxesAndSources();
        }
        else
        {
            ElementResidualVector residual(fvGeometry_.numScv());
            residual = 0.0;

            if (!elementIsGhost_)
            {
                // add the terms associated with previous stages
                for (std::size_t k = 0; k < stageParams_->size()-1; ++k)
                    addStageTerms_(residual, k,
                                   makeLocalContext(fvGeometry_, prevElemVariables_[k]));

                // add the terms of the current stage
                addStageTerms_(residual, stageParams_->size()-1,
                               makeLocalContext(fvGeometry_, elemVariables_));
            }

            return residual;
        }
    }

    //! Return true if a stationary problem is assembled
    bool isStationary() const
    { return !stageParams_; }

protected:

    //! add the terms of a stage to the current element residual
    void addStageTerms_(ElementResidualVector& r,
                        std::size_t stageIdx,
                        const LocalContext& context) const
    {
        LocalOperator localOperator(context);
        if (!stageParams_->skipTemporal(stageIdx))
            r.axpy(stageParams_->temporalWeight(stageIdx),
                   localOperator.evalStorage());
        if (!stageParams_->skipSpatial(stageIdx))
            r.axpy(stageParams_->spatialWeight(stageIdx),
                   localOperator.evalFluxesAndSources());
    }

    /*!
     * \brief Computes the derivatives with respect to the dofs of the given
     *        element and adds them to the global matrix.
     * \return The element residual at the current solution.
     */
    template<class JacobianMatrix, class PartialReassembler = DefaultPartialReassembler>
    NumEqVector assembleJacobianAndResidual_(JacobianMatrix& A,
                                             const PartialReassembler* partialReassembler = nullptr)
    {
        // TODO: DO we need this to be constexpr?
        if (diffMethod_ == DiffMethod::numeric)
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
    template<class JacobianMatrix, class PartialReassembler = DefaultPartialReassembler>
    NumEqVector assembleJacobianAndResidualNumeric_(JacobianMatrix& A,
                                                    const PartialReassembler* partialReassembler = nullptr)
    {
        // alias for the variables of the current stage
        auto& curVariables = elemVariables_;
        auto& curElemVolVars = curVariables.elemVolVars();
        auto& curElemFluxVarsCache = curVariables.elemFluxVarsCache();
        const auto& x = curVariables.gridVariables().dofs();

        // get stencil informations
        const auto& element = fvGeometry_.element();
        const auto& gridGeometry = fvGeometry_.gridGeometry();
        const auto& connectivityMap = gridGeometry.connectivityMap();
        const auto globalI = gridGeometry.elementMapper().index(element);
        const auto numNeighbors = connectivityMap[globalI].size();

        // deduce the problem type
        const auto& problem = curElemVolVars.gridVolVars().problem();
        using Problem = std::decay_t<decltype(problem)>;

        // assemble the undeflected residual
        using Residuals = ReservedBlockVector<NumEqVector, maxElementStencilSize>;
        Residuals origResiduals(numNeighbors + 1); origResiduals = 0.0;
        origResiduals[0] = evalLocalResidual()[0];

        // lambda for convenient evaluation of the fluxes across scvfs in the neighbors
        auto evalNeighborFlux = [&] (const auto& neighbor, const auto& scvf)
        {
            const auto context = makeLocalContext(fvGeometry_, elemVariables_);
            return LocalOperator(context).computeFlux(neighbor, scvf);
        };

        // get the elements in which we need to evaluate the fluxes
        // and calculate these in the undeflected state
        Dune::ReservedVector<Element, maxElementStencilSize> neighborElements;
        neighborElements.resize(numNeighbors);

        unsigned int j = 1;
        for (const auto& dataJ : connectivityMap[globalI])
        {
            neighborElements[j-1] = gridGeometry.element(dataJ.globalJ);

            if (neighborElements[j-1].partitionType() != Dune::GhostEntity)
                for (const auto scvfIdx : dataJ.scvfsJ)
                    origResiduals[j] += evalNeighborFlux(neighborElements[j-1],
                                                         fvGeometry_.scvf(scvfIdx));

            ++j;
        }

        // reference to the element's scv (needed later) and corresponding vol vars
        const auto& scv = fvGeometry_.scv(globalI);
        auto& curVolVars = getVolVarAcess_(curElemVolVars, gridVariables_.gridVolVars(), scv);

        // save a copy of the original privars and vol vars in order
        // to restore the original solution after deflection
        const auto origPriVars = x[globalI];
        const auto origVolVars = curVolVars;

        // element solution to be deflected
        const Experimental::SolutionStateView solStateView(curVariables.gridVariables());
        const auto origElemSol = elementSolution(element, solStateView, fvGeometry_.gridGeometry());
        auto elemSol = origElemSol;

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
                curVolVars.update(elemSol, problem, element, scv);
                curElemFluxVarsCache.update(element, fvGeometry_, curElemVolVars);
                if (enableGridFluxVarsCache)
                    gridVariables_.gridFluxVarsCache().updateElement(element, fvGeometry_, curElemVolVars);

                // calculate the residual with the deflected primary variables
                partialDerivsTmp[0] = evalLocalResidual()[0];

                // calculate the fluxes in the neighbors with the deflected primary variables
                for (std::size_t k = 0; k < numNeighbors; ++k)
                    if (neighborElements[k].partitionType() != Dune::GhostEntity)
                        for (auto scvfIdx : connectivityMap[globalI][k].scvfsJ)
                            partialDerivsTmp[k+1] += evalNeighborFlux(neighborElements[k],
                                                                      fvGeometry_.scvf(scvfIdx));

                return partialDerivsTmp;
            };

            // derive the residuals numerically
            static const NumericEpsilon<Scalar, numEq> numEps{problem.paramGroup()};
            static const int numDiffMethod = getParamFromGroup<int>(problem.paramGroup(), "Assembly.NumericDifferenceMethod");

            const auto eps = numEps(elemSol[0][pvIdx], pvIdx);
            NumericDifferentiation::partialDerivative(evalResiduals, elemSol[0][pvIdx],
                                                      partialDerivs, origResiduals,
                                                      eps, numDiffMethod);

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
                const auto internalDirichletConstraintsOwnElement = problem.hasInternalDirichletConstraint(element, scv);
                const auto dirichletValues = problem.internalDirichlet(element, scv);

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
                    const auto& neighborElement = gridGeometry.element(dataJ.globalJ);
                    const auto& neighborScv = fvGeometry_.scv(dataJ.globalJ);
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
        curElemFluxVarsCache.update(element, fvGeometry_, curElemVolVars);
        if (enableGridFluxVarsCache)
            gridVariables_.gridFluxVarsCache().updateElement(element, fvGeometry_, curElemVolVars);

        // return the original residual
        return origResiduals[0];
    }

    //! Returns the volume variables of the local view (case of no caching)
    template<class ElemVolVars, class GridVolVars, std::enable_if_t<!GridVolVars::cachingEnabled, int> = 0>
    auto& getVolVarAcess_(ElemVolVars& elemVolVars, GridVolVars& gridVolVars, const SubControlVolume& scv)
    { return elemVolVars[scv]; }

    //! Returns the volume variables of the grid vol vars (case of global caching)
    template<class ElemVolVars, class GridVolVars, std::enable_if_t<GridVolVars::cachingEnabled, int> = 0>
    auto& getVolVarAcess_(ElemVolVars& elemVolVars, GridVolVars& gridVolVars, const SubControlVolume& scv)
    { return gridVolVars.volVars(scv); }

    DiffMethod diffMethod_;                           //!< the type of differentiation method
    GridVariables& gridVariables_;                    //!< reference to the grid variables
    FVElementGeometry fvGeometry_;                    //!< element-local finite volume geometry
    ElementVariables elemVariables_;                  //!< element variables of the current stage
    std::vector<ElementVariables> prevElemVariables_; //!< element variables of prior stages

    bool elementIsGhost_;
    std::shared_ptr<const StageParams> stageParams_;
};

} // end namespace Dumux::Experimental

#endif
