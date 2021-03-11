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

#include <dumux/common/numeqvector.hh>
#include <dumux/common/typetraits/problem.hh>

#include <dumux/discretization/method.hh>
#include <dumux/discretization/solutionstate.hh>
#include <dumux/timestepping/multistagetimestepper.hh>


namespace Dumux::Experimental {

/*!
 * \ingroup Assembly
 * \brief An assembler for Jacobian and residual contribution
 *        per element for the box scheme.
 * \tparam LO The element-local operator
 */
template<class LO>
class BoxLocalAssembler
{
    using LocalContext = typename LO::LocalContext;
    using ElementVariables = typename LocalContext::ElementVariables;
    using ElementGridGeometry = typename LocalContext::ElementGridGeometry;

    using GG = typename ElementGridGeometry::GridGeometry;
    using GV = typename ElementVariables::GridVariables;

    using FVElementGeometry = typename GG::LocalView;
    using SubControlVolume = typename GG::SubControlVolume;
    using Element = typename GG::GridView::template Codim<0>::Entity;

    using PrimaryVariables = typename GV::PrimaryVariables;
    using Scalar = typename GV::Scalar;

    static constexpr int numEq = NumEqVectorTraits<PrimaryVariables>::numEq;

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
    explicit BoxLocalAssembler(const Element& element,
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
    explicit BoxLocalAssembler(const Element& element,
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
        const auto eIdxGlobal = fvGeometry_.gridGeometry().elementMapper().index(element);

        if (partialReassembler && partialReassembler->elementColor(eIdxGlobal) == EntityColor::green)
        {
            const auto residual = evalLocalResidual();
            for (const auto& scv : scvs(fvGeometry_))
                res[scv.dofIndex()] += residual[scv.localDofIndex()];
        }
        else if (!elementIsGhost_)
        {
            const auto residual = assembleJacobianAndResidual_(jac, partialReassembler);
            for (const auto& scv : scvs(fvGeometry_))
                res[scv.dofIndex()] += residual[scv.localDofIndex()];
        }
        else
        {
            static constexpr int dim = GG::GridView::dimension;

            for (const auto& scv : scvs(fvGeometry_))
            {
                const auto vIdxLocal = scv.indexInElement();
                const auto& v = fvGeometry_.element().template subEntity<dim>(vIdxLocal);

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
            const auto& priVars = elemVariables_.elemVolVars()[scvI].priVars();
            res[scvI.dofIndex()][eqIdx] = priVars[pvIdx] - dirichletValues[pvIdx];

            auto& row = jac[scvI.dofIndex()];
            for (auto col = row.begin(); col != row.end(); ++col)
                row[col.index()][eqIdx] = 0.0;

            jac[scvI.dofIndex()][scvI.dofIndex()][eqIdx][pvIdx] = 1.0;

            // if a periodic dof has Dirichlet values also apply the same Dirichlet values to the other dof
            if (fvGeometry_.gridGeometry().dofOnPeriodicBoundary(scvI.dofIndex()))
            {
                const auto periodicDof = fvGeometry_.gridGeometry().periodicallyMappedDof(scvI.dofIndex());
                res[periodicDof][eqIdx] = priVars[pvIdx] - dirichletValues[pvIdx];
                const auto end = jac[periodicDof].end();
                for (auto it = jac[periodicDof].begin(); it != end; ++it)
                    (*it) = periodicDof != it.index() ? 0.0 : 1.0;
            }
        };

        enforceDirichletConstraints(applyDirichlet);
    }

    /*!
     * \brief Computes the derivatives with respect to the dofs of the given
     *        element and adds them to the given jacobian matrix.
     */
    template<class JacobianMatrix>
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

    //! Assemble the residual only
    template<class ResidualVector>
    void assembleResidual(ResidualVector& res)
    {
        const auto residual = evalLocalResidual();
        for (const auto& scv : scvs(fvGeometry_))
            res[scv.dofIndex()] += residual[scv.localDofIndex()];

        auto applyDirichlet = [&] (const auto& scvI,
                                   const auto& dirichletValues,
                                   const auto eqIdx,
                                   const auto pvIdx)
        {
            const auto& priVars = elemVariables_.elemVolVars()[scvI].priVars();
            res[scvI.dofIndex()][eqIdx] = priVars[pvIdx] - dirichletValues[pvIdx];
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

    //! Evaluates Dirichlet boundaries
    template< typename ApplyDirichletFunctionType >
    void evalDirichletBoundaries(ApplyDirichletFunctionType applyDirichlet)
    {
        const auto& element = fvGeometry_.element();
        const auto& problem = gridVariables_.gridVolVars().problem();

        for (const auto& scvI : scvs(fvGeometry_))
        {
            const auto bcTypes = problem.boundaryTypes(element, scvI);
            if (bcTypes.hasDirichlet())
            {
                const auto dirichletValues = getDirichletValues_(problem, element, scvI,
                                                                 gridVariables_.timeLevel());

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
    ElementResidualVector assembleJacobianAndResidual_(JacobianMatrix& A,
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
    ElementResidualVector assembleJacobianAndResidualNumeric_(JacobianMatrix& A,
                                                              const PartialReassembler* partialReassembler = nullptr)
    {
        // alias for the variables of the current stage
        auto& curVariables = elemVariables_;
        auto& curElemVolVars = curVariables.elemVolVars();
        const auto& element = fvGeometry_.element();
        const auto& problem = curElemVolVars.gridVolVars().problem();

        const auto origResiduals = evalLocalResidual();
        auto elemSolState = makeElementSolutionState(element, curVariables.gridVariables());
        auto& elemSol = elemSolState.elementSolution();
        const auto origElemSol = elemSol;

        //////////////////////////////////////////////////////////////////////////////////////////////
        // Calculate derivatives of the residual of all dofs in element with respect to themselves. //
        //////////////////////////////////////////////////////////////////////////////////////////////

        ElementResidualVector partialDerivs(fvGeometry_.numScv());
        for (const auto& scvI : scvs(fvGeometry_))
        {
            // dof index and corresponding actual pri vars
            const auto globalI = scvI.dofIndex();
            const auto localI = scvI.localDofIndex();

            const auto origCurVolVars = curElemVolVars[scvI];
            auto& curVolVars = getVolVarAcess_(curElemVolVars, gridVariables_.gridVolVars(), scvI);

            // calculate derivatives w.r.t to the privars at the dof at hand
            for (int pvIdx = 0; pvIdx < numEq; pvIdx++)
            {
                partialDerivs = 0.0;
                auto evalResiduals = [&](Scalar priVar)
                {
                    // update the volume variables and compute element residual
                    elemSol[scvI.localDofIndex()][pvIdx] = priVar;
                    curVolVars.update(elemSolState, problem, element, scvI);
                    return evalLocalResidual();
                };

                static const NumericEpsilon<Scalar, numEq> numEps{problem.paramGroup()};
                static const int numDiffMethod = getParamFromGroup<int>(problem.paramGroup(),
                                                                        "Assembly.NumericDifferenceMethod");

                // epsilon used for numeric differentiation
                const auto eps = numEps(elemSol[localI][pvIdx], pvIdx);

                // derive the residuals numerically
                NumericDifferentiation::partialDerivative(evalResiduals, elemSol[localI][pvIdx], partialDerivs,
                                                          origResiduals, eps, numDiffMethod);

                // TODO: Distinguish between implicit/explicit here. For explicit schemes,
                //       no entries between different scvs of an element are reserved. Thus,
                //       we currently get an error when using explicit schemes.
                // TODO: Doesn't this mean we only have to evaluate the residual for a single
                //       scv instead of calling evalLocalResidual()? That computes the residuals
                //       and derivs for all other scvs of the element, too, which are never used.
                //       Note: this is the same in the current implementation of master.
                //       Should we try to optimize this for explicit schemes? Or adjust the Jacobian pattern?
                // update the global stiffness matrix with the current partial derivatives
                for (const auto& scvJ : scvs(fvGeometry_))
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

    //! Returns the volume variables of the grid vol vars (case of global caching)
    template<class ElemVolVars, class GridVolVars>
    auto& getVolVarAcess_(ElemVolVars& elemVolVars,
                          GridVolVars& gridVolVars,
                          const SubControlVolume& scv)
    {
        if constexpr (GridVolVars::cachingEnabled)
            return gridVolVars.volVars(scv);
        else
            return elemVolVars[scv];
    }

    //! get the user-defined Dirichlet boundary conditions
    template<class Problem, class TimeLevel>
    auto getDirichletValues_(const Problem& problem,
                             const Element& element,
                             const SubControlVolume& scv,
                             const TimeLevel& timeLevel)
    {
        if constexpr(ProblemTraits<Problem>::hasTransientDirichletInterface)
            return problem.dirichlet(element, scv, timeLevel);
        else
            return problem.dirichlet(element, scv);
    }

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
