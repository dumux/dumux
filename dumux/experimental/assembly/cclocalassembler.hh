// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup Experimental
 * \ingroup Assembly
 * \ingroup CCDiscretization
 * \brief An assembler for Jacobian and residual contribution per element (cell-centered methods)
 */
#ifndef DUMUX_EXPERIMENTAL_CC_LOCAL_ASSEMBLER_HH
#define DUMUX_EXPERIMENTAL_CC_LOCAL_ASSEMBLER_HH

#include <dune/common/reservedvector.hh>
#include <dune/grid/common/gridenums.hh> // for GhostEntity
#include <dune/istl/matrixindexset.hh>

#include <dumux/common/typetraits/typetraits.hh>
#include <dumux/common/reservedblockvector.hh>
#include <dumux/common/properties.hh>
#include <dumux/common/parameters.hh>
#include <dumux/common/numericdifferentiation.hh>
#include <dumux/common/numeqvector.hh>
#include <dumux/assembly/numericepsilon.hh>
#include <dumux/assembly/diffmethod.hh>
#include <dumux/assembly/entitycolor.hh>
#include <dumux/assembly/partialreassembler.hh>
#include <dumux/discretization/fluxstencil.hh>
#include <dumux/discretization/cellcentered/elementsolution.hh>

#include <dumux/experimental/assembly/fvlocalassemblerbase.hh>

namespace Dumux::Experimental {

/*!
 * \ingroup Experimental
 * \ingroup Assembly
 * \ingroup CCDiscretization
 * \brief A base class for all local cell-centered assemblers
 * \tparam TypeTag The TypeTag
 * \tparam Assembler The assembler type
 * \tparam Implementation The actual implementation
 */
template<class TypeTag, class Assembler, class Implementation>
class CCLocalAssemblerBase : public Experimental::FVLocalAssemblerBase<TypeTag, Assembler, Implementation>
{
    using ParentType = Experimental::FVLocalAssemblerBase<TypeTag, Assembler, Implementation>;
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using FVElementGeometry =  typename GridGeometry::LocalView;
    using Element =  typename FVElementGeometry::Element;
    using SubControlVolumeFace =  typename GridGeometry::SubControlVolumeFace;
    using GridView = typename GridGeometry::GridView;
    using JacobianMatrix = GetPropType<TypeTag, Properties::JacobianMatrix>;
    using GridVariables = GetPropType<TypeTag, Properties::GridVariables>;
    using Problem = GetPropType<TypeTag, Properties::Problem>;
    using ElementVolumeVariables = typename GetPropType<TypeTag, Properties::GridVolumeVariables>::LocalView;
    using NumEqVector = Dumux::NumEqVector<GetPropType<TypeTag, Properties::PrimaryVariables>>;

public:

    using ParentType::ParentType;

    template <class ResidualVector, class StageParams, class PartialReassembler = DefaultPartialReassembler, class CouplingFunction = Noop>
    void assembleJacobianAndResidual(JacobianMatrix& jac, ResidualVector& res, GridVariables& gridVariables,
                                     const StageParams& stageParams, ResidualVector& temporal, ResidualVector& spatial,
                                     ResidualVector& constrainedDofs,
                                     const PartialReassembler* partialReassembler = nullptr,
                                     const CouplingFunction& maybeAssembleCouplingBlocks = noop)
    {
        this->asImp_().bindLocalViews();
        const auto globalI = this->fvGeometry().gridGeometry().elementMapper().index(this->element());

        this->localResidual().spatialWeight(1.0);
        this->localResidual().temporalWeight(1.0);

        spatial[globalI] = this->evalLocalFluxAndSourceResidual(this->curElemVolVars())[0];
        temporal[globalI] = this->localResidual().evalStorage(this->fvGeometry(), this->curElemVolVars())[0];
        res[globalI] = spatial[globalI]*stageParams.spatialWeight(stageParams.size()-1)
                       + temporal[globalI]*stageParams.temporalWeight(stageParams.size()-1);

        this->localResidual().spatialWeight(stageParams.spatialWeight(stageParams.size()-1));
        this->localResidual().temporalWeight(stageParams.temporalWeight(stageParams.size()-1));


        if (partialReassembler
            && partialReassembler->elementColor(globalI) == EntityColor::green)
        {
            // assemble the coupling blocks for coupled models (does nothing if not coupled)
            maybeAssembleCouplingBlocks(res[globalI]);
        }
        else
        {
            this->asImp_().assembleJacobian(jac, gridVariables, res[globalI]); // forward to the internal implementation

            // assemble the coupling blocks for coupled models (does nothing if not coupled)
            maybeAssembleCouplingBlocks(res[globalI]);
        }
    }

    /*!
     * \brief Assemble the residual only
     */
    template <class ResidualVector>
    void assembleResidual(ResidualVector& res)
    {
        this->asImp_().bindLocalViews();
        const auto globalI = this->assembler().gridGeometry().elementMapper().index(this->element());
        res[globalI] = this->asImp_().evalLocalResidual()[0]; // forward to the internal implementation

        using Problem = GetPropType<TypeTag, Properties::Problem>;
        if constexpr (Problem::enableInternalDirichletConstraints())
        {
            const auto applyDirichlet = [&] (const auto& scvI,
                                             const auto& dirichletValues,
                                             const auto eqIdx,
                                             const auto pvIdx)
            {
                res[scvI.dofIndex()][eqIdx] = this->curElemVolVars()[scvI].priVars()[pvIdx] - dirichletValues[pvIdx];
            };

            this->asImp_().enforceInternalDirichletConstraints(applyDirichlet);
        }
    }

    /*!
     * \brief Evaluates the fluxes (element can potentially be a neighbor)
     */
    NumEqVector evalFlux(const Element& neighbor,
                         const SubControlVolumeFace& scvf) const
    {
        return this->localResidual().evalFlux(
            this->asImp_().problem(), neighbor,
            this->fvGeometry(), this->curElemVolVars(),
            this->elemFluxVarsCache(), scvf
        );
    }

    /*!
     * \brief Assemble the residual only
     */
    template<class SubResidualVector>
    void assembleCurrentResidual(SubResidualVector& spatialRes, SubResidualVector& temporalRes)
    {
        this->asImp_().bindLocalViews();
        const auto globalI = this->fvGeometry().gridGeometry().elementMapper().index(this->element());
        spatialRes[globalI] = this->evalLocalFluxAndSourceResidual(this->curElemVolVars())[0];
        temporalRes[globalI] = this->localResidual().evalStorage(this->fvGeometry(), this->curElemVolVars())[0];

        if constexpr (Problem::enableInternalDirichletConstraints())
            DUNE_THROW(Dune::NotImplemented, "Not implemented");
    }

    /*!
     * \brief Update the coupling context for coupled models.
     * \note This does nothing per default (not a coupled model).
     */
    template<class... Args>
    void maybeUpdateCouplingContext(Args&&...) {}

    /*!
     * \brief Update the additional domain derivatives for coupled models.
     * \note This does nothing per default (not a coupled model).
     */
    template<class... Args>
    void maybeEvalAdditionalDomainDerivatives(Args&&...) {}
};

/*!
 * \ingroup Experimental
 * \ingroup Assembly
 * \ingroup CCDiscretization
 * \brief An assembler for Jacobian and residual contribution per element (cell-centered methods)
 * \tparam TypeTag The TypeTag
 * \tparam diffMethod The differentiation method to residual compute derivatives
 * \tparam Implementation The actual implementation via CRTP
 */
template<class TypeTag, class Assembler, DiffMethod diffMethod = DiffMethod::numeric, class Implementation = void>
class CCLocalAssembler;

/*!
 * \ingroup Experimental
 * \ingroup Assembly
 * \ingroup CCDiscretization
 * \brief Cell-centered scheme local assembler using numeric differentiation
 */
template<class TypeTag, class Assembler, class Implementation>
class CCLocalAssembler<TypeTag, Assembler, DiffMethod::numeric, Implementation>
: public CCLocalAssemblerBase<TypeTag, Assembler,
                              NonVoidOr<CCLocalAssembler<TypeTag, Assembler, DiffMethod::numeric, Implementation>, Implementation>>
{
    using ThisType = CCLocalAssembler<TypeTag, Assembler, DiffMethod::numeric, Implementation>;
    using ParentType = CCLocalAssemblerBase<TypeTag, Assembler, NonVoidOr<ThisType, Implementation>>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using NumEqVector = Dumux::NumEqVector<GetPropType<TypeTag, Properties::PrimaryVariables>>;
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using FVElementGeometry = typename GridGeometry::LocalView;
    using Element = typename FVElementGeometry::Element;
    using GridVariables = GetPropType<TypeTag, Properties::GridVariables>;
    using JacobianMatrix = GetPropType<TypeTag, Properties::JacobianMatrix>;
    using Problem = typename GridVariables::GridVolumeVariables::Problem;

    static constexpr int numEq = GetPropType<TypeTag, Properties::ModelTraits>::numEq();
    static constexpr int dim = GetPropType<TypeTag, Properties::GridGeometry>::GridView::dimension;

    using FluxStencil = Dumux::FluxStencil<FVElementGeometry>;
    static constexpr int maxElementStencilSize = GridGeometry::maxElementStencilSize;
    static constexpr bool enableGridFluxVarsCache = GetPropType<TypeTag, Properties::GridVariables>::GridFluxVariablesCache::cachingEnabled;

public:
    using ParentType::ParentType;

    using LocalResidual = typename ParentType::LocalResidual;
    using ElementResidualVector = typename LocalResidual::ElementResidualVector;

    /*!
     * \brief Computes the derivatives with respect to the given element and adds them
     *        to the global matrix. Calculates the element residual at the current solution.
     */
    void assembleJacobian(JacobianMatrix& A, GridVariables& gridVariables, const NumEqVector& origResidual)
    {
        if (this->isImplicit())
            assembleJacobianImplicit_(A, gridVariables, origResidual);
        else
            assembleJacobianExplicit_(A, gridVariables, origResidual);
    }

private:

    /*!
     * \brief Computes the derivatives with respect to the given element and adds them to the global matrix.
     *        Calculates the element residual at the current solution.
     */
    void assembleJacobianImplicit_(JacobianMatrix& A, GridVariables& gridVariables, const NumEqVector& origResidual)
    {
        //////////////////////////////////////////////////////////////////////////////////////////////////
        // Calculate derivatives of all dofs in stencil with respect to the dofs in the element. In the //
        // neighboring elements we do so by computing the derivatives of the fluxes which depend on the //
        // actual element. In the actual element we evaluate the derivative of the entire residual.     //
        //////////////////////////////////////////////////////////////////////////////////////////////////

        // get some aliases for convenience
        const auto& element = this->element();
        const auto& fvGeometry = this->fvGeometry();
        const auto& gridGeometry = this->fvGeometry().gridGeometry();
        auto&& curElemVolVars = this->curElemVolVars();
        auto&& elemFluxVarsCache = this->elemFluxVarsCache();

        // get stencil information
        const auto globalI = gridGeometry.elementMapper().index(element);
        const auto& connectivityMap = gridGeometry.connectivityMap();
        const auto numNeighbors = connectivityMap[globalI].size();

        // container to store the neighboring elements
        Dune::ReservedVector<Element, maxElementStencilSize> neighborElements;
        neighborElements.resize(numNeighbors);

        // assemble the undeflected residual
        using Residuals = ReservedBlockVector<NumEqVector, maxElementStencilSize>;
        Residuals origResiduals(numNeighbors + 1); origResiduals = 0.0;
        origResiduals[0] = origResidual;

        // lambda for convenient evaluation of the fluxes across scvfs in the neighbors
        // if the neighbor is a ghost we don't want to add anything to their residual
        // so we return 0 and omit computing the flux
        auto evalNeighborFlux = [&] (const auto& neighbor, const auto& scvf)
        {
            if (neighbor.partitionType() == Dune::GhostEntity)
                return NumEqVector(0.0);
            else
                return this->evalFlux(neighbor, scvf);
        };

        // get the elements in which we need to evaluate the fluxes
        // and calculate these in the undeflected state
        unsigned int j = 1;
        for (const auto& dataJ : connectivityMap[globalI])
        {
            neighborElements[j-1] = gridGeometry.element(dataJ.globalJ);
            for (const auto scvfIdx : dataJ.scvfsJ)
                origResiduals[j] += evalNeighborFlux(neighborElements[j-1], fvGeometry.scvf(scvfIdx));

            ++j;
        }

        // reference to the element's scv (needed later) and corresponding vol vars
        const auto& scv = fvGeometry.scv(globalI);
        auto& curVolVars = ParentType::getVolVarAccess(gridVariables.curGridVolVars(), curElemVolVars, scv);

        // save a copy of the original privars and vol vars in order
        // to restore the original solution after deflection
        const auto& curSol = this->asImp_().curSol();
        const auto origPriVars = curSol[globalI];
        const auto origVolVars = curVolVars;

        // element solution container to be deflected
        auto elemSol = elementSolution<FVElementGeometry>(origPriVars);

        // derivatives in the neighbors with respect to the current elements
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
                this->asImp_().maybeUpdateCouplingContext(scv, elemSol, pvIdx);
                curVolVars.update(elemSol, this->asImp_().problem(), element, scv);
                elemFluxVarsCache.update(element, fvGeometry, curElemVolVars);
                if (enableGridFluxVarsCache)
                    gridVariables.gridFluxVarsCache().updateElement(element, fvGeometry, curElemVolVars);

                // calculate the residual with the deflected primary variables
                partialDerivsTmp[0] = this->evalLocalResidual()[0];

                // calculate the fluxes in the neighbors with the deflected primary variables
                for (std::size_t k = 0; k < numNeighbors; ++k)
                    for (auto scvfIdx : connectivityMap[globalI][k].scvfsJ)
                        partialDerivsTmp[k+1] += evalNeighborFlux(neighborElements[k], fvGeometry.scvf(scvfIdx));

                return partialDerivsTmp;
            };

            // derive the residuals numerically
            static const NumericEpsilon<Scalar, numEq> eps_{this->asImp_().problem().paramGroup()};
            static const int numDiffMethod = getParamFromGroup<int>(this->asImp_().problem().paramGroup(), "Assembly.NumericDifferenceMethod");
            NumericDifferentiation::partialDerivative(evalResiduals, elemSol[0][pvIdx], partialDerivs, origResiduals,
                                                      eps_(elemSol[0][pvIdx], pvIdx), numDiffMethod);

            // Correct derivative for ghost elements, i.e. set a 1 for the derivative w.r.t. the
            // current primary variable and a 0 elsewhere. As we always solve for a delta of the
            // solution with respect to the initial one, this results in a delta of 0 for ghosts.
            if (this->elementIsGhost())
            {
                partialDerivs[0] = 0.0;
                partialDerivs[0][pvIdx] = 1.0;
            }

            // restore the original state of the scv's volume variables
            curVolVars = origVolVars;

            // restore the current element solution
            elemSol[0][pvIdx] = origPriVars[pvIdx];

            // restore the undeflected state of the coupling context
            this->asImp_().maybeUpdateCouplingContext(scv, elemSol, pvIdx);

            // add the current partial derivatives to the global jacobian matrix
            // no special treatment is needed if globalJ is a ghost because then derivatives have been assembled to 0 above
            if constexpr (Problem::enableInternalDirichletConstraints())
            {
                // check if own element has internal Dirichlet constraint
                const auto internalDirichletConstraintsOwnElement = this->asImp_().problem().hasInternalDirichletConstraint(this->element(), scv);
                const auto dirichletValues = this->asImp_().problem().internalDirichlet(this->element(), scv);

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
                    const auto& neighborScv = fvGeometry.scv(dataJ.globalJ);
                    const auto internalDirichletConstraintsNeighbor = this->asImp_().problem().hasInternalDirichletConstraint(neighborElement, neighborScv);

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
        }

        // restore original state of the flux vars cache in case of global caching.
        // This has to be done in order to guarantee that everything is in an undeflected
        // state before the assembly of another element is called. In the case of local caching
        // this is obsolete because the elemFluxVarsCache used here goes out of scope after this.
        // We only have to do this for the last primary variable, for all others the flux var cache
        // is updated with the correct element volume variables before residual evaluations
        if (enableGridFluxVarsCache)
            gridVariables.gridFluxVarsCache().updateElement(element, fvGeometry, curElemVolVars);

        // evaluate additional derivatives that might arise from the coupling (no-op if not coupled)
        ElementResidualVector orig{origResidual};
        this->asImp_().maybeEvalAdditionalDomainDerivatives(orig, A, gridVariables);
    }

    /*!
     * \brief Computes the derivatives with respect to the given element and adds them
     *        to the global matrix. Calculates the element residual at the current solution.
     */
    void assembleJacobianExplicit_(JacobianMatrix& A, GridVariables& gridVariables, const NumEqVector& origResidual)
    {
        if (this->assembler().isStationaryProblem())
            DUNE_THROW(Dune::InvalidStateException, "Using explicit jacobian assembler with stationary local residual");

        // assemble the undeflected residual
        auto storageResidual = origResidual;

        //////////////////////////////////////////////////////////////////////////////////////////////////
        // Calculate derivatives of all dofs in stencil with respect to the dofs in the element. In the //
        // neighboring elements all derivatives are zero. For the assembled element only the storage    //
        // derivatives are non-zero.                                                                    //
        //////////////////////////////////////////////////////////////////////////////////////////////////

        // get some aliases for convenience
        const auto& element = this->element();
        const auto& fvGeometry = this->fvGeometry();
        const auto& gridGeometry = this->fvGeometry().gridGeometry();
        auto&& curElemVolVars = this->curElemVolVars();

        // reference to the element's scv (needed later) and corresponding vol vars
        const auto globalI = gridGeometry.elementMapper().index(element);
        const auto& scv = fvGeometry.scv(globalI);
        auto& curVolVars = ParentType::getVolVarAccess(gridVariables.curGridVolVars(), curElemVolVars, scv);

        // save a copy of the original privars and vol vars in order
        // to restore the original solution after deflection
        const auto& curSol = this->asImp_().curSol();
        const auto origPriVars = curSol[globalI];
        const auto origVolVars = curVolVars;

        // element solution container to be deflected
        auto elemSol = elementSolution<FVElementGeometry>(origPriVars);

        // derivatives in the neighbors with respect to the current elements
        NumEqVector partialDeriv;
        for (int pvIdx = 0; pvIdx < numEq; ++pvIdx)
        {
            // reset derivatives of element dof with respect to itself
            partialDeriv = 0.0;

            auto evalStorage = [&](Scalar priVar)
            {
                // update the volume variables and calculate
                // the residual with the deflected primary variables
                elemSol[0][pvIdx] = priVar;
                curVolVars.update(elemSol, this->asImp_().problem(), element, scv);
                return this->evalStorage()[0];
            };

            // for non-ghosts compute the derivative numerically
            if (!this->elementIsGhost())
            {
                static const NumericEpsilon<Scalar, numEq> eps_{this->asImp_().problem().paramGroup()};
                static const int numDiffMethod = getParamFromGroup<int>(this->asImp_().problem().paramGroup(), "Assembly.NumericDifferenceMethod");
                NumericDifferentiation::partialDerivative(evalStorage, elemSol[0][pvIdx], partialDeriv, storageResidual,
                                                          eps_(elemSol[0][pvIdx], pvIdx), numDiffMethod);
            }

            // for ghost elements we assemble a 1.0 where the primary variable and zero everywhere else
            // as we always solve for a delta of the solution with respect to the initial solution this
            // results in a delta of zero for ghosts
            else partialDeriv[pvIdx] = 1.0;

            // restore the original state of the scv's volume variables
            curVolVars = origVolVars;

            // restore the current element solution
            elemSol[0][pvIdx] = origPriVars[pvIdx];

            // add the current partial derivatives to the global jacobian matrix
            // only diagonal entries for explicit jacobians
            if constexpr (Problem::enableInternalDirichletConstraints())
            {
                // check if own element has internal Dirichlet constraint
                const auto internalDirichletConstraints = this->asImp_().problem().hasInternalDirichletConstraint(this->element(), scv);
                const auto dirichletValues = this->asImp_().problem().internalDirichlet(this->element(), scv);

                for (int eqIdx = 0; eqIdx < numEq; ++eqIdx)
                {
                    // TODO: we need to set constrained dof flags for these dofs and not modify the residual here
                    // this function only takes care of the Jacobian
                    if (internalDirichletConstraints[eqIdx])
                    {
                        storageResidual[eqIdx] = origVolVars.priVars()[eqIdx] - dirichletValues[eqIdx];
                        A[globalI][globalI][eqIdx][pvIdx] = (eqIdx == pvIdx) ? 1.0 : 0.0;
                    }
                    else
                        A[globalI][globalI][eqIdx][pvIdx] += partialDeriv[eqIdx];
                }
            }
            else
            {
                for (int eqIdx = 0; eqIdx < numEq; eqIdx++)
                    A[globalI][globalI][eqIdx][pvIdx] += partialDeriv[eqIdx];
            }
        }
    }
};

} // end namespace Dumux

#endif
