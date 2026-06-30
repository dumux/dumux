// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup Assembly
 * \ingroup CVFEDiscretization
 * \brief An assembler for Jacobian and residual contribution per element (CVFE methods)
 */
#ifndef DUMUX_CVFE_LOCAL_ASSEMBLER__HH
#define DUMUX_CVFE_LOCAL_ASSEMBLER__HH

#include <cassert>
#include <type_traits>
#include <utility>

#include <dune/common/hybridutilities.hh>
#include <dune/grid/common/gridenums.hh>

#include <dumux/common/properties.hh>
#include <dumux/common/parameters.hh>
#include <dumux/common/numericdifferentiation.hh>
#include <dumux/common/multimapperview.hh>
#include <dumux/common/typetraits/localdofs_.hh>

#include <dumux/assembly/numericepsilon.hh>
#include <dumux/assembly/diffmethod.hh>
#include <dumux/assembly/localassemblerbase.hh>
#include <dumux/assembly/partialreassembler.hh>
#include <dumux/assembly/entitycolor.hh>

#include <dumux/discretization/cvfe/elementsolution.hh>
#include <dumux/discretization/cvfe/localdof.hh>

#include "cvfevolvarsdeflectionpolicy_.hh"

namespace Dumux::Experimental {

#ifndef DOXYGEN
namespace Detail::CVFE {

struct NoOperator
{
    template<class... Args>
    constexpr void operator()(Args&&...) const {}
};

template<class X, class Y>
using Impl = std::conditional_t<!std::is_same_v<X, void>, X, Y>;

} // end namespace Detail
#endif // DOXYGEN

/*!
 * \ingroup Assembly
 * \ingroup CVFEDiscretization
 * \brief A base class for all local CVFE assemblers
 * \tparam TypeTag The TypeTag
 * \tparam Assembler The assembler type
 * \tparam Implementation The actual implementation
 * \tparam implicit Specifies whether the time discretization is implicit or not (i.e. explicit)
 */
template<class TypeTag, class Assembler, class Implementation, bool implicit>
class CVFELocalAssemblerBase : public LocalAssemblerBase<TypeTag, Assembler, Implementation, implicit>
{
    using ParentType = LocalAssemblerBase<TypeTag, Assembler, Implementation, implicit>;
    using JacobianMatrix = GetPropType<TypeTag, Properties::JacobianMatrix>;
    using GridVariables = GetPropType<TypeTag, Properties::GridVariables>;
    using GridGeometry =  GetPropType<TypeTag, Properties::GridGeometry>;

    static constexpr int numEq = GetPropType<TypeTag, Properties::ModelTraits>::numEq();
    static constexpr int dim = GridGeometry::GridView::dimension;

public:

    using ParentType::ParentType;
    using LocalResidual = typename ParentType::LocalResidual;
    using ElementResidualVector = typename LocalResidual::ElementResidualVector;

    void bindLocalViews()
    {
        ParentType::bindLocalViews();
    }

    /*!
     * \brief Computes the derivatives with respect to the given element and adds them
     *        to the global matrix. The element residual is written into the right hand side.
     */
    template <class ResidualVector, class PartialReassembler = DefaultPartialReassembler, class CouplingFunction = Detail::CVFE::NoOperator>
    void assembleJacobianAndResidual(JacobianMatrix& jac, ResidualVector& res, GridVariables& gridVariables,
                                     const PartialReassembler* partialReassembler = nullptr,
                                     const CouplingFunction& maybeAssembleCouplingBlocks = {})
    {
        this->asImp_().bindLocalViews();
        const auto eIdxGlobal = this->asImp_().problem().gridGeometry().elementMapper().index(this->element());
        if (partialReassembler
            && partialReassembler->elementColor(eIdxGlobal) == EntityColor::green)
        {
            const auto residual = this->asImp_().evalLocalResidual(); // forward to the internal implementation

            for (const auto& localDof : localDofs(this->fvGeometry()))
                res[localDof.dofIndex()] += residual[localDof.index()];

            // assemble the coupling blocks for coupled models (does nothing if not coupled)
            maybeAssembleCouplingBlocks(residual);
        }
        else if (!this->elementIsGhost())
        {
            const auto residual = this->asImp_().assembleJacobianAndResidualImpl(jac, gridVariables, partialReassembler); // forward to the internal implementation

            for (const auto& localDof : localDofs(this->fvGeometry()))
                res[localDof.dofIndex()] += residual[localDof.index()];

            // assemble the coupling blocks for coupled models (does nothing if not coupled)
            maybeAssembleCouplingBlocks(residual);
        }
        else
        {
            // Treatment of ghost elements
            assert(this->elementIsGhost());

            // handle dofs per codimension
            const auto& gridGeometry = this->asImp_().problem().gridGeometry();
            Dune::Hybrid::forEach(std::make_integer_sequence<int, dim+1>{}, [&](auto d)
            {
                constexpr int codim = dim - d;
                const auto& localCoeffs = gridGeometry.feCache().get(this->element().type()).localCoefficients();
                for (int idx = 0; idx < localCoeffs.size(); ++idx)
                {
                    const auto& localKey = localCoeffs.localKey(idx);

                    // skip if we are not handling this codim right now
                    if (localKey.codim() != codim)
                        continue;

                    // do not change the non-ghost entities
                    auto entity = this->element().template subEntity<codim>(localKey.subEntity());
                    if (entity.partitionType() == Dune::InteriorEntity || entity.partitionType() == Dune::BorderEntity)
                        continue;

                    // Set identity rows for ALL DOFs of this ghost entity.
                    // Entities with multiple DOFs (e.g. PQ3 edge interior DOFs with 2 per edge)
                    // require iterating over all DOF indices via asMultiMapper(dofMapper()).indices(entity).
                    using BlockType = typename JacobianMatrix::block_type;
                    for (const auto dofIndex : asMultiMapper(gridGeometry.dofMapper()).indices(entity))
                    {
                        BlockType &J = jac[dofIndex][dofIndex];
                        for (int j = 0; j < BlockType::rows; ++j)
                            J[j][j] = 1.0;
                        res[dofIndex] = 0;
                    }
                }
            });
        }
    }

    /*!
     * \brief Multi-stage assembly: assembles Jacobian and residual while separately
     *        accumulating the current stage's temporal and spatial operator evaluations.
     *
     * \param jac            The Jacobian matrix to add to
     * \param res            The residual to add the stage-weighted contribution to
     * \param gridVariables  The grid variables for this subdomain
     * \param stageParams    Multi-stage parameters (Butcher tableau weights for this stage)
     * \param temporal       Accumulator for the temporal (storage) operator at this stage
     * \param spatial        Accumulator for the spatial (flux+source) operator at this stage
     * \param constrainedDofs Vector tracking which dofs have Dirichlet constraints
     * \param maybeAssembleCouplingBlocks Optional functor for coupling Jacobian blocks
     */
    template<class ResidualVector, class StageParams, class CouplingFunction = Detail::CVFE::NoOperator>
    void assembleJacobianAndResidual(JacobianMatrix& jac, ResidualVector& res, GridVariables& gridVariables,
                                     const StageParams& stageParams,
                                     ResidualVector& temporal, ResidualVector& spatial,
                                     ResidualVector& constrainedDofs,
                                     const CouplingFunction& maybeAssembleCouplingBlocks = {})
    {
        this->asImp_().bindLocalViews();

        const auto sWeight = stageParams.spatialWeight(stageParams.size()-1);
        const auto tWeight = stageParams.temporalWeight(stageParams.size()-1);

        if (!this->elementIsGhost())
        {
            // evaluate current stage spatial and temporal contributions separately
            const auto spatialResidual = this->evalLocalFluxAndSourceResidual(this->curElemVars());
            const auto storageResidual = this->localResidual().evalStorageCurrentLevel(
                this->element(), this->fvGeometry(), this->curElemVars());

            // accumulate into stage vectors and form weighted stage residual
            ElementResidualVector origResidual(spatialResidual.size());
            origResidual = 0.0;
            for (const auto& localDof : localDofs(this->fvGeometry()))
            {
                const auto li = localDof.index();
                const auto di = localDof.dofIndex();
                spatial[di] += spatialResidual[li];
                temporal[di] += storageResidual[li];
                origResidual[li] += spatialResidual[li]*sWeight + storageResidual[li]*tWeight;
                res[di] += origResidual[li];
            }

            // assemble the Jacobian differentiated w.r.t. the stage-weighted residual
            this->asImp_().assembleJacobianAndResidualImpl(jac, gridVariables, tWeight, sWeight);

            maybeAssembleCouplingBlocks(origResidual);
        }
    }

    /*!
     * \brief Convenience method for assembling only the current residual (no Jacobian)
     *        used during stage 0 preparation in multi-stage assembly.
     */
    template<class ResidualVector>
    void assembleCurrentResidual(ResidualVector& temporal, ResidualVector& spatial)
    {
        this->asImp_().bindLocalViews();

        if (!this->elementIsGhost())
        {
            const auto spatialResidual = this->evalLocalFluxAndSourceResidual(this->curElemVars());
            const auto storageResidual = this->localResidual().evalStorageCurrentLevel(
                this->element(), this->fvGeometry(), this->curElemVars());

            for (const auto& localDof : localDofs(this->fvGeometry()))
            {
                spatial[localDof.dofIndex()] += spatialResidual[localDof.index()];
                temporal[localDof.dofIndex()] += storageResidual[localDof.index()];
            }
        }
    }

    /*!
     * \brief Computes the derivatives with respect to the given element and adds them
     *        to the global matrix.
     */
    void assembleJacobian(JacobianMatrix& jac, GridVariables& gridVariables)
    {
        this->asImp_().bindLocalViews();
        this->asImp_().assembleJacobianAndResidualImpl(jac, gridVariables); // forward to the internal implementation
    }

    /*!
     * \brief Assemble the residual only
     */
    template <class ResidualVector>
    void assembleResidual(ResidualVector& res)
    {
        this->asImp_().bindLocalViews();
        const auto residual = this->evalLocalResidual();

        for (const auto& localDof : localDofs(this->fvGeometry()))
            res[localDof.dofIndex()] += residual[localDof.index()];
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
 * \ingroup Assembly
 * \ingroup CVFEDiscretization
 * \brief An assembler for Jacobian and residual contribution per element (CVFE methods)
 * \tparam TypeTag The TypeTag
 * \tparam diffMethod The differentiation method used to compute derivatives of the residual
 * \tparam implicit Specifies whether the time discretization is implicit or not (i.e. explicit)
 * \tparam Implementation via CRTP
 */
template<class TypeTag, class Assembler, DiffMethod diffMethod = DiffMethod::numeric, bool implicit = true, class Implementation = void>
class CVFELocalAssembler;

/*!
 * \ingroup Assembly
 * \ingroup CVFEDiscretization
 * \brief Control volume finite element local assembler using numeric differentiation and implicit time discretization
 */
template<class TypeTag, class Assembler, class Implementation>
class CVFELocalAssembler<TypeTag, Assembler, DiffMethod::numeric, /*implicit=*/true, Implementation>
: public CVFELocalAssemblerBase<TypeTag, Assembler,
                                Detail::CVFE::Impl<Implementation, CVFELocalAssembler<TypeTag, Assembler, DiffMethod::numeric, true, Implementation>>,
                                true>
{
    using ThisType = CVFELocalAssembler<TypeTag, Assembler, DiffMethod::numeric, true, Implementation>;
    using ParentType = CVFELocalAssemblerBase<TypeTag, Assembler, Detail::CVFE::Impl<Implementation, ThisType>, true>;
    using PrimaryVariable = typename GetPropType<TypeTag, Properties::PrimaryVariables>::value_type;
    using GridVariables = GetPropType<TypeTag, Properties::GridVariables>;
    using JacobianMatrix = GetPropType<TypeTag, Properties::JacobianMatrix>;

    static constexpr int numEq = GetPropType<TypeTag, Properties::ModelTraits>::numEq();
    static constexpr int dim = GetPropType<TypeTag, Properties::GridGeometry>::GridView::dimension;

public:

    using LocalResidual = typename ParentType::LocalResidual;
    using ElementResidualVector = typename LocalResidual::ElementResidualVector;
    using ParentType::ParentType;

    /*!
     * \brief Computes the derivatives with respect to the given element and adds them
     *        to the global matrix.
     *
     * \return The element residual at the current solution.
     */
    template <class PartialReassembler = DefaultPartialReassembler>
    ElementResidualVector assembleJacobianAndResidualImpl(JacobianMatrix& A, GridVariables& gridVariables,
                                                          const PartialReassembler* partialReassembler = nullptr)
    {
        // get some aliases for convenience
        const auto& element = this->element();
        const auto& fvGeometry = this->fvGeometry();
        const auto& curSol = this->asImp_().curSol();

        auto&& curElemVars = this->curElemVars();

        // get the vector of the actual element residuals
        const auto origResiduals = this->evalLocalResidual();

        //////////////////////////////////////////////////////////////////////////////////////////////////
        //                                                                                              //
        // Calculate derivatives of all dofs in stencil with respect to the dofs in the element. In the //
        // neighboring elements we do so by computing the derivatives of the fluxes which depend on the //
        // actual element. In the actual element we evaluate the derivative of the entire residual.     //
        //                                                                                              //
        //////////////////////////////////////////////////////////////////////////////////////////////////

        // if all variables in the stencil have to be updated or if it's enough to only update the
        // variables whose associated dof has been deflected
        static const bool updateAllVars = getParamFromGroup<bool>(
            this->asImp_().problem().paramGroup(), "Assembly.VarsDependOnAllElementDofs", false
        );

        // create the element solution
        auto elemSol = elementSolution(element, curSol, fvGeometry.gridGeometry());

        // create the vector storing the partial derivatives
        ElementResidualVector partialDerivs(Dumux::Detail::LocalDofs::numLocalDofs(fvGeometry));

        auto deflectionPolicy = Dumux::Detail::CVFE::makeVariablesDeflectionPolicy(
            gridVariables.curGridVars(),
            curElemVars,
            fvGeometry,
            updateAllVars
        );

        auto assembleDerivative = [&, this](const auto& localDof)
        {
            // dof index and corresponding actual privars
            const auto dofIdx = localDof.dofIndex();
            const auto localIdx = localDof.index();
            deflectionPolicy.store(localDof);

            // calculate derivatives w.r.t to the privars at the dof at hand
            for (int pvIdx = 0; pvIdx < numEq; pvIdx++)
            {
                partialDerivs = 0.0;

                auto evalResiduals = [&](PrimaryVariable priVar)
                {
                    // update the variables and compute element residual
                    elemSol[localIdx][pvIdx] = priVar;
                    deflectionPolicy.update(elemSol, localDof, this->asImp_().problem());
                    this->asImp_().maybeUpdateCouplingContext(localDof, elemSol, pvIdx);
                    return this->evalLocalResidual();
                };

                // derive the residuals numerically
                static const NumericEpsilon<PrimaryVariable, numEq> eps_{this->asImp_().problem().paramGroup()};
                static const int numDiffMethod = getParamFromGroup<int>(this->asImp_().problem().paramGroup(), "Assembly.NumericDifferenceMethod");
                NumericDifferentiation::partialDerivative(evalResiduals, elemSol[localIdx][pvIdx], partialDerivs, origResiduals,
                                                          eps_(elemSol[localIdx][pvIdx], pvIdx), numDiffMethod);

                // update the global stiffness matrix with the current partial derivatives
                for (const auto& localDofJ : localDofs(fvGeometry))
                {
                    // don't add derivatives for green dofs
                    if (!partialReassembler
                        || partialReassembler->dofColor(localDofJ.dofIndex()) != EntityColor::green)
                    {
                        for (int eqIdx = 0; eqIdx < numEq; eqIdx++)
                        {
                            // A[i][col][eqIdx][pvIdx] is the rate of change of
                            // the residual of equation 'eqIdx' at dof 'i'
                            // depending on the primary variable 'pvIdx' at dof
                            // 'col'.
                            A[localDofJ.dofIndex()][dofIdx][eqIdx][pvIdx] += partialDerivs[localDofJ.index()][eqIdx];
                        }
                    }
                }

                // restore the original state of the localDof's variables
                deflectionPolicy.restore(localDof);

                // restore the original element solution
                elemSol[localIdx][pvIdx] = curSol[localDof.dofIndex()][pvIdx];
                this->asImp_().maybeUpdateCouplingContext(localDof, elemSol, pvIdx);
            }
        };

        // calculation of the derivatives
        for (const auto& localDof : localDofs(fvGeometry))
            assembleDerivative(localDof);

        // evaluate additional derivatives that might arise from the coupling (no-op if not coupled)
        this->asImp_().maybeEvalAdditionalDomainDerivatives(origResiduals, A, gridVariables);

        return origResiduals;
    }

    /*!
     * \brief Multi-stage variant: computes the Jacobian of the stage-weighted residual.
     *
     * Differentiates `tWeight * storage(u) + sWeight * (flux+source)(u)` w.r.t. dof values.
     * Used by the multi-stage overload of `assembleJacobianAndResidual`.
     */
    template<class Scalar>
    ElementResidualVector assembleJacobianAndResidualImpl(JacobianMatrix& A, GridVariables& gridVariables,
                                                           Scalar tWeight, Scalar sWeight)
    {
        const auto& element = this->element();
        const auto& fvGeometry = this->fvGeometry();
        const auto& curSol = this->asImp_().curSol();

        auto&& curElemVars = this->curElemVars();

        const auto origResiduals = this->evalLocalResidualForStage(curElemVars, tWeight, sWeight);

        static const bool updateAllVars = getParamFromGroup<bool>(
            this->asImp_().problem().paramGroup(), "Assembly.VarsDependOnAllElementDofs", false
        );

        auto elemSol = elementSolution(element, curSol, fvGeometry.gridGeometry());
        ElementResidualVector partialDerivs(Dumux::Detail::LocalDofs::numLocalDofs(fvGeometry));

        auto deflectionPolicy = Dumux::Detail::CVFE::makeVariablesDeflectionPolicy(
            gridVariables.curGridVars(), curElemVars, fvGeometry, updateAllVars
        );

        auto assembleDerivative = [&, this](const auto& localDof)
        {
            const auto dofIdx = localDof.dofIndex();
            const auto localIdx = localDof.index();
            deflectionPolicy.store(localDof);

            for (int pvIdx = 0; pvIdx < numEq; pvIdx++)
            {
                partialDerivs = 0.0;

                auto evalResiduals = [&](PrimaryVariable priVar)
                {
                    elemSol[localIdx][pvIdx] = priVar;
                    deflectionPolicy.update(elemSol, localDof, this->asImp_().problem());
                    this->asImp_().maybeUpdateCouplingContext(localDof, elemSol, pvIdx);
                    return this->evalLocalResidualForStage(curElemVars, tWeight, sWeight);
                };

                static const NumericEpsilon<PrimaryVariable, numEq> eps_{this->asImp_().problem().paramGroup()};
                static const int numDiffMethod = getParamFromGroup<int>(this->asImp_().problem().paramGroup(), "Assembly.NumericDifferenceMethod");
                NumericDifferentiation::partialDerivative(evalResiduals, elemSol[localIdx][pvIdx], partialDerivs, origResiduals,
                                                          eps_(elemSol[localIdx][pvIdx], pvIdx), numDiffMethod);

                for (const auto& localDofJ : localDofs(fvGeometry))
                    for (int eqIdx = 0; eqIdx < numEq; eqIdx++)
                        A[localDofJ.dofIndex()][dofIdx][eqIdx][pvIdx] += partialDerivs[localDofJ.index()][eqIdx];

                deflectionPolicy.restore(localDof);
                elemSol[localIdx][pvIdx] = curSol[localDof.dofIndex()][pvIdx];
                this->asImp_().maybeUpdateCouplingContext(localDof, elemSol, pvIdx);
            }
        };

        for (const auto& localDof : localDofs(fvGeometry))
            assembleDerivative(localDof);

        this->asImp_().maybeEvalAdditionalDomainDerivatives(origResiduals, A, gridVariables);

        return origResiduals;
    }

}; // implicit CVFEAssembler with numeric Jacobian

} // end namespace Dumux

#endif
