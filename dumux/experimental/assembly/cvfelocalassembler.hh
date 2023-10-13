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
 * \ingroup CVFEDiscretization
 * \brief An assembler for Jacobian and residual contribution per element (CVFE methods)
 */
#ifndef DUMUX_EXPERIMENTAL_CVFE_LOCAL_ASSEMBLER_HH
#define DUMUX_EXPERIMENTAL_CVFE_LOCAL_ASSEMBLER_HH

#include <dune/common/exceptions.hh>
#include <dune/common/hybridutilities.hh>
#include <dune/common/reservedvector.hh>
#include <dune/grid/common/gridenums.hh>
#include <dune/istl/matrixindexset.hh>
#include <dune/istl/bvector.hh>

#include <dumux/common/typetraits/typetraits.hh>
#include <dumux/common/properties.hh>
#include <dumux/common/parameters.hh>
#include <dumux/common/numericdifferentiation.hh>

#include <dumux/assembly/numericepsilon.hh>
#include <dumux/assembly/diffmethod.hh>
#include <dumux/experimental/assembly/fvlocalassemblerbase.hh>
#include <dumux/assembly/partialreassembler.hh>
#include <dumux/assembly/entitycolor.hh>

#include <dumux/discretization/cvfe/elementsolution.hh>

#include <dumux/assembly/volvardeflectionhelper_.hh>

namespace Dumux::Experimental {

/*!
 * \ingroup Experimental
 * \ingroup Assembly
 * \ingroup CVFEDiscretization
 * \brief A base class for all local CVFE assemblers
 * \tparam TypeTag The TypeTag
 * \tparam Assembler The assembler type
 * \tparam Implementation The actual implementation
 */
template<class TypeTag, class Assembler, class Implementation>
class CVFELocalAssemblerBase : public Experimental::FVLocalAssemblerBase<TypeTag, Assembler, Implementation>
{
    using ParentType = Experimental::FVLocalAssemblerBase<TypeTag, Assembler, Implementation>;
    using JacobianMatrix = GetPropType<TypeTag, Properties::JacobianMatrix>;
    using GridVariables = GetPropType<TypeTag, Properties::GridVariables>;
    using Problem = GetPropType<TypeTag, Properties::Problem>;
    using ElementVolumeVariables = typename GridVariables::GridVolumeVariables::LocalView;
    using SolutionVector = typename Assembler::SolutionVector;

    static constexpr int numEq = GetPropType<TypeTag, Properties::ModelTraits>::numEq();
    static constexpr int dim = GetPropType<TypeTag, Properties::GridGeometry>::GridView::dimension;

public:

    using ParentType::ParentType;

    using LocalResidual = typename ParentType::LocalResidual;
    using ElementResidualVector = typename LocalResidual::ElementResidualVector;


    void bindLocalViews()
    {
        ParentType::bindLocalViews();
        this->elemBcTypes().update(this->asImp_().problem(), this->element(), this->fvGeometry());
    }

    /*!
     * \brief Computes the derivatives with respect to the given element and adds them
     *        to the global matrix. The element residual is written into the right hand side.
     */
    template <class ResidualVector, class StageParams, class PartialReassembler = DefaultPartialReassembler, class CouplingFunction = Noop>
    void assembleJacobianAndResidual(JacobianMatrix& jac, ResidualVector& res, GridVariables& gridVariables,
                                     const StageParams& stageParams, ResidualVector& temporal, ResidualVector& spatial,
                                     ResidualVector& constrainedDofs,
                                     const PartialReassembler* partialReassembler = nullptr,
                                     const CouplingFunction& maybeAssembleCouplingBlocks = noop)
    {
        this->asImp_().bindLocalViews();
        const auto eIdxGlobal = this->asImp_().problem().gridGeometry().elementMapper().index(this->element());

        this->localResidual().spatialWeight(1.0);
        this->localResidual().temporalWeight(1.0);

        const auto sWeight = stageParams.spatialWeight(stageParams.size()-1);
        const auto tWeight = stageParams.temporalWeight(stageParams.size()-1);

        const auto flux = this->evalLocalFluxAndSourceResidual(this->curElemVolVars());
        const auto storage = this->localResidual().evalStorage(this->fvGeometry(), this->curElemVolVars());
        ElementResidualVector origResidual(flux.size());
        origResidual = 0.0;
        for (const auto& scv : scvs(this->fvGeometry()))
        {
            spatial[scv.dofIndex()] += flux[scv.localDofIndex()];
            temporal[scv.dofIndex()] += storage[scv.localDofIndex()];
            origResidual[scv.localDofIndex()] += flux[scv.localDofIndex()]*sWeight + storage[scv.localDofIndex()]*tWeight;
            res[scv.dofIndex()] += origResidual[scv.localDofIndex()];
        }

        this->localResidual().spatialWeight(sWeight);
        this->localResidual().temporalWeight(tWeight);

        if (partialReassembler && partialReassembler->elementColor(eIdxGlobal) == EntityColor::green)
        {
            // assemble the coupling blocks for coupled models (does nothing if not coupled)
            maybeAssembleCouplingBlocks(origResidual);
        }
        else if (!this->elementIsGhost())
        {
            this->asImp_().assembleJacobian(jac, gridVariables, origResidual, partialReassembler); // forward to the internal implementation

            // assemble the coupling blocks for coupled models (does nothing if not coupled)
            maybeAssembleCouplingBlocks(origResidual);
        }
        else
        {
            // Treatment of ghost elements
            assert(this->elementIsGhost());

            // handle dofs per codimension
            const auto& gridGeometry = this->asImp_().problem().gridGeometry();
            Dune::Hybrid::forEach(std::make_integer_sequence<int, dim>{}, [&](auto d)
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

                    // WARNING: this only works if the mapping from codim+subEntity to
                    // global dofIndex is unique (on dof per entity of this codim).
                    // For more general mappings, we should use a proper local-global mapping here.
                    // For example through dune-functions.
                    const auto dofIndex = gridGeometry.dofMapper().index(entity);

                    // this might be a vector-valued dof
                    using BlockType = typename JacobianMatrix::block_type;
                    BlockType &J = jac[dofIndex][dofIndex];
                    for (int j = 0; j < BlockType::rows; ++j)
                        J[j][j] = 1.0;

                    // set residual for the ghost dof
                    res[dofIndex] = 0;
                    constrainedDofs[dofIndex] = 1;
                }
            });
        }

        auto applyDirichlet = [&] (const auto& scvI,
                                   const auto& dirichletValues,
                                   const auto eqIdx,
                                   const auto pvIdx)
        {
            res[scvI.dofIndex()][eqIdx] = this->curElemVolVars()[scvI].priVars()[pvIdx] - dirichletValues[pvIdx];
            constrainedDofs[scvI.dofIndex()][eqIdx] = 1;

            auto& row = jac[scvI.dofIndex()];
            for (auto col = row.begin(); col != row.end(); ++col)
                row[col.index()][eqIdx] = 0.0;

            jac[scvI.dofIndex()][scvI.dofIndex()][eqIdx][pvIdx] = 1.0;

            // if a periodic dof has Dirichlet values also apply the same Dirichlet values to the other dof
            if (this->asImp_().problem().gridGeometry().dofOnPeriodicBoundary(scvI.dofIndex()))
            {
                const auto periodicDof = this->asImp_().problem().gridGeometry().periodicallyMappedDof(scvI.dofIndex());
                res[periodicDof][eqIdx] = this->curElemVolVars()[scvI].priVars()[pvIdx] - dirichletValues[pvIdx];
                constrainedDofs[periodicDof][eqIdx] = 1;
                const auto end = jac[periodicDof].end();
                for (auto it = jac[periodicDof].begin(); it != end; ++it)
                    (*it) = periodicDof != it.index() ? 0.0 : 1.0;
            }
        };

        this->asImp_().enforceDirichletConstraints(applyDirichlet);
    }

    /*!
     * \brief Computes the derivatives with respect to the given element and adds them
     *        to the global matrix.
     */
    void assembleJacobian(JacobianMatrix& jac, GridVariables& gridVariables)
    {
        this->asImp_().bindLocalViews();
        this->asImp_().assembleJacobian(jac, gridVariables); // forward to the internal implementation

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

        this->asImp_().enforceDirichletConstraints(applyDirichlet);
    }

    /*!
     * \brief Assemble the residual only
     */
    template <class ResidualVector>
    void assembleResidual(ResidualVector& res)
    {
        this->asImp_().bindLocalViews();
        const auto residual = this->evalLocalResidual();

        for (const auto& scv : scvs(this->fvGeometry()))
            res[scv.dofIndex()] += residual[scv.localDofIndex()];

        auto applyDirichlet = [&] (const auto& scvI,
                                   const auto& dirichletValues,
                                   const auto eqIdx,
                                   const auto pvIdx)
        {
            res[scvI.dofIndex()][eqIdx] = this->curElemVolVars()[scvI].priVars()[pvIdx] - dirichletValues[pvIdx];
        };

        this->asImp_().enforceDirichletConstraints(applyDirichlet);
    }

    /*!
     * \brief Assemble the residual only
     */
    template<class ResidualVector>
    void assembleCurrentResidual(ResidualVector& spatialRes, ResidualVector& temporalRes)
    {
        this->asImp_().bindLocalViews();
        const auto flux = this->evalLocalFluxAndSourceResidual(this->curElemVolVars());
        const auto storage = this->localResidual().evalStorage(this->fvGeometry(), this->curElemVolVars());
        for (const auto& scv : scvs(this->fvGeometry()))
        {
            spatialRes[scv.dofIndex()] += flux[scv.localDofIndex()];
            temporalRes[scv.dofIndex()] += storage[scv.localDofIndex()];
        }
    }

    //! Enforce Dirichlet constraints
    template<typename ApplyFunction>
    void enforceDirichletConstraints(const ApplyFunction& applyDirichlet)
    {
        // enforce Dirichlet boundary conditions
        this->asImp_().evalDirichletBoundaries(applyDirichlet);
        // take care of internal Dirichlet constraints (if enabled)
        this->asImp_().enforceInternalDirichletConstraints(applyDirichlet);
    }

    /*!
     * \brief Evaluates Dirichlet boundaries
     */
    template< typename ApplyDirichletFunctionType >
    void evalDirichletBoundaries(ApplyDirichletFunctionType applyDirichlet)
    {
        // enforce Dirichlet boundaries by overwriting partial derivatives with 1 or 0
        // and set the residual to (privar - dirichletvalue)
        if (this->elemBcTypes().hasDirichlet())
        {
            for (const auto& scvI : scvs(this->fvGeometry()))
            {
                const auto bcTypes = this->elemBcTypes().get(this->fvGeometry(), scvI);
                if (bcTypes.hasDirichlet())
                {
                    const auto dirichletValues = this->asImp_().problem().dirichlet(this->element(), scvI);

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
 * \ingroup CVFEDiscretization
 * \brief An assembler for Jacobian and residual contribution per element (CVFE methods)
 * \tparam TypeTag The TypeTag
 * \tparam diffMethod The differentiation method to residual compute derivatives
 * \tparam Implementation via CRTP
 */
template<class TypeTag, class Assembler, DiffMethod diffMethod = DiffMethod::numeric, class Implementation = void>
class CVFELocalAssembler;

/*!
 * \ingroup Experimental
 * \ingroup Assembly
 * \ingroup CVFEDiscretization
 * \brief Control volume finite element local assembler using numeric differentiation
 */
template<class TypeTag, class Assembler, class Implementation>
class CVFELocalAssembler<TypeTag, Assembler, DiffMethod::numeric, Implementation>
: public CVFELocalAssemblerBase<TypeTag, Assembler,
                                NonVoidOr<CVFELocalAssembler<TypeTag, Assembler, DiffMethod::numeric, Implementation>, Implementation>>
{
    using ThisType = CVFELocalAssembler<TypeTag, Assembler, DiffMethod::numeric, Implementation>;
    using ParentType = CVFELocalAssemblerBase<TypeTag, Assembler, NonVoidOr<ThisType, Implementation>>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using GridVariables = GetPropType<TypeTag, Properties::GridVariables>;
    using VolumeVariables = GetPropType<TypeTag, Properties::VolumeVariables>;
    using JacobianMatrix = GetPropType<TypeTag, Properties::JacobianMatrix>;

    static constexpr int numEq = GetPropType<TypeTag, Properties::ModelTraits>::numEq();
    static constexpr int dim = GetPropType<TypeTag, Properties::GridGeometry>::GridView::dimension;

    static constexpr bool enableGridFluxVarsCache
        = GridVariables::GridFluxVariablesCache::cachingEnabled;
    static constexpr bool solutionDependentFluxVarsCache
        = GridVariables::GridFluxVariablesCache::FluxVariablesCache::isSolDependent;

public:

    using LocalResidual = typename ParentType::LocalResidual;
    using ElementResidualVector = typename LocalResidual::ElementResidualVector;
    using ParentType::ParentType;

    template <class PartialReassembler = DefaultPartialReassembler>
    void assembleJacobian(JacobianMatrix& A, GridVariables& gridVariables,
                          const ElementResidualVector& origResiduals,
                          const PartialReassembler* partialReassembler = nullptr)
    {
        if (this->isImplicit())
            assembleJacobianImplicit_(A, gridVariables, origResiduals, partialReassembler);
        else
            assembleJacobianExplicit_(A, gridVariables, origResiduals, partialReassembler);
    }

private:
    /*!
     * \brief Computes the derivatives with respect to the given element and adds them
     *        to the global matrix.
     *
     * \return The element residual at the current solution.
     */
    template <class PartialReassembler = DefaultPartialReassembler>
    void assembleJacobianImplicit_(JacobianMatrix& A, GridVariables& gridVariables,
                                   const ElementResidualVector& origResiduals,
                                   const PartialReassembler* partialReassembler = nullptr)
    {
        // get some aliases for convenience
        const auto& element = this->element();
        const auto& fvGeometry = this->fvGeometry();
        const auto& curSol = this->asImp_().curSol();

        auto&& curElemVolVars = this->curElemVolVars();
        auto&& elemFluxVarsCache = this->elemFluxVarsCache();

        //////////////////////////////////////////////////////////////////////////////////////////////////
        //                                                                                              //
        // Calculate derivatives of all dofs in stencil with respect to the dofs in the element. In the //
        // neighboring elements we do so by computing the derivatives of the fluxes which depend on the //
        // actual element. In the actual element we evaluate the derivative of the entire residual.     //
        //                                                                                              //
        //////////////////////////////////////////////////////////////////////////////////////////////////

        // if all volvars in the stencil have to be updated or if it's enough to only update the
        // volVars for the scv whose associated dof has been deflected
        static const bool updateAllVolVars = getParamFromGroup<bool>(
            this->asImp_().problem().paramGroup(), "Assembly.BoxVolVarsDependOnAllElementDofs", false
        );

        // create the element solution
        auto elemSol = elementSolution(element, curSol, fvGeometry.gridGeometry());

        // create the vector storing the partial derivatives
        ElementResidualVector partialDerivs(fvGeometry.numScv());

        Dumux::Detail::VolVarsDeflectionHelper deflectionHelper(
            [&] (const auto& scv) -> VolumeVariables& {
                return this->getVolVarAccess(gridVariables.curGridVolVars(), curElemVolVars, scv);
            },
            fvGeometry,
            updateAllVolVars
        );

        // calculation of the derivatives
        for (const auto& scv : scvs(fvGeometry))
        {
            // dof index and corresponding actual pri vars
            const auto dofIdx = scv.dofIndex();
            deflectionHelper.setCurrent(scv);

            // calculate derivatives w.r.t to the privars at the dof at hand
            for (int pvIdx = 0; pvIdx < numEq; pvIdx++)
            {
                partialDerivs = 0.0;

                auto evalResiduals = [&](Scalar priVar)
                {
                    // update the volume variables and compute element residual
                    elemSol[scv.localDofIndex()][pvIdx] = priVar;
                    deflectionHelper.deflect(elemSol, scv, this->asImp_().problem());
                    if constexpr (solutionDependentFluxVarsCache)
                    {
                        elemFluxVarsCache.update(element, fvGeometry, curElemVolVars);
                        if constexpr (enableGridFluxVarsCache)
                            gridVariables.gridFluxVarsCache().updateElement(element, fvGeometry, curElemVolVars);
                    }
                    this->asImp_().maybeUpdateCouplingContext(scv, elemSol, pvIdx);
                    return this->evalLocalResidual();
                };

                // derive the residuals numerically
                static const NumericEpsilon<Scalar, numEq> eps_{this->asImp_().problem().paramGroup()};
                static const int numDiffMethod = getParamFromGroup<int>(this->asImp_().problem().paramGroup(), "Assembly.NumericDifferenceMethod");
                NumericDifferentiation::partialDerivative(evalResiduals, elemSol[scv.localDofIndex()][pvIdx], partialDerivs, origResiduals,
                                                          eps_(elemSol[scv.localDofIndex()][pvIdx], pvIdx), numDiffMethod);

                // update the global stiffness matrix with the current partial derivatives
                for (const auto& scvJ : scvs(fvGeometry))
                {
                    // don't add derivatives for green dofs
                    if (!partialReassembler
                        || partialReassembler->dofColor(scvJ.dofIndex()) != EntityColor::green)
                    {
                        for (int eqIdx = 0; eqIdx < numEq; eqIdx++)
                        {
                            // A[i][col][eqIdx][pvIdx] is the rate of change of
                            // the residual of equation 'eqIdx' at dof 'i'
                            // depending on the primary variable 'pvIdx' at dof
                            // 'col'.
                            A[scvJ.dofIndex()][dofIdx][eqIdx][pvIdx] += partialDerivs[scvJ.localDofIndex()][eqIdx];
                        }
                    }
                }

                // restore the original state of the scv's volume variables
                deflectionHelper.restore(scv);

                // restore the original element solution
                elemSol[scv.localDofIndex()][pvIdx] = curSol[scv.dofIndex()][pvIdx];
                this->asImp_().maybeUpdateCouplingContext(scv, elemSol, pvIdx);
            }
        }

        // restore original state of the flux vars cache in case of global caching.
        // In the case of local caching this is obsolete because the elemFluxVarsCache used here goes out of scope after this.
        if constexpr (enableGridFluxVarsCache)
            gridVariables.gridFluxVarsCache().updateElement(element, fvGeometry, curElemVolVars);

        // evaluate additional derivatives that might arise from the coupling (no-op if not coupled)
        this->asImp_().maybeEvalAdditionalDomainDerivatives(origResiduals, A, gridVariables);
    }

    /*!
     * \brief Computes the derivatives with respect to the given element and adds them
     *        to the global matrix.
     * \return The element residual at the current solution.
     */
    template <class PartialReassembler = DefaultPartialReassembler>
    void assembleJacobianExplicit_(JacobianMatrix& A, GridVariables& gridVariables,
                                   const ElementResidualVector& origResiduals,
                                   const PartialReassembler* partialReassembler = nullptr)
    {
        if (partialReassembler)
            DUNE_THROW(Dune::NotImplemented, "partial reassembly for explicit time discretization");

        // get some aliases for convenience
        const auto& element = this->element();
        const auto& fvGeometry = this->fvGeometry();
        const auto& curSol = this->asImp_().curSol();
        auto&& curElemVolVars = this->curElemVolVars();

        // create the element solution
        auto elemSol = elementSolution(element, curSol, fvGeometry.gridGeometry());

        // create the vector storing the partial derivatives
        ElementResidualVector partialDerivs(fvGeometry.numScv());

        // calculation of the derivatives
        for (const auto& scv : scvs(fvGeometry))
        {
            // dof index and corresponding actual pri vars
            const auto dofIdx = scv.dofIndex();
            auto& curVolVars = this->getVolVarAccess(gridVariables.curGridVolVars(), curElemVolVars, scv);
            const VolumeVariables origVolVars(curVolVars);

            // calculate derivatives w.r.t to the privars at the dof at hand
            for (int pvIdx = 0; pvIdx < numEq; pvIdx++)
            {
                partialDerivs = 0.0;

                auto evalStorage = [&](Scalar priVar)
                {
                    elemSol[scv.localDofIndex()][pvIdx] = priVar;
                    curVolVars.update(elemSol, this->asImp_().problem(), element, scv);
                    return this->evalStorage();
                };

                // derive the residuals numerically
                static const NumericEpsilon<Scalar, numEq> eps_{this->asImp_().problem().paramGroup()};
                static const int numDiffMethod = getParamFromGroup<int>(this->asImp_().problem().paramGroup(), "Assembly.NumericDifferenceMethod");
                NumericDifferentiation::partialDerivative(evalStorage, elemSol[scv.localDofIndex()][pvIdx], partialDerivs, origResiduals,
                                                          eps_(elemSol[scv.localDofIndex()][pvIdx], pvIdx), numDiffMethod);

                // update the global stiffness matrix with the current partial derivatives
                for (int eqIdx = 0; eqIdx < numEq; eqIdx++)
                {
                    // A[i][col][eqIdx][pvIdx] is the rate of change of
                    // the residual of equation 'eqIdx' at dof 'i'
                    // depending on the primary variable 'pvIdx' at dof
                    // 'col'.
                    A[dofIdx][dofIdx][eqIdx][pvIdx] += partialDerivs[scv.localDofIndex()][eqIdx];
                }

                // restore the original state of the scv's volume variables
                curVolVars = origVolVars;

                // restore the original element solution
                elemSol[scv.localDofIndex()][pvIdx] = curSol[scv.dofIndex()][pvIdx];
            }
        }
    }
};

} // end namespace Dumux

#endif
