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
#ifndef DUMUX_CVFE_LOCAL_ASSEMBLER_HH
#define DUMUX_CVFE_LOCAL_ASSEMBLER_HH

#include <dune/common/exceptions.hh>
#include <dune/common/hybridutilities.hh>
#include <dune/common/reservedvector.hh>
#include <dune/grid/common/gridenums.hh>
#include <dune/istl/matrixindexset.hh>
#include <dune/istl/bvector.hh>

#include <dumux/common/properties.hh>
#include <dumux/common/parameters.hh>
#include <dumux/common/numericdifferentiation.hh>
#include <dumux/common/typetraits/localdofs_.hh>
#include <dumux/common/typetraits/boundary_.hh>

#include <dumux/assembly/numericepsilon.hh>
#include <dumux/assembly/diffmethod.hh>
#include <dumux/assembly/fvlocalassemblerbase.hh>
#include <dumux/assembly/partialreassembler.hh>
#include <dumux/assembly/entitycolor.hh>

#include <dumux/discretization/cvfe/elementsolution.hh>
#include <dumux/discretization/cvfe/localdof.hh>

#include "cvfevolvarsdeflectionpolicy_.hh"

namespace Dumux {

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
 * \tparam implicit Specifies whether the time discretization is implicit or not not (i.e. explicit)
 */
template<class TypeTag, class Assembler, class Implementation, bool implicit>
class CVFELocalAssemblerBase : public FVLocalAssemblerBase<TypeTag, Assembler, Implementation, implicit>
{
    using ParentType = FVLocalAssemblerBase<TypeTag, Assembler, Implementation, implicit>;
    using Problem = GetPropType<TypeTag, Properties::Problem>;
    using JacobianMatrix = GetPropType<TypeTag, Properties::JacobianMatrix>;
    using GridVariables = GetPropType<TypeTag, Properties::GridVariables>;
    using GridVolumeVariables = typename GridVariables::GridVolumeVariables;
    using ElementVolumeVariables = typename GridVariables::GridVolumeVariables::LocalView;
    using VolumeVariables = GetPropType<TypeTag, Properties::VolumeVariables>;
    using SolutionVector = typename Assembler::SolutionVector;
    using GridGeometry =  GetPropType<TypeTag, Properties::GridGeometry>;
    using FVElementGeometry = typename GridGeometry::LocalView;

    static constexpr int numEq = GetPropType<TypeTag, Properties::ModelTraits>::numEq();
    static constexpr int dim = GridGeometry::GridView::dimension;

public:

    using ParentType::ParentType;

    void bindLocalViews()
    {
        ParentType::bindLocalViews();
        this->elemBcTypes().update(this->asImp_().problem(), this->element(), this->fvGeometry());
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
                }
            });
        }

        auto applyDirichlet = [&] (const auto& scvOrLocalDofI,
                                   const auto& dirichletValues,
                                   const auto eqIdx,
                                   const auto pvIdx)
        {
            res[scvOrLocalDofI.dofIndex()][eqIdx] = this->curElemVolVars()[scvOrLocalDofI].priVars()[pvIdx] - dirichletValues[pvIdx];

            auto& row = jac[scvOrLocalDofI.dofIndex()];
            for (auto col = row.begin(); col != row.end(); ++col)
                row[col.index()][eqIdx] = 0.0;

            jac[scvOrLocalDofI.dofIndex()][scvOrLocalDofI.dofIndex()][eqIdx][pvIdx] = 1.0;

            // if a periodic dof has Dirichlet values also apply the same Dirichlet values to the other dof
            if (this->asImp_().problem().gridGeometry().dofOnPeriodicBoundary(scvOrLocalDofI.dofIndex()))
            {
                const auto periodicDof = this->asImp_().problem().gridGeometry().periodicallyMappedDof(scvOrLocalDofI.dofIndex());
                res[periodicDof][eqIdx] = this->asImp_().curSol()[periodicDof][pvIdx] - dirichletValues[pvIdx];

                auto& rowP = jac[periodicDof];
                for (auto col = rowP.begin(); col != rowP.end(); ++col)
                    rowP[col.index()][eqIdx] = 0.0;

                rowP[periodicDof][eqIdx][pvIdx] = 1.0;
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
        this->asImp_().assembleJacobianAndResidualImpl(jac, gridVariables); // forward to the internal implementation

        auto applyDirichlet = [&] (const auto& scvOrLocalDofI,
                                   const auto& dirichletValues,
                                   const auto eqIdx,
                                   const auto pvIdx)
        {
            auto& row = jac[scvOrLocalDofI.dofIndex()];
            for (auto col = row.begin(); col != row.end(); ++col)
                row[col.index()][eqIdx] = 0.0;

            jac[scvOrLocalDofI.dofIndex()][scvOrLocalDofI.dofIndex()][eqIdx][pvIdx] = 1.0;
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

        for (const auto& localDof : localDofs(this->fvGeometry()))
            res[localDof.dofIndex()] += residual[localDof.index()];

        auto applyDirichlet = [&] (const auto& scvOrLocalDofI,
                                   const auto& dirichletValues,
                                   const auto eqIdx,
                                   const auto pvIdx)
        {
            res[scvOrLocalDofI.dofIndex()][eqIdx] = this->curElemVolVars()[scvOrLocalDofI].priVars()[pvIdx] - dirichletValues[pvIdx];
        };

        this->asImp_().enforceDirichletConstraints(applyDirichlet);
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
        // when having the new boundary interface Dirichlet conditions are incorporated via constraints
        if constexpr (!Detail::hasProblemBoundaryTypesForIntersectionFunction<Problem, FVElementGeometry, typename GridGeometry::GridView::Intersection>())
        {
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
 * \tparam diffMethod The differentiation method to residual compute derivatives
 * \tparam implicit Specifies whether the time discretization is implicit or not not (i.e. explicit)
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
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using GridVariables = GetPropType<TypeTag, Properties::GridVariables>;
    using GridVolumeVariables = GetPropType<TypeTag, Properties::GridVolumeVariables>;
    using ElementVolumeVariables = typename GridVariables::GridVolumeVariables::LocalView;
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

        auto&& curElemVolVars = this->curElemVolVars();
        auto&& elemFluxVarsCache = this->elemFluxVarsCache();

        // get the vector of the actual element residuals
        const auto origResiduals = this->evalLocalResidual();

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
        ElementResidualVector partialDerivs(Detail::LocalDofs::numLocalDofs(fvGeometry));

        auto deflectionPolicy = Detail::CVFE::makeVariablesDeflectionPolicy(
            gridVariables.curGridVolVars(),
            curElemVolVars,
            fvGeometry,
            updateAllVolVars
        );

        auto assembleDerivative = [&, this](const auto& localDof)
        {
            // dof index and corresponding actual pri vars
            const auto dofIdx = localDof.dofIndex();
            const auto localIdx = localDof.index();
            deflectionPolicy.store(localDof);

            // calculate derivatives w.r.t to the privars at the dof at hand
            for (int pvIdx = 0; pvIdx < numEq; pvIdx++)
            {
                partialDerivs = 0.0;

                auto evalResiduals = [&](Scalar priVar)
                {
                    // update the volume variables and compute element residual
                    elemSol[localIdx][pvIdx] = priVar;
                    deflectionPolicy.update(elemSol, localDof, this->asImp_().problem());
                    if constexpr (solutionDependentFluxVarsCache)
                    {
                        elemFluxVarsCache.update(element, fvGeometry, curElemVolVars);
                        if constexpr (enableGridFluxVarsCache)
                            gridVariables.gridFluxVarsCache().updateElement(element, fvGeometry, curElemVolVars);
                    }
                    this->asImp_().maybeUpdateCouplingContext(localDof, elemSol, pvIdx);
                    return this->evalLocalResidual();
                };

                // derive the residuals numerically
                static const NumericEpsilon<Scalar, numEq> eps_{this->asImp_().problem().paramGroup()};
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

                // restore the original state of the scv's volume variables
                deflectionPolicy.restore(localDof);

                // restore the original element solution
                elemSol[localIdx][pvIdx] = curSol[localDof.dofIndex()][pvIdx];
                this->asImp_().maybeUpdateCouplingContext(localDof, elemSol, pvIdx);
            }
        };

        // calculation of the derivatives
        for (const auto& localDof : localDofs(fvGeometry))
            assembleDerivative(localDof);

        // restore original state of the flux vars cache in case of global caching.
        // In the case of local caching this is obsolete because the elemFluxVarsCache used here goes out of scope after this.
        if constexpr (enableGridFluxVarsCache)
            gridVariables.gridFluxVarsCache().updateElement(element, fvGeometry, curElemVolVars);

        // evaluate additional derivatives that might arise from the coupling (no-op if not coupled)
        this->asImp_().maybeEvalAdditionalDomainDerivatives(origResiduals, A, gridVariables);

        return origResiduals;
    }

}; // implicit CVFEAssembler with numeric Jacobian

/*!
 * \ingroup Assembly
 * \ingroup CVFEDiscretization
 * \brief Control volume finite element local assembler using numeric differentiation and explicit time discretization
 */
template<class TypeTag, class Assembler, class Implementation>
class CVFELocalAssembler<TypeTag, Assembler, DiffMethod::numeric, /*implicit=*/false, Implementation>
: public CVFELocalAssemblerBase<TypeTag, Assembler,
                                Detail::CVFE::Impl<Implementation, CVFELocalAssembler<TypeTag, Assembler, DiffMethod::numeric, false, Implementation>>,
                                false>
{
    using ThisType = CVFELocalAssembler<TypeTag, Assembler, DiffMethod::numeric, false, Implementation>;
    using ParentType = CVFELocalAssemblerBase<TypeTag, Assembler, Detail::CVFE::Impl<Implementation, ThisType>, false>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using GridVariables = GetPropType<TypeTag, Properties::GridVariables>;
    using VolumeVariables = GetPropType<TypeTag, Properties::VolumeVariables>;
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
        if (partialReassembler)
            DUNE_THROW(Dune::NotImplemented, "partial reassembly for explicit time discretization");

        // get some aliases for convenience
        const auto& element = this->element();
        const auto& fvGeometry = this->fvGeometry();
        const auto& curSol = this->asImp_().curSol();
        auto&& curElemVolVars = this->curElemVolVars();

        // get the vecor of the actual element residuals
        const auto origResiduals = this->evalLocalResidual();
        const auto origStorageResiduals = this->evalLocalStorageResidual();

        //////////////////////////////////////////////////////////////////////////////////////////////////
        //                                                                                              //
        // Calculate derivatives of all dofs in stencil with respect to the dofs in the element. In the //
        // neighboring elements we do so by computing the derivatives of the fluxes which depend on the //
        // actual element. In the actual element we evaluate the derivative of the entire residual.     //
        //                                                                                              //
        //////////////////////////////////////////////////////////////////////////////////////////////////

        // create the element solution
        auto elemSol = elementSolution(element, curSol, fvGeometry.gridGeometry());

        // create the vector storing the partial derivatives
        ElementResidualVector partialDerivs(Detail::LocalDofs::numLocalDofs(fvGeometry));

        // calculation of the derivatives
        for (const auto& scv : scvs(fvGeometry))
        {
            // dof index and corresponding actual pri vars
            const auto localIdx = scv.localDofIndex();
            const auto dofIdx = scv.dofIndex();
            auto& curVolVars = this->getVolVarAccess(gridVariables.curGridVolVars(), curElemVolVars, scv);
            const VolumeVariables origVolVars(curVolVars);

            // calculate derivatives w.r.t to the privars at the dof at hand
            for (int pvIdx = 0; pvIdx < numEq; pvIdx++)
            {
                partialDerivs = 0.0;

                auto evalStorage = [&](Scalar priVar)
                {
                    // auto partialDerivsTmp = partialDerivs;
                    elemSol[localIdx][pvIdx] = priVar;
                    curVolVars.update(elemSol, this->asImp_().problem(), element, scv);
                    return this->evalLocalStorageResidual();
                };

                // derive the residuals numerically
                static const NumericEpsilon<Scalar, numEq> eps_{this->asImp_().problem().paramGroup()};
                static const int numDiffMethod = getParamFromGroup<int>(this->asImp_().problem().paramGroup(), "Assembly.NumericDifferenceMethod");
                NumericDifferentiation::partialDerivative(evalStorage, elemSol[localIdx][pvIdx], partialDerivs, origStorageResiduals,
                                                          eps_(elemSol[localIdx][pvIdx], pvIdx), numDiffMethod);

                // update the global stiffness matrix with the current partial derivatives
                for (int eqIdx = 0; eqIdx < numEq; eqIdx++)
                {
                    // A[i][col][eqIdx][pvIdx] is the rate of change of
                    // the residual of equation 'eqIdx' at dof 'i'
                    // depending on the primary variable 'pvIdx' at dof
                    // 'col'.
                    A[dofIdx][dofIdx][eqIdx][pvIdx] += partialDerivs[localIdx][eqIdx];
                }

                // restore the original state of the scv's volume variables
                curVolVars = origVolVars;

                // restore the original element solution
                elemSol[localIdx][pvIdx] = curSol[dofIdx][pvIdx];
                // TODO additional dof dependencies
            }
        }

        return origResiduals;
    }
}; // explicit CVFEAssembler with numeric Jacobian

/*!
 * \ingroup Assembly
 * \ingroup CVFEDiscretization
 * \brief Control volume finite element local assembler using analytic differentiation and implicit time discretization
 */
template<class TypeTag, class Assembler, class Implementation>
class CVFELocalAssembler<TypeTag, Assembler, DiffMethod::analytic, /*implicit=*/true, Implementation>
: public CVFELocalAssemblerBase<TypeTag, Assembler,
                                Detail::CVFE::Impl<Implementation, CVFELocalAssembler<TypeTag, Assembler, DiffMethod::analytic, true, Implementation>>,
                                true>
{
    using ThisType = CVFELocalAssembler<TypeTag, Assembler, DiffMethod::analytic, true, Implementation>;
    using ParentType = CVFELocalAssemblerBase<TypeTag, Assembler, Detail::CVFE::Impl<Implementation, ThisType>, true>;
    using GridVariables = GetPropType<TypeTag, Properties::GridVariables>;
    using JacobianMatrix = GetPropType<TypeTag, Properties::JacobianMatrix>;

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
        if (partialReassembler)
            DUNE_THROW(Dune::NotImplemented, "partial reassembly for analytic differentiation");

        // get some aliases for convenience
        const auto& element = this->element();
        const auto& fvGeometry = this->fvGeometry();
        const auto& problem = this->asImp_().problem();
        const auto& curElemVolVars = this->curElemVolVars();
        const auto& elemFluxVarsCache = this->elemFluxVarsCache();

        // get the vecor of the actual element residuals
        const auto origResiduals = this->evalLocalResidual();

        //////////////////////////////////////////////////////////////////////////////////////////////////
        //                                                                                              //
        // Calculate derivatives of all dofs in stencil with respect to the dofs in the element. In the //
        // neighboring elements we do so by computing the derivatives of the fluxes which depend on the //
        // actual element. In the actual element we evaluate the derivative of the entire residual.     //
        //                                                                                              //
        //////////////////////////////////////////////////////////////////////////////////////////////////

        // calculation of the source and storage derivatives
        for (const auto& scv : scvs(fvGeometry))
        {
            // dof index and corresponding actual pri vars
            const auto dofIdx = scv.dofIndex();
            const auto& volVars = curElemVolVars[scv];

            // derivative of this scv residual w.r.t the d.o.f. of the same scv (because of mass lumping)
            // only if the problem is instationary we add derivative of storage term
            // TODO if e.g. porosity depends on all dofs in the element, we would have off-diagonal matrix entries!?
            if (!this->assembler().isStationaryProblem())
                this->localResidual().addStorageDerivatives(A[dofIdx][dofIdx],
                                                            problem,
                                                            element,
                                                            fvGeometry,
                                                            volVars,
                                                            scv);

            // derivative of this scv residual w.r.t the d.o.f. of the same scv (because of mass lumping)
            // add source term derivatives
            this->localResidual().addSourceDerivatives(A[dofIdx][dofIdx],
                                                       problem,
                                                       element,
                                                       fvGeometry,
                                                       volVars,
                                                       scv);
        }

        // localJacobian[scvIdx][otherScvIdx][eqIdx][priVarIdx] of the fluxes
        for (const auto& scvf : scvfs(fvGeometry))
        {
            if (!scvf.boundary())
            {
                // add flux term derivatives
                this->localResidual().addFluxDerivatives(A,
                                                         problem,
                                                         element,
                                                         fvGeometry,
                                                         curElemVolVars,
                                                         elemFluxVarsCache,
                                                         scvf);
            }

            // the boundary gets special treatment to simplify
            // for the user
            else
            {
                const auto& insideScv = fvGeometry.scv(scvf.insideScvIdx());
                if (this->elemBcTypes().get(fvGeometry, insideScv).hasNeumann())
                {
                    // add flux term derivatives
                    this->localResidual().addRobinFluxDerivatives(A[insideScv.dofIndex()],
                                                                  problem,
                                                                  element,
                                                                  fvGeometry,
                                                                  curElemVolVars,
                                                                  elemFluxVarsCache,
                                                                  scvf);
                }
            }
        }

        return origResiduals;
    }

}; // implicit CVFEAssembler with analytic Jacobian

/*!
 * \ingroup Assembly
 * \ingroup CVFEDiscretization
 * \brief Control volume finite element local assembler using analytic differentiation and explicit time discretization
 */
template<class TypeTag, class Assembler, class Implementation>
class CVFELocalAssembler<TypeTag, Assembler, DiffMethod::analytic, /*implicit=*/false, Implementation>
: public CVFELocalAssemblerBase<TypeTag, Assembler,
                                Detail::CVFE::Impl<Implementation, CVFELocalAssembler<TypeTag, Assembler, DiffMethod::analytic, false, Implementation>>,
                                false>
{
    using ThisType = CVFELocalAssembler<TypeTag, Assembler, DiffMethod::analytic, false, Implementation>;
    using ParentType = CVFELocalAssemblerBase<TypeTag, Assembler, Detail::CVFE::Impl<Implementation, ThisType>, false>;
    using GridVariables = GetPropType<TypeTag, Properties::GridVariables>;
    using JacobianMatrix = GetPropType<TypeTag, Properties::JacobianMatrix>;

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
        if (partialReassembler)
            DUNE_THROW(Dune::NotImplemented, "partial reassembly for explicit time discretization");

        // get some aliases for convenience
        const auto& element = this->element();
        const auto& fvGeometry = this->fvGeometry();
        const auto& problem = this->asImp_().problem();
        const auto& curElemVolVars = this->curElemVolVars();

        // get the vecor of the actual element residuals
        const auto origResiduals = this->evalLocalResidual();

        //////////////////////////////////////////////////////////////////////////////////////////////////
        //                                                                                              //
        // Calculate derivatives of all dofs in stencil with respect to the dofs in the element. In the //
        // neighboring elements we do so by computing the derivatives of the fluxes which depend on the //
        // actual element. In the actual element we evaluate the derivative of the entire residual.     //
        //                                                                                              //
        //////////////////////////////////////////////////////////////////////////////////////////////////

        // calculation of the source and storage derivatives
        for (const auto& scv : scvs(fvGeometry))
        {
            // dof index and corresponding actual pri vars
            const auto dofIdx = scv.dofIndex();
            const auto& volVars = curElemVolVars[scv];

            // derivative of this scv residual w.r.t the d.o.f. of the same scv (because of mass lumping)
            // only if the problem is instationary we add derivative of storage term
            this->localResidual().addStorageDerivatives(A[dofIdx][dofIdx],
                                                        problem,
                                                        element,
                                                        fvGeometry,
                                                        volVars,
                                                        scv);
        }

        return origResiduals;
    }

}; // explicit CVFEAssembler with analytic Jacobian

} // end namespace Dumux

#endif
