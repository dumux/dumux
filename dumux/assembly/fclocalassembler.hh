// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup Assembly
 * \ingroup FaceCenteredStaggeredDiscretization
 * \brief An assembler for Jacobian and residual contribution per element (face-centered staggered methods)
 */
#ifndef DUMUX_FC_LOCAL_ASSEMBLER_HH
#define DUMUX_FC_LOCAL_ASSEMBLER_HH

#include <dune/grid/common/gridenums.hh>

#include <dumux/common/properties.hh>
#include <dumux/common/parameters.hh>
#include <dumux/common/numericdifferentiation.hh>
#include <dumux/assembly/numericepsilon.hh>
#include <dumux/assembly/diffmethod.hh>
#include <dumux/assembly/fvlocalassemblerbase.hh>
#include <dumux/assembly/entitycolor.hh>
#include <dumux/assembly/partialreassembler.hh>
#include <dumux/discretization/facecentered/staggered/elementsolution.hh>

namespace Dumux {

namespace Detail {

struct NoOpFunctor
{
    template<class... Args>
    constexpr void operator()(Args&&...) const {}
};

template<class T, class Default>
using NonVoidOrDefault_t = std::conditional_t<!std::is_same_v<T, void>, T, Default>;

} // end namespace Detail

/*!
 * \ingroup Assembly
 * \ingroup FaceCenteredStaggeredDiscretization
 * \brief A base class for all local cell-centered assemblers
 * \tparam TypeTag The TypeTag
 * \tparam Assembler The assembler type
 * \tparam Implementation The actual implementation
 * \tparam implicit Specifies whether the time discretization is implicit or not (i.e. explicit)
 */
template<class TypeTag, class Assembler, class Implementation, bool implicit>
class FaceCenteredLocalAssemblerBase : public FVLocalAssemblerBase<TypeTag, Assembler, Implementation, implicit>
{
    using ParentType = FVLocalAssemblerBase<TypeTag, Assembler, Implementation, implicit>;
    using GridView = typename GetPropType<TypeTag, Properties::GridGeometry>::GridView;
    using JacobianMatrix = GetPropType<TypeTag, Properties::JacobianMatrix>;
    using GridVariables = GetPropType<TypeTag, Properties::GridVariables>;
    using PrimaryVariables = GetPropType<TypeTag, Properties::PrimaryVariables>;

    static constexpr auto numEq = GetPropType<TypeTag, Properties::ModelTraits>::numEq();

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
    template <class ResidualVector, class PartialReassembler = DefaultPartialReassembler, class CouplingFunction = Detail::NoOpFunctor>
    void assembleJacobianAndResidual(JacobianMatrix& jac, ResidualVector& res, GridVariables& gridVariables,
                                     const PartialReassembler* partialReassembler,
                                     const CouplingFunction& maybeAssembleCouplingBlocks = CouplingFunction{})
    {
        static_assert(!std::decay_t<decltype(this->asImp_().problem())>::enableInternalDirichletConstraints(),
            "Internal Dirichlet constraints are currently not implemented for face-centered staggered models!");

        this->asImp_().bindLocalViews();
        const auto& gridGeometry = this->asImp_().problem().gridGeometry();
        const auto eIdxGlobal = gridGeometry.elementMapper().index(this->element());
        if (partialReassembler
            && partialReassembler->elementColor(eIdxGlobal) == EntityColor::green)
        {
            const auto residual = this->asImp_().evalLocalResidual(); // forward to the internal implementation
            for (const auto& scv : scvs(this->fvGeometry()))
                res[scv.dofIndex()] += residual[scv.localDofIndex()];

            // assemble the coupling blocks for coupled models (does nothing if not coupled)
            maybeAssembleCouplingBlocks(residual);
        }
        else if (!this->elementIsGhost())
        {
            const auto residual = this->asImp_().assembleJacobianAndResidualImpl(jac, gridVariables, partialReassembler); // forward to the internal implementation

            if (this->element().partitionType() == Dune::InteriorEntity)
            {
                for (const auto& scv : scvs(this->fvGeometry()))
                    res[scv.dofIndex()] += residual[scv.localDofIndex()];
            }
            else
            {
                // handle residual and matrix entries for parallel runs
                for (const auto& scv : scvs(this->fvGeometry()))
                {
                    const auto& facet = this->element().template subEntity <1> (scv.indexInElement());
                    // make sure that the residual at border entities is consistent by adding the
                    // the contribution from the neighboring overlap element's scv
                    if (facet.partitionType() == Dune::BorderEntity)
                        res[scv.dofIndex()] += residual[scv.localDofIndex()];

                    // set the matrix entries of all DOFs within the overlap region (except the border DOF)
                    // to 1.0 and the residual entries to 0.0
                    else
                    {
                        const auto idx = scv.dofIndex();
                        jac[idx][idx] = 0.0;
                        for (int i = 0; i < jac[idx][idx].size(); ++i)
                            jac[idx][idx][i][i] = 1.0;
                        res[idx] = 0;
                    }
                }
            }

            // assemble the coupling blocks for coupled models (does nothing if not coupled)
            maybeAssembleCouplingBlocks(residual);
        }
        else
            DUNE_THROW(Dune::NotImplemented, "Ghost elements not supported");


        auto applyDirichlet = [&] (const auto& scvI,
                                   const auto& dirichletValues,
                                   const auto eqIdx,
                                   const auto pvIdx)
        {
            res[scvI.dofIndex()][eqIdx] = this->curElemVolVars()[scvI].priVars()[pvIdx] - dirichletValues[pvIdx];

            auto& row = jac[scvI.dofIndex()];
            for (auto col = row.begin(); col != row.end(); ++col)
                row[col.index()][eqIdx] = 0.0;

            jac[scvI.dofIndex()][scvI.dofIndex()][eqIdx][pvIdx] = 1.0;

            // if a periodic dof has Dirichlet values also apply the same Dirichlet values to the other dof TODO periodic
            // if (this->assembler().gridGeometry().dofOnPeriodicBoundary(scvI.dofIndex()))
            // {
            //     const auto periodicDof = this->assembler().gridGeometry().periodicallyMappedDof(scvI.dofIndex());
            //     res[periodicDof][eqIdx] = this->curElemVolVars()[scvI].priVars()[pvIdx] - dirichletValues[pvIdx];
            //     const auto end = jac[periodicDof].end();
            //     for (auto it = jac[periodicDof].begin(); it != end; ++it)
            //         (*it) = periodicDof != it.index() ? 0.0 : 1.0;
            // }
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
    template<class ResidualVector>
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
            res[scvI.dofIndex()][eqIdx] = this->curElemVolVars()[scvI].priVars()[eqIdx] - dirichletValues[pvIdx];
        };

        this->asImp_().enforceDirichletConstraints(applyDirichlet);
    }

    /*!
     * \brief Enforce Dirichlet constraints
     */
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
            for (const auto& scvf : scvfs(this->fvGeometry()))
            {
                if (scvf.isFrontal() && scvf.boundary())
                {
                    const auto bcTypes = this->elemBcTypes()[scvf.localIndex()];
                    if (bcTypes.hasDirichlet())
                    {
                        const auto& scv = this->fvGeometry().scv(scvf.insideScvIdx());
                        const auto dirichletValues = this->asImp_().problem().dirichlet(this->element(), scvf);

                        // set the Dirichlet conditions in residual and jacobian
                        for (int eqIdx = 0; eqIdx < numEq; ++eqIdx)
                        {
                            for (int pvIdx = 0; pvIdx < GridView::dimension; ++pvIdx)
                            {
                                if (bcTypes.isDirichlet(pvIdx) && pvIdx == scv.dofAxis()) // TODO?
                                    applyDirichlet(scv, dirichletValues, eqIdx, pvIdx);
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
 * \ingroup FaceCenteredStaggeredDiscretization
 * \brief An assembler for Jacobian and residual contribution per element (Face-centered methods)
 * \tparam TypeTag The TypeTag
 * \tparam diffMethod The differentiation method to residual compute derivatives
 * \tparam implicit Specifies whether the time discretization is implicit or not (i.e. explicit)
 * \tparam Implementation The actual implementation, if void this class is the actual implementation
 */
template<class TypeTag, class Assembler, DiffMethod diffMethod = DiffMethod::numeric, bool implicit = true, class Implementation = void>
class FaceCenteredLocalAssembler;

/*!
 * \ingroup Assembly
 * \ingroup FaceCenteredStaggeredDiscretization
 * \brief Face-centered scheme local assembler using numeric differentiation and implicit time discretization
 */
template<class TypeTag, class Assembler, class Implementation>
class FaceCenteredLocalAssembler<TypeTag, Assembler, DiffMethod::numeric, /*implicit=*/true, Implementation>
: public FaceCenteredLocalAssemblerBase<
    TypeTag, Assembler,
    Detail::NonVoidOrDefault_t<Implementation, FaceCenteredLocalAssembler<TypeTag, Assembler, DiffMethod::numeric, true, Implementation>>,
    /*implicit=*/true
>
{
    using ThisType = FaceCenteredLocalAssembler<TypeTag, Assembler, DiffMethod::numeric, true, Implementation>;
    using ParentType = FaceCenteredLocalAssemblerBase<TypeTag, Assembler, Detail::NonVoidOrDefault_t<Implementation, ThisType>, true>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using Element = typename GetPropType<TypeTag, Properties::GridGeometry>::GridView::template Codim<0>::Entity;
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using FVElementGeometry = typename GridGeometry::LocalView;
    using SubControlVolume = typename FVElementGeometry::SubControlVolume;
    using SubControlVolumeFace = typename FVElementGeometry::SubControlVolumeFace;
    using GridVariables = GetPropType<TypeTag, Properties::GridVariables>;
    using JacobianMatrix = GetPropType<TypeTag, Properties::JacobianMatrix>;
    using VolumeVariables = GetPropType<TypeTag, Properties::VolumeVariables>;

    static constexpr auto numEq = GetPropType<TypeTag, Properties::ModelTraits>::numEq();
    static constexpr bool enableGridFluxVarsCache = GetPropType<TypeTag, Properties::GridVariables>::GridFluxVariablesCache::cachingEnabled;

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
        const auto& problem = this->asImp_().problem();
        const auto& element = this->element();
        const auto& fvGeometry = this->fvGeometry();
        const auto& curSol = this->asImp_().curSol();
        auto&& curElemVolVars = this->curElemVolVars();

        // get the vector of the actual element residuals
        const auto origResiduals = this->evalLocalResidual();

        //////////////////////////////////////////////////////////////////////////////////////////////////
        //                                                                                              //
        // Calculate derivatives of all dofs in stencil with respect to the dofs in the element. In the //
        // neighboring elements we do so by computing the derivatives of the fluxes which depend on the //
        // actual element. In the actual element we evaluate the derivative of the entire residual.     //
        //                                                                                              //
        //////////////////////////////////////////////////////////////////////////////////////////////////

        // one residual per element facet
        const auto numElementResiduals = fvGeometry.numScv();

        // create the vector storing the partial derivatives
        ElementResidualVector partialDerivs(numElementResiduals);

        const auto evalSource = [&](ElementResidualVector& residual, const SubControlVolume& scv)
        {
            this->localResidual().evalSource(residual, problem, element, fvGeometry, curElemVolVars, scv);
        };

        const auto evalStorage = [&](ElementResidualVector& residual, const SubControlVolume& scv)
        {
            this->localResidual().evalStorage(residual, problem, element, fvGeometry, this->prevElemVolVars(), curElemVolVars, scv);
        };

        const auto evalFlux = [&](ElementResidualVector& residual, const SubControlVolumeFace& scvf)
        {
            if (!scvf.processorBoundary())
                this->localResidual().evalFlux(residual, problem, element, fvGeometry, curElemVolVars, this->elemBcTypes(), this->elemFluxVarsCache(), scvf);
        };

        const auto evalDerivative = [&] (const auto& scvI, const auto& scvJ)
        {
            // derivative w.r.t. own DOF
            for (int pvIdx = 0; pvIdx < numEq; pvIdx++)
            {
                partialDerivs = 0.0;
                const auto& otherElement = fvGeometry.gridGeometry().element(scvJ.elementIndex());
                auto otherElemSol = elementSolution(otherElement, curSol, fvGeometry.gridGeometry()); // TODO allow selective creation of elemsol (for one scv)
                auto& curOtherVolVars = this->getVolVarAccess(gridVariables.curGridVolVars(), curElemVolVars, scvJ);
                const VolumeVariables origOtherVolVars(curOtherVolVars);

                auto evalResiduals = [&](Scalar priVar)
                {
                    // update the volume variables and compute element residual
                    otherElemSol[scvJ.localDofIndex()][pvIdx] = priVar;
                    curOtherVolVars.update(otherElemSol, problem, otherElement, scvJ);
                    this->asImp_().maybeUpdateCouplingContext(scvJ, otherElemSol, pvIdx);

                    ElementResidualVector residual(numElementResiduals);
                    residual = 0;

                    evalSource(residual, scvI);

                    if (!this->assembler().isStationaryProblem())
                        evalStorage(residual, scvI);

                    for (const auto& scvf : scvfs(fvGeometry, scvI))
                        evalFlux(residual, scvf);

                    return residual;
                };

                // derive the residuals numerically
                static const NumericEpsilon<Scalar, numEq> eps_{this->asImp_().problem().paramGroup()};
                static const int numDiffMethod = getParamFromGroup<int>(this->asImp_().problem().paramGroup(), "Assembly.NumericDifferenceMethod");
                NumericDifferentiation::partialDerivative(evalResiduals, otherElemSol[scvJ.localDofIndex()][pvIdx], partialDerivs, origResiduals,
                                                        eps_(otherElemSol[scvJ.localDofIndex()][pvIdx], pvIdx), numDiffMethod);

                const auto updateJacobian = [&]()
                {
                    for (int eqIdx = 0; eqIdx < numEq; eqIdx++)
                    {
                        // A[i][col][eqIdx][pvIdx] is the rate of change of
                        // the residual of equation 'eqIdx' at dof 'i'
                        // depending on the primary variable 'pvIdx' at dof
                        // 'col'.
                        A[scvI.dofIndex()][scvJ.dofIndex()][eqIdx][pvIdx] += partialDerivs[scvI.localDofIndex()][eqIdx];
                    }
                };

                using GeometryHelper = typename std::decay_t<decltype(fvGeometry.gridGeometry())>::GeometryHelper;
                using LocalIntersectionMapper = typename std::decay_t<decltype(fvGeometry.gridGeometry())>::LocalIntersectionMapper;
                LocalIntersectionMapper localIsMapper;

                const bool isParallel = fvGeometry.gridGeometry().gridView().comm().size() > 1;
                if (isParallel)
                    localIsMapper.update(fvGeometry.gridGeometry().gridView(), element);

                if (element.partitionType() == Dune::InteriorEntity)
                    updateJacobian();
                else
                {
                    const auto localIdxI = scvI.indexInElement();
                    const auto localIdxJ = scvJ.indexInElement();

                    const auto& facetI = GeometryHelper::facet(localIsMapper.refToRealIdx(localIdxI), element);
                    // add contribution of opposite scv lying within the overlap/ghost zone
                    if (facetI.partitionType() == Dune::BorderEntity &&
                        (localIdxJ == GeometryHelper::localOppositeIdx(localIdxI) || scvJ.dofIndex() == scvI.dofIndex()))
                        updateJacobian();
                }

                if (isParallel && element.partitionType() == Dune::InteriorEntity)
                {
                    const auto localIdxI = scvI.indexInElement();
                    const auto& facetI = GeometryHelper::facet(localIsMapper.refToRealIdx(localIdxI), element);
                    if (facetI.partitionType() == Dune::BorderEntity)
                    {
                        for (const auto& scvf : scvfs(fvGeometry, scvI))
                        {
                            if (scvf.isFrontal() || scvf.boundary())
                                continue;

                            // parallel scvs TODO drawing
                            if (scvf.outsideScvIdx() == scvJ.index())
                                updateJacobian();
                            else
                            {
                                // normal scvs
                                const auto& orthogonalScvf = fvGeometry.lateralOrthogonalScvf(scvf);
                                if (orthogonalScvf.boundary())
                                    continue;

                                if (orthogonalScvf.insideScvIdx() == scvJ.index() || orthogonalScvf.outsideScvIdx() == scvJ.index())
                                    updateJacobian();
                            }
                        }
                    }
                }

                // restore the original state of the scv's volume variables
                curOtherVolVars = origOtherVolVars;

                // restore the original element solution
                otherElemSol[scvJ.localDofIndex()][pvIdx] = curSol[scvJ.dofIndex()][pvIdx];
                this->asImp_().maybeUpdateCouplingContext(scvJ, otherElemSol, pvIdx);
                // TODO additional dof dependencies
            }
        };

        // calculation of the derivatives
        for (const auto& scvI : scvs(fvGeometry))
        {
            // derivative w.r.t. own DOFs
            evalDerivative(scvI, scvI);

            // derivative w.r.t. other DOFs
            const auto& otherScvIndices = fvGeometry.gridGeometry().connectivityMap()[scvI.index()];
            for (const auto globalJ : otherScvIndices)
                evalDerivative(scvI, fvGeometry.scv(globalJ));
        }

        // evaluate additional derivatives that might arise from the coupling
        this->asImp_().maybeEvalAdditionalDomainDerivatives(origResiduals, A, gridVariables);

        return origResiduals;
    }
};


/*!
 * \ingroup Assembly
 * \ingroup FaceCenteredStaggeredDiscretization
 * \brief TODO docme
 */
template<class TypeTag, class Assembler, class Implementation>
class FaceCenteredLocalAssembler<TypeTag, Assembler, DiffMethod::numeric, /*implicit=*/false, Implementation>
: public FaceCenteredLocalAssemblerBase<
    TypeTag, Assembler,
    Detail::NonVoidOrDefault_t<Implementation, FaceCenteredLocalAssembler<TypeTag, Assembler, DiffMethod::numeric, false, Implementation>>,
    /*implicit=*/false
>
{
    using ThisType = FaceCenteredLocalAssembler<TypeTag, Assembler, DiffMethod::numeric, false, Implementation>;
    using ParentType = FaceCenteredLocalAssemblerBase<TypeTag, Assembler, Detail::NonVoidOrDefault_t<Implementation, ThisType>, false>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using Element = typename GetPropType<TypeTag, Properties::GridGeometry>::GridView::template Codim<0>::Entity;
    using GridVariables = GetPropType<TypeTag, Properties::GridVariables>;
    using JacobianMatrix = GetPropType<TypeTag, Properties::JacobianMatrix>;
    using VolumeVariables = GetPropType<TypeTag, Properties::VolumeVariables>;

    enum { numEq = GetPropType<TypeTag, Properties::ModelTraits>::numEq() };

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
        const auto& problem = this->asImp_().problem();
        const auto& element = this->element();
        const auto& fvGeometry = this->fvGeometry();
        const auto& curSol = this->asImp_().curSol();
        auto&& curElemVolVars = this->curElemVolVars();

        // get the vector of the actual element residuals
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
                    // auto partialDerivsTmp = partialDerivs;
                    elemSol[scv.localDofIndex()][pvIdx] = priVar;
                    curVolVars.update(elemSol, problem, element, scv);
                    return this->evalLocalStorageResidual();
                };

                // derive the residuals numerically
                static const NumericEpsilon<Scalar, numEq> eps_{problem.paramGroup()};
                static const int numDiffMethod = getParamFromGroup<int>(problem.paramGroup(), "Assembly.NumericDifferenceMethod");
                NumericDifferentiation::partialDerivative(evalStorage, elemSol[scv.localDofIndex()][pvIdx], partialDerivs, origStorageResiduals,
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
        return origResiduals;
    }
};

/*!
 * \ingroup Assembly
 * \ingroup FaceCenteredStaggeredDiscretization
 * \brief TODO docme
 */
template<class TypeTag, class Assembler, class Implementation>
class FaceCenteredLocalAssembler<TypeTag, Assembler, DiffMethod::analytic, /*implicit=*/true, Implementation>
: public FaceCenteredLocalAssemblerBase<
    TypeTag, Assembler,
    FaceCenteredLocalAssembler<TypeTag, Assembler, DiffMethod::analytic, true, Implementation>,
    /*implicit=*/true
>
{
    using ThisType = FaceCenteredLocalAssembler<TypeTag, Assembler, DiffMethod::analytic, true, Implementation>;
    using ParentType = FaceCenteredLocalAssemblerBase<TypeTag, Assembler, Detail::NonVoidOrDefault_t<Implementation, ThisType>, true>;
    using JacobianMatrix = GetPropType<TypeTag, Properties::JacobianMatrix>;
    using GridVariables = GetPropType<TypeTag, Properties::GridVariables>;
    using VolumeVariables = GetPropType<TypeTag, Properties::VolumeVariables>;

    enum { numEq = GetPropType<TypeTag, Properties::ModelTraits>::numEq() };

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
                if (this->elemBcTypes()[insideScv.localDofIndex()].hasNeumann())
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
};

/*!
 * \ingroup Assembly
 * \ingroup FaceCenteredStaggeredDiscretization
 * \brief TODO docme
 */
template<class TypeTag, class Assembler, class Implementation>
class FaceCenteredLocalAssembler<TypeTag, Assembler, DiffMethod::analytic, /*implicit=*/false, Implementation>
: public FaceCenteredLocalAssemblerBase<
    TypeTag, Assembler,
    FaceCenteredLocalAssembler<TypeTag, Assembler, DiffMethod::analytic, false, Implementation>,
    /*implicit=*/false
>
{
    using ThisType = FaceCenteredLocalAssembler<TypeTag, Assembler, DiffMethod::analytic, false, Implementation>;
    using ParentType = FaceCenteredLocalAssemblerBase<TypeTag, Assembler, Detail::NonVoidOrDefault_t<Implementation, ThisType>, false>;
    using JacobianMatrix = GetPropType<TypeTag, Properties::JacobianMatrix>;
    using GridVariables = GetPropType<TypeTag, Properties::GridVariables>;
    using VolumeVariables = GetPropType<TypeTag, Properties::VolumeVariables>;

    enum { numEq = GetPropType<TypeTag, Properties::ModelTraits>::numEq() };

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

        // get the vector of the actual element residuals
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
};

} // end namespace Dumux

#endif
