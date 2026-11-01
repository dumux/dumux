// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup Experimental
 * \ingroup Assembly
 * \ingroup FaceCenteredStaggeredDiscretization
 * \brief A local assembler for the weighted operators of time integration methods
 *        (face-centered staggered methods)
 */
#ifndef DUMUX_EXPERIMENTAL_FC_LOCAL_ASSEMBLER_HH
#define DUMUX_EXPERIMENTAL_FC_LOCAL_ASSEMBLER_HH

#include <dune/common/exceptions.hh>
#include <dune/grid/common/gridenums.hh>

#include <dumux/common/typetraits/typetraits.hh>
#include <dumux/common/properties.hh>
#include <dumux/common/parameters.hh>
#include <dumux/common/numericdifferentiation.hh>
#include <dumux/common/numeqvector.hh>
#include <dumux/assembly/numericepsilon.hh>
#include <dumux/assembly/diffmethod.hh>
#include <dumux/discretization/facecentered/staggered/elementsolution.hh>

#include <dumux/experimental/assembly/fvlocalassemblerbase.hh>
#include <dumux/experimental/timestepping/implicitsolve.hh>

namespace Dumux::Experimental {

/*!
 * \ingroup Experimental
 * \ingroup Assembly
 * \ingroup FaceCenteredStaggeredDiscretization
 * \brief A base class for the local face-centered staggered assemblers of the implicit solves of
 *        time integration methods
 * \tparam TypeTag The TypeTag
 * \tparam Assembler The assembler type
 * \tparam Implementation The actual implementation
 */
template<class TypeTag, class Assembler, class Implementation>
class FaceCenteredLocalAssemblerBase : public Experimental::FVLocalAssemblerBase<TypeTag, Assembler, Implementation>
{
    using ParentType = Experimental::FVLocalAssemblerBase<TypeTag, Assembler, Implementation>;
    using JacobianMatrix = GetPropType<TypeTag, Properties::JacobianMatrix>;
    using GridVariables = GetPropType<TypeTag, Properties::GridVariables>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using Problem = GetPropType<TypeTag, Properties::Problem>;

    static constexpr int numEq = GetPropType<TypeTag, Properties::ModelTraits>::numEq();

public:
    //! The Dirichlet boundary conditions are enforced by constraint rows applied after the assembly
    static constexpr bool enforcesElementConstraints = true;

    using ParentType::ParentType;

    void bindLocalViews()
    {
        ParentType::bindLocalViews();
        this->elemBcTypes().update(this->asImp_().problem(), this->element(), this->fvGeometry());
    }

    /*!
     * \brief Computes the derivatives of the weighted residual with respect to the given element
     *        and adds them to the global matrix. The element residual is added to the right hand side.
     */
    template <class ResidualVector, class Scalar_>
    void assembleJacobianAndResidual(JacobianMatrix& jac, ResidualVector& res, GridVariables& gridVariables,
                                     const OperatorWeights<Scalar_>& weights,
                                     ResidualVector& constrainedDofs)
    {
        checkSupported_();

        this->localResidual().spatialWeight(weights.spatial, weights.spatial == 0.0);
        this->localResidual().temporalWeight(weights.storage, weights.storage == 0.0);
        this->localResidual().inertiaWeight(weights.inertia, weights.inertia == 0.0);

        this->asImp_().bindLocalViews();

        const auto residual = this->asImp_().assembleJacobianAndResidualImpl(jac, gridVariables);
        addToGlobal_(res, residual);

        // the equations of the facets in the overlap are left to their owners
        if (this->element().partitionType() != Dune::InteriorEntity)
        {
            for (const auto& scv : scvs(this->fvGeometry()))
            {
                if (this->element().template subEntity<1>(scv.indexInElement()).partitionType() == Dune::BorderEntity)
                    continue;

                const auto idx = scv.dofIndex();
                jac[idx][idx] = 0.0;
                for (int i = 0; i < numEq; ++i)
                    jac[idx][idx][i][i] = 1.0;
                res[idx] = 0.0;
            }
        }
    }

    /*!
     * \brief Set the unweighted operators of this element in the requested global vectors
     * \note The vectors have to be zero as the operators of the elements sharing a face add up
     */
    template<class SubResidualVector>
    void assembleOperators(OperatorEvaluations<SubResidualVector> evaluations)
    {
        checkSupported_();

        this->localResidual().spatialWeight(1.0, !evaluations.spatial);
        this->asImp_().bindLocalViews();

        if (evaluations.spatial)
            addToGlobal_(*evaluations.spatial, this->evalLocalFluxAndSourceResidual(this->curElemVolVars()));
        if (evaluations.storage)
            addToGlobal_(*evaluations.storage, this->localResidual().evalStorageOperator(this->fvGeometry(), this->curElemVolVars()));
        if (evaluations.inertia)
            addToGlobal_(*evaluations.inertia, this->localResidual().evalInertiaOperator(this->fvGeometry(), this->curElemVolVars()));
    }

    /*!
     * \brief Replace the equations of the degrees of freedom with a Dirichlet boundary
     *        condition by the constraint x - g = 0
     * \param curSol The current solution
     * \param res The residual, in which the constraint rows are set
     * \param jac The Jacobian, in which the constraint rows are set (nullptr if only the residual is assembled)
     * \param constrainedDofs Marks the constrained equations
     */
    template<class SolutionVector, class ResidualVector>
    void enforceElementConstraints(const SolutionVector& curSol,
                                   ResidualVector& res,
                                   JacobianMatrix* jac,
                                   ResidualVector& constrainedDofs)
    {
        this->fvGeometry().bind(this->element());
        this->elemBcTypes().update(this->asImp_().problem(), this->element(), this->fvGeometry());

        if (!this->elemBcTypes().hasDirichlet())
            return;

        for (const auto& scvf : scvfs(this->fvGeometry()))
        {
            if (!(scvf.isFrontal() && scvf.boundary()))
                continue;

            const auto bcTypes = this->elemBcTypes()[scvf.localIndex()];
            if (!bcTypes.hasDirichlet())
                continue;

            const auto& scv = this->fvGeometry().scv(scvf.insideScvIdx());
            const int componentIdx = scv.dofAxis();
            if (!bcTypes.isDirichlet(componentIdx))
                continue;

            static_assert(numEq == 1, "Not yet implemented for more than one vector-valued primary variable");
            const auto dirichletValues = this->asImp_().problem().dirichlet(this->element(), scvf);
            const auto dofIdx = scv.dofIndex();

            res[dofIdx][0] = curSol[dofIdx][0] - dirichletValues[componentIdx];
            constrainedDofs[dofIdx][0] = 1.0;

            if (jac)
            {
                auto& row = (*jac)[dofIdx];
                for (auto col = row.begin(); col != row.end(); ++col)
                    row[col.index()][0] = 0.0;
                (*jac)[dofIdx][dofIdx][0][0] = 1.0;
            }
        }
    }

    /*!
     * \brief The weighted residual of the equation of the facet scv at the current state of the local views
     * \note Only the entry of scv is nonzero. It equals the entry of scv in evalLocalResidual().
     */
    template<class SubControlVolume>
    auto evalScvResidual(const SubControlVolume& scv) const
    {
        return this->localResidual().evalScvResidual(this->fvGeometry(), this->curElemVolVars(),
                                                     this->elemFluxVarsCache(), this->elemBcTypes(), scv);
    }

    /*!
     * \brief Update the coupling context for coupled models.
     * \note This does nothing per default (not a coupled model).
     */
    template<class... Args>
    void maybeUpdateCouplingContext(Args&&...) {}

private:
    void checkSupported_() const
    {
        if constexpr (Problem::enableInternalDirichletConstraints())
            DUNE_THROW(Dune::NotImplemented, "Internal Dirichlet constraints in the weighted operator assembly of face-centered staggered schemes");

        if (this->elementIsGhost())
            DUNE_THROW(Dune::NotImplemented, "Ghost elements not supported");
    }

    //! Add the element contributions of the facets to the global vector (in parallel runs, overlap facets are left to their owners)
    template<class GlobalVector, class ElementVector>
    void addToGlobal_(GlobalVector& global, const ElementVector& local) const
    {
        for (const auto& scv : scvs(this->fvGeometry()))
        {
            if (this->element().partitionType() != Dune::InteriorEntity
                && this->element().template subEntity<1>(scv.indexInElement()).partitionType() != Dune::BorderEntity)
                continue;

            global[scv.dofIndex()] += local[scv.localDofIndex()];
        }
    }
};

/*!
 * \ingroup Experimental
 * \ingroup Assembly
 * \ingroup FaceCenteredStaggeredDiscretization
 * \brief Face-centered staggered local assembler for the implicit solves of time integration methods
 * \tparam TypeTag The TypeTag
 * \tparam Assembler The assembler type
 * \tparam diffMethod The differentiation method to compute derivatives of the residual
 * \tparam Implementation The actual implementation via CRTP
 */
template<class TypeTag, class Assembler, DiffMethod diffMethod = DiffMethod::numeric, class Implementation = void>
class FaceCenteredLocalAssembler;

/*!
 * \ingroup Experimental
 * \ingroup Assembly
 * \ingroup FaceCenteredStaggeredDiscretization
 * \brief Face-centered staggered local assembler using numeric differentiation
 * \note The fluxes over the faces of a facet enter only the equation of that facet, so each
 *       derivative is that of a single facet equation with respect to one degree of freedom in its stencil.
 */
template<class TypeTag, class Assembler, class Implementation>
class FaceCenteredLocalAssembler<TypeTag, Assembler, DiffMethod::numeric, Implementation>
: public FaceCenteredLocalAssemblerBase<TypeTag, Assembler,
                                        NonVoidOr<FaceCenteredLocalAssembler<TypeTag, Assembler, DiffMethod::numeric, Implementation>, Implementation>>
{
    using ThisType = FaceCenteredLocalAssembler<TypeTag, Assembler, DiffMethod::numeric, Implementation>;
    using ParentType = FaceCenteredLocalAssemblerBase<TypeTag, Assembler, NonVoidOr<ThisType, Implementation>>;
    using PrimaryVariable = typename GetPropType<TypeTag, Properties::PrimaryVariables>::value_type;
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using FVElementGeometry = typename GridGeometry::LocalView;
    using SubControlVolume = typename FVElementGeometry::SubControlVolume;
    using GridVariables = GetPropType<TypeTag, Properties::GridVariables>;
    using JacobianMatrix = GetPropType<TypeTag, Properties::JacobianMatrix>;
    using VolumeVariables = GetPropType<TypeTag, Properties::VolumeVariables>;

    static constexpr int numEq = GetPropType<TypeTag, Properties::ModelTraits>::numEq();

public:
    using ParentType::ParentType;

    using LocalResidual = typename ParentType::LocalResidual;
    using ElementResidualVector = typename LocalResidual::ElementResidualVector;

    /*!
     * \brief Computes the derivatives with respect to the given element and adds them
     *        to the global matrix.
     * \return The element residual at the current solution.
     */
    ElementResidualVector assembleJacobianAndResidualImpl(JacobianMatrix& A, GridVariables& gridVariables)
    {
        if (!this->isImplicit())
            DUNE_THROW(Dune::NotImplemented, "Explicit time integration with face-centered staggered schemes");

        const auto& problem = this->asImp_().problem();
        const auto& element = this->element();
        const auto& fvGeometry = this->fvGeometry();
        const auto& gridGeometry = fvGeometry.gridGeometry();
        const auto& curSol = this->asImp_().curSol();
        auto&& curElemVolVars = this->curElemVolVars();

        const auto origResiduals = this->evalLocalResidual();

        const auto numElementResiduals = fvGeometry.numScv();
        ElementResidualVector partialDerivs(numElementResiduals);

        const auto evalDerivatives = [&](const SubControlVolume& scvI, const SubControlVolume& scvJ)
        {
            const auto& otherElement = gridGeometry.element(scvJ.elementIndex());
            auto otherElemSol = elementSolution(otherElement, curSol, gridGeometry);
            auto& curOtherVolVars = this->getVolVarAccess(gridVariables.curGridVolVars(), curElemVolVars, scvJ);
            const VolumeVariables origOtherVolVars(curOtherVolVars);

            for (int pvIdx = 0; pvIdx < numEq; ++pvIdx)
            {
                partialDerivs = 0.0;

                auto evalResiduals = [&](PrimaryVariable priVar)
                {
                    otherElemSol[scvJ.localDofIndex()][pvIdx] = priVar;
                    curOtherVolVars.update(otherElemSol, problem, otherElement, scvJ);
                    this->asImp_().maybeUpdateCouplingContext(scvJ, otherElemSol, pvIdx);
                    return this->evalScvResidual(scvI);
                };

                static const NumericEpsilon<PrimaryVariable, numEq> eps_{problem.paramGroup()};
                static const int numDiffMethod = getParamFromGroup<int>(problem.paramGroup(), "Assembly.NumericDifferenceMethod");
                NumericDifferentiation::partialDerivative(evalResiduals, otherElemSol[scvJ.localDofIndex()][pvIdx], partialDerivs, origResiduals,
                                                          eps_(otherElemSol[scvJ.localDofIndex()][pvIdx], pvIdx), numDiffMethod);

                for (int eqIdx = 0; eqIdx < numEq; ++eqIdx)
                    A[scvI.dofIndex()][scvJ.dofIndex()][eqIdx][pvIdx] += partialDerivs[scvI.localDofIndex()][eqIdx];

                curOtherVolVars = origOtherVolVars;
                otherElemSol[scvJ.localDofIndex()][pvIdx] = curSol[scvJ.dofIndex()][pvIdx];
                this->asImp_().maybeUpdateCouplingContext(scvJ, otherElemSol, pvIdx);
            }
        };

        const auto& connectivityMap = gridGeometry.connectivityMap();
        for (const auto& scvI : scvs(fvGeometry))
        {
            // the equations of the facets in the overlap are left to their owners
            if (element.partitionType() != Dune::InteriorEntity
                && element.template subEntity<1>(scvI.indexInElement()).partitionType() != Dune::BorderEntity)
                continue;

            evalDerivatives(scvI, scvI);
            for (const auto globalJ : connectivityMap[scvI.index()])
                evalDerivatives(scvI, fvGeometry.scv(globalJ));
        }

        return origResiduals;
    }
};

} // end namespace Dumux::Experimental

#endif
