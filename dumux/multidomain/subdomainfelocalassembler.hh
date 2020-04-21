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
 * \ingroup FEMDiscretization
 * \ingroup MultiDomain
 * \brief An assembler for Jacobian and residual contribution per element
 *        (finite element methods) for multidomain problems
 */
#ifndef DUMUX_MULTIDOMAIN_FE_LOCAL_ASSEMBLER_HH
#define DUMUX_MULTIDOMAIN_FE_LOCAL_ASSEMBLER_HH

#include <dune/common/indices.hh>
#include <dune/common/hybridutilities.hh>
#include <dune/geometry/referenceelements.hh>

#include <dumux/common/properties.hh>
#include <dumux/common/parameters.hh>
#include <dumux/common/numericdifferentiation.hh>

#include <dumux/assembly/numericepsilon.hh>
#include <dumux/assembly/diffmethod.hh>
#include <dumux/assembly/felocalassembler.hh>

namespace Dumux {

/*!
 * \ingroup Assembly
 * \ingroup FEMDiscretization
 * \ingroup MultiDomain
 * \brief A base class for all fe local assemblers
 * \tparam id the id of the sub domain
 * \tparam TypeTag the TypeTag
 * \tparam Assembler the assembler type
 * \tparam Implementation the actual implementation type
 * \tparam implicit Specifies whether the time discretization is implicit or not not (i.e. explicit)
 */
template<std::size_t id,
         class TypeTag,
         class Assembler,
         class Implementation,
         bool useImplicitAssembly>
class SubDomainFELocalAssemblerBase
: public FELocalAssemblerBase<TypeTag, Assembler, Implementation, useImplicitAssembly>
{
    using ParentType = FELocalAssemblerBase<TypeTag, Assembler, Implementation, useImplicitAssembly>;

    using Problem = GetPropType<TypeTag, Properties::Problem>;
    using NumEqVector = GetPropType<TypeTag, Properties::NumEqVector>;

    using PrimaryVariables = GetPropType<TypeTag, Properties::PrimaryVariables>;
    using GridVariables = GetPropType<TypeTag, Properties::GridVariables>;
    using Scalar = typename GridVariables::Scalar;

    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using FEElementGeometry = typename GridGeometry::LocalView;
    using ElementSolution = FEElementSolution<FEElementGeometry, PrimaryVariables>;

    using GridView = typename GridGeometry::GridView;
    using Element = typename GridView::template Codim<0>::Entity;

    using SolutionVector = typename Assembler::SolutionVector;
    using JacobianMatrix = typename Assembler::JacobianMatrix;
    using SubSolutionVector = GetPropType<TypeTag, Properties::SolutionVector>;
    using ElementBoundaryTypes = GetPropType<TypeTag, Properties::ElementBoundaryTypes>;
    using ReferenceElements = typename Dune::ReferenceElements<Scalar, GridView::dimension>;

    using CouplingManager = typename Assembler::CouplingManager;
    static constexpr auto numEq = GetPropType<TypeTag, Properties::ModelTraits>::numEq();

public:
    //! pull up public export of base class
    using typename ParentType::ElementResidualVector;

    //! Export the domain id of this sub-domain
    static constexpr auto domainId = typename Dune::index_constant<id>();

    //! The constructor
    explicit SubDomainFELocalAssemblerBase(const Assembler& assembler,
                                           const Element& element,
                                           const SolutionVector& curSol,
                                           CouplingManager& couplingManager)
    : ParentType(assembler,
                 element,
                 curSol,
                 localView(assembler.gridGeometry(domainId)),
                 ElementSolution(),
                 ElementSolution(),
                 assembler.localResidual(domainId),
                 (element.partitionType() == Dune::GhostEntity))
    , couplingManager_(couplingManager)
    {}

    /*!
     * \brief Prepares all local views necessary for local assembly.
     */
    void bindLocalViews()
    {
        // get some references for convenience
        const auto& element = this->element();
        const auto& curSol = this->curSol()[domainId];
        const auto& prevSol = this->assembler().prevSol()[domainId];
        const auto& gridGeometry = this->assembler().gridGeometry(domainId);

        // first bind the coupling context
        couplingManager_.bindCouplingContext(domainId, element, this->assembler());

        // then bind the local views
        this->feGeometry().bind(element);
        this->elemBcTypes().update(this->problem(), element, this->feGeometry());

        if (ParentType::isImplicit())
        {
            this->curElemSol() = elementSolution(element, curSol, gridGeometry);
            if (!this->assembler().isStationaryProblem())
                this->prevElemSol() = elementSolution(element, prevSol, gridGeometry);
        }
        else
        {
            this->curElemSol() = elementSolution(element, curSol, gridGeometry);
            this->prevElemSol() = elementSolution(element, prevSol, gridGeometry);
        }
    }

    /*!
     * \brief Computes the derivatives with respect to the given element and adds them
     *        to the global matrix. The element residual is written into the right hand side.
     */
    template<class JacobianMatrixRow, class GridVariablesTuple>
    void assembleJacobianAndResidual(JacobianMatrixRow& jacRow, SubSolutionVector& res, GridVariablesTuple& gridVariables)
    {
        this->asImp_().bindLocalViews();

        // for the diagonal jacobian block
        // forward to the internal implementation
        const auto residual = this->asImp_().assembleJacobianAndResidualImpl(jacRow[domainId], *std::get<domainId>(gridVariables));

        // update the residual vector
        const auto& tsLocalView = this->trialSpaceBasisLocalView();
        for (unsigned int i = 0; i < tsLocalView.size(); ++i)
            res[tsLocalView.index(i)] += residual[i];

        // assemble the coupling blocks
        using namespace Dune::Hybrid;
        forEach(integralRange(Dune::Hybrid::size(jacRow)), [&](auto&& i)
        {
            if (i != id)
                this->assembleJacobianCoupling(i, jacRow, residual, gridVariables);
        });

        // lambda for the incorporation of Dirichlet Bcs
        auto applyDirichlet = [&] (const auto& dirichletValues,
                                   const auto localDofIdx,
                                   const auto dofIdx,
                                   const auto eqIdx,
                                   const auto pvIdx)
        {
            res[dofIdx][eqIdx] = this->curElemSol()[localDofIdx][pvIdx] - dirichletValues[pvIdx];

            auto& row = jacRow[domainId][dofIdx];
            for (auto col = row.begin(); col != row.end(); ++col)
                row[col.index()][eqIdx] = 0.0;

            jacRow[domainId][dofIdx][dofIdx][eqIdx][pvIdx] = 1.0;
        };

        // incorporate Dirichlet BCs
        this->asImp_().enforceDirichletConstraints(applyDirichlet);
    }

    /*!
     * \brief Assemble the entries in a coupling block of the jacobian.
     *        There is no coupling block between a domain and itself.
     */
    template<std::size_t otherId, class JacRow, class GridVariables,
             typename std::enable_if_t<(otherId == id), int> = 0>
    void assembleJacobianCoupling(Dune::index_constant<otherId> domainJ, JacRow& jacRow,
                                  const ElementResidualVector& res, GridVariables& gridVariables)
    {}

    /*!
     * \brief Assemble the entries in a coupling block of the jacobian.
     */
    template<std::size_t otherId, class JacRow, class GridVariables,
             typename std::enable_if_t<(otherId != id), int> = 0>
    void assembleJacobianCoupling(Dune::index_constant<otherId> domainJ, JacRow& jacRow,
                                  const ElementResidualVector& res, GridVariables& gridVariables)
    {
        this->asImp_().assembleJacobianCoupling(domainJ, jacRow[domainJ], res, *std::get<domainId>(gridVariables));
    }

    /*!
     * \brief Assemble the residual vector entries only
     */
    void assembleResidual(SubSolutionVector& res)
    {
        this->asImp_().bindLocalViews();

        const auto residual = this->evalLocalResidual();
        const auto& tsLocalView = this->trialSpaceBasisLocalView();
        for (unsigned int i = 0; i < tsLocalView.size(); ++i)
            res[tsLocalView.index(i)] += residual[i];

        auto applyDirichlet = [&] (const auto& dirichletValues,
                                   const auto localDofIdx,
                                   const auto dofIdx,
                                   const auto eqIdx,
                                   const auto pvIdx)
        {
            res[dofIdx][eqIdx] = this->curElemSol()[localDofIdx][pvIdx] - dirichletValues[pvIdx];
        };

        this->asImp_().enforceDirichletConstraints(applyDirichlet);
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
            const auto& tsLocalView = this->trialSpaceBasisLocalView();
            const auto& element = this->element();
            const auto& fe = tsLocalView.tree().finiteElement();

            for (unsigned int localDofIdx = 0; localDofIdx < tsLocalView.size(); localDofIdx++)
            {
                const auto& bcTypes = this->elemBcTypes()[localDofIdx];

                if (!bcTypes.hasDirichlet())
                    continue;

                const auto dofIdx = tsLocalView.index(localDofIdx);
                const auto& localKey = fe.localCoefficients().localKey(localDofIdx);
                const auto subEntity = localKey.subEntity();
                const auto codim = localKey.codim();

                // values of dirichlet BCs
                PrimaryVariables dirichletValues;
                if (codim == 1) dirichletValues = this->template getDirichletValues_<1>(problem(), element, subEntity);
                else if (codim == 2) dirichletValues = this->template getDirichletValues_<2>(problem(), element, subEntity);
                else if (codim == 3) dirichletValues = this->template getDirichletValues_<3>(problem(), element, subEntity);

                for (int eqIdx = 0; eqIdx < numEq; ++eqIdx)
                {
                    if (bcTypes.isDirichlet(eqIdx))
                    {
                        const auto pvIdx = bcTypes.eqToDirichletIndex(eqIdx);
                        assert(0 <= pvIdx && pvIdx < numEq);
                        applyDirichlet(dirichletValues, localDofIdx, dofIdx, eqIdx, pvIdx);
                    }
                }
            }
        }
    }

    //! return reference to the underlying problem
    const Problem& problem() const
    { return this->assembler().problem(domainId); }

    //! return reference to the coupling manager
    CouplingManager& couplingManager()
    { return couplingManager_; }

private:
    CouplingManager& couplingManager_; //!< the coupling manager
};

/*!
 * \ingroup Assembly
 * \ingroup FEMDiscretization
 * \ingroup MultiDomain
 * \brief The finite element schemes' multidomain local assembler
 * \tparam id the id of the sub domain
 * \tparam TypeTag the TypeTag
 * \tparam Assembler the assembler type
 * \tparam DM the numeric differentiation method
 * \tparam implicit whether the assembler is explicit or implicit in time
 */
template<std::size_t id, class TypeTag, class Assembler, DiffMethod DM = DiffMethod::numeric, bool implicit = true>
class SubDomainFELocalAssembler;

/*!
 * \ingroup Assembly
 * \ingroup FEMDiscretization
 * \ingroup MultiDomain
 * \brief Finite element schemes' multidomain local assembler using numeric differentiation and implicit time discretization
 */
template<std::size_t id, class TypeTag, class Assembler>
class SubDomainFELocalAssembler<id, TypeTag, Assembler, DiffMethod::numeric, /*implicit=*/true>
: public SubDomainFELocalAssemblerBase<id, TypeTag, Assembler,
             SubDomainFELocalAssembler<id, TypeTag, Assembler, DiffMethod::numeric, true>, true >
{
    using ThisType = SubDomainFELocalAssembler<id, TypeTag, Assembler, DiffMethod::numeric, /*implicit=*/true>;
    using ParentType = SubDomainFELocalAssemblerBase<id, TypeTag, Assembler, ThisType, /*implicit=*/true>;

    using GridVariables = GetPropType<TypeTag, Properties::GridVariables>;
    using Scalar = typename GridVariables::Scalar;

    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using GridView = typename GridGeometry::GridView;

    static constexpr int numEq = GetPropType<TypeTag, Properties::ModelTraits>::numEq();
    static constexpr int dim = GridView::dimension;
    static constexpr auto domainI = Dune::index_constant<id>();

public:
    //! pull up public export of base class
    using typename ParentType::ElementResidualVector;

    //! Pull up base class constructor
    using ParentType::ParentType;

    /*!
     * \brief Computes the derivatives with respect to the given element and adds them
     *        to the global matrix.
     *
     * \return The element residual at the current solution.
     */
    template<class JacobianMatrixDiagBlock, class GridVariables>
    ElementResidualVector assembleJacobianAndResidualImpl(JacobianMatrixDiagBlock& A, GridVariables& gridVariables)
    {
        // get some aliases for convenience
        const auto& tsLocalView = this->trialSpaceBasisLocalView();
        const auto& curSol = this->curSol()[domainI];

        // get the vector of the actual element residuals
        const auto origResiduals = this->evalLocalResidual();

        //////////////////////////////////////////////////////////////////////////////////////////////
        // Calculate derivatives of the residual of all dofs in element with respect to themselves. //
        //////////////////////////////////////////////////////////////////////////////////////////////

        // create the element solution
        auto& elemSol = this->curElemSol();

        // create the vector storing the partial derivatives
        ElementResidualVector partialDerivs(tsLocalView.size());

        // calculation of the derivatives
        for (unsigned int localI = 0; localI < tsLocalView.size(); ++localI)
        {
            // dof index and corresponding actual pri vars
            const auto globalI = tsLocalView.index(localI);

            // calculate derivatives w.r.t to the privars at the dof at hand
            for (int pvIdx = 0; pvIdx < numEq; pvIdx++)
            {
                partialDerivs = 0.0;

                auto evalResiduals = [&](Scalar priVar)
                {
                    // update the volume variables and compute element residual
                    elemSol[localI][pvIdx] = priVar;
                    this->couplingManager().updateCouplingContext(domainI, *this, domainI, globalI, elemSol[localI], pvIdx);
                    return this->evalLocalResidual();
                };

                // derive the residuals numerically
                static const NumericEpsilon<Scalar, numEq> eps_{this->problem().paramGroup()};
                static const int numDiffMethod = getParamFromGroup<int>(this->problem().paramGroup(), "Assembly.NumericDifferenceMethod");
                NumericDifferentiation::partialDerivative(evalResiduals,
                                                          elemSol[localI][pvIdx],
                                                          partialDerivs,
                                                          origResiduals,
                                                          eps_(elemSol[localI][pvIdx], pvIdx),
                                                          numDiffMethod);

                // update the global stiffness matrix with the current partial derivatives
                for (unsigned int localJ = 0; localJ < tsLocalView.size(); ++localJ)
                {
                    const auto globalJ = tsLocalView.index(localJ);

                    // don't add derivatives for green entities
                    for (int eqIdx = 0; eqIdx < numEq; eqIdx++)
                    {
                        // A[i][col][eqIdx][pvIdx] is the rate of change of
                        // the residual of equation 'eqIdx' at dof 'i'
                        // depending on the primary variable 'pvIdx' at dof
                        // 'col'.
                        A[globalJ][globalI][eqIdx][pvIdx] += partialDerivs[localJ][eqIdx];
                    }
                }

                // restore the original element solution
                elemSol[localI][pvIdx] = curSol[globalI][pvIdx];
                this->couplingManager().updateCouplingContext(domainI, *this, domainI, globalI, elemSol[localI], pvIdx);
            }
        }

        // evaluate additional derivatives that might arise from the coupling
        this->couplingManager().evalAdditionalDomainDerivatives(domainI, *this, origResiduals, A, gridVariables);

        return origResiduals;
    }

    /*!
     * \brief Computes the derivatives with respect to the given element and adds them
     *        to the global matrix.
     *
     * \return The element residual at the current solution.
     */
    template<std::size_t otherId, class JacobianBlock, class GridVariables>
    void assembleJacobianCoupling(Dune::index_constant<otherId> domainJ, JacobianBlock& A,
                                  const ElementResidualVector& res, GridVariables& gridVariables)
    {
        // get some references for convenience
        const auto& element = this->element();
        const auto& stencil = this->couplingManager().couplingStencil(domainI, element, domainJ);

        const auto& curSolJ = this->curSol()[domainJ];
        for (const auto globalJ : stencil)
        {
            // undeflected privars and privars to be deflected
            const auto origPriVarsJ = curSolJ[globalJ];
            auto priVarsJ = origPriVarsJ;

            // the undeflected coupling residual
            const auto origResidual = this->couplingManager().evalCouplingResidual(domainI, *this, domainJ, globalJ);

            for (int pvIdx = 0; pvIdx < JacobianBlock::block_type::cols; ++pvIdx)
            {
                auto evalCouplingResidual = [&](Scalar priVar)
                {
                    priVarsJ[pvIdx] = priVar;
                    this->couplingManager().updateCouplingContext(domainI, *this, domainJ, globalJ, priVarsJ, pvIdx);
                    return this->couplingManager().evalCouplingResidual(domainI, *this, domainJ, globalJ);
                };

                // derive the residuals numerically
                ElementResidualVector partialDerivs(element.subEntities(dim));

                const auto& paramGroup = this->assembler().problem(domainJ).paramGroup();
                static const int numDiffMethod = getParamFromGroup<int>(paramGroup, "Assembly.NumericDifferenceMethod");
                static const auto epsCoupl = this->couplingManager().numericEpsilon(domainJ, paramGroup);

                NumericDifferentiation::partialDerivative(evalCouplingResidual, origPriVarsJ[pvIdx], partialDerivs, origResidual,
                                                          epsCoupl(origPriVarsJ[pvIdx], pvIdx), numDiffMethod);

                // update the global stiffness matrix with the current partial derivatives
                const auto tsLocalView = this->trialSpaceBasisLocalView();
                for (unsigned int localI = 0; localI < tsLocalView.size(); ++localI)
                {
                    const auto globalI = tsLocalView.index(localI);

                    // don't add derivatives for green entities
                    for (int eqIdx = 0; eqIdx < numEq; eqIdx++)
                    {
                        // A[i][col][eqIdx][pvIdx] is the rate of change of
                        // the residual of equation 'eqIdx' at dof 'i'
                        // depending on the primary variable 'pvIdx' at dof
                        // 'col'.

                        // If the dof is coupled by a Dirichlet condition,
                        // set the derived value only once (i.e. overwrite existing values).
                        // For other dofs, add the contribution of the partial derivative.
                        const auto bcTypes = this->elemBcTypes()[localI];
                        if (bcTypes.isCouplingDirichlet(eqIdx))
                            A[globalI][globalJ][eqIdx][pvIdx] = partialDerivs[localI][eqIdx];
                        else if (bcTypes.isDirichlet(eqIdx))
                            A[globalI][globalJ][eqIdx][pvIdx] = 0.0;
                        else
                            A[globalI][globalJ][eqIdx][pvIdx] += partialDerivs[localI][eqIdx];
                    }
                }

                // restore the current element solution
                priVarsJ[pvIdx] = origPriVarsJ[pvIdx];

                // restore the undeflected state of the coupling context
                this->couplingManager().updateCouplingContext(domainI, *this, domainJ, globalJ, priVarsJ, pvIdx);
            }
        }
    }
};

} // end namespace Dumux

#endif
