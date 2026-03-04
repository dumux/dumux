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
 * \ingroup MultiDomain
 * \brief An assembler for Jacobian and residual contribution per element (CVFE methods) for multidomain problems
 */
#ifndef DUMUX_MULTIDOMAIN_SUBDOMAIN_CVFE_LOCAL_ASSEMBLER__HH
#define DUMUX_MULTIDOMAIN_SUBDOMAIN_CVFE_LOCAL_ASSEMBLER__HH

#include <type_traits>

#include <dune/common/reservedvector.hh>
#include <dune/common/indices.hh>
#include <dune/common/hybridutilities.hh>
#include <dune/grid/common/gridenums.hh> // for GhostEntity
#include <dune/istl/matrixindexset.hh>

#include <dumux/common/typetraits/localdofs_.hh>
#include <dumux/common/typetraits/boundary_.hh>
#include <dumux/common/reservedblockvector.hh>
#include <dumux/common/properties.hh>
#include <dumux/common/parameters.hh>
#include <dumux/common/numericdifferentiation.hh>
#include <dumux/assembly/numericepsilon.hh>
#include <dumux/assembly/diffmethod.hh>
#include <dumux/assembly/cvfelocalassembler_.hh>
#include <dumux/discretization/cvfe/localdof.hh>

namespace Dumux::Experimental {

/*!
 * \ingroup Assembly
 * \ingroup CVFEDiscretization
 * \ingroup MultiDomain
 * \brief A base class for all CVFE subdomain local assemblers
 * \tparam id the id of the sub domain
 * \tparam TypeTag the TypeTag
 * \tparam Assembler the assembler type
 * \tparam Implementation the actual implementation type
 * \tparam implicit Specifies whether the time discretization is implicit or not (i.e. explicit)
 */
template<std::size_t id, class TypeTag, class Assembler, class Implementation, DiffMethod dm, bool implicit>
class SubDomainCVFELocalAssemblerBase : public CVFELocalAssembler<TypeTag, Assembler, dm, implicit, Implementation>
{
    using ParentType = CVFELocalAssembler<TypeTag, Assembler, dm, implicit, Implementation>;

    using Problem = GetPropType<TypeTag, Properties::Problem>;
    using SolutionVector = typename Assembler::SolutionVector;

    using GridVariables = GetPropType<TypeTag, Properties::GridVariables>;
    using GridVariablesCache = typename GridVariables::GridVariablesCache;
    using ElementVariables = typename GridVariablesCache::LocalView;
    using Scalar = typename GridVariables::Scalar;

    using GridGeometry = typename GridVariables::GridGeometry;
    using FVElementGeometry = typename GridGeometry::LocalView;
    using GridView = typename GridGeometry::GridView;
    using Element = typename GridView::template Codim<0>::Entity;

    using CouplingManager = typename Assembler::CouplingManager;

    static constexpr auto numEq = GetPropType<TypeTag, Properties::ModelTraits>::numEq();

public:
    //! export the domain id of this sub-domain
    static constexpr auto domainId = typename Dune::index_constant<id>();
    //! pull up constructor of parent class
    using ParentType::ParentType;
    //! export element residual vector type
    using ElementResidualVector = typename ParentType::LocalResidual::ElementResidualVector;

    // the constructor
    explicit SubDomainCVFELocalAssemblerBase(
        const Assembler& assembler,
        const Element& element,
        const SolutionVector& curSol,
        CouplingManager& couplingManager
    )
    : ParentType(assembler,
                 element,
                 curSol,
                 localView(assembler.gridGeometry(domainId)),
                 localView(assembler.gridVariables(domainId).curGridVars()),
                 localView(assembler.gridVariables(domainId).prevGridVars()),
                 assembler.localResidual(domainId),
                 (element.partitionType() == Dune::GhostEntity))
    , couplingManager_(couplingManager)
    {}

    /*!
     * \brief Computes the derivatives with respect to the given element and adds them
     *        to the global matrix. The element residual is written into the right hand side.
     */
    template<class JacobianMatrixRow, class SubResidualVector, class GridVariablesTuple>
    void assembleJacobianAndResidual(JacobianMatrixRow& jacRow, SubResidualVector& res, GridVariablesTuple& gridVariables)
    {
        auto assembleCouplingBlocks = [&](const auto& residual)
        {
            // assemble the coupling blocks
            using namespace Dune::Hybrid;
            forEach(integralRange(Dune::Hybrid::size(jacRow)), [&](auto&& i)
            {
                if constexpr (std::decay_t<decltype(i)>{} != id)
                    this->assembleJacobianCoupling(i, jacRow, residual, gridVariables);
            });
        };

        // the coupled model does not support partial reassembly yet
        const DefaultPartialReassembler* noReassembler = nullptr;
        ParentType::assembleJacobianAndResidual(jacRow[domainId], res, *std::get<domainId>(gridVariables), noReassembler, assembleCouplingBlocks);
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
     * \brief Evaluates the local source term for an element and given element variables
     */
    ElementResidualVector evalLocalSourceResidual(const Element& element, const ElementVariables& elemVars) const
    {
        static_assert(!Dumux::Detail::LocalDofs::hasNonCVLocalDofsInterface<FVElementGeometry>(), "Separate source calculation not implemented for hybrid schemes.");

        // initialize the residual vector for all scvs in this element
        ElementResidualVector residual(Dumux::Detail::LocalDofs::numLocalDofs(this->fvGeometry()));

        // evaluate the source term
        // forward to the local residual specialized for the discretization methods
        for (const auto& scv : scvs(this->fvGeometry()))
        {
            residual[scv.localDofIndex()] = this->localResidual().sourceIntegral(this->fvGeometry(), elemVars, scv);
        }

        return residual;
    }

    /*!
     * \brief Evaluates the local source term depending on time discretization scheme
     */
    ElementResidualVector evalLocalSourceResidual(const Element& neighbor) const
    { return this->evalLocalSourceResidual(neighbor, implicit ? this->curElemVars() : this->prevElemVars()); }

    /*!
     * \brief Prepares all local views necessary for local assembly.
     */
    void bindLocalViews()
    {
        // get some references for convenience
        const auto& element = this->element();
        const auto& curSol = this->curSol(domainId);
        auto&& fvGeometry = this->fvGeometry();
        auto&& curElemVars = this->curElemVars();

        // bind the caches
        couplingManager_.bindCouplingContext(domainId, element, this->assembler());
        fvGeometry.bind(element);

        if constexpr (implicit)
        {
            curElemVars.bind(element, fvGeometry, curSol);
            if (!this->assembler().isStationaryProblem())
                this->prevElemVars().bindElement(element, fvGeometry, this->assembler().prevSol()[domainId]);
        }
        else
        {
            auto& prevElemVars = this->prevElemVars();
            const auto& prevSol = this->assembler().prevSol()[domainId];

            curElemVars.bindElement(element, fvGeometry, curSol);
            prevElemVars.bind(element, fvGeometry, prevSol);
        }

        this->elemBcTypes().update(problem(), this->element(), this->fvGeometry());
    }

    //! return reference to the underlying problem
    template<std::size_t i = domainId>
    const Problem& problem(Dune::index_constant<i> dId = domainId) const
    { return this->assembler().problem(dId); }

    //! return reference to the underlying problem
    template<std::size_t i = domainId>
    const auto& curSol(Dune::index_constant<i> dId = domainId) const
    { return ParentType::curSol()[dId]; }

    //! return reference to the coupling manager
    CouplingManager& couplingManager()
    { return couplingManager_; }

private:
    CouplingManager& couplingManager_; //!< the coupling manager
};

/*!
 * \ingroup Assembly
 * \ingroup CVFEDiscretization
 * \ingroup MultiDomain
 * \brief The CVFE scheme multidomain local assembler
 * \tparam id the id of the sub domain
 * \tparam TypeTag the TypeTag
 * \tparam Assembler the assembler type
 * \tparam DM the numeric differentiation method
 * \tparam implicit whether the assembler is explicit or implicit in time
 */
template<std::size_t id, class TypeTag, class Assembler, DiffMethod DM = DiffMethod::numeric, bool implicit = true>
class SubDomainCVFELocalAssembler;

/*!
 * \ingroup Assembly
 * \ingroup CVFEDiscretization
 * \ingroup MultiDomain
 * \brief CVFE scheme multi domain local assembler using numeric differentiation and implicit time discretization
 */
template<std::size_t id, class TypeTag, class Assembler>
class SubDomainCVFELocalAssembler<id, TypeTag, Assembler, DiffMethod::numeric, /*implicit=*/true>
: public SubDomainCVFELocalAssemblerBase<id, TypeTag, Assembler,
             SubDomainCVFELocalAssembler<id, TypeTag, Assembler, DiffMethod::numeric, true>, DiffMethod::numeric, true >
{
    using ThisType = SubDomainCVFELocalAssembler<id, TypeTag, Assembler, DiffMethod::numeric, /*implicit=*/true>;
    using ParentType = SubDomainCVFELocalAssemblerBase<id, TypeTag, Assembler, ThisType, DiffMethod::numeric, /*implicit=*/true>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;

    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using GridView = typename GridGeometry::GridView;
    using FVElementGeometry = typename GridGeometry::LocalView;
    using Element = typename GridView::template Codim<0>::Entity;
    using Problem = GetPropType<TypeTag, Properties::Problem>;

    static constexpr int numEq = GetPropType<TypeTag, Properties::ModelTraits>::numEq();
    static constexpr bool enableGridVarsCache = GetPropType<TypeTag, Properties::GridVariables>::GridVariablesCache::cachingEnabled;
    static constexpr auto domainI = Dune::index_constant<id>();

public:
    using ParentType::ParentType;
    //! export element residual vector type
    using ElementResidualVector = typename ParentType::LocalResidual::ElementResidualVector;

    /*!
     * \brief Update the coupling context for coupled models.
     */
    template<class ScvOrLocalDof, class ElemSol>
    void maybeUpdateCouplingContext(const ScvOrLocalDof& scvOrLocalDof, ElemSol& elemSol, const int pvIdx)
    {
        this->couplingManager().updateCouplingContext(domainI, *this, domainI, scvOrLocalDof.dofIndex(), elemSol[Dumux::Detail::LocalDofs::index(scvOrLocalDof)], pvIdx);
    }

    /*!
     * \brief Update the additional domain derivatives for coupled models.
     */
    template<class JacobianMatrixDiagBlock, class GridVariables>
    void maybeEvalAdditionalDomainDerivatives(const ElementResidualVector& origResiduals, JacobianMatrixDiagBlock& A, GridVariables& gridVariables)
    {
        this->couplingManager().evalAdditionalDomainDerivatives(domainI, *this, origResiduals, A, gridVariables);
    }

    /*!
     * \brief Computes the derivatives with respect to the given element and adds them
     *        to the global matrix.
     */
    template<std::size_t otherId, class JacobianBlock, class GridVariables>
    void assembleJacobianCoupling(Dune::index_constant<otherId> domainJ, JacobianBlock& A,
                                  const ElementResidualVector& res, GridVariables& gridVariables)
    {
        ////////////////////////////////////////////////////////////////////////////////////////////////////////
        // Calculate derivatives of all dofs in the element with respect to all dofs in the coupling stencil. //
        ////////////////////////////////////////////////////////////////////////////////////////////////////////

        // get some aliases for convenience
        const auto& element = this->element();
        const auto& fvGeometry = this->fvGeometry();
        auto&& curElemVars = this->curElemVars();

        // convenience lambda for call to update self
        auto updateCoupledVariables = [&] ()
        {
            // Update ourself after the context has been modified. Depending on the
            // type of caching, other objects might have to be updated. All ifs can be optimized away.
            if constexpr (enableGridVarsCache)
                this->couplingManager().updateCoupledVariables(domainI, *this, gridVariables.curGridVars());
            else
                this->couplingManager().updateCoupledVariables(domainI, *this, curElemVars);
        };

        // get element stencil information
        const auto& stencil = this->couplingManager().couplingStencil(domainI, element, domainJ);
        const auto& curSolJ = this->curSol(domainJ);
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
                    updateCoupledVariables();
                    return this->couplingManager().evalCouplingResidual(domainI, *this, domainJ, globalJ);
                };

                // derive the residuals numerically
                ElementResidualVector partialDerivs(Dumux::Detail::LocalDofs::numLocalDofs(fvGeometry));

                const auto& paramGroup = this->assembler().problem(domainJ).paramGroup();
                static const int numDiffMethod = getParamFromGroup<int>(paramGroup, "Assembly.NumericDifferenceMethod");
                static const auto epsCoupl = this->couplingManager().numericEpsilon(domainJ, paramGroup);

                NumericDifferentiation::partialDerivative(evalCouplingResidual, origPriVarsJ[pvIdx], partialDerivs, origResidual,
                                                          epsCoupl(origPriVarsJ[pvIdx], pvIdx), numDiffMethod);

                // update the global stiffness matrix with the current partial derivatives
                // Note: For the new interface Dirichlet constraints are incorporated later
                for (const auto& localDof : localDofs(fvGeometry))
                {
                    for (int eqIdx = 0; eqIdx < numEq; eqIdx++)
                    {
                        // A[i][col][eqIdx][pvIdx] is the rate of change of
                        // the residual of equation 'eqIdx' at dof 'i'
                        // depending on the primary variable 'pvIdx' at dof
                        // 'col'.
                        A[localDof.dofIndex()][globalJ][eqIdx][pvIdx] += partialDerivs[localDof.index()][eqIdx];
                    }
                }

                // restore the current element solution
                priVarsJ[pvIdx] = origPriVarsJ[pvIdx];

                // restore the undeflected state of the coupling context
                this->couplingManager().updateCouplingContext(domainI, *this, domainJ, globalJ, priVarsJ, pvIdx);
            }

            // Restore original state of the element variables.
            // This has to be done in case they depend on variables of domainJ before
            // we continue with the numeric derivative w.r.t the next globalJ. Otherwise,
            // the next "origResidual" will be incorrect.
            updateCoupledVariables();
        }
    }
};

/*!
 * \ingroup Assembly
 * \ingroup CVFEDiscretization
 * \ingroup MultiDomain
 * \brief CVFE scheme multi domain local assembler using numeric differentiation and explicit time discretization
 */
template<std::size_t id, class TypeTag, class Assembler>
class SubDomainCVFELocalAssembler<id, TypeTag, Assembler, DiffMethod::numeric, /*implicit=*/false>
: public SubDomainCVFELocalAssemblerBase<id, TypeTag, Assembler,
             SubDomainCVFELocalAssembler<id, TypeTag, Assembler, DiffMethod::numeric, false>, DiffMethod::numeric, false >
{
    using ThisType = SubDomainCVFELocalAssembler<id, TypeTag, Assembler, DiffMethod::numeric, /*implicit=*/false>;
    using ParentType = SubDomainCVFELocalAssemblerBase<id, TypeTag, Assembler, ThisType, DiffMethod::numeric, /*implicit=*/false>;

public:
    using ParentType::ParentType;
    //! export element residual vector type
    using ElementResidualVector = typename ParentType::LocalResidual::ElementResidualVector;

    /*!
     * \brief Computes the coupling derivatives with respect to the given element and adds them
     *        to the global matrix.
     * \note Since the coupling can only enter sources or fluxes and these are evaluated on
     *       the old time level (explicit scheme), the coupling blocks are empty.
     */
    template<std::size_t otherId, class JacobianBlock, class GridVariables>
    void assembleJacobianCoupling(Dune::index_constant<otherId> domainJ, JacobianBlock& A,
                                  const ElementResidualVector& res, GridVariables& gridVariables)
    {}
};

} // end namespace Dumux::Experimental

#endif
