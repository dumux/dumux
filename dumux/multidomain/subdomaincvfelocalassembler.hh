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
#ifndef DUMUX_MULTIDOMAIN_SUBDOMAIN_CVFE_LOCAL_ASSEMBLER_HH
#define DUMUX_MULTIDOMAIN_SUBDOMAIN_CVFE_LOCAL_ASSEMBLER_HH

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
#include <dumux/common/numeqvector.hh>
#include <dumux/assembly/numericepsilon.hh>
#include <dumux/assembly/diffmethod.hh>
#include <dumux/assembly/cvfelocalassembler.hh>
#include <dumux/discretization/extrusion.hh>

namespace Dumux {

/*!
 * \ingroup Assembly
 * \ingroup CVFEDiscretization
 * \ingroup MultiDomain
 * \brief A base class for all CVFE subdomain local assemblers
 * \tparam id the id of the sub domain
 * \tparam TypeTag the TypeTag
 * \tparam Assembler the assembler type
 * \tparam Implementation the actual implementation type
 * \tparam implicit Specifies whether the time discretization is implicit or not not (i.e. explicit)
 */
template<std::size_t id, class TypeTag, class Assembler, class Implementation, DiffMethod dm, bool implicit>
class SubDomainCVFELocalAssemblerBase : public CVFELocalAssembler<TypeTag, Assembler, dm, implicit, Implementation>
{
    using ParentType = CVFELocalAssembler<TypeTag, Assembler, dm, implicit, Implementation>;

    using Problem = GetPropType<TypeTag, Properties::Problem>;
    using LocalResidualValues = Dumux::NumEqVector<GetPropType<TypeTag, Properties::PrimaryVariables>>;
    using JacobianMatrix = GetPropType<TypeTag, Properties::JacobianMatrix>;
    using SolutionVector = typename Assembler::SolutionVector;
    using ElementBoundaryTypes = GetPropType<TypeTag, Properties::ElementBoundaryTypes>;

    using GridVariables = GetPropType<TypeTag, Properties::GridVariables>;
    using GridVolumeVariables = typename GridVariables::GridVolumeVariables;
    using ElementVolumeVariables = typename GridVolumeVariables::LocalView;
    using ElementFluxVariablesCache = typename GridVariables::GridFluxVariablesCache::LocalView;
    using Scalar = typename GridVariables::Scalar;

    using GridGeometry = typename GridVariables::GridGeometry;
    using Extrusion = Extrusion_t<GridGeometry>;
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
                 localView(assembler.gridVariables(domainId).curGridVolVars()),
                 localView(assembler.gridVariables(domainId).prevGridVolVars()),
                 localView(assembler.gridVariables(domainId).gridFluxVarsCache()),
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
     * \brief Evaluates the local source term for an element and given element volume variables
     */
    ElementResidualVector evalLocalSourceResidual(const Element& element, const ElementVolumeVariables& elemVolVars) const
    {
        // initialize the residual vector for all scvs in this element
        ElementResidualVector residual(Detail::LocalDofs::numLocalDofs(this->fvGeometry()));

        // evaluate the volume terms (storage + source terms)
        // forward to the local residual specialized for the discretization methods
        for (const auto& scv : scvs(this->fvGeometry()))
        {
            const auto& curVolVars = elemVolVars[scv];
            auto source = this->localResidual().computeSource(problem(), element, this->fvGeometry(), elemVolVars, scv);
            source *= -Extrusion::volume(this->fvGeometry(), scv)*curVolVars.extrusionFactor();
            residual[scv.localDofIndex()] = std::move(source);
        }

        return residual;
    }

    /*!
     * \brief Evaluates the local source term depending on time discretization scheme
     */
    ElementResidualVector evalLocalSourceResidual(const Element& neighbor) const
    { return this->evalLocalSourceResidual(neighbor, implicit ? this->curElemVolVars() : this->prevElemVolVars()); }

    /*!
     * \brief Prepares all local views necessary for local assembly.
     */
    void bindLocalViews()
    {
        // get some references for convenience
        const auto& element = this->element();
        const auto& curSol = this->curSol(domainId);
        auto&& fvGeometry = this->fvGeometry();
        auto&& curElemVolVars = this->curElemVolVars();
        auto&& elemFluxVarsCache = this->elemFluxVarsCache();

        // bind the caches
        couplingManager_.bindCouplingContext(domainId, element, this->assembler());
        fvGeometry.bind(element);

        if constexpr (implicit)
        {
            curElemVolVars.bind(element, fvGeometry, curSol);
            elemFluxVarsCache.bind(element, fvGeometry, curElemVolVars);
            if (!this->assembler().isStationaryProblem())
                this->prevElemVolVars().bindElement(element, fvGeometry, this->assembler().prevSol()[domainId]);
        }
        else
        {
            auto& prevElemVolVars = this->prevElemVolVars();
            const auto& prevSol = this->assembler().prevSol()[domainId];

            curElemVolVars.bindElement(element, fvGeometry, curSol);
            prevElemVolVars.bind(element, fvGeometry, prevSol);
            elemFluxVarsCache.bind(element, fvGeometry, prevElemVolVars);
        }

        this->elemBcTypes().update(problem(), this->element(), this->fvGeometry());
    }

    //! return reference to the underlying problem
    template<std::size_t i = domainId>
    const Problem& problem(Dune::index_constant<i> dId = domainId) const
    { return this->assembler().problem(domainId); }

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
    static constexpr int dim = GridView::dimension;

    static constexpr bool enableGridFluxVarsCache = GetPropType<TypeTag, Properties::GridVariables>::GridFluxVariablesCache::cachingEnabled;
    static constexpr bool enableGridVolVarsCache = GetPropType<TypeTag, Properties::GridVariables>::GridVolumeVariables::cachingEnabled;
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
        this->couplingManager().updateCouplingContext(domainI, *this, domainI, scvOrLocalDof.dofIndex(), elemSol[Detail::LocalDofs::index(scvOrLocalDof)], pvIdx);
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
        auto&& curElemVolVars = this->curElemVolVars();
        auto&& elemFluxVarsCache = this->elemFluxVarsCache();

        // convenience lambda for call to update self
        auto updateCoupledVariables = [&] ()
        {
            // Update ourself after the context has been modified. Depending on the
            // type of caching, other objects might have to be updated. All ifs can be optimized away.
            if constexpr (enableGridFluxVarsCache)
            {
                if constexpr (enableGridVolVarsCache)
                    this->couplingManager().updateCoupledVariables(domainI, *this, gridVariables.curGridVolVars(), gridVariables.gridFluxVarsCache());
                else
                    this->couplingManager().updateCoupledVariables(domainI, *this, curElemVolVars, gridVariables.gridFluxVarsCache());
            }
            else
            {
                if constexpr (enableGridVolVarsCache)
                    this->couplingManager().updateCoupledVariables(domainI, *this, gridVariables.curGridVolVars(), elemFluxVarsCache);
                else
                    this->couplingManager().updateCoupledVariables(domainI, *this, curElemVolVars, elemFluxVarsCache);
            }
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
                ElementResidualVector partialDerivs(Detail::LocalDofs::numLocalDofs(fvGeometry));

                const auto& paramGroup = this->assembler().problem(domainJ).paramGroup();
                static const int numDiffMethod = getParamFromGroup<int>(paramGroup, "Assembly.NumericDifferenceMethod");
                static const auto epsCoupl = this->couplingManager().numericEpsilon(domainJ, paramGroup);

                NumericDifferentiation::partialDerivative(evalCouplingResidual, origPriVarsJ[pvIdx], partialDerivs, origResidual,
                                                          epsCoupl(origPriVarsJ[pvIdx], pvIdx), numDiffMethod);

                // update the global stiffness matrix with the current partial derivatives
                for (const auto& scv : scvs(fvGeometry))
                {
                    for (int eqIdx = 0; eqIdx < numEq; eqIdx++)
                    {
                        // A[i][col][eqIdx][pvIdx] is the rate of change of
                        // the residual of equation 'eqIdx' at dof 'i'
                        // depending on the primary variable 'pvIdx' at dof
                        // 'col'.
                        A[scv.dofIndex()][globalJ][eqIdx][pvIdx] += partialDerivs[scv.localDofIndex()][eqIdx];

                        // For the old boundary interface we set them for each scv
                        // for the new interface we enforce them globally by constraints, i.e. do nothing here
                        if constexpr (!Detail::hasProblemBoundaryTypesForIntersectionFunction<Problem, FVElementGeometry, typename GridGeometry::GridView::Intersection>())
                        {
                            // If the dof is coupled by a Dirichlet condition,
                            // set the derived value only once (i.e. overwrite existing values).
                            if (this->elemBcTypes().hasDirichlet())
                            {
                                const auto bcTypes = this->elemBcTypes().get(fvGeometry, scv);
                                if (bcTypes.isCouplingDirichlet(eqIdx))
                                    A[scv.dofIndex()][globalJ][eqIdx][pvIdx] = partialDerivs[scv.localDofIndex()][eqIdx];
                                else if (bcTypes.isDirichlet(eqIdx))
                                    A[scv.dofIndex()][globalJ][eqIdx][pvIdx] = 0.0;
                            }
                        }

                        // enforce internal Dirichlet constraints
                        if constexpr (Problem::enableInternalDirichletConstraints())
                        {
                            const auto internalDirichletConstraints = this->problem().hasInternalDirichletConstraint(this->element(), scv);
                            if (internalDirichletConstraints[eqIdx])
                                A[scv.dofIndex()][globalJ][eqIdx][pvIdx] = 0.0;
                        }

                    }
                }

                // restore the current element solution
                priVarsJ[pvIdx] = origPriVarsJ[pvIdx];

                // restore the undeflected state of the coupling context
                this->couplingManager().updateCouplingContext(domainI, *this, domainJ, globalJ, priVarsJ, pvIdx);
            }

            // Restore original state of the flux vars cache and/or vol vars.
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

} // end namespace Dumux

#endif
