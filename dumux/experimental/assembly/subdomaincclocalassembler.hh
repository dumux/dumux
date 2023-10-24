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
 * \ingroup MultiDomain
 * \brief A multidomain local assembler for Jacobian and residual contribution per element (cell-centered methods)
 */
#ifndef DUMUX_EXPERIMENTAL_SUBDOMAIN_CC_LOCAL_ASSEMBLER_HH
#define DUMUX_EXPERIMENTAL_SUBDOMAIN_CC_LOCAL_ASSEMBLER_HH

#include <dune/common/indices.hh>
#include <dune/common/reservedvector.hh>
#include <dune/common/hybridutilities.hh>
#include <dune/grid/common/gridenums.hh> // for GhostEntity
#include <dune/istl/matrixindexset.hh>

#include <dumux/common/reservedblockvector.hh>
#include <dumux/common/properties.hh>
#include <dumux/common/parameters.hh>
#include <dumux/common/numericdifferentiation.hh>
#include <dumux/common/numeqvector.hh>
#include <dumux/discretization/elementsolution.hh>
#include <dumux/discretization/extrusion.hh>
#include <dumux/assembly/numericepsilon.hh>
#include <dumux/assembly/diffmethod.hh>

#include <dumux/experimental/assembly/cclocalassembler.hh>

namespace Dumux::Experimental {

/*!
 * \ingroup Experimental
 * \ingroup Assembly
 * \ingroup CCDiscretization
 * \ingroup MultiDomain
 * \brief A base class for all multidomain local assemblers
 * \tparam id the id of the sub domain
 * \tparam TypeTag the TypeTag
 * \tparam Assembler the assembler type
 * \tparam Implementation the actual assembler implementation
 */
template<std::size_t id, class TypeTag, class Assembler, class Implementation, DiffMethod dm>
class SubDomainCCLocalAssemblerBase : public Experimental::CCLocalAssembler<TypeTag, Assembler, dm, Implementation>
{
    using ParentType = Experimental::CCLocalAssembler<TypeTag, Assembler, dm, Implementation>;

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
    using FVElementGeometry = typename GridGeometry::LocalView;
    using SubControlVolume = typename GridGeometry::SubControlVolume;
    using SubControlVolumeFace = typename GridGeometry::SubControlVolumeFace;
    using Extrusion = Extrusion_t<GridGeometry>;
    using GridView = typename GridGeometry::GridView;
    using Element = typename GridView::template Codim<0>::Entity;

    using CouplingManager = typename Assembler::CouplingManager;
    using ElementResidualVector = typename ParentType::LocalResidual::ElementResidualVector;

public:
    //! export the domain id of this sub-domain
    static constexpr auto domainId = typename Dune::index_constant<id>();
    //! the local residual type of this domain
    using LocalResidual = GetPropType<TypeTag, Properties::LocalResidual>;
    //! pull up constructor of parent class
    using ParentType::ParentType;

    //! the constructor
    explicit SubDomainCCLocalAssemblerBase(const Assembler& assembler,
                                           const Element& element,
                                           const SolutionVector& curSol,
                                           CouplingManager& couplingManager)
    : ParentType(assembler,
                 element,
                 curSol,
                 localView(assembler.gridGeometry(domainId)),
                 localView(assembler.gridVariables(domainId).curGridVolVars()),
                 localView(assembler.gridVariables(domainId).gridFluxVarsCache()),
                 assembler.localResidual(domainId),
                 (element.partitionType() == Dune::GhostEntity),
                 assembler.isImplicit())
    , couplingManager_(couplingManager)
    {}

    /*!
     * \brief Computes the derivatives with respect to the given element and adds them
     *        to the global matrix. The element residual is written into the right hand side.
     */
    template<class JacobianMatrixRow, class SubResidualVector, class GridVariablesTuple, class StageParams>
    void assembleJacobianAndResidual(JacobianMatrixRow& jacRow, SubResidualVector& res, GridVariablesTuple& gridVariables,
                                     const StageParams& stageParams, SubResidualVector& temporal, SubResidualVector& spatial,
                                     SubResidualVector& constrainedDofs)
    {
        auto assembleCouplingBlocks = [&](const auto& residual)
        {
            // for the coupling blocks
            using namespace Dune::Hybrid;
            forEach(integralRange(Dune::Hybrid::size(jacRow)), [&](auto&& i)
            {
                if (std::decay_t<decltype(i)>{} != id)
                    this->assembleJacobianCoupling(i, jacRow, residual, gridVariables);
            });
        };

        // the coupled model does not support partial reassembly yet
        const DefaultPartialReassembler* noReassembler = nullptr;
        ParentType::assembleJacobianAndResidual(
            jacRow[domainId], res,
            *std::get<domainId>(gridVariables),
            stageParams, temporal, spatial,
            constrainedDofs,
            noReassembler,
            assembleCouplingBlocks
        );
    }

    /*!
     * \brief Assemble the entries in a coupling block of the jacobian.
     *        There is no coupling block between a domain and itself.
     */
    template<std::size_t otherId, class JacRow, class GridVariables>
    void assembleJacobianCoupling(Dune::index_constant<otherId> domainJ, JacRow& jacRow,
                                  const LocalResidualValues& res, GridVariables& gridVariables)
    {
        if constexpr (id != otherId)
            this->asImp_().assembleJacobianCoupling(domainJ, jacRow[domainJ], res, *std::get<domainId>(gridVariables));
    }

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
        curElemVolVars.bind(element, fvGeometry, curSol);
        elemFluxVarsCache.bind(element, fvGeometry, curElemVolVars);
    }

    //! return reference to the subdomain problem
    template<std::size_t i = domainId>
    const Problem& problem(Dune::index_constant<i> dId = domainId) const
    { return this->assembler().problem(domainId); }

    //! return reference to the subdomain solution
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
 * \ingroup Experimental
 * \ingroup Assembly
 * \ingroup CCDiscretization
 * \ingroup MultiDomain
 * \brief The cell-centered scheme multidomain local assembler
 * \tparam id the id of the sub domain
 * \tparam TypeTag the TypeTag
 * \tparam Assembler the assembler type
 * \tparam DM the numeric differentiation method
 */
template<std::size_t id, class TypeTag, class Assembler, DiffMethod DM = DiffMethod::numeric>
class SubDomainCCLocalAssembler;

/*!
 * \ingroup Experimental
 * \ingroup Assembly
 * \ingroup CCDiscretization
 * \ingroup MultiDomain
 * \brief Cell-centered scheme multidomain local assembler using numeric differentiation
 */
template<std::size_t id, class TypeTag, class Assembler>
class SubDomainCCLocalAssembler<id, TypeTag, Assembler, DiffMethod::numeric>
: public SubDomainCCLocalAssemblerBase<id, TypeTag, Assembler,
            SubDomainCCLocalAssembler<id, TypeTag, Assembler, DiffMethod::numeric>, DiffMethod::numeric>
{
    using ThisType = SubDomainCCLocalAssembler<id, TypeTag, Assembler, DiffMethod::numeric>;
    using ParentType = SubDomainCCLocalAssemblerBase<id, TypeTag, Assembler, ThisType, DiffMethod::numeric>;
    using Problem = GetPropType<TypeTag, Properties::Problem>;

    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using LocalResidualValues = Dumux::NumEqVector<GetPropType<TypeTag, Properties::PrimaryVariables>>;

    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using FVElementGeometry = typename GridGeometry::LocalView;
    using SubControlVolume = typename GridGeometry::SubControlVolume;
    using GridView = typename GridGeometry::GridView;
    using Element = typename GridView::template Codim<0>::Entity;

    enum { numEq = GetPropType<TypeTag, Properties::ModelTraits>::numEq() };
    enum { dim = GridView::dimension };

    static constexpr bool enableGridFluxVarsCache = GetPropType<TypeTag, Properties::GridVariables>::GridFluxVariablesCache::cachingEnabled;
    static constexpr bool enableGridVolVarsCache = GetPropType<TypeTag, Properties::GridVariables>::GridVolumeVariables::cachingEnabled;
    static constexpr int maxElementStencilSize = GridGeometry::maxElementStencilSize;
    static constexpr auto domainI = Dune::index_constant<id>();

public:
    using ParentType::ParentType;
    using ElementResidualVector = typename ParentType::LocalResidual::ElementResidualVector;

    /*!
     * \brief Update the coupling context for coupled models.
     */
    template<class ElemSol>
    void maybeUpdateCouplingContext(const SubControlVolume& scv, ElemSol& elemSol, const int pvIdx)
    {
        if (this->assembler().isImplicit())
            this->couplingManager().updateCouplingContext(domainI, *this, domainI, scv.dofIndex(), elemSol[0], pvIdx);
    }

    /*!
     * \brief Update the additional domain derivatives for coupled models.
     */
    template<class JacobianMatrixDiagBlock, class GridVariables>
    void maybeEvalAdditionalDomainDerivatives(const ElementResidualVector& origResiduals, JacobianMatrixDiagBlock& A, GridVariables& gridVariables)
    {
        if (this->assembler().isImplicit())
            this->couplingManager().evalAdditionalDomainDerivatives(domainI, *this, origResiduals, A, gridVariables);
    }

    /*!
     * \brief Computes the derivatives with respect to the given element and adds them
     *        to the global matrix.
     */
    template<std::size_t otherId, class JacobianBlock, class GridVariables>
    void assembleJacobianCoupling(Dune::index_constant<otherId> domainJ, JacobianBlock& A,
                                  const LocalResidualValues& res, GridVariables& gridVariables)
    {
        if (!this->assembler().isImplicit())
            return;

        ////////////////////////////////////////////////////////////////////////////////////////////////////////
        // Calculate derivatives of all dofs in the element with respect to all dofs in the coupling stencil. //
        ////////////////////////////////////////////////////////////////////////////////////////////////////////

        // get some aliases for convenience
        const auto& element = this->element();
        const auto& fvGeometry = this->fvGeometry();
        const auto& gridGeometry = fvGeometry.gridGeometry();
        auto&& curElemVolVars = this->curElemVolVars();
        auto&& elemFluxVarsCache = this->elemFluxVarsCache();

        // get stencil information
        const auto globalI = gridGeometry.elementMapper().index(element);
        const auto& stencil = this->couplingManager().couplingStencil(domainI, element, domainJ);
        const auto& curSolJ = this->curSol(domainJ);

        // make sure the flux variables cache does not contain any artifacts from prior differentiation
        elemFluxVarsCache.update(element, fvGeometry, curElemVolVars);

        // convenience lambda for call to update self
        auto updateCoupledVariables = [&] ()
        {
            // Update ourself after the context has been modified. Depending on the
            // type of caching, other objects might have to be updated. All ifs can be optimized away.
            if (enableGridFluxVarsCache)
            {
                if (enableGridVolVarsCache)
                {
                    this->couplingManager().updateCoupledVariables(domainI, *this, gridVariables.curGridVolVars(), gridVariables.gridFluxVarsCache());
                    this->couplingManager().updateCoupledVariables(domainI, *this, gridVariables.curGridVolVars(), elemFluxVarsCache);
                }
                else
                {
                    this->couplingManager().updateCoupledVariables(domainI, *this, curElemVolVars, gridVariables.gridFluxVarsCache());
                    this->couplingManager().updateCoupledVariables(domainI, *this, curElemVolVars, elemFluxVarsCache);
                }
            }
            else
            {
                if (enableGridVolVarsCache)
                    this->couplingManager().updateCoupledVariables(domainI, *this, gridVariables.curGridVolVars(), elemFluxVarsCache);
                else
                    this->couplingManager().updateCoupledVariables(domainI, *this, curElemVolVars, elemFluxVarsCache);
            }
        };

        for (const auto& globalJ : stencil)
        {
            // undeflected privars and privars to be deflected
            const auto origPriVarsJ = curSolJ[globalJ];
            auto priVarsJ = origPriVarsJ;

            const auto origResidual = this->couplingManager().evalCouplingResidual(domainI, *this, domainJ, globalJ)[0];

            for (int pvIdx = 0; pvIdx < JacobianBlock::block_type::cols; ++pvIdx)
            {
                auto evalCouplingResidual = [&](Scalar priVar)
                {
                    // update the volume variables and the flux var cache
                    priVarsJ[pvIdx] = priVar;
                    this->couplingManager().updateCouplingContext(domainI, *this, domainJ, globalJ, priVarsJ, pvIdx);
                    updateCoupledVariables();
                    return this->couplingManager().evalCouplingResidual(domainI, *this, domainJ, globalJ)[0];
                };

                // derive the residuals numerically
                LocalResidualValues partialDeriv(0.0);
                const auto& paramGroup = this->assembler().problem(domainJ).paramGroup();
                static const int numDiffMethod = getParamFromGroup<int>(paramGroup, "Assembly.NumericDifferenceMethod");
                static const auto epsCoupl = this->couplingManager().numericEpsilon(domainJ, paramGroup);
                NumericDifferentiation::partialDerivative(evalCouplingResidual, priVarsJ[pvIdx], partialDeriv, origResidual,
                                                          epsCoupl(priVarsJ[pvIdx], pvIdx), numDiffMethod);

                // add the current partial derivatives to the global jacobian matrix
                if constexpr (Problem::enableInternalDirichletConstraints())
                {
                    const auto& scv = fvGeometry.scv(globalI);
                    const auto internalDirichletConstraints = this->problem().hasInternalDirichletConstraint(this->element(), scv);
                    if (internalDirichletConstraints.none())
                    {
                        for (int eqIdx = 0; eqIdx < numEq; eqIdx++)
                            A[globalI][globalJ][eqIdx][pvIdx] += partialDeriv[eqIdx];
                    }
                    else
                    {
                        for (int eqIdx = 0; eqIdx < numEq; ++eqIdx)
                        {
                            if (internalDirichletConstraints[eqIdx])
                                A[globalI][globalJ][eqIdx][pvIdx] = 0.0;
                            else
                                A[globalI][globalJ][eqIdx][pvIdx] += partialDeriv[eqIdx];
                        }
                    }
                }
                else
                {
                    for (int eqIdx = 0; eqIdx < numEq; eqIdx++)
                        A[globalI][globalJ][eqIdx][pvIdx] += partialDeriv[eqIdx];
                }

                // restore the current element solution
                priVarsJ[pvIdx] = origPriVarsJ[pvIdx];

                // restore the undeflected state of the coupling context
                this->couplingManager().updateCouplingContext(domainI, *this, domainJ, globalJ, priVarsJ, pvIdx);
            }

            // restore original state of the flux vars cache and/or vol vars in case of global caching.
            // This has to be done in order to guarantee that the undeflected residual computation done
            // above is correct when jumping to the next coupled dof of domainJ.
            updateCoupledVariables();
        }
    }
};

} // end namespace Dumux

#endif
