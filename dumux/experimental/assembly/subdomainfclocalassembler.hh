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
 * \ingroup MultiDomain
 * \brief A multidomain local assembler for the weighted operators of time integration methods
 *        (face-centered staggered methods)
 */
#ifndef DUMUX_EXPERIMENTAL_SUBDOMAIN_FC_LOCAL_ASSEMBLER_HH
#define DUMUX_EXPERIMENTAL_SUBDOMAIN_FC_LOCAL_ASSEMBLER_HH

#include <dune/common/exceptions.hh>
#include <dune/common/indices.hh>
#include <dune/common/hybridutilities.hh>
#include <dune/grid/common/gridenums.hh>

#include <dumux/common/properties.hh>
#include <dumux/common/parameters.hh>
#include <dumux/common/numericdifferentiation.hh>
#include <dumux/assembly/diffmethod.hh>

#include <dumux/experimental/assembly/fclocalassembler.hh>
#include <dumux/experimental/timestepping/implicitsolve.hh>

namespace Dumux::Experimental {

/*!
 * \ingroup Experimental
 * \ingroup Assembly
 * \ingroup FaceCenteredStaggeredDiscretization
 * \ingroup MultiDomain
 * \brief A base class for the face-centered staggered multidomain local assemblers of time integration methods
 * \tparam id the id of the sub domain
 * \tparam TypeTag the TypeTag
 * \tparam Assembler the assembler type
 * \tparam Implementation the actual assembler implementation
 * \tparam dm the differentiation method
 */
template<std::size_t id, class TypeTag, class Assembler, class Implementation, DiffMethod dm>
class SubDomainFaceCenteredLocalAssemblerBase
: public Experimental::FaceCenteredLocalAssembler<TypeTag, Assembler, dm, Implementation>
{
    using ParentType = Experimental::FaceCenteredLocalAssembler<TypeTag, Assembler, dm, Implementation>;

    using Problem = GetPropType<TypeTag, Properties::Problem>;
    using SolutionVector = typename Assembler::SolutionVector;
    using GridView = typename GetPropType<TypeTag, Properties::GridGeometry>::GridView;
    using Element = typename GridView::template Codim<0>::Entity;
    using CouplingManager = typename Assembler::CouplingManager;

public:
    //! export the domain id of this sub-domain
    static constexpr auto domainId = typename Dune::index_constant<id>();
    //! pull up constructor of parent class
    using ParentType::ParentType;

    explicit SubDomainFaceCenteredLocalAssemblerBase(const Assembler& assembler,
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
     * \brief Computes the derivatives of the weighted residual with respect to the given element
     *        and adds them to the global matrix row. The element residual is added to the right hand side.
     */
    template<class JacobianMatrixRow, class SubResidualVector, class GridVariablesTuple, class Scalar>
    void assembleJacobianAndResidual(JacobianMatrixRow& jacRow, SubResidualVector& res, GridVariablesTuple& gridVariables,
                                     const OperatorWeights<Scalar>& weights,
                                     SubResidualVector& constrainedDofs)
    {
        ParentType::assembleJacobianAndResidual(jacRow[domainId], res, *std::get<domainId>(gridVariables), weights, constrainedDofs);

        using namespace Dune::Hybrid;
        forEach(integralRange(Dune::Hybrid::size(jacRow)), [&](const auto domainJ)
        {
            if constexpr (domainJ != id)
                this->asImp_().assembleJacobianCoupling(domainJ, jacRow[domainJ], *std::get<domainId>(gridVariables));
        });
    }

    /*!
     * \brief Replace the equations of the degrees of freedom with a Dirichlet boundary condition
     *        by the constraint x - g = 0, also in the coupling blocks of the matrix row
     * \note Dirichlet conditions that are set by a coupling are not supported
     */
    template<class SubSolutionVector, class SubResidualVector, class JacobianMatrixRow>
    void enforceElementConstraints(const SubSolutionVector& curSol,
                                   SubResidualVector& res,
                                   JacobianMatrixRow& jacRow,
                                   SubResidualVector& constrainedDofs)
    {
        ParentType::enforceElementConstraints(curSol, res, &jacRow[domainId], constrainedDofs);

        if (!this->elemBcTypes().hasDirichlet())
            return;

        for (const auto& scvf : scvfs(this->fvGeometry()))
        {
            if (!(scvf.isFrontal() && scvf.boundary()))
                continue;

            const auto bcTypes = this->elemBcTypes()[scvf.localIndex()];
            const auto& scv = this->fvGeometry().scv(scvf.insideScvIdx());
            const int componentIdx = scv.dofAxis();
            if (!bcTypes.isDirichlet(componentIdx))
                continue;

            if (bcTypes.isCouplingDirichlet(componentIdx))
                DUNE_THROW(Dune::NotImplemented, "Coupling Dirichlet conditions with face-centered staggered schemes in the weighted operator assembly");

            using namespace Dune::Hybrid;
            forEach(integralRange(Dune::Hybrid::size(jacRow)), [&](const auto domainJ)
            {
                if constexpr (domainJ != id)
                {
                    auto& row = jacRow[domainJ][scv.dofIndex()];
                    for (auto col = row.begin(); col != row.end(); ++col)
                        row[col.index()][componentIdx] = 0.0;
                }
            });
        }
    }

    /*!
     * \brief Prepares all local views necessary for local assembly.
     */
    void bindLocalViews()
    {
        const auto& element = this->element();
        couplingManager_.bindCouplingContext(domainId, element, this->assembler());
        this->fvGeometry().bind(element);
        this->curElemVolVars().bind(element, this->fvGeometry(), this->curSol(domainId));
        this->elemFluxVarsCache().bind(element, this->fvGeometry(), this->curElemVolVars());
        this->elemBcTypes().update(problem(), element, this->fvGeometry());
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
    CouplingManager& couplingManager_;
};

/*!
 * \ingroup Experimental
 * \ingroup Assembly
 * \ingroup FaceCenteredStaggeredDiscretization
 * \ingroup MultiDomain
 * \brief The face-centered staggered multidomain local assembler of time integration methods
 * \tparam id the id of the sub domain
 * \tparam TypeTag the TypeTag
 * \tparam Assembler the assembler type
 * \tparam DM the numeric differentiation method
 */
template<std::size_t id, class TypeTag, class Assembler, DiffMethod DM = DiffMethod::numeric>
class SubDomainFaceCenteredLocalAssembler;

/*!
 * \ingroup Experimental
 * \ingroup Assembly
 * \ingroup FaceCenteredStaggeredDiscretization
 * \ingroup MultiDomain
 * \brief The face-centered staggered multidomain local assembler using numeric differentiation
 */
template<std::size_t id, class TypeTag, class Assembler>
class SubDomainFaceCenteredLocalAssembler<id, TypeTag, Assembler, DiffMethod::numeric>
: public SubDomainFaceCenteredLocalAssemblerBase<id, TypeTag, Assembler,
            SubDomainFaceCenteredLocalAssembler<id, TypeTag, Assembler, DiffMethod::numeric>, DiffMethod::numeric>
{
    using ThisType = SubDomainFaceCenteredLocalAssembler<id, TypeTag, Assembler, DiffMethod::numeric>;
    using ParentType = SubDomainFaceCenteredLocalAssemblerBase<id, TypeTag, Assembler, ThisType, DiffMethod::numeric>;
    using Problem = GetPropType<TypeTag, Properties::Problem>;
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using SubControlVolume = typename GridGeometry::LocalView::SubControlVolume;

    static constexpr int numEq = GetPropType<TypeTag, Properties::ModelTraits>::numEq();
    static constexpr bool enableGridFluxVarsCache = GetPropType<TypeTag, Properties::GridVariables>::GridFluxVariablesCache::cachingEnabled;
    static constexpr bool enableGridVolVarsCache = GetPropType<TypeTag, Properties::GridVariables>::GridVolumeVariables::cachingEnabled;
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
        this->couplingManager().updateCouplingContext(domainI, *this, domainI, scv.dofIndex(), elemSol[scv.localDofIndex()], pvIdx);
    }

    /*!
     * \brief Computes the derivatives of the weighted residual of the facets in the given element
     *        with respect to the primary variables of domainJ and adds them to the global matrix.
     */
    template<std::size_t otherId, class JacobianBlock, class GridVariables>
    void assembleJacobianCoupling(Dune::index_constant<otherId> domainJ, JacobianBlock& A, GridVariables& gridVariables)
    {
        static_assert(otherId != id);

        const auto& element = this->element();
        const auto& fvGeometry = this->fvGeometry();
        auto&& curElemVolVars = this->curElemVolVars();
        auto&& elemFluxVarsCache = this->elemFluxVarsCache();
        const auto& curSolJ = this->curSol(domainJ);

        auto updateCoupledVariables = [&] ()
        {
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

        const auto origResidual = this->evalLocalResidual();
        ElementResidualVector partialDerivs(fvGeometry.numScv());

        const bool interiorElement = element.partitionType() == Dune::InteriorEntity;
        for (const auto& scv : scvs(fvGeometry))
        {
            // of an overlap element only the border facets have assembled equations
            if (!interiorElement && element.template subEntity<1>(scv.indexInElement()).partitionType() != Dune::BorderEntity)
                continue;

            for (const auto globalJ : this->couplingManager().couplingStencil(domainI, element, scv, domainJ))
            {
                const auto origPriVarsJ = curSolJ[globalJ];
                auto priVarsJ = origPriVarsJ;

                for (int pvIdx = 0; pvIdx < JacobianBlock::block_type::cols; ++pvIdx)
                {
                    auto evalResidual = [&](auto priVar)
                    {
                        priVarsJ[pvIdx] = priVar;
                        this->couplingManager().updateCouplingContext(domainI, *this, domainJ, globalJ, priVarsJ, pvIdx);
                        updateCoupledVariables();
                        return this->evalScvResidual(scv);
                    };

                    const auto& paramGroup = this->assembler().problem(domainJ).paramGroup();
                    static const int numDiffMethod = getParamFromGroup<int>(paramGroup, "Assembly.NumericDifferenceMethod");
                    static const auto epsCoupl = this->couplingManager().numericEpsilon(domainJ, paramGroup);

                    partialDerivs = 0.0;
                    NumericDifferentiation::partialDerivative(evalResidual, origPriVarsJ[pvIdx], partialDerivs, origResidual,
                                                              epsCoupl(origPriVarsJ[pvIdx], pvIdx), numDiffMethod);

                    for (int eqIdx = 0; eqIdx < numEq; ++eqIdx)
                        A[scv.dofIndex()][globalJ][eqIdx][pvIdx] += partialDerivs[scv.localDofIndex()][eqIdx];

                    priVarsJ[pvIdx] = origPriVarsJ[pvIdx];
                    this->couplingManager().updateCouplingContext(domainI, *this, domainJ, globalJ, priVarsJ, pvIdx);
                }

                // the undeflected state is needed for the next coupled dof
                updateCoupledVariables();
            }
        }
    }
};

} // end namespace Dumux::Experimental

#endif
