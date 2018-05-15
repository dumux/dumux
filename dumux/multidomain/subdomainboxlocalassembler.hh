// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   See the file COPYING for full copying permissions.                      *
 *                                                                           *
 *   This program is free software: you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation, either version 2 of the License, or       *
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
 * \ingroup BoxDiscretization
 * \brief An assembler for Jacobian and residual contribution per element (box methods)
 * \tparam TypeTag the TypeTag
 * \tparam DM the differentiation method to residual compute derivatives
 * \tparam implicit if to use an implicit or explicit time discretization
 */
#ifndef DUMUX_MULTIDOMAIN_BOX_LOCAL_ASSEMBLER_HH
#define DUMUX_MULTIDOMAIN_BOX_LOCAL_ASSEMBLER_HH

#include <dune/common/reservedvector.hh>
#include <dune/grid/common/gridenums.hh> // for GhostEntity
#include <dune/istl/matrixindexset.hh>

#include <dumux/common/reservedblockvector.hh>
#include <dumux/common/properties.hh>
#include <dumux/common/parameters.hh>
#include <dumux/common/numericdifferentiation.hh>
#include <dumux/assembly/diffmethod.hh>
#include <dumux/assembly/fvlocalassemblerbase.hh>

namespace Dumux {

/*!
 * \ingroup Assembly
 * \ingroup CCDiscretization
 * \brief A base class for all local assemblers
 * \tparam id the id of the sub domain
 * \tparam TypeTag the TypeTag
 * \tparam Assembler the assembler type
 */
template<std::size_t id, class TypeTag, class Assembler, class Implementation>
class SubDomainBoxLocalAssemblerBase : public FVLocalAssemblerBase<TypeTag, Assembler,Implementation, true>
{
    using ParentType = FVLocalAssemblerBase<TypeTag, Assembler,Implementation, true>;

    using Problem = typename GET_PROP_TYPE(TypeTag, Problem);
    using LocalResidualValues = typename GET_PROP_TYPE(TypeTag, NumEqVector);
    using LocalResidual = typename GET_PROP_TYPE(TypeTag, LocalResidual);
    using ElementResidualVector = typename LocalResidual::ElementResidualVector;
    using JacobianMatrix = typename GET_PROP_TYPE(TypeTag, JacobianMatrix);
    using SolutionVector = typename Assembler::SolutionVector;
    using SubSolutionVector = typename GET_PROP_TYPE(TypeTag, SolutionVector);
    using ElementBoundaryTypes = typename GET_PROP_TYPE(TypeTag, ElementBoundaryTypes);

    using GridVariables = typename GET_PROP_TYPE(TypeTag, GridVariables);
    using GridVolumeVariables = typename GridVariables::GridVolumeVariables;
    using ElementVolumeVariables = typename GridVolumeVariables::LocalView;
    using ElementFluxVariablesCache = typename GridVariables::GridFluxVariablesCache::LocalView;
    using Scalar = typename GridVariables::Scalar;

    using FVGridGeometry = typename GridVariables::GridGeometry;
    using FVElementGeometry = typename FVGridGeometry::LocalView;
    using SubControlVolume = typename FVGridGeometry::SubControlVolume;
    using SubControlVolumeFace = typename FVGridGeometry::SubControlVolumeFace;
    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
    using Element = typename GridView::template Codim<0>::Entity;

    using CouplingManager = typename Assembler::CouplingManager;

    static constexpr auto numEq = GET_PROP_TYPE(TypeTag, ModelTraits)::numEq();

public:
    //! export the domain id of this sub-domain
    static constexpr auto domainId = typename Dune::index_constant<id>();
    //! pull up constructor of parent class
    using ParentType::ParentType;

    // the constructor
    explicit SubDomainBoxLocalAssemblerBase(const Assembler& assembler,
                                            const Element& element,
                                            const SolutionVector& curSol,
                                            CouplingManager& couplingManager)
    : ParentType(assembler,
                 element,
                 curSol,
                 localView(assembler.fvGridGeometry(domainId)),
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
    template<class JacobianMatrixRow, class GridVariablesTuple>
    void assembleJacobianAndResidual(JacobianMatrixRow& jacRow, SubSolutionVector& res, GridVariablesTuple& gridVariables)
    {
        this->asImp_().bindLocalViews();
        this->elemBcTypes().update(problem(), this->element(), this->fvGeometry());

        // for the diagonal jacobian block
        // forward to the internal implementation
        const auto residual = this->asImp_().assembleJacobianAndResidualImpl(jacRow[domainId], *std::get<domainId>(gridVariables));

        for (const auto& scv : scvs(this->fvGeometry()))
            res[scv.dofIndex()] += residual[scv.indexInElement()];

        auto applyDirichlet = [&] (const auto& scvI,
                                   const auto& dirichletValues,
                                   const auto eqIdx,
                                   const auto pvIdx)
        {
            res[scvI.dofIndex()][eqIdx] = this->curElemVolVars()[scvI].priVars()[pvIdx] - dirichletValues[pvIdx];
            for (const auto& scvJ : scvs(this->fvGeometry()))
            {
                jacRow[domainId][scvI.dofIndex()][scvJ.dofIndex()][eqIdx] = 0.0;
                if (scvI.indexInElement() == scvJ.indexInElement())
                    jacRow[domainId][scvI.dofIndex()][scvI.dofIndex()][eqIdx][pvIdx] = 1.0;
            }
        };

        // for the coupling blocks
        using namespace Dune::Hybrid;
        forEach(integralRange(Dune::Hybrid::size(jacRow)), [&, domainId = domainId](auto&& i)
        {
            if (i != domainId)
                this->assembleJacobianCoupling(i, jacRow, residual, gridVariables);
        });

        this->asImp_().evalDirichletBoundaries(applyDirichlet);
    }

    template<std::size_t otherId, class JacRow, class GridVariables,
             typename std::enable_if_t<(otherId == id), int> = 0>
    void assembleJacobianCoupling(Dune::index_constant<otherId> domainJ, JacRow& jacRow,
                                  const ElementResidualVector& res, GridVariables& gridVariables)
    {}

    template<std::size_t otherId, class JacRow, class GridVariables,
             typename std::enable_if_t<(otherId != id), int> = 0>
    void assembleJacobianCoupling(Dune::index_constant<otherId> domainJ, JacRow& jacRow,
                                  const ElementResidualVector& res, GridVariables& gridVariables)
    {
        this->asImp_().assembleJacobianCoupling(domainJ, jacRow[domainJ], res, *std::get<otherId>(gridVariables));
    }

    /*!
     * \brief Assemble the residual only
     */
    void assembleResidual(SubSolutionVector& res)
    {
        this->asImp_().bindLocalViews();
        this->elemBcTypes().update(problem(), this->element(), this->fvGeometry());

        const auto residual = this->asImp_().assembleResidualImpl();

        for (const auto& scv : scvs(this->fvGeometry()))
            res[scv.dofIndex()] += residual[scv.indexInElement()];
    }


    ElementResidualVector evalLocalSourceResidual(const Element& neighbor, const ElementVolumeVariables& elemVolVars, const SubControlVolume& scv) const
    {
        const auto& curVolVars = elemVolVars[scv];
        auto source = this->localResidual().computeSource(problem(), neighbor, this->fvGeometry(), elemVolVars, scv);
        source *= scv.volume()*curVolVars.extrusionFactor();
        return source;
    }

    ElementResidualVector evalFluxResidual(const Element& neighbor,
                                         const SubControlVolumeFace& scvf) const
    {
        return this->localResidual().evalFlux(problem(), neighbor, this->fvGeometry(), this->curElemVolVars(), this->elemFluxVarsCache(), scvf);
    }

    const Problem& problem() const
    { return this->assembler().problem(domainId); }

    CouplingManager& couplingManager()
    { return couplingManager_; }

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
                const auto bcTypes = this->elemBcTypes()[scvI.indexInElement()];
                if (bcTypes.hasDirichlet())
                {
                    const auto dirichletValues = this->problem().dirichlet(this->element(), scvI);

                    // set the dirichlet conditions in residual and jacobian
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

private:
    CouplingManager& couplingManager_; //!< the coupling manager
};

/*!
 * \ingroup Assembly
 * \ingroup CCDiscretization
 * \brief A base class for all implicit local assemblers
 * \tparam TypeTag the TypeTag
 * \tparam Assembler the assembler type
 */
template<std::size_t id, class TypeTag, class Assembler, class Implementation>
class SubDomainBoxLocalAssemblerImplicitBase : public SubDomainBoxLocalAssemblerBase<id, TypeTag, Assembler, Implementation>
{
    using ParentType = SubDomainBoxLocalAssemblerBase<id, TypeTag, Assembler, Implementation>;
    using LocalResidual = typename GET_PROP_TYPE(TypeTag, LocalResidual);
    using ElementResidualVector = typename LocalResidual::ElementResidualVector;
    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
    using FVGridGeometry = typename GET_PROP_TYPE(TypeTag, FVGridGeometry);
    using SubControlVolume = typename FVGridGeometry::SubControlVolume;
    using Element = typename GridView::template Codim<0>::Entity;
    static constexpr auto domainId = Dune::index_constant<id>();
public:
    using ParentType::ParentType;

    void bindLocalViews()
    {
        // get some references for convenience
        auto& couplingManager = this->couplingManager();
        const auto& element = this->element();
        const auto& curSol = this->curSol()[domainId];
        auto&& fvGeometry = this->fvGeometry();
        auto&& curElemVolVars = this->curElemVolVars();
        auto&& elemFluxVarsCache = this->elemFluxVarsCache();

        // bind the caches
        couplingManager.bindCouplingContext(domainId, element, this->assembler());
        fvGeometry.bind(element);
        curElemVolVars.bind(element, fvGeometry, curSol);
        elemFluxVarsCache.bind(element, fvGeometry, curElemVolVars);
        if (!this->assembler().isStationaryProblem())
            this->prevElemVolVars().bindElement(element, fvGeometry, this->assembler().prevSol()[domainId]);
    }


    using ParentType::evalLocalFluxAndSourceResidual;
    ElementResidualVector evalLocalFluxAndSourceResidual() const
    { return this->evalLocalFluxAndSourceResidual(this->curElemVolVars()); }

    using ParentType::evalLocalSourceResidual;
    ElementResidualVector evalLocalSourceResidual(const Element& neighbor, const SubControlVolume& scv) const
    { return this->evalLocalSourceResidual(neighbor, this->curElemVolVars(), scv); }

    /*!
     * \brief Computes the residual
     * \return The element residual at the current solution.
     */
    ElementResidualVector assembleResidualImpl()
    { return this->evalLocalResidual(); }
};

/*!
 * \ingroup Assembly
 * \ingroup CCDiscretization
 * \brief An assembler for Jacobian and residual contribution per element (cell-centered methods)
 * \tparam TypeTag the TypeTag
 * \tparam DM the differentiation method to residual compute derivatives
 * \tparam implicit if to use an implicit or explicit time discretization
 */
template<std::size_t id, class TypeTag, class Assembler, DiffMethod DM = DiffMethod::numeric, bool implicit = true>
class SubDomainBoxLocalAssembler;

/*!
 * \ingroup Assembly
 * \ingroup CCDiscretization
 * \brief Cell-centered scheme local assembler using numeric differentiation and implicit time discretization
 */
template<std::size_t id, class TypeTag, class Assembler>
class SubDomainBoxLocalAssembler<id, TypeTag, Assembler, DiffMethod::numeric, /*implicit=*/true>
: public SubDomainBoxLocalAssemblerImplicitBase<id, TypeTag, Assembler,
            SubDomainBoxLocalAssembler<id, TypeTag, Assembler, DiffMethod::numeric, true> >
{
    using ThisType = SubDomainBoxLocalAssembler<id, TypeTag, Assembler, DiffMethod::numeric, /*implicit=*/true>;
    using ParentType = SubDomainBoxLocalAssemblerImplicitBase<id, TypeTag, Assembler, ThisType>;
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using VolumeVariables = typename GET_PROP_TYPE(TypeTag, VolumeVariables);
    using LocalResidual = typename GET_PROP_TYPE(TypeTag, LocalResidual);
    using ElementResidualVector = typename LocalResidual::ElementResidualVector;

    using FVGridGeometry = typename GET_PROP_TYPE(TypeTag, FVGridGeometry);
    using FVElementGeometry = typename FVGridGeometry::LocalView;
    using GridView = typename FVGridGeometry::GridView;
    using Element = typename GridView::template Codim<0>::Entity;

    enum { numEq = GET_PROP_TYPE(TypeTag, ModelTraits)::numEq() };
    enum { dim = GET_PROP_TYPE(TypeTag, GridView)::dimension };

    static constexpr bool enableGridFluxVarsCache = GET_PROP_VALUE(TypeTag, EnableGridFluxVariablesCache);
    static constexpr int maxNeighbors = 4*(2*dim);
    static constexpr auto domainI = Dune::index_constant<id>();

public:
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
        //////////////////////////////////////////////////////////////////////////////////////////////////
        // Calculate derivatives of all dofs in stencil with respect to the dofs in the element. In the //
        // neighboring elements we do so by computing the derivatives of the fluxes which depend on the //
        // actual element. In the actual element we evaluate the derivative of the entire residual.     //
        //////////////////////////////////////////////////////////////////////////////////////////////////

        // get some aliases for convenience
        const auto& element = this->element();
        const auto& fvGeometry = this->fvGeometry();
        const auto& curSol = this->curSol()[domainI];

        // get the vecor of the acutal element residuals
        const auto origResiduals = this->evalLocalResidual();

        // compute the derivatives of this element and add them to the Jacobian
        computeAndAddDerivatives(element, fvGeometry, origResiduals, curSol, A, gridVariables);

        // compute additional derivatives of this element with respect to other elements
        const auto eIdxI = fvGeometry.fvGridGeometry().elementMapper().index(element);
        for (const auto eIdxJ : this->couplingManager().getAdditionalDofDependencies(domainI, eIdxI))
        {
            const auto elementJ = fvGeometry.fvGridGeometry().element(eIdxJ);
            auto fvGeometryJ = localView(fvGeometry.fvGridGeometry());
            fvGeometryJ.bindElement(elementJ);
            computeAndAddDerivatives(elementJ, fvGeometryJ, origResiduals, curSol, A, gridVariables);
        }

        return origResiduals;
    }

    //////////////////////////////////////////////////////////////////////////////////////////////////
    //                                                                                              //
    // Calculate derivatives of all dofs in stencil with respect to the dofs in the element. In the //
    // neighboring elements we do so by computing the derivatives of the fluxes which depend on the //
    // actual element. In the actual element we evaluate the derivative of the entire residual.     //
    //                                                                                              //
    //////////////////////////////////////////////////////////////////////////////////////////////////
    template<class ElementResidual, class SolutionVector,
             class JacobianMatrixDiagBlock, class GridVariables>
    void computeAndAddDerivatives(const Element& element,
                                  const FVElementGeometry& fvGeometry,
                                  const ElementResidual& origResiduals,
                                  const SolutionVector& curSol,
                                  JacobianMatrixDiagBlock& A,
                                  GridVariables& gridVariables)
    {
        // create the element solution
        auto elemSol = elementSolution(element, curSol, fvGeometry.fvGridGeometry());

        auto partialDerivs = origResiduals;
        partialDerivs = 0.0;

        // calculation of the derivatives
        for (auto&& scv : scvs(fvGeometry))
        {
            // dof index and corresponding actual pri vars
            const auto dofIdx = scv.dofIndex();
            auto&& curElemVolVars = this->curElemVolVars();
            auto& curVolVars = this->getVolVarAccess(gridVariables.curGridVolVars(), curElemVolVars, scv);
            const VolumeVariables origVolVars(curVolVars);

            // calculate derivatives w.r.t to the privars at the dof at hand
            for (int pvIdx = 0; pvIdx < numEq; pvIdx++)
            {
                partialDerivs = 0.0;

                auto evalResiduals = [&](Scalar priVar)
                {
                    // update the volume variables and compute element residual
                    const auto localDofIndex = scv.indexInElement();
                    elemSol[localDofIndex][pvIdx] = priVar;
                    curVolVars.update(elemSol, this->problem(), element, scv);
                    this->couplingManager().updateCouplingContext(domainI, domainI, element, elemSol, localDofIndex, pvIdx, this->assembler());
                    return this->evalLocalResidual();
                };

                // derive the residuals numerically
                static const int numDiffMethod = getParam<int>("Assembly.NumericDifferenceMethod");
                NumericDifferentiation::partialDerivative(evalResiduals, elemSol[scv.indexInElement()][pvIdx], partialDerivs, origResiduals, numDiffMethod);

                // update the global stiffness matrix with the current partial derivatives
                for (auto&& scvJ : scvs(fvGeometry))
                {
                      for (int eqIdx = 0; eqIdx < numEq; eqIdx++)
                      {
                          // A[i][col][eqIdx][pvIdx] is the rate of change of
                          // the residual of equation 'eqIdx' at dof 'i'
                          // depending on the primary variable 'pvIdx' at dof
                          // 'col'.
                          A[scvJ.dofIndex()][dofIdx][eqIdx][pvIdx] += partialDerivs[scvJ.indexInElement()][eqIdx];
                      }
                }

                // restore the original state of the scv's volume variables
                curVolVars = origVolVars;

                // restore the original element solution
                elemSol[scv.indexInElement()][pvIdx] = curSol[scv.dofIndex()][pvIdx];
            }
        }
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
        ////////////////////////////////////////////////////////////////////////////////////////////////////////
        // Calculate derivatives of all dofs in the element with respect to all dofs in the coupling stencil. //
        ////////////////////////////////////////////////////////////////////////////////////////////////////////

        // get some aliases for convenience
        const auto& element = this->element();
        const auto& fvGeometry = this->fvGeometry();
        const auto& curElemVolVars = this->curElemVolVars();
        const auto& elemFluxVarsCache = this->elemFluxVarsCache();

        const auto& gridGeometryJ = this->assembler().fvGridGeometry(domainJ);
        const auto& curSolJ = this->curSol()[domainJ];

        // get element stencil informations
        const auto& stencil = this->couplingManager().couplingElementStencil(element, domainI, domainJ);
        for (const auto globalJ : stencil)
        {
            const auto& elementJ = gridGeometryJ.element(globalJ);
            auto elemSolJ = elementSolution(elementJ, curSolJ, gridGeometryJ);
            const auto origResidual = this->couplingManager().evalCouplingResidual(domainI, element, fvGeometry, curElemVolVars, this->elemBcTypes(), elemFluxVarsCache, this->assembler().localResidual(domainI),
                                                                                   domainJ, elementJ);

            // compute derivatives w.r.t. each dof in the coupled element
            for (const auto& dofData : this->couplingManager().coupledElementDofData(domainI, element, domainJ, globalJ))
            {
                // the element-local index of the current dof
                const auto localDofIdx = dofData.localIndex;

                for (int pvIdx = 0; pvIdx < JacobianBlock::block_type::cols; ++pvIdx)
                {
                    // store undeflected privars for reset
                    const auto origPriVarsJ = elemSolJ[localDofIdx];

                    auto evalCouplingResidual = [&](Scalar priVar)
                    {
                        elemSolJ[localDofIdx][pvIdx] = priVar;
                        this->couplingManager().updateCouplingContext(domainI, domainJ, elementJ, elemSolJ, localDofIdx, pvIdx, this->assembler());
                        return this->couplingManager().evalCouplingResidual(domainI, element, fvGeometry, curElemVolVars, this->elemBcTypes(), elemFluxVarsCache,
                                                                            this->assembler().localResidual(domainI), domainJ, elementJ);
                    };

                    // derive the residuals numerically
                    ElementResidualVector partialDerivs(element.subEntities(dim));
                    NumericDifferentiation::partialDerivative(evalCouplingResidual, origPriVarsJ[pvIdx], partialDerivs, origResidual);

                    // update the global stiffness matrix with the current partial derivatives
                    for (auto&& scv : scvs(fvGeometry))
                    {
                        for (int eqIdx = 0; eqIdx < numEq; eqIdx++)
                        {
                            // A[i][col][eqIdx][pvIdx] is the rate of change of
                            // the residual of equation 'eqIdx' at dof 'i'
                            // depending on the primary variable 'pvIdx' at dof
                            // 'col'.

                            // If the dof is coupled by a Dirichlet condition,
                            // set the derived value only once (i.e. overwrite existing values).
                            // For other dofs, add the contribution of the partial derivative.
                            const auto bcTypes = this->elemBcTypes()[scv.localDofIndex()];
                            if (bcTypes.isDirichlet(eqIdx))
                                A[scv.dofIndex()][dofData.index][eqIdx][pvIdx] = partialDerivs[scv.localDofIndex()][eqIdx];
                            else
                                A[scv.dofIndex()][dofData.index][eqIdx][pvIdx] += partialDerivs[scv.localDofIndex()][eqIdx];
                        }
                    }

                    // restore the current element solution
                    elemSolJ[localDofIdx][pvIdx] = origPriVarsJ[pvIdx];

                    // restore the undeflected state of the coupling context
                    this->couplingManager().updateCouplingContext(domainI, domainJ, elementJ, elemSolJ, localDofIdx, pvIdx, this->assembler());
                }
            }
        }
    }
};

} // end namespace Dumux

#endif
