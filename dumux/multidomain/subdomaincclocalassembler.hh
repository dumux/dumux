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
 * \ingroup CCDiscretization
 * \ingroup MultiDomain
 * \brief A multidomain local assembler for Jacobian and residual contribution per element (cell-centered methods)
 */
#ifndef DUMUX_MULTIDOMAIN_CC_LOCAL_ASSEMBLER_HH
#define DUMUX_MULTIDOMAIN_CC_LOCAL_ASSEMBLER_HH

#include <dune/common/indices.hh>
#include <dune/common/reservedvector.hh>
#include <dune/common/hybridutilities.hh>
#include <dune/grid/common/gridenums.hh> // for GhostEntity
#include <dune/istl/matrixindexset.hh>

#include <dumux/common/reservedblockvector.hh>
#include <dumux/common/properties.hh>
#include <dumux/common/parameters.hh>
#include <dumux/common/numericdifferentiation.hh>
#include <dumux/discretization/elementsolution.hh>
#include <dumux/discretization/extrusion.hh>
#include <dumux/assembly/numericepsilon.hh>
#include <dumux/assembly/diffmethod.hh>
#include <dumux/assembly/fvlocalassemblerbase.hh>

namespace Dumux {

/*!
 * \ingroup Assembly
 * \ingroup CCDiscretization
 * \ingroup MultiDomain
 * \brief A base class for all multidomain local assemblers
 * \tparam id the id of the sub domain
 * \tparam TypeTag the TypeTag
 * \tparam Assembler the assembler type
 * \tparam Implementation the actual assembler implementation
 * \tparam implicit Specifies whether the time discretization is implicit or not not (i.e. explicit)
 */
template<std::size_t id, class TypeTag, class Assembler, class Implementation, bool implicit>
class SubDomainCCLocalAssemblerBase : public FVLocalAssemblerBase<TypeTag, Assembler,Implementation, implicit>
{
    using ParentType = FVLocalAssemblerBase<TypeTag, Assembler,Implementation, implicit>;

    using Problem = GetPropType<TypeTag, Properties::Problem>;
    using LocalResidualValues = GetPropType<TypeTag, Properties::NumEqVector>;
    using JacobianMatrix = GetPropType<TypeTag, Properties::JacobianMatrix>;
    using SolutionVector = typename Assembler::SolutionVector;
    using SubSolutionVector = GetPropType<TypeTag, Properties::SolutionVector>;
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

        // for the diagonal jacobian block
        const auto globalI = this->fvGeometry().gridGeometry().elementMapper().index(this->element());
        res[globalI] = this->asImp_().assembleJacobianAndResidualImplInverse(jacRow[domainId], *std::get<domainId>(gridVariables));

        // for the coupling blocks
        using namespace Dune::Hybrid;
        forEach(integralRange(Dune::Hybrid::size(jacRow)), [&](auto&& i)
        {
            if (i != id)
                this->assembleJacobianCoupling(i, jacRow, res[globalI], gridVariables);
        });
    }

    /*!
     * \brief Assemble the entries in a coupling block of the jacobian.
     *        There is no coupling block between a domain and itself.
     */
    template<std::size_t otherId, class JacRow, class GridVariables,
             typename std::enable_if_t<(otherId == id), int> = 0>
    void assembleJacobianCoupling(Dune::index_constant<otherId> domainJ, JacRow& jacRow,
                                  const LocalResidualValues& res, GridVariables& gridVariables)
    {}

    /*!
     * \brief Assemble the entries in a coupling block of the jacobian.
     */
    template<std::size_t otherId, class JacRow, class GridVariables,
             typename std::enable_if_t<(otherId != id), int> = 0>
    void assembleJacobianCoupling(Dune::index_constant<otherId> domainJ, JacRow& jacRow,
                                  const LocalResidualValues& res, GridVariables& gridVariables)
    {
        this->asImp_().assembleJacobianCoupling(domainJ, jacRow[domainJ], res, *std::get<domainId>(gridVariables));
    }

    /*!
     * \brief Assemble the residual only
     */
    void assembleResidual(SubSolutionVector& res)
    {
        this->asImp_().bindLocalViews();
        const auto globalI = this->fvGeometry().gridGeometry().elementMapper().index(this->element());
        res[globalI] = this->evalLocalResidual()[0]; // forward to the internal implementation
    }

    /*!
     * \brief Evaluates the local source term for an element and given element volume variables
     */
    ElementResidualVector evalLocalSourceResidual(const Element& element, const ElementVolumeVariables& elemVolVars) const
    {
        // initialize the residual vector for all scvs in this element
        ElementResidualVector residual(this->fvGeometry().numScv());

        // evaluate the volume terms (storage + source terms)
        // forward to the local residual specialized for the discretization methods
        for (auto&& scv : scvs(this->fvGeometry()))
        {
            const auto& curVolVars = elemVolVars[scv];
            auto source = this->localResidual().computeSource(problem(), element, this->fvGeometry(), elemVolVars, scv);
            source *= -Extrusion::volume(scv)*curVolVars.extrusionFactor();
            residual[scv.indexInElement()] = std::move(source);
        }

        return residual;
    }

    /*!
     * \brief Evaluates the local source term depending on time discretization scheme
     */
    ElementResidualVector evalLocalSourceResidual(const Element& neighbor) const
    { return this->evalLocalSourceResidual(neighbor, implicit ? this->curElemVolVars() : this->prevElemVolVars()); }

    /*!
     * \brief Evaluates the storage terms within the element
     */
    LocalResidualValues evalLocalStorageResidual() const
    {
        return this->localResidual().evalStorage(this->element(), this->fvGeometry(), this->prevElemVolVars(), this->curElemVolVars())[0];
    }

    /*!
     * \brief Evaluates the fluxes depending on the chose time discretization scheme
     */
    LocalResidualValues evalFluxResidual(const Element& neighbor,
                                         const SubControlVolumeFace& scvf) const
    {
        const auto& elemVolVars = implicit ? this->curElemVolVars() : this->prevElemVolVars();
        return this->localResidual().evalFlux(problem(), neighbor, this->fvGeometry(), elemVolVars, this->elemFluxVarsCache(), scvf);
    }

    /*!
     * \brief Prepares all local views necessary for local assembly.
     */
    void bindLocalViews()
    {
        // get some references for convenience
        const auto& element = this->element();
        const auto& curSol = this->curSol()[domainId];
        auto&& fvGeometry = this->fvGeometry();
        auto&& curElemVolVars = this->curElemVolVars();
        auto&& elemFluxVarsCache = this->elemFluxVarsCache();

        // bind the caches
        couplingManager_.bindCouplingContext(domainId, element, this->assembler());
        fvGeometry.bind(element);

        if (implicit)
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
 * \ingroup CCDiscretization
 * \ingroup MultiDomain
 * \brief The cell-centered scheme multidomain local assembler
 * \tparam id the id of the sub domain
 * \tparam TypeTag the TypeTag
 * \tparam Assembler the assembler type
 * \tparam DM the numeric differentiation method
 * \tparam implicit whether the assembler is explicit or implicit in time
 */
template<std::size_t id, class TypeTag, class Assembler, DiffMethod DM = DiffMethod::numeric, bool implicit = true>
class SubDomainCCLocalAssembler;

/*!
 * \ingroup Assembly
 * \ingroup CCDiscretization
 * \ingroup MultiDomain
 * \brief Cell-centered scheme multidomain local assembler using numeric differentiation and implicit time discretization
 */
template<std::size_t id, class TypeTag, class Assembler>
class SubDomainCCLocalAssembler<id, TypeTag, Assembler, DiffMethod::numeric, /*implicit=*/true>
: public SubDomainCCLocalAssemblerBase<id, TypeTag, Assembler,
            SubDomainCCLocalAssembler<id, TypeTag, Assembler, DiffMethod::numeric, true>, true>
{
    using ThisType = SubDomainCCLocalAssembler<id, TypeTag, Assembler, DiffMethod::numeric, /*implicit=*/true>;
    using ParentType = SubDomainCCLocalAssemblerBase<id, TypeTag, Assembler, ThisType, /*implicit=*/true>;
    using Problem = GetPropType<TypeTag, Properties::Problem>;

    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using LocalResidualValues = GetPropType<TypeTag, Properties::NumEqVector>;

    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using FVElementGeometry = typename GridGeometry::LocalView;
    using GridView = typename GridGeometry::GridView;
    using Element = typename GridView::template Codim<0>::Entity;

    enum { numEq = GetPropType<TypeTag, Properties::ModelTraits>::numEq() };
    enum { dim = GridView::dimension };

    static constexpr bool enableGridFluxVarsCache = getPropValue<TypeTag, Properties::EnableGridFluxVariablesCache>();
    static constexpr bool enableGridVolVarsCache = getPropValue<TypeTag, Properties::EnableGridVolumeVariablesCache>();
    static constexpr int maxElementStencilSize = GridGeometry::maxElementStencilSize;
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
    LocalResidualValues assembleJacobianAndResidualImplInverse(JacobianMatrixDiagBlock& A, GridVariables& gridVariables)
    {
        //////////////////////////////////////////////////////////////////////////////////////////////////
        // Calculate derivatives of all dofs in stencil with respect to the dofs in the element. In the //
        // neighboring elements we do so by computing the derivatives of the fluxes which depend on the //
        // actual element. In the actual element we evaluate the derivative of the entire residual.     //
        //////////////////////////////////////////////////////////////////////////////////////////////////

        // get some aliases for convenience
        const auto& element = this->element();
        const auto& fvGeometry = this->fvGeometry();
        const auto& gridGeometry = fvGeometry.gridGeometry();
        auto&& curElemVolVars = this->curElemVolVars();
        auto&& elemFluxVarsCache = this->elemFluxVarsCache();

        // get stencil informations
        const auto globalI = gridGeometry.elementMapper().index(element);
        const auto& connectivityMap = gridGeometry.connectivityMap();
        const auto numNeighbors = connectivityMap[globalI].size();

        // container to store the neighboring elements
        Dune::ReservedVector<Element, maxElementStencilSize> neighborElements;
        neighborElements.resize(numNeighbors);

        // assemble the undeflected residual
        using Residuals = ReservedBlockVector<LocalResidualValues, maxElementStencilSize>;
        Residuals origResiduals(numNeighbors + 1); origResiduals = 0.0;
        origResiduals[0] = this->evalLocalResidual()[0];

        // lambda for convenient evaluation of the fluxes across scvfs in the neighbors
        // if the neighbor is a ghost we don't want to add anything to their residual
        // so we return 0 and omit computing the flux
        const auto evalNeighborFlux = [&] (const auto& neighbor, const auto& scvf)
        {
            if (neighbor.partitionType() == Dune::GhostEntity)
                return LocalResidualValues(0.0);
            else
                return this->evalFluxResidual(neighbor, scvf);
        };

        // get the elements in which we need to evaluate the fluxes
        // and calculate these in the undeflected state
        unsigned int j = 1;
        for (const auto& dataJ : connectivityMap[globalI])
        {
            neighborElements[j-1] = gridGeometry.element(dataJ.globalJ);
            for (const auto scvfIdx : dataJ.scvfsJ)
                origResiduals[j] += evalNeighborFlux(neighborElements[j-1], fvGeometry.scvf(scvfIdx));

            ++j;
        }

        // reference to the element's scv (needed later) and corresponding vol vars
        const auto& scv = fvGeometry.scv(globalI);
        auto& curVolVars = ParentType::getVolVarAccess(gridVariables.curGridVolVars(), curElemVolVars, scv);

        // save a copy of the original privars and vol vars in order
        // to restore the original solution after deflection
        const auto& curSol = this->curSol()[domainI];
        const auto origPriVars = curSol[globalI];
        const auto origVolVars = curVolVars;

        // element solution container to be deflected
        auto elemSol = elementSolution<FVElementGeometry>(origPriVars);

        // derivatives in the neighbors with repect to the current elements
        // in index 0 we save the derivative of the element residual with respect to it's own dofs
        Residuals partialDerivs(numNeighbors + 1);

        for (int pvIdx = 0; pvIdx < numEq; ++pvIdx)
        {
            partialDerivs = 0.0;

            auto evalResiduals = [&](Scalar priVar)
            {
                Residuals partialDerivsTmp(numNeighbors + 1);
                partialDerivsTmp = 0.0;
                // update the volume variables and the flux var cache
                elemSol[0][pvIdx] = priVar;
                this->couplingManager().updateCouplingContext(domainI, *this, domainI, globalI, elemSol[0], pvIdx);
                curVolVars.update(elemSol, this->problem(), element, scv);
                elemFluxVarsCache.update(element, fvGeometry, curElemVolVars);
                if (enableGridFluxVarsCache)
                    gridVariables.gridFluxVarsCache().updateElement(element, fvGeometry, curElemVolVars);

                // calculate the residual with the deflected primary variables
                partialDerivsTmp[0] = this->evalLocalResidual()[0];

                // calculate the fluxes in the neighbors with the deflected primary variables
                for (std::size_t k = 0; k < connectivityMap[globalI].size(); ++k)
                    for (auto scvfIdx : connectivityMap[globalI][k].scvfsJ)
                        partialDerivsTmp[k+1] += evalNeighborFlux(neighborElements[k], fvGeometry.scvf(scvfIdx));

                return partialDerivsTmp;
            };

            // derive the residuals numerically
            static const int numDiffMethod = getParamFromGroup<int>(this->problem().paramGroup(), "Assembly.NumericDifferenceMethod");
            static const NumericEpsilon<Scalar, numEq> eps_{this->problem().paramGroup()};
            NumericDifferentiation::partialDerivative(evalResiduals, elemSol[0][pvIdx], partialDerivs, origResiduals,
                                                      eps_(elemSol[0][pvIdx], pvIdx), numDiffMethod);

            // Correct derivative for ghost elements, i.e. set a 1 for the derivative w.r.t. the
            // current primary variable and a 0 elsewhere. As we always solve for a delta of the
            // solution with repect to the initial one, this results in a delta of 0 for ghosts.
            if (this->elementIsGhost())
            {
                partialDerivs[0] = 0.0;
                partialDerivs[0][pvIdx] = 1.0;
            }

            // add the current partial derivatives to the global jacobian matrix
            if constexpr (Problem::enableInternalDirichletConstraints())
            {
                // check if own element has internal Dirichlet constraint
                const auto internalDirichletConstraintsOwnElement = this->problem().hasInternalDirichletConstraint(this->element(), scv);
                const auto dirichletValues = this->problem().internalDirichlet(this->element(), scv);

                for (int eqIdx = 0; eqIdx < numEq; ++eqIdx)
                {
                    if (internalDirichletConstraintsOwnElement[eqIdx])
                    {
                        origResiduals[0][eqIdx] = origVolVars.priVars()[eqIdx] - dirichletValues[eqIdx];
                        A[globalI][globalI][eqIdx][pvIdx] = (eqIdx == pvIdx) ? 1.0 : 0.0;
                    }
                    else
                        A[globalI][globalI][eqIdx][pvIdx] += partialDerivs[0][eqIdx];
                }

                // off-diagonal entries
                j = 1;
                for (const auto& dataJ : connectivityMap[globalI])
                {
                    const auto& neighborElement = neighborElements[j-1];
                    const auto& neighborScv = fvGeometry.scv(dataJ.globalJ);
                    const auto internalDirichletConstraintsNeighbor = this->problem().hasInternalDirichletConstraint(neighborElement, neighborScv);

                    for (int eqIdx = 0; eqIdx < numEq; ++eqIdx)
                    {
                        if (internalDirichletConstraintsNeighbor[eqIdx])
                            A[dataJ.globalJ][globalI][eqIdx][pvIdx] = 0.0;
                        else
                            A[dataJ.globalJ][globalI][eqIdx][pvIdx] += partialDerivs[j][eqIdx];
                    }

                    ++j;
                }
            }
            else // no internal Dirichlet constraints specified
            {
                for (int eqIdx = 0; eqIdx < numEq; eqIdx++)
                {
                    // the diagonal entries
                    A[globalI][globalI][eqIdx][pvIdx] += partialDerivs[0][eqIdx];

                    // off-diagonal entries
                    j = 1;
                    for (const auto& dataJ : connectivityMap[globalI])
                        A[dataJ.globalJ][globalI][eqIdx][pvIdx] += partialDerivs[j++][eqIdx];
                }
            }

            // restore the original state of the scv's volume variables
            curVolVars = origVolVars;

            // restore the current element solution
            elemSol[0][pvIdx] = origPriVars[pvIdx];

            // restore the undeflected state of the coupling context
            this->couplingManager().updateCouplingContext(domainI, *this, domainI, globalI, elemSol[0], pvIdx);
        }

        // restore original state of the flux vars cache in case of global caching.
        // This has to be done in order to guarantee that everything is in an undeflected
        // state before the assembly of another element is called. In the case of local caching
        // this is obsolete because the elemFluxVarsCache used here goes out of scope after this.
        // We only have to do this for the last primary variable, for all others the flux var cache
        // is updated with the correct element volume variables before residual evaluations
        if (enableGridFluxVarsCache)
            gridVariables.gridFluxVarsCache().updateElement(element, fvGeometry, curElemVolVars);

        // compute potential additional derivatives causes by the coupling
        typename ParentType::LocalResidual::ElementResidualVector origElementResidual({origResiduals[0]});
        this->couplingManager().evalAdditionalDomainDerivatives(domainI, *this, origElementResidual, A, gridVariables);

        // return the original residual
        return origResiduals[0];
    }

    /*!
     * \brief Computes the derivatives with respect to the given element and adds them
     *        to the global matrix.
     */
    template<std::size_t otherId, class JacobianBlock, class GridVariables>
    void assembleJacobianCoupling(Dune::index_constant<otherId> domainJ, JacobianBlock& A,
                                  const LocalResidualValues& res, GridVariables& gridVariables)
    {
        ////////////////////////////////////////////////////////////////////////////////////////////////////////
        // Calculate derivatives of all dofs in the element with respect to all dofs in the coupling stencil. //
        ////////////////////////////////////////////////////////////////////////////////////////////////////////

        // get some aliases for convenience
        const auto& element = this->element();
        const auto& fvGeometry = this->fvGeometry();
        const auto& gridGeometry = fvGeometry.gridGeometry();
        auto&& curElemVolVars = this->curElemVolVars();
        auto&& elemFluxVarsCache = this->elemFluxVarsCache();

        // get stencil informations
        const auto globalI = gridGeometry.elementMapper().index(element);
        const auto& stencil = this->couplingManager().couplingStencil(domainI, element, domainJ);
        const auto& curSolJ = this->curSol()[domainJ];

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

/*!
 * \ingroup Assembly
 * \ingroup CCDiscretization
 * \ingroup MultiDomain
 * \brief Cell-centered scheme multidomain local assembler using numeric differentiation and explicit time discretization
 */
template<std::size_t id, class TypeTag, class Assembler>
class SubDomainCCLocalAssembler<id, TypeTag, Assembler, DiffMethod::numeric, /*implicit=*/false>
: public SubDomainCCLocalAssemblerBase<id, TypeTag, Assembler,
            SubDomainCCLocalAssembler<id, TypeTag, Assembler, DiffMethod::numeric, false>, false>
{
    using ThisType = SubDomainCCLocalAssembler<id, TypeTag, Assembler, DiffMethod::numeric, /*implicit=*/false>;
    using ParentType = SubDomainCCLocalAssemblerBase<id, TypeTag, Assembler, ThisType, /*implicit=*/false>;

    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using LocalResidualValues = GetPropType<TypeTag, Properties::NumEqVector>;
    using Problem = GetPropType<TypeTag, Properties::Problem>;

    static constexpr int numEq = GetPropType<TypeTag, Properties::ModelTraits>::numEq();
    static constexpr auto domainI = Dune::index_constant<id>();

public:
    using ParentType::ParentType;

    /*!
     * \brief Computes the derivatives with respect to the given element and adds them
     *        to the global matrix.
     *
     * \note In an explicit scheme, only the storage terms need to be differentiated.
     *       Thus, this can be done as in the uncoupled case as the coupling can only
     *       enter sources or fluxes.
     *
     * \return The element residual at the current solution.
     */
    template<class JacobianMatrixDiagBlock, class GridVariables>
    LocalResidualValues assembleJacobianAndResidualImplInverse(JacobianMatrixDiagBlock& A, GridVariables& gridVariables)
    {
        if (this->assembler().isStationaryProblem())
            DUNE_THROW(Dune::InvalidStateException, "Using explicit jacobian assembler with stationary local residual");

        // assemble the undeflected residual
        const auto residual = this->evalLocalResidual()[0];
        const auto storageResidual = this->evalLocalStorageResidual();

        //////////////////////////////////////////////////////////////////////////////////////////////////
        // Calculate derivatives of all dofs in stencil with respect to the dofs in the element. In the //
        // neighboring elements all derivatives are zero. For the assembled element only the storage    //
        // derivatives are non-zero.                                                                    //
        //////////////////////////////////////////////////////////////////////////////////////////////////

        // get some aliases for convenience
        const auto& element = this->element();
        const auto& fvGeometry = this->fvGeometry();
        const auto& gridGeometry = fvGeometry.gridGeometry();
        auto&& curElemVolVars = this->curElemVolVars();

        // reference to the element's scv (needed later) and corresponding vol vars
        const auto globalI = gridGeometry.elementMapper().index(element);
        const auto& scv = fvGeometry.scv(globalI);
        auto& curVolVars = ParentType::getVolVarAccess(gridVariables.curGridVolVars(), curElemVolVars, scv);

        // save a copy of the original privars and vol vars in order
        // to restore the original solution after deflection
        const auto& curSol = this->curSol()[domainI];
        const auto origPriVars = curSol[globalI];
        const auto origVolVars = curVolVars;

        // element solution container to be deflected
        auto elemSol = elementSolution(element, curSol, gridGeometry);

        // derivatives in the neighbors with repect to the current elements
        LocalResidualValues partialDeriv;
        for (int pvIdx = 0; pvIdx < numEq; ++pvIdx)
        {
            // reset derivatives of element dof with respect to itself
            partialDeriv = 0.0;

            auto evalStorage = [&](Scalar priVar)
            {
                // update the volume variables and calculate
                // the residual with the deflected primary variables
                elemSol[0][pvIdx] = priVar;
                curVolVars.update(elemSol, this->problem(), element, scv);
                return this->evalLocalStorageResidual();
            };

            // for non-ghosts compute the derivative numerically
            if (!this->elementIsGhost())
            {
                static const NumericEpsilon<Scalar, numEq> eps_{this->problem().paramGroup()};
                static const int numDiffMethod = getParamFromGroup<int>(this->problem().paramGroup(), "Assembly.NumericDifferenceMethod");
                NumericDifferentiation::partialDerivative(evalStorage, elemSol[0][pvIdx], partialDeriv, storageResidual,
                                                          eps_(elemSol[0][pvIdx], pvIdx), numDiffMethod);
            }

            // for ghost elements we assemble a 1.0 where the primary variable and zero everywhere else
            // as we always solve for a delta of the solution with repect to the initial solution this
            // results in a delta of zero for ghosts
            else partialDeriv[pvIdx] = 1.0;

            // add the current partial derivatives to the global jacobian matrix
            // only diagonal entries for explicit jacobians
            if constexpr (Problem::enableInternalDirichletConstraints())
            {
                // check if own element has internal Dirichlet constraint
                const auto internalDirichletConstraints = this->problem().hasInternalDirichletConstraint(this->element(), scv);
                const auto dirichletValues = this->problem().internalDirichlet(this->element(), scv);

                for (int eqIdx = 0; eqIdx < numEq; ++eqIdx)
                {
                    if (internalDirichletConstraints[eqIdx])
                    {
                        residual[eqIdx] = origVolVars.priVars()[eqIdx] - dirichletValues[eqIdx];
                        A[globalI][globalI][eqIdx][pvIdx] = (eqIdx == pvIdx) ? 1.0 : 0.0;
                    }
                    else
                        A[globalI][globalI][eqIdx][pvIdx] += partialDeriv[eqIdx];
                }
            }
            else
            {
                for (int eqIdx = 0; eqIdx < numEq; eqIdx++)
                    A[globalI][globalI][eqIdx][pvIdx] += partialDeriv[eqIdx];
            }

            // restore the original state of the scv's volume variables
            curVolVars = origVolVars;

            // restore the current element solution
            elemSol[0][pvIdx] = origPriVars[pvIdx];
        }

        // return the original residual
        return residual;
    }

    /*!
     * \brief Computes the coupling derivatives with respect to the given element and adds them
     *        to the global matrix.
     * \note Since the coupling can only enter sources or fluxes and these are evaluated on
     *       the old time level (explicit scheme), the coupling blocks are empty.
     */
    template<std::size_t otherId, class JacobianBlock, class GridVariables>
    void assembleJacobianCoupling(Dune::index_constant<otherId> domainJ, JacobianBlock& A,
                                  const LocalResidualValues& res, GridVariables& gridVariables)
    {}
};

/*!
 * \ingroup Assembly
 * \ingroup CCDiscretization
 * \brief Cell-centered scheme local assembler using analytic differentiation and implicit time discretization
 */
template<std::size_t id, class TypeTag, class Assembler>
class SubDomainCCLocalAssembler<id, TypeTag, Assembler, DiffMethod::analytic, /*implicit=*/true>
: public SubDomainCCLocalAssemblerBase<id, TypeTag, Assembler,
            SubDomainCCLocalAssembler<id, TypeTag, Assembler, DiffMethod::analytic, true>, true>
{
    using ThisType = SubDomainCCLocalAssembler<id, TypeTag, Assembler, DiffMethod::analytic, /*implicit=*/true>;
    using ParentType = SubDomainCCLocalAssemblerBase<id, TypeTag, Assembler, ThisType, /*implicit=*/true>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using LocalResidualValues = GetPropType<TypeTag, Properties::NumEqVector>;
    using GridView = typename GetPropType<TypeTag, Properties::GridGeometry>::GridView;
    using Element = typename GridView::template Codim<0>::Entity;
    using Problem = GetPropType<TypeTag, Properties::Problem>;

    enum { numEq = GetPropType<TypeTag, Properties::ModelTraits>::numEq() };
    enum { dim = GridView::dimension };

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
    LocalResidualValues assembleJacobianAndResidualImpl(JacobianMatrixDiagBlock& A, GridVariables& gridVariables)
    {
        // treat ghost separately, we always want zero update for ghosts
        if (this->elementIsGhost())
        {
            const auto globalI = this->assembler().gridGeometry(domainI).elementMapper().index(this->element());
            for (int pvIdx = 0; pvIdx < numEq; ++pvIdx)
                A[globalI][globalI][pvIdx][pvIdx] = 1.0;

            // return zero residual
            return LocalResidualValues(0.0);
        }

        // get some aliases for convenience
        const auto& problem = this->problem();
        const auto& element = this->element();
        const auto& fvGeometry = this->fvGeometry();
        const auto& curElemVolVars = this->curElemVolVars();
        const auto& elemFluxVarsCache = this->elemFluxVarsCache();

        // get reference to the element's current vol vars
        const auto globalI = this->assembler().gridGeometry(domainI).elementMapper().index(element);
        const auto& scv = fvGeometry.scv(globalI);
        const auto& volVars = curElemVolVars[scv];

        // if the problem is instationary, add derivative of storage term
        if (!this->assembler().isStationaryProblem())
            this->localResidual().addStorageDerivatives(A[globalI][globalI], problem, element, fvGeometry, volVars, scv);

        // add source term derivatives
        this->localResidual().addSourceDerivatives(A[globalI][globalI], problem, element, fvGeometry, volVars, scv);

        // add flux derivatives for each scvf
        for (const auto& scvf : scvfs(fvGeometry))
        {
            // inner faces
            if (!scvf.boundary())
                this->localResidual().addFluxDerivatives(A[globalI], problem, element, fvGeometry, curElemVolVars, elemFluxVarsCache, scvf);

            // boundary faces
            else
            {
                const auto& bcTypes = problem.boundaryTypes(element, scvf);

                // add Dirichlet boundary flux derivatives
                if (bcTypes.hasDirichlet() && !bcTypes.hasNeumann())
                    this->localResidual().addCCDirichletFluxDerivatives(A[globalI], problem, element, fvGeometry, curElemVolVars, elemFluxVarsCache, scvf);

                // add Robin ("solution dependent Neumann") boundary flux derivatives
                else if (bcTypes.hasNeumann() && !bcTypes.hasDirichlet())
                    this->localResidual().addRobinFluxDerivatives(A[globalI], problem, element, fvGeometry, curElemVolVars, elemFluxVarsCache, scvf);

                else
                    DUNE_THROW(Dune::NotImplemented, "Mixed boundary conditions. Use pure boundary conditions by converting Dirichlet BCs to Robin BCs");
            }
        }

        if constexpr (Problem::enableInternalDirichletConstraints())
        {
            // check if own element has internal Dirichlet constraint
            const auto internalDirichletConstraints = this->problem().hasInternalDirichletConstraint(this->element(), scv);
            const auto dirichletValues = this->problem().internalDirichlet(this->element(), scv);

            auto residual = this->evalLocalResidual()[0];

            for (int pvIdx = 0; pvIdx < numEq; ++pvIdx)
            {
                for (int eqIdx = 0; eqIdx < numEq; ++eqIdx)
                {
                    if (internalDirichletConstraints[eqIdx])
                    {
                        residual[eqIdx] = volVars.priVars()[eqIdx] - dirichletValues[eqIdx];
                        A[globalI][globalI][eqIdx][pvIdx] = (eqIdx == pvIdx) ? 1.0 : 0.0;

                        // inner faces
                        for (const auto& scvf : scvfs(fvGeometry))
                            if (!scvf.boundary())
                                A[globalI][fvGeometry.scv(scvf.outsideScvIdx()).dofIndex()][eqIdx][pvIdx] = 0.0;
                    }
                }
            }
            return residual;
        }
        else
            // return element residual
            return this->evalLocalResidual()[0];
    }

    /*!
     * \brief Computes the derivatives with respect to the given element and adds them
     *        to the global matrix.
     *
     * \return The element residual at the current solution.
     */
    template<std::size_t otherId, class JacobianBlock, class GridVariables>
    void assembleJacobianCoupling(Dune::index_constant<otherId> domainJ, JacobianBlock& A,
                                  const LocalResidualValues& res, GridVariables& gridVariables)
    {
        ////////////////////////////////////////////////////////////////////////////////////////////////////////
        // Calculate derivatives of all dofs in the element with respect to all dofs in the coupling stencil. //
        ////////////////////////////////////////////////////////////////////////////////////////////////////////

        // get some aliases for convenience
        const auto& element = this->element();
        const auto& fvGeometry = this->fvGeometry();
        const auto& gridGeometry = fvGeometry.gridGeometry();
        auto&& curElemVolVars = this->curElemVolVars();
        // auto&& elemFluxVarsCache = this->elemFluxVarsCache();

        // get stencil informations
        const auto globalI = gridGeometry.elementMapper().index(element);
        const auto& stencil = this->couplingManager().couplingStencil(domainI, element, domainJ);

        for (const auto globalJ : stencil)
        {
            const auto& elementJ = this->assembler().gridGeometry(domainJ).element(globalJ);
            this->couplingManager().addCouplingDerivatives(A[globalI][globalJ], domainI, element, fvGeometry, curElemVolVars, domainJ, elementJ);

            // add the current partial derivatives to the global jacobian matrix
            if constexpr (Problem::enableInternalDirichletConstraints())
            {
                const auto& scv = fvGeometry.scv(globalI);
                const auto internalDirichletConstraints = this->problem().hasInternalDirichletConstraint(this->element(), scv);
                if (internalDirichletConstraints.any())
                {
                    for (int pvIdx = 0; pvIdx < numEq; ++pvIdx)
                        for (int eqIdx = 0; eqIdx < numEq; ++eqIdx)
                            if (internalDirichletConstraints[eqIdx])
                                A[globalI][globalJ][eqIdx][pvIdx] = 0.0;
                }
            }
        }
    }
};

} // end namespace Dumux

#endif
