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
 * \brief An assembler for Jacobian and residual contribution per element (cell-centered methods)
 */
#ifndef DUMUX_CC_LOCAL_ASSEMBLER_HH
#define DUMUX_CC_LOCAL_ASSEMBLER_HH

#include <dune/common/reservedvector.hh>
#include <dune/grid/common/gridenums.hh> // for GhostEntity
#include <dune/istl/matrixindexset.hh>

#include <dumux/common/reservedblockvector.hh>
#include <dumux/common/properties.hh>
#include <dumux/common/parameters.hh>
#include <dumux/common/numericdifferentiation.hh>
#include <dumux/assembly/numericepsilon.hh>
#include <dumux/assembly/diffmethod.hh>
#include <dumux/assembly/fvlocalassemblerbase.hh>
#include <dumux/assembly/entitycolor.hh>
#include <dumux/assembly/partialreassembler.hh>
#include <dumux/discretization/fluxstencil.hh>
#include <dumux/discretization/cellcentered/elementsolution.hh>

namespace Dumux {

/*!
 * \ingroup Assembly
 * \ingroup CCDiscretization
 * \brief A base class for all local cell-centered assemblers
 * \tparam TypeTag The TypeTag
 * \tparam Assembler The assembler type
 * \tparam Implementation The actual implementation
 * \tparam implicit Specifies whether the time discretization is implicit or not (i.e. explicit)
 */
template<class TypeTag, class Assembler, class Implementation, bool implicit>
class CCLocalAssemblerBase : public FVLocalAssemblerBase<TypeTag, Assembler, Implementation, implicit>
{
    using ParentType = FVLocalAssemblerBase<TypeTag, Assembler, Implementation, implicit>;
    using GridView = typename GetPropType<TypeTag, Properties::GridGeometry>::GridView;
    using JacobianMatrix = GetPropType<TypeTag, Properties::JacobianMatrix>;
    using GridVariables = GetPropType<TypeTag, Properties::GridVariables>;
    using SolutionVector = GetPropType<TypeTag, Properties::SolutionVector>;
    using ElementVolumeVariables = typename GetPropType<TypeTag, Properties::GridVolumeVariables>::LocalView;
    using NumEqVector = GetPropType<TypeTag, Properties::NumEqVector>;

public:

    using ParentType::ParentType;

    /*!
     * \brief Computes the derivatives with respect to the given element and adds them
     *        to the global matrix. The element residual is written into the right hand side.
     */
    template <class PartialReassembler = DefaultPartialReassembler>
    void assembleJacobianAndResidual(JacobianMatrix& jac, SolutionVector& res, GridVariables& gridVariables,
                                     const PartialReassembler* partialReassembler)
    {
        this->asImp_().bindLocalViews();
        const auto globalI = this->assembler().gridGeometry().elementMapper().index(this->element());
        if (partialReassembler
            && partialReassembler->elementColor(globalI) == EntityColor::green)
        {
            res[globalI] = this->asImp_().evalLocalResidual()[0]; // forward to the internal implementation
        }
        else
        {
            res[globalI] = this->asImp_().assembleJacobianAndResidualImpl(jac, gridVariables); // forward to the internal implementation
        }
    }

    /*!
     * \brief Computes the derivatives with respect to the given element and adds them
     *        to the global matrix.
     */
    void assembleJacobian(JacobianMatrix& jac, GridVariables& gridVariables)
    {
        this->asImp_().bindLocalViews();
        this->asImp_().assembleJacobianAndResidualImpl(jac, gridVariables); // forward to the internal implementation
    }

    /*!
     * \brief Assemble the residual only
     */
    void assembleResidual(SolutionVector& res)
    {
        this->asImp_().bindLocalViews();
        const auto globalI = this->assembler().gridGeometry().elementMapper().index(this->element());
        res[globalI] = this->asImp_().evalLocalResidual()[0]; // forward to the internal implementation
    }
};

/*!
 * \ingroup Assembly
 * \ingroup CCDiscretization
 * \brief An assembler for Jacobian and residual contribution per element (cell-centered methods)
 * \tparam TypeTag The TypeTag
 * \tparam diffMethod The differentiation method to residual compute derivatives
 * \tparam implicit Specifies whether the time discretization is implicit or not not (i.e. explicit)
 */
template<class TypeTag, class Assembler, DiffMethod diffMethod = DiffMethod::numeric, bool implicit = true>
class CCLocalAssembler;

/*!
 * \ingroup Assembly
 * \ingroup CCDiscretization
 * \brief Cell-centered scheme local assembler using numeric differentiation and implicit time discretization
 */
template<class TypeTag, class Assembler>
class CCLocalAssembler<TypeTag, Assembler, DiffMethod::numeric, /*implicit=*/true>
: public CCLocalAssemblerBase<TypeTag, Assembler,
                              CCLocalAssembler<TypeTag, Assembler, DiffMethod::numeric, true>, true >
{
    using ThisType = CCLocalAssembler<TypeTag, Assembler, DiffMethod::numeric, true>;
    using ParentType = CCLocalAssemblerBase<TypeTag, Assembler, ThisType, true>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using NumEqVector = GetPropType<TypeTag, Properties::NumEqVector>;
    using Element = typename GetPropType<TypeTag, Properties::GridGeometry>::GridView::template Codim<0>::Entity;
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using FVElementGeometry = typename GridGeometry::LocalView;
    using GridVariables = GetPropType<TypeTag, Properties::GridVariables>;
    using JacobianMatrix = GetPropType<TypeTag, Properties::JacobianMatrix>;
    using Problem = typename GridVariables::GridVolumeVariables::Problem;

    enum { numEq = GetPropType<TypeTag, Properties::ModelTraits>::numEq() };
    enum { dim = GetPropType<TypeTag, Properties::GridGeometry>::GridView::dimension };

    using FluxStencil = Dumux::FluxStencil<FVElementGeometry>;
    static constexpr int maxElementStencilSize = GridGeometry::maxElementStencilSize;
    static constexpr bool enableGridFluxVarsCache = getPropValue<TypeTag, Properties::EnableGridFluxVariablesCache>();

public:

    using ParentType::ParentType;

    /*!
     * \brief Computes the derivatives with respect to the given element and adds them
     *        to the global matrix.
     *
     * \return The element residual at the current solution.
     */
    NumEqVector assembleJacobianAndResidualImpl(JacobianMatrix& A, GridVariables& gridVariables)
    {
        //////////////////////////////////////////////////////////////////////////////////////////////////
        // Calculate derivatives of all dofs in stencil with respect to the dofs in the element. In the //
        // neighboring elements we do so by computing the derivatives of the fluxes which depend on the //
        // actual element. In the actual element we evaluate the derivative of the entire residual.     //
        //////////////////////////////////////////////////////////////////////////////////////////////////

        // get some aliases for convenience
        const auto& element = this->element();
        const auto& fvGeometry = this->fvGeometry();
        const auto& gridGeometry = this->assembler().gridGeometry();
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
        using Residuals = ReservedBlockVector<NumEqVector, maxElementStencilSize>;
        Residuals origResiduals(numNeighbors + 1); origResiduals = 0.0;
        origResiduals[0] = this->evalLocalResidual()[0];

        // lambda for convenient evaluation of the fluxes across scvfs in the neighbors
        // if the neighbor is a ghost we don't want to add anything to their residual
        // so we return 0 and omit computing the flux
        auto evalNeighborFlux = [&] (const auto& neighbor, const auto& scvf)
        {
            if (neighbor.partitionType() == Dune::GhostEntity)
                return NumEqVector(0.0);
            else
                return this->localResidual().evalFlux(this->problem(),
                                                      neighbor,
                                                      this->fvGeometry(),
                                                      this->curElemVolVars(),
                                                      this->elemFluxVarsCache(), scvf);
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
        const auto& curSol = this->curSol();
        const auto origPriVars = curSol[globalI];
        const auto origVolVars = curVolVars;

        // element solution container to be deflected
        auto elemSol = elementSolution(element, curSol, gridGeometry);

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
                curVolVars.update(elemSol, this->problem(), element, scv);
                elemFluxVarsCache.update(element, fvGeometry, curElemVolVars);
                if (enableGridFluxVarsCache)
                    gridVariables.gridFluxVarsCache().updateElement(element, fvGeometry, curElemVolVars);

                // calculate the residual with the deflected primary variables
                partialDerivsTmp[0] = this->evalLocalResidual()[0];

                // calculate the fluxes in the neighbors with the deflected primary variables
                for (std::size_t k = 0; k < numNeighbors; ++k)
                    for (auto scvfIdx : connectivityMap[globalI][k].scvfsJ)
                        partialDerivsTmp[k+1] += evalNeighborFlux(neighborElements[k], fvGeometry.scvf(scvfIdx));

                return partialDerivsTmp;
            };

            // derive the residuals numerically
            static const NumericEpsilon<Scalar, numEq> eps_{this->problem().paramGroup()};
            static const int numDiffMethod = getParamFromGroup<int>(this->problem().paramGroup(), "Assembly.NumericDifferenceMethod");
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
            // no special treatment is needed if globalJ is a ghost because then derivatives have been assembled to 0 above
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
        }

        // restore original state of the flux vars cache in case of global caching.
        // This has to be done in order to guarantee that everything is in an undeflected
        // state before the assembly of another element is called. In the case of local caching
        // this is obsolete because the elemFluxVarsCache used here goes out of scope after this.
        // We only have to do this for the last primary variable, for all others the flux var cache
        // is updated with the correct element volume variables before residual evaluations
        if (enableGridFluxVarsCache)
            gridVariables.gridFluxVarsCache().updateElement(element, fvGeometry, curElemVolVars);

        // return the original residual
        return origResiduals[0];
    }
};


/*!
 * \ingroup Assembly
 * \ingroup CCDiscretization
 * \brief Cell-centered scheme local assembler using numeric differentiation and explicit time discretization
 */
template<class TypeTag, class Assembler>
class CCLocalAssembler<TypeTag, Assembler, DiffMethod::numeric, /*implicit=*/false>
: public CCLocalAssemblerBase<TypeTag, Assembler,
            CCLocalAssembler<TypeTag, Assembler, DiffMethod::numeric, false>, false>
{
    using ThisType = CCLocalAssembler<TypeTag, Assembler, DiffMethod::numeric, false>;
    using ParentType = CCLocalAssemblerBase<TypeTag, Assembler, ThisType, false>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using NumEqVector = GetPropType<TypeTag, Properties::NumEqVector>;
    using Element = typename GetPropType<TypeTag, Properties::GridGeometry>::GridView::template Codim<0>::Entity;
    using GridVariables = GetPropType<TypeTag, Properties::GridVariables>;
    using JacobianMatrix = GetPropType<TypeTag, Properties::JacobianMatrix>;
    using Problem = typename GridVariables::GridVolumeVariables::Problem;

    enum { numEq = GetPropType<TypeTag, Properties::ModelTraits>::numEq() };

public:
    using ParentType::ParentType;

    /*!
     * \brief Computes the derivatives with respect to the given element and adds them
     *        to the global matrix.
     *
     * \return The element residual at the current solution.
     */
    NumEqVector assembleJacobianAndResidualImpl(JacobianMatrix& A, GridVariables& gridVariables)
    {
        if (this->assembler().isStationaryProblem())
            DUNE_THROW(Dune::InvalidStateException, "Using explicit jacobian assembler with stationary local residual");

        // assemble the undeflected residual
        auto residual = this->evalLocalResidual()[0];
        const auto storageResidual = this->evalLocalStorageResidual()[0];

        //////////////////////////////////////////////////////////////////////////////////////////////////
        // Calculate derivatives of all dofs in stencil with respect to the dofs in the element. In the //
        // neighboring elements all derivatives are zero. For the assembled element only the storage    //
        // derivatives are non-zero.                                                                    //
        //////////////////////////////////////////////////////////////////////////////////////////////////

        // get some aliases for convenience
        const auto& element = this->element();
        const auto& fvGeometry = this->fvGeometry();
        const auto& gridGeometry = this->assembler().gridGeometry();
        auto&& curElemVolVars = this->curElemVolVars();

        // reference to the element's scv (needed later) and corresponding vol vars
        const auto globalI = gridGeometry.elementMapper().index(element);
        const auto& scv = fvGeometry.scv(globalI);
        auto& curVolVars = ParentType::getVolVarAccess(gridVariables.curGridVolVars(), curElemVolVars, scv);

        // save a copy of the original privars and vol vars in order
        // to restore the original solution after deflection
        const auto& curSol = this->curSol();
        const auto origPriVars = curSol[globalI];
        const auto origVolVars = curVolVars;

        // element solution container to be deflected
        auto elemSol = elementSolution(element, curSol, gridGeometry);

        NumEqVector partialDeriv;

        // derivatives in the neighbors with repect to the current elements
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
                return this->evalLocalStorageResidual()[0];
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
};

/*!
 * \ingroup Assembly
 * \ingroup CCDiscretization
 * \brief Cell-centered scheme local assembler using analytic (hand-coded) differentiation and implicit time discretization
 */
template<class TypeTag, class Assembler>
class CCLocalAssembler<TypeTag, Assembler, DiffMethod::analytic, /*implicit=*/true>
: public CCLocalAssemblerBase<TypeTag, Assembler,
            CCLocalAssembler<TypeTag, Assembler, DiffMethod::analytic, true>, true>
{
    using ThisType = CCLocalAssembler<TypeTag, Assembler, DiffMethod::analytic, true>;
    using ParentType = CCLocalAssemblerBase<TypeTag, Assembler, ThisType, true>;
    using NumEqVector = GetPropType<TypeTag, Properties::NumEqVector>;
    using JacobianMatrix = GetPropType<TypeTag, Properties::JacobianMatrix>;
    using GridVariables = GetPropType<TypeTag, Properties::GridVariables>;
    using Problem = typename GridVariables::GridVolumeVariables::Problem;

    enum { numEq = GetPropType<TypeTag, Properties::ModelTraits>::numEq() };

public:
    using ParentType::ParentType;

    /*!
     * \brief Computes the derivatives with respect to the given element and adds them
     *        to the global matrix.
     *
     * \return The element residual at the current solution.
     */
    NumEqVector assembleJacobianAndResidualImpl(JacobianMatrix& A, const GridVariables& gridVariables)
    {
        // treat ghost separately, we always want zero update for ghosts
        if (this->elementIsGhost())
        {
            const auto globalI = this->assembler().gridGeometry().elementMapper().index(this->element());
            for (int pvIdx = 0; pvIdx < numEq; ++pvIdx)
                A[globalI][globalI][pvIdx][pvIdx] = 1.0;

            // return zero residual
            return NumEqVector(0.0);
        }

        // assemble the undeflected residual
        auto residual = this->evalLocalResidual()[0];

        // get some aliases for convenience
        const auto& problem = this->problem();
        const auto& element = this->element();
        const auto& fvGeometry = this->fvGeometry();
        const auto& curElemVolVars = this->curElemVolVars();
        const auto& elemFluxVarsCache = this->elemFluxVarsCache();

        // get reference to the element's current vol vars
        const auto globalI = this->assembler().gridGeometry().elementMapper().index(element);
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
        }

        // return element residual
        return residual;
    }
};

/*!
 * \ingroup Assembly
 * \ingroup CCDiscretization
 * \brief Cell-centered scheme local assembler using analytic (hand-coded) differentiation and explicit time discretization
 */
template<class TypeTag, class Assembler>
class CCLocalAssembler<TypeTag, Assembler, DiffMethod::analytic, /*implicit=*/false>
: public CCLocalAssemblerBase<TypeTag, Assembler,
            CCLocalAssembler<TypeTag, Assembler, DiffMethod::analytic, false>, false>
{
    using ThisType = CCLocalAssembler<TypeTag, Assembler, DiffMethod::analytic, false>;
    using ParentType = CCLocalAssemblerBase<TypeTag, Assembler, ThisType, false>;
    using NumEqVector = GetPropType<TypeTag, Properties::NumEqVector>;
    using JacobianMatrix = GetPropType<TypeTag, Properties::JacobianMatrix>;
    using GridVariables = GetPropType<TypeTag, Properties::GridVariables>;
    using Problem = typename GridVariables::GridVolumeVariables::Problem;

    enum { numEq = GetPropType<TypeTag, Properties::ModelTraits>::numEq() };

public:
    using ParentType::ParentType;

    /*!
     * \brief Computes the derivatives with respect to the given element and adds them
     *        to the global matrix.
     *
     * \return The element residual at the current solution.
     */
    NumEqVector assembleJacobianAndResidualImpl(JacobianMatrix& A, const GridVariables& gridVariables)
    {
        // treat ghost separately, we always want zero update for ghosts
        if (this->elementIsGhost())
        {
            const auto globalI = this->assembler().gridGeometry().elementMapper().index(this->element());
            for (int pvIdx = 0; pvIdx < numEq; ++pvIdx)
                A[globalI][globalI][pvIdx][pvIdx] = 1.0;

            // return zero residual
            return NumEqVector(0.0);
        }

        // assemble the undeflected residual
        const auto residual = this->evalLocalResidual()[0];

        // get reference to the element's current vol vars
        const auto globalI = this->assembler().gridGeometry().elementMapper().index(this->element());
        const auto& scv = this->fvGeometry().scv(globalI);
        const auto& volVars = this->curElemVolVars()[scv];

        // add hand-code derivative of storage term
        this->localResidual().addStorageDerivatives(A[globalI][globalI], this->problem(), this->element(), this->fvGeometry(), volVars, scv);

        if constexpr (Problem::enableInternalDirichletConstraints())
        {
            // check if own element has internal Dirichlet constraint
            const auto internalDirichletConstraints = this->problem().hasInternalDirichletConstraint(this->element(), scv);
            const auto dirichletValues = this->problem().internalDirichlet(this->element(), scv);

            for (int pvIdx = 0; pvIdx < numEq; ++pvIdx)
            {
                for (int eqIdx = 0; eqIdx < numEq; ++eqIdx)
                {
                    if (internalDirichletConstraints[eqIdx])
                    {
                        residual[eqIdx] = volVars.priVars()[eqIdx] - dirichletValues[eqIdx];
                        A[globalI][globalI][eqIdx][pvIdx] = (eqIdx == pvIdx) ? 1.0 : 0.0;
                    }
                }
            }
        }

        // return the original residual
        return residual;
    }
};

} // end namespace Dumux

#endif
