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
    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
    using JacobianMatrix = typename GET_PROP_TYPE(TypeTag, JacobianMatrix);
    using GridVariables = typename GET_PROP_TYPE(TypeTag, GridVariables);
    using SolutionVector = typename GET_PROP_TYPE(TypeTag, SolutionVector);
    using ElementVolumeVariables = typename GET_PROP_TYPE(TypeTag, ElementVolumeVariables);
    using NumEqVector = typename GET_PROP_TYPE(TypeTag, NumEqVector);

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
        const auto globalI = this->assembler().fvGridGeometry().elementMapper().index(this->element());
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
        const auto globalI = this->assembler().fvGridGeometry().elementMapper().index(this->element());
        res[globalI] = this->asImp_().evalLocalResidual()[0]; // forward to the internal implementation
    }
};

/*!
 * \ingroup Assembly
 * \ingroup CCDiscretization
 * \brief An assembler for Jacobian and residual contribution per element (cell-centered methods)
 * \tparam TypeTag The TypeTag
 * \tparam DM The differentiation method to residual compute derivatives
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
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using NumEqVector = typename GET_PROP_TYPE(TypeTag, NumEqVector);
    using Element = typename GET_PROP_TYPE(TypeTag, GridView)::template Codim<0>::Entity;
    using FVGridGeometry = typename GET_PROP_TYPE(TypeTag, FVGridGeometry);
    using FVElementGeometry = typename FVGridGeometry::LocalView;
    using GridVariables = typename GET_PROP_TYPE(TypeTag, GridVariables);
    using JacobianMatrix = typename GET_PROP_TYPE(TypeTag, JacobianMatrix);

    enum { numEq = GET_PROP_TYPE(TypeTag, ModelTraits)::numEq() };
    enum { dim = GET_PROP_TYPE(TypeTag, GridView)::dimension };

    using FluxStencil = Dumux::FluxStencil<FVElementGeometry>;
    static constexpr int maxElementStencilSize = FVGridGeometry::maxElementStencilSize;
    static constexpr bool enableGridFluxVarsCache = GET_PROP_VALUE(TypeTag, EnableGridFluxVariablesCache);

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
        const auto& fvGridGeometry = this->assembler().fvGridGeometry();
        auto&& curElemVolVars = this->curElemVolVars();
        auto&& elemFluxVarsCache = this->elemFluxVarsCache();

        // get stencil informations
        const auto globalI = fvGridGeometry.elementMapper().index(element);
        const auto& connectivityMap = fvGridGeometry.connectivityMap();
        const auto numNeighbors = connectivityMap[globalI].size();

        // container to store the neighboring elements
        Dune::ReservedVector<Element, maxElementStencilSize> neighborElements;
        neighborElements.resize(numNeighbors);

        // assemble the undeflected residual
        using Residuals = ReservedBlockVector<NumEqVector, maxElementStencilSize>;
        Residuals origResiduals(numNeighbors + 1); origResiduals = 0.0;
        origResiduals[0] = this->evalLocalResidual()[0];

        // lambda for convenient evaluation of the fluxes across scvfs in the neighbors
        auto evalNeighborFlux = [&] (const auto& neighbor, const auto& scvf)
        {
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
            neighborElements[j-1] = fvGridGeometry.element(dataJ.globalJ);
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
        auto elemSol = elementSolution(element, curSol, fvGridGeometry);

        // derivatives in the neighbors with repect to the current elements
        // in index 0 we save the derivative of the element residual with respect to it's own dofs
        Residuals partialDerivs(numNeighbors + 1);

        for (int pvIdx = 0; pvIdx < numEq; ++pvIdx)
        {

            // for ghost elements we assemble a 1.0 where the primary variable and zero everywhere else
            // as we always solve for a delta of the solution with repect to the initial solution this
            // results in a delta of zero for ghosts, we still need to do the neighbor derivatives though
            // so we are not done yet here.
            partialDerivs = 0.0;
            if (this->elementIsGhost()) partialDerivs[0][pvIdx] = 1.0;

            auto evalResiduals = [&](Scalar priVar)
            {
                Residuals partialDerivsTmp(numNeighbors + 1);
                partialDerivsTmp = 0.0;
                // update the volume variables and the flux var cache
                elemSol[0][pvIdx] = priVar;
                curVolVars.update(elemSol, this->problem(), element, scv);
                if (enableGridFluxVarsCache)
                    gridVariables.gridFluxVarsCache().updateElement(element, fvGeometry, curElemVolVars);
                else
                    elemFluxVarsCache.update(element, fvGeometry, curElemVolVars);

                // calculate the residual with the deflected primary variables
                if (!this->elementIsGhost()) partialDerivsTmp[0] = this->evalLocalResidual()[0];

                // calculate the fluxes in the neighbors with the deflected primary variables
                for (std::size_t k = 0; k < numNeighbors; ++k)
                    for (auto scvfIdx : connectivityMap[globalI][k].scvfsJ)
                        partialDerivsTmp[k+1] += evalNeighborFlux(neighborElements[k], fvGeometry.scvf(scvfIdx));

                return partialDerivsTmp;
            };

            // derive the residuals numerically
            static const NumericEpsilon<Scalar, numEq> eps_{GET_PROP_VALUE(TypeTag, ModelParameterGroup)};
            static const int numDiffMethod = getParamFromGroup<int>(GET_PROP_VALUE(TypeTag, ModelParameterGroup), "Assembly.NumericDifferenceMethod");
            NumericDifferentiation::partialDerivative(evalResiduals, elemSol[0][pvIdx], partialDerivs, origResiduals,
                                                      eps_(elemSol[0][pvIdx], pvIdx), numDiffMethod);

            // add the current partial derivatives to the global jacobian matrix
            for (int eqIdx = 0; eqIdx < numEq; eqIdx++)
            {
                // the diagonal entries
                A[globalI][globalI][eqIdx][pvIdx] += partialDerivs[0][eqIdx];

                // off-diagonal entries
                j = 1;
                for (const auto& dataJ : connectivityMap[globalI])
                    A[dataJ.globalJ][globalI][eqIdx][pvIdx] += partialDerivs[j++][eqIdx];
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
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using NumEqVector = typename GET_PROP_TYPE(TypeTag, NumEqVector);
    using Element = typename GET_PROP_TYPE(TypeTag, GridView)::template Codim<0>::Entity;
    using GridVariables = typename GET_PROP_TYPE(TypeTag, GridVariables);
    using JacobianMatrix = typename GET_PROP_TYPE(TypeTag, JacobianMatrix);

    enum { numEq = GET_PROP_TYPE(TypeTag, ModelTraits)::numEq() };

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
        const auto residual = this->evalLocalResidual()[0];

        //////////////////////////////////////////////////////////////////////////////////////////////////
        // Calculate derivatives of all dofs in stencil with respect to the dofs in the element. In the //
        // neighboring elements all derivatives are zero. For the assembled element only the storage    //
        // derivatives are non-zero.                                                                    //
        //////////////////////////////////////////////////////////////////////////////////////////////////

        // get some aliases for convenience
        const auto& element = this->element();
        const auto& fvGeometry = this->fvGeometry();
        const auto& fvGridGeometry = this->assembler().fvGridGeometry();
        auto&& curElemVolVars = this->curElemVolVars();

        // reference to the element's scv (needed later) and corresponding vol vars
        const auto globalI = fvGridGeometry.elementMapper().index(element);
        const auto& scv = fvGeometry.scv(globalI);
        auto& curVolVars = ParentType::getVolVarAccess(gridVariables.curGridVolVars(), curElemVolVars, scv);

        // save a copy of the original privars and vol vars in order
        // to restore the original solution after deflection
        const auto& curSol = this->curSol();
        const auto origPriVars = curSol[globalI];
        const auto origVolVars = curVolVars;

        // element solution container to be deflected
        auto elemSol = elementSolution(element, curSol, fvGridGeometry);

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
                static const NumericEpsilon<Scalar, numEq> eps_{GET_PROP_VALUE(TypeTag, ModelParameterGroup)};
                static const int numDiffMethod = getParamFromGroup<int>(GET_PROP_VALUE(TypeTag, ModelParameterGroup), "Assembly.NumericDifferenceMethod");
                NumericDifferentiation::partialDerivative(evalStorage, elemSol[0][pvIdx], partialDeriv, residual,
                                                          eps_(elemSol[0][pvIdx], pvIdx), numDiffMethod);
            }

            // for ghost elements we assemble a 1.0 where the primary variable and zero everywhere else
            // as we always solve for a delta of the solution with repect to the initial solution this
            // results in a delta of zero for ghosts
            else partialDeriv[pvIdx] = 1.0;

            // add the current partial derivatives to the global jacobian matrix
            // only diagonal entries for explicit jacobians
            for (int eqIdx = 0; eqIdx < numEq; eqIdx++)
                A[globalI][globalI][eqIdx][pvIdx] += partialDeriv[eqIdx];

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
    using NumEqVector = typename GET_PROP_TYPE(TypeTag, NumEqVector);
    using JacobianMatrix = typename GET_PROP_TYPE(TypeTag, JacobianMatrix);
    using GridVariables = typename GET_PROP_TYPE(TypeTag, GridVariables);

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
        // assemble the undeflected residual
        const auto residual = this->evalLocalResidual()[0];

        // get some aliases for convenience
        const auto& problem = this->problem();
        const auto& element = this->element();
        const auto& fvGeometry = this->fvGeometry();
        const auto& curElemVolVars = this->curElemVolVars();
        const auto& elemFluxVarsCache = this->elemFluxVarsCache();

        // get reference to the element's current vol vars
        const auto globalI = this->assembler().fvGridGeometry().elementMapper().index(element);
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
    using NumEqVector = typename GET_PROP_TYPE(TypeTag, NumEqVector);
    using JacobianMatrix = typename GET_PROP_TYPE(TypeTag, JacobianMatrix);
    using GridVariables = typename GET_PROP_TYPE(TypeTag, GridVariables);

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
        // assemble the undeflected residual
        const auto residual = this->evalLocalResidual()[0];

        // get reference to the element's current vol vars
        const auto globalI = this->assembler().fvGridGeometry().elementMapper().index(this->element());
        const auto& scv = this->fvGeometry().scv(globalI);
        const auto& volVars = this->curElemVolVars()[scv];

        // add hand-code derivative of storage term
        this->localResidual().addStorageDerivatives(A[globalI][globalI], this->problem(), this->element(), this->fvGeometry(), volVars, scv);

        // return the original residual
        return residual;
    }
};

} // end namespace Dumux

#endif
