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
 * \brief An assembler for the global linear system for fully implicit models
 *        and cell-centered discretization schemes using Newton's method.
 */
#ifndef DUMUX_CC_IMPLICIT_NEWTON_ASSEMBLER_HH
#define DUMUX_CC_IMPLICIT_NEWTON_ASSEMBLER_HH

#include <dune/istl/matrixindexset.hh>

#include <dumux/implicit/properties.hh>

namespace Dumux {

/*!
 * \ingroup ImplicitModel
 * \brief An assembler for the global linear system for fully implicit models
 *        and cell-centered discretization schemes using Newton's method.
 */
template<class TypeTag>
class CCImplicitNewtonAssembler
{
    using NumEqVector = typename GET_PROP_TYPE(TypeTag, NumEqVector);
    using ElementBoundaryTypes = typename GET_PROP_TYPE(TypeTag, ElementBoundaryTypes);

public:
    template<class Assembler>
    static void assemble(Assembler& assembler)
    {
        // get some references for convenience
        const auto& problem = assembler.problem_();
        const auto& gridFvGeometry = assembler.gridFvGeometry_();
        const auto& assemblyMap = assembler.assemblyMap_();
        auto& localResidual = assembler.localResidual_();
        auto& gridVariables = assembler.gridVariables_();
        auto& curSol = gridVariables.curSol();
        auto& prevSol = gridVariables.prevSol();
        auto& A = assembler.matrix_();
        auto& rhs = assembler.rhs_();

        // loop over elements and insert element contributions into global system
        for (const auto element : elements(gridFvGeometry.gridView()))
        {
            const bool isGhost = (element.partitionType() == Dune::GhostEntity);

            // prepare the local views
            auto fvGeometry = localView(globalFvGeometry());
            fvGeometry.bind(element);

            auto curElemVolVars = localView(gridVariables.curGlobalVolVars());
            curElemVolVars.bind(element, fvGeometry, curSol);

            auto prevElemVolVars = localView(gridVariables.prevGlobalVolVars());
            prevElemVolVars.bindElement(element, fvGeometry, prevSol);

            auto elemFluxVarsCache = localView(gridVariables.globalFluxVarsCache());
            elemFluxVarsCache.bind(element, fvGeometry, curElemVolVars);

            // the global dof of the actual element
            const auto globalI = gridFvGeometry.elementMapper().index(element);

            // check for boundaries on the element
            ElementBoundaryTypes elemBcTypes;
            elemBcTypes.update(problem, element, fvGeometry);

            // the actual element's current residual
            const NumEqVector residual = isGhost ? 0.0 : localResidual.eval(element,
                                                                            fvGeometry,
                                                                            prevElemVolVars,
                                                                            curElemVolVars,
                                                                            elemBcTypes,
                                                                            elemFluxVarsCache);

            // store this also in the global container
            rhs[globalI] = residual;

            // TODO Do we really need this??????????
            // this->model_().updatePVWeights(fvGeometry);

            //////////////////////////////////////////////////////////////////////////////////////////////////
            //                                                                                              //
            // Calculate derivatives of all dofs in stencil with respect to the dofs in the element. In the //
            // neighboring elements we do so by computing the derivatives of the fluxes which depend on the //
            // actual element. In the actual element we evaluate the derivative of the entire residual.     //
            //                                                                                              //
            //////////////////////////////////////////////////////////////////////////////////////////////////

            static const int numericDifferenceMethod = GET_PARAM_FROM_GROUP(TypeTag, int, Implicit, NumericDifferenceMethod);

            // get stencil informations
            const auto numNeighbors = assemblyMap[globalI].size();

            // container to store the neighboring elements
            std::vector<Element> neighborElements;
            neighborElements.reserve(numNeighbors);

            // get the elements in which we need to evaluate the fluxes
            // and calculate these in the undeflected state
            Dune::BlockVector<NumEqVector> origFlux(numNeighbors);
            origFlux = 0.0;
            unsigned int j = 0;
            for (const auto& dataJ : assemblyMap[globalI])
            {
                neighborElements.emplace_back(gridFvGeometry().element(dataJ.globalJ));
                for (const auto scvfIdx : dataJ.scvfsJ)
                    origFlux[j] += localResidual.evalFlux_(neighborElements.back(),
                                                           fvGeometry,
                                                           curElemVolVars,
                                                           fvGeometry.scvf(scvfIdx),
                                                           elemFluxVarsCache);
                // increment neighbor counter
                ++j;
            }

            // reference to the element's scv (needed later) and corresponding vol vars
            const auto& scv = fvGeometry.scv(globalI);
            auto& curVolVars = curElemVolVars[scv];

            // save a copy of the original privars and vol vars in order
            // to restore the original solution after deflection
            const auto origPriVars = curSol[globalI];
            const auto origVolVars = curVolVars;

            // derivatives in the neighbors with repect to the current elements
            Dune::BlockVector<NumEqVector> neighborDeriv(numNeighbors);
            for (int pvIdx = 0; pvIdx < numEq; pvIdx++)
            {
                // reset derivatives of element dof with respect to itself
                // as well as neighbor derivatives
                NumEqVector partialDeriv(0.0);
                neighborDeriv = 0.0;

                if (isGhost)
                    partialDeriv[pvIdx] = 1.0;

                Scalar eps = assembler.numericEpsilon(curVolVars.priVar(pvIdx));
                Scalar delta = 0;

                if (numericDifferenceMethod >= 0)
                {
                    // we are not using backward differences, i.e. we need to
                    // calculate f(x + \epsilon)

                    // deflect primary variables
                    curSol[globalI][pvIdx] += eps;
                    delta += eps;

                    // update the volume variables and the flux var cache
                    curVolVars.update(gridVariables.elementSolution(element, curSol), problem, element, scv);
                    elemFluxVarsCache.update(element, fvGeometry, curElemVolVars);

                    // calculate the residual with the deflected primary variables
                    if (!isGhost)
                        partialDeriv = localResidual.eval(element,
                                                          fvGeometry,
                                                          prevElemVolVars,
                                                          curElemVolVars,
                                                          elemBcTypes,
                                                          elemFluxVarsCache);

                    // calculate the fluxes in the neighbors with the deflected primary variables
                    for (std::size_t k = 0; k < numNeighbors; ++k)
                        for (auto scvfIdx : assemblyMap[globalI][k].scvfsJ)
                            neighborDeriv[k] += localResidual.evalFlux_(neighborElements[k],
                                                                        fvGeometry,
                                                                        curElemVolVars,
                                                                        fvGeometry.scvf(scvfIdx),
                                                                        elemFluxVarsCache);
                }
                else
                {
                    // we are using backward differences, i.e. we don't need
                    // to calculate f(x + \epsilon) and we can recycle the
                    // (already calculated) residual f(x)
                    if (!isGhost)
                        partialDeriv = residual;
                    neighborDeriv = origFlux;
                }

                if (numericDifferenceMethod_ <= 0)
                {
                    // we are not using forward differences, i.e. we
                    // need to calculate f(x - \epsilon)

                    // deflect the primary variables
                    curSol[globalI][pvIdx] -= delta + eps;
                    delta += eps;

                    // update the volume variables and the flux var cache
                    curVolVars.update(gridVariables.elementSolution(element, curSol), problem, element, scv);
                    elemFluxVarsCache.update(element, fvGeometry, curElemVolVars);

                    // calculate the residual with the deflected primary variables and subtract it
                    if (!isGhost)
                        partialDeriv -= localResidual.eval(element,
                                                           fvGeometry,
                                                           prevElemVolVars,
                                                           curElemVolVars,
                                                           elemBcTypes,
                                                           elemFluxVarsCache);

                    // calculate the fluxes into element with the deflected primary variables
                    for (std::size_t k = 0; k < numNeighbors; ++k)
                        for (auto scvfIdx : assemblyMap_[globalI][k].scvfsJ)
                            neighborDeriv[k] -= this->localResidual().evalFlux_(neighborElements[k],
                                                                                fvGeometry,
                                                                                curElemVolVars,
                                                                                fvGeometry.scvf(scvfIdx),
                                                                                elemFluxVarsCache);
                }
                else
                {
                    // we are using forward differences, i.e. we don't need to
                    // calculate f(x - \epsilon) and we can recycle the
                    // (already calculated) residual f(x)
                    if (!isGhost)
                        partialDeriv -= residual;
                    neighborDeriv -= origFlux;
                }

                // divide difference in residuals by the magnitude of the
                // deflections between the two function evaluation
                if (!isGhost)
                    partialDeriv /= delta;
                neighborDeriv /= delta;

                // restore the original state of the scv's volume variables
                curVolVars = origVolVars;

                // restore the current element solution
                curSol[globalI_] = origPriVars;

                // add the current partial derivatives to the global jacobian matrix
                for (int eqIdx = 0; eqIdx < numEq; eqIdx++)
                {
                    // the diagonal entries
                    A[globalI][globalI][eqIdx][pvIdx] += partialDeriv[eqIdx];

                    // off-diagonal entries
                    j = 0;
                    for (const auto& dataJ : assemblyMap_[globalI])
                        A[dataJ.globalJ][globalI][eqIdx][pvIdx] += neighborDeriv[j++][pvIdx];
                }
            }

            //////////////////////////////////////////////////////////////////////////////////////////////
            //                                                                                          //
            // Calculate derivatives of the dofs in the element with respect to user-defined additional //
            // dof dependencies. We do so by evaluating the change in the source term of the current    //
            // element with respect to the primary variables at the given additional dofs.              //
            //                                                                                          //
            //////////////////////////////////////////////////////////////////////////////////////////////

            const auto& additionalDofDepedencies = problem.getAdditionalDofDependencies(globalI);
            if (!additionalDofDepedencies.empty() && !isGhost)
            {
                // compute the source in the undeflected state
                auto source = localResidual.computeSource(element, fvGeometry, curElemVolVars, scv);
                source *= -scv.volume()*curVolVarsI.extrusionFactor();

                // deflect solution at given dofs and recalculate the source
                for (auto globalJ : additionalDofDependencies)
                {
                    const auto& scvJ = fvGeometry.scv(globalJ);
                    auto& curVolVarsJ = curElemVolVars[scv];
                    const auto& elementJ = gridFvGeometry.element(globalJ);

                    // save a copy of the original privars and volvars
                    // to restore original solution after deflection
                    const auto origPriVars = curSol[globalJ];
                    const auto origVolVarsJ = curVolVarsJ;

                    // derivatives with repect to the additional DOF we depend on
                    for (int pvIdx = 0; pvIdx < numEq; pvIdx++)
                    {
                        // derivatives of element dof with respect to itself
                        NumEqVector partialDeriv(0.0);
                        const auto eps = assembler.numericEpsilon(curVolVarsJ.priVar(pvIdx));
                        Scalar delta = 0;

                        if (numericDifferenceMethod_ >= 0)
                        {
                            // we are not using backward differences, i.e. we need to
                            // calculate f(x + \epsilon)

                            // deflect primary variables
                            curSol[globalJ][pvIdx] += eps;
                            delta += eps;

                            // update the volume variables and the flux var cache
                            curVolVarsJ.update(gridVariables.elementSolution(elementJ, curSol), problem, elementJ, scvJ);

                            // calculate the source with the deflected primary variables
                            auto deflSource = localResidual.computeSource(element, fvGeometry, curElemVolVars, scv);
                            deflSource *= -scv.volume()*curVolVarsI.extrusionFactor();
                            partialDeriv = std::move(deflSource);
                        }
                        else
                        {
                            // we are using backward differences, i.e. we don't need
                            // to calculate f(x + \epsilon) and we can recycle the
                            // (already calculated) source f(x)
                            partialDeriv = source;
                        }

                        if (numericDifferenceMethod_ <= 0)
                        {
                            // we are not using forward differences, i.e. we
                            // need to calculate f(x - \epsilon)

                            // deflect the primary variables
                            curSol[globalJ][pvIdx] -= delta + eps;
                            delta += eps;

                            // update the volume variables and the flux var cache
                            curVolVarsJ.update(gridVariables.elementSolution(elementJ, curSol), problem, elementJ, scvJ);

                            // calculate the source with the deflected primary variables and subtract
                            auto deflSource = localResidual.computeSource(element, fvGeometry, curElemVolVars, scv);
                            deflSource *= -scv.volume()*curVolVarsI.extrusionFactor();
                            partialDeriv -= std::move(deflSource);
                        }
                        else
                        {
                            // we are using forward differences, i.e. we don't need to
                            // calculate f(x - \epsilon) and we can recycle the
                            // (already calculated) source f(x)
                            partialDeriv -= source;
                        }

                        // divide difference in residuals by the magnitude of the
                        // deflections between the two function evaluation
                        partialDeriv /= delta;

                        // restore the original state of the dofs privars and the volume variables
                        curSol[globalJ] = origPriVars;
                        curVolVarsJ = origVolVarsJ;

                        // add the current partial derivatives to the global jacobian matrix
                        for (int eqIdx = 0; eqIdx < numEq; eqIdx++)
                            A[globalI][globalJ][eqIdx][pvIdx] += partialDeriv[eqIdx];
                    }
                }
            }
        }
    }
};

}