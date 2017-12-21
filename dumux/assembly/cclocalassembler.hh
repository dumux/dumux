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
 * \tparam TypeTag the TypeTag
 * \tparam DM the differentiation method to residual compute derivatives
 * \tparam implicit if to use an implicit or explicit time discretization
 */
#ifndef DUMUX_CC_LOCAL_ASSEMBLER_HH
#define DUMUX_CC_LOCAL_ASSEMBLER_HH

#include <dune/grid/common/gridenums.hh> // for GhostEntity
#include <dune/istl/matrixindexset.hh>
#include <dune/istl/bvector.hh>

#include <dumux/common/properties.hh>
#include <dumux/common/parameters.hh>
#include <dumux/assembly/diffmethod.hh>

namespace Dumux {

/*!
 * \ingroup Assembly
 * \ingroup CCDiscretization
 * \brief An assembler for Jacobian and residual contribution per element (cell-centered methods)
 * \tparam TypeTag the TypeTag
 * \tparam DM the differentiation method to residual compute derivatives
 * \tparam implicit if to use an implicit or explicit time discretization
 */
template<class TypeTag,
         DiffMethod DM = DiffMethod::numeric,
         bool implicit = true>
class CCLocalAssembler;

/*!
 * \ingroup Assembly
 * \ingroup CCDiscretization
 * \brief Cell-centered scheme local assembler using numeric differentiation and implicit time discretization
 */
template<class TypeTag>
class CCLocalAssembler<TypeTag,
                       DiffMethod::numeric,
                       /*implicit=*/true>
{
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using NumEqVector = typename GET_PROP_TYPE(TypeTag, NumEqVector);
    using ElementBoundaryTypes = typename GET_PROP_TYPE(TypeTag, ElementBoundaryTypes);
    using Element = typename GET_PROP_TYPE(TypeTag, GridView)::template Codim<0>::Entity;
    using SolutionVector = typename GET_PROP_TYPE(TypeTag, SolutionVector);
    using ElementSolutionVector = typename GET_PROP_TYPE(TypeTag, ElementSolutionVector);
    using ElementVolumeVariables = typename GET_PROP_TYPE(TypeTag, ElementVolumeVariables);
    using GridVolumeVariables = typename GET_PROP_TYPE(TypeTag, GridVolumeVariables);
    using VolumeVariables = typename GET_PROP_TYPE(TypeTag, VolumeVariables);
    using SubControlVolume = typename GET_PROP_TYPE(TypeTag, SubControlVolume);
    using JacobianMatrix = typename GET_PROP_TYPE(TypeTag, JacobianMatrix);

    enum { numEq = GET_PROP_VALUE(TypeTag, NumEq) };

    static constexpr bool enableGridFluxVarsCache = GET_PROP_VALUE(TypeTag, EnableGridFluxVariablesCache);

public:

    /*!
     * \brief Computes the derivatives with respect to the given element and adds them
     *        to the global matrix. The element residual is written into the right hand side.
     */
    template<class Assembler>
    static void assemble(Assembler& assembler, JacobianMatrix& jac, SolutionVector& res,
                         const Element& element, const SolutionVector& curSol)
    {
        const auto globalI = assembler.fvGridGeometry().elementMapper().index(element);
        res[globalI] = assemble_(assembler, jac, element, curSol);
    }

    /*!
     * \brief Computes the derivatives with respect to the given element and adds them
     *        to the global matrix.
     */
    template<class Assembler>
    static void assemble(Assembler& assembler, JacobianMatrix& jac,
                         const Element& element, const SolutionVector& curSol)
    {
        assemble_(assembler, jac, element, curSol);
    }

    /*!
     * \brief Assemble the residual only
     */
    template<class Assembler>
    static void assemble(Assembler& assembler, SolutionVector& res,
                         const Element& element, const SolutionVector& curSol)
    {
        const auto globalI = assembler.fvGridGeometry().elementMapper().index(element);
        res[globalI] = assemble_(assembler, element, curSol);
    }

    /*!
     * \brief Computes the epsilon used for numeric differentiation
     *        for a given value of a primary variable.
     *
     * \param priVar The value of the primary variable
     */
    static Scalar numericEpsilon(const Scalar priVar)
    {
        // define the base epsilon as the geometric mean of 1 and the
        // resolution of the scalar type. E.g. for standard 64 bit
        // floating point values, the resolution is about 10^-16 and
        // the base epsilon is thus approximately 10^-8.
        /*
        static const Scalar baseEps
            = Dumux::geometricMean<Scalar>(std::numeric_limits<Scalar>::epsilon(), 1.0);
        */
        static const Scalar baseEps = 1e-10;
        assert(std::numeric_limits<Scalar>::epsilon()*1e4 < baseEps);
        // the epsilon value used for the numeric differentiation is
        // now scaled by the absolute value of the primary variable...
        return baseEps*(std::abs(priVar) + 1.0);
    }

private:
    /*!
     * \brief Computes the residual
     *
     * \return The element residual at the current solution.
     */
    template<class Assembler>
    static NumEqVector assemble_(Assembler& assembler,
                                 const Element& element, const SolutionVector& curSol)
    {
        // is the actual element a ghost element?
        const bool isGhost = (element.partitionType() == Dune::GhostEntity);
        if (isGhost) return NumEqVector(0.0);

        // get some references for convenience
        const auto& problem = assembler.problem();
        auto& localResidual = assembler.localResidual();
        auto& gridVariables = assembler.gridVariables();

        // prepare the local views
        auto fvGeometry = localView(assembler.fvGridGeometry());
        fvGeometry.bind(element);

        auto curElemVolVars = localView(gridVariables.curGridVolVars());
        curElemVolVars.bind(element, fvGeometry, curSol);

        auto elemFluxVarsCache = localView(gridVariables.gridFluxVarsCache());
        elemFluxVarsCache.bind(element, fvGeometry, curElemVolVars);

        auto prevElemVolVars = localView(gridVariables.prevGridVolVars());

        // for compatibility with box models
        ElementBoundaryTypes elemBcTypes;

        // the actual element's current residual
        NumEqVector residual(0.0);
        if (localResidual.isStationary())
        {
            residual = localResidual.eval(problem,
                                          element,
                                          fvGeometry,
                                          curElemVolVars,
                                          elemBcTypes,
                                          elemFluxVarsCache)[0];
        }
        else
        {
            prevElemVolVars.bindElement(element, fvGeometry, localResidual.prevSol());
            residual = localResidual.eval(problem,
                                          element,
                                          fvGeometry,
                                          prevElemVolVars,
                                          curElemVolVars,
                                          elemBcTypes,
                                          elemFluxVarsCache)[0];
        }

        return residual;
    }

    /*!
     * \brief Computes the derivatives with respect to the given element and adds them
     *        to the global matrix.
     *
     * \return The element residual at the current solution.
     */
    template<class Assembler>
    static NumEqVector assemble_(Assembler& assembler, JacobianMatrix& A,
                                 const Element& element, const SolutionVector& curSol)
    {
        // get some references for convenience
        const auto& problem = assembler.problem();
        const auto& fvGridGeometry = assembler.fvGridGeometry();
        const auto& connectivityMap = fvGridGeometry.connectivityMap();
        auto& localResidual = assembler.localResidual();
        auto& gridVariables = assembler.gridVariables();

        // prepare the local views
        auto fvGeometry = localView(assembler.fvGridGeometry());
        fvGeometry.bind(element);

        auto curElemVolVars = localView(gridVariables.curGridVolVars());
        curElemVolVars.bind(element, fvGeometry, curSol);

        auto elemFluxVarsCache = localView(gridVariables.gridFluxVarsCache());
        elemFluxVarsCache.bind(element, fvGeometry, curElemVolVars);

        auto prevElemVolVars = localView(gridVariables.prevGridVolVars());

        // the global dof of the actual element
        const auto globalI = fvGridGeometry.elementMapper().index(element);

        // check for boundaries on the element
        // TODO Do we need them for cell-centered models?
        ElementBoundaryTypes elemBcTypes;
        elemBcTypes.update(problem, element, fvGeometry);

        // is the actual element a ghost element?
        const bool isGhost = (element.partitionType() == Dune::GhostEntity);
        // is the local residual stationary?
        const bool isStationary = localResidual.isStationary();

        // the actual element's current residual
        NumEqVector residual(0.0);
        if (!isGhost)
        {
            if (isStationary)
            {
                residual = localResidual.eval(problem,
                                              element,
                                              fvGeometry,
                                              curElemVolVars,
                                              elemBcTypes,
                                              elemFluxVarsCache)[0];
            }
            else
            {
                prevElemVolVars.bindElement(element, fvGeometry, localResidual.prevSol());
                residual = localResidual.eval(problem,
                                              element,
                                              fvGeometry,
                                              prevElemVolVars,
                                              curElemVolVars,
                                              elemBcTypes,
                                              elemFluxVarsCache)[0];
            }
        }


        // TODO Do we really need this??????????
        // this->model_().updatePVWeights(fvGeometry);

        //////////////////////////////////////////////////////////////////////////////////////////////////
        //                                                                                              //
        // Calculate derivatives of all dofs in stencil with respect to the dofs in the element. In the //
        // neighboring elements we do so by computing the derivatives of the fluxes which depend on the //
        // actual element. In the actual element we evaluate the derivative of the entire residual.     //
        //                                                                                              //
        //////////////////////////////////////////////////////////////////////////////////////////////////

        static const std::string group = GET_PROP_VALUE(TypeTag, ModelParameterGroup);
        static const int numericDifferenceMethod = getParamFromGroup<int>(group, "Implicit.NumericDifferenceMethod");

        // get stencil informations
        const auto numNeighbors = connectivityMap[globalI].size();

        // container to store the neighboring elements
        std::vector<Element> neighborElements;
        neighborElements.reserve(numNeighbors);

        // get the elements in which we need to evaluate the fluxes
        // and calculate these in the undeflected state
        Dune::BlockVector<NumEqVector> origFlux(numNeighbors);
        origFlux = 0.0;
        unsigned int j = 0;
        for (const auto& dataJ : connectivityMap[globalI])
        {
            neighborElements.emplace_back(fvGridGeometry.element(dataJ.globalJ));
            for (const auto scvfIdx : dataJ.scvfsJ)
            {
                origFlux[j] += localResidual.evalFlux(problem,
                                                      neighborElements.back(),
                                                      fvGeometry,
                                                      curElemVolVars,
                                                      elemFluxVarsCache,
                                                      fvGeometry.scvf(scvfIdx));
            }
            // increment neighbor counter
            ++j;
        }

        // reference to the element's scv (needed later) and corresponding vol vars
        const auto& scv = fvGeometry.scv(globalI);
        auto& curVolVars = getVolVarAccess(gridVariables.curGridVolVars(), curElemVolVars, scv);

        // save a copy of the original privars and vol vars in order
        // to restore the original solution after deflection
        const auto origPriVars = curSol[globalI];
        const auto origVolVars = curVolVars;

        // element solution container to be deflected
        ElementSolutionVector elemSol(origPriVars);

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

            Scalar eps = numericEpsilon(curVolVars.priVar(pvIdx));
            Scalar delta = 0;

            if (numericDifferenceMethod >= 0)
            {
                // we are not using backward differences, i.e. we need to
                // calculate f(x + \epsilon)

                // deflect primary variables
                elemSol[0][pvIdx] += eps;
                delta += eps;

                // update the volume variables and the flux var cache
                curVolVars.update(elemSol, problem, element, scv);
                if (enableGridFluxVarsCache)
                    gridVariables.gridFluxVarsCache().updateElement(element, fvGeometry, curElemVolVars);
                else
                    elemFluxVarsCache.update(element, fvGeometry, curElemVolVars);

                // calculate the residual with the deflected primary variables
                if (!isGhost)
                {
                    if (isStationary)
                    {
                        partialDeriv = localResidual.eval(problem,
                                                          element,
                                                          fvGeometry,
                                                          curElemVolVars,
                                                          elemBcTypes,
                                                          elemFluxVarsCache)[0];
                    }
                    else
                    {
                        partialDeriv = localResidual.eval(problem,
                                                          element,
                                                          fvGeometry,
                                                          prevElemVolVars,
                                                          curElemVolVars,
                                                          elemBcTypes,
                                                          elemFluxVarsCache)[0];
                    }
                }

                // calculate the fluxes in the neighbors with the deflected primary variables
                for (std::size_t k = 0; k < numNeighbors; ++k)
                    for (auto scvfIdx : connectivityMap[globalI][k].scvfsJ)
                    {
                        neighborDeriv[k] += localResidual.evalFlux(problem,
                                                                   neighborElements[k],
                                                                   fvGeometry,
                                                                   curElemVolVars,
                                                                   elemFluxVarsCache,
                                                                   fvGeometry.scvf(scvfIdx));
                    }
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

            if (numericDifferenceMethod <= 0)
            {
                // we are not using forward differences, i.e. we
                // need to calculate f(x - \epsilon)

                // deflect the primary variables
                elemSol[0][pvIdx] -= delta + eps;
                delta += eps;

                // update the volume variables and the flux var cache
                curVolVars.update(elemSol, problem, element, scv);
                if (enableGridFluxVarsCache)
                    gridVariables.gridFluxVarsCache().updateElement(element, fvGeometry, curElemVolVars);
                else
                    elemFluxVarsCache.update(element, fvGeometry, curElemVolVars);

                // calculate the residual with the deflected primary variables and subtract it
                if (!isGhost)
                {
                    if (isStationary)
                    {
                        partialDeriv -= localResidual.eval(problem,
                                                           element,
                                                           fvGeometry,
                                                           curElemVolVars,
                                                           elemBcTypes,
                                                           elemFluxVarsCache)[0];
                    }
                    else
                    {
                        partialDeriv -= localResidual.eval(problem,
                                                           element,
                                                           fvGeometry,
                                                           prevElemVolVars,
                                                           curElemVolVars,
                                                           elemBcTypes,
                                                           elemFluxVarsCache)[0];
                    }
                }

                // calculate the fluxes into element with the deflected primary variables
                for (std::size_t k = 0; k < numNeighbors; ++k)
                    for (auto scvfIdx : connectivityMap[globalI][k].scvfsJ)
                    {
                        neighborDeriv[k] -= localResidual.evalFlux(problem,
                                                                   neighborElements[k],
                                                                   fvGeometry,
                                                                   curElemVolVars,
                                                                   elemFluxVarsCache,
                                                                   fvGeometry.scvf(scvfIdx));
                    }
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
            elemSol[0][pvIdx] = origPriVars[pvIdx];

            // add the current partial derivatives to the global jacobian matrix
            for (int eqIdx = 0; eqIdx < numEq; eqIdx++)
            {
                // the diagonal entries
                A[globalI][globalI][eqIdx][pvIdx] += partialDeriv[eqIdx];

                // off-diagonal entries
                j = 0;
                for (const auto& dataJ : connectivityMap[globalI])
                    A[dataJ.globalJ][globalI][eqIdx][pvIdx] += neighborDeriv[j++][eqIdx];
            }
        }

        // Restore original state of the flux vars cache in case of global caching.
        // This has to be done in order to guarantee that everything is in an undeflected
        // state before the assembly of another element is called. In the case of local caching
        // this is obsolete because the elemFluxVarsCache used here goes out of scope after this.
        if (enableGridFluxVarsCache)
            gridVariables.gridFluxVarsCache().updateElement(element, fvGeometry, curElemVolVars);

        //////////////////////////////////////////////////////////////////////////////////////////////
        //                                                                                          //
        // Calculate derivatives of the dofs in the element with respect to user-defined additional //
        // dof dependencies. We do so by evaluating the change in the source term of the current    //
        // element with respect to the primary variables at the given additional dofs.              //
        //                                                                                          //
        //////////////////////////////////////////////////////////////////////////////////////////////

        // const auto& additionalDofDepedencies = problem.getAdditionalDofDependencies(globalI);
        // if (!additionalDofDepedencies.empty() && !isGhost)
        // {
        //     // compute the source in the undeflected state
        //     auto source = localResidual.computeSource(element, fvGeometry, curElemVolVars, scv);
        //     source *= -scv.volume()*curVolVarsI.extrusionFactor();

        //     // deflect solution at given dofs and recalculate the source
        //     for (auto globalJ : additionalDofDependencies)
        //     {
        //         const auto& scvJ = fvGeometry.scv(globalJ);
        //         auto& curVolVarsJ = curElemVolVars[scv];
        //         const auto& elementJ = fvGridGeometry.element(globalJ);

        //         // save a copy of the original privars and volvars
        //         // to restore original solution after deflection
        //         const auto origPriVars = curSol[globalJ];
        //         const auto origVolVarsJ = curVolVarsJ;

        //         // derivatives with repect to the additional DOF we depend on
        //         for (int pvIdx = 0; pvIdx < numEq; pvIdx++)
        //         {
        //             // derivatives of element dof with respect to itself
        //             NumEqVector partialDeriv(0.0);
        //             const auto eps = numericEpsilon(curVolVarsJ.priVar(pvIdx));
        //             Scalar delta = 0;

        //             if (numericDifferenceMethod >= 0)
        //             {
        //                 // we are not using backward differences, i.e. we need to
        //                 // calculate f(x + \epsilon)

        //                 // deflect primary variables
        //                 curSol[globalJ][pvIdx] += eps;
        //                 delta += eps;

        //                 // update the volume variables and the flux var cache
        //                 curVolVarsJ.update(gridVariables.elementSolution(elementJ, curSol), problem, elementJ, scvJ);

        //                 // calculate the source with the deflected primary variables
        //                 auto deflSource = localResidual.computeSource(element, fvGeometry, curElemVolVars, scv);
        //                 deflSource *= -scv.volume()*curVolVarsI.extrusionFactor();
        //                 partialDeriv = std::move(deflSource);
        //             }
        //             else
        //             {
        //                 // we are using backward differences, i.e. we don't need
        //                 // to calculate f(x + \epsilon) and we can recycle the
        //                 // (already calculated) source f(x)
        //                 partialDeriv = source;
        //             }

        //             if (numericDifferenceMethod <= 0)
        //             {
        //                 // we are not using forward differences, i.e. we
        //                 // need to calculate f(x - \epsilon)

        //                 // deflect the primary variables
        //                 curSol[globalJ][pvIdx] -= delta + eps;
        //                 delta += eps;

        //                 // update the volume variables and the flux var cache
        //                 curVolVarsJ.update(gridVariables.elementSolution(elementJ, curSol), problem, elementJ, scvJ);

        //                 // calculate the source with the deflected primary variables and subtract
        //                 auto deflSource = localResidual.computeSource(element, fvGeometry, curElemVolVars, scv);
        //                 deflSource *= -scv.volume()*curVolVarsI.extrusionFactor();
        //                 partialDeriv -= std::move(deflSource);
        //             }
        //             else
        //             {
        //                 // we are using forward differences, i.e. we don't need to
        //                 // calculate f(x - \epsilon) and we can recycle the
        //                 // (already calculated) source f(x)
        //                 partialDeriv -= source;
        //             }

        //             // divide difference in residuals by the magnitude of the
        //             // deflections between the two function evaluation
        //             partialDeriv /= delta;

        //             // restore the original state of the dofs privars and the volume variables
        //             curSol[globalJ] = origPriVars;
        //             curVolVarsJ = origVolVarsJ;

        //             // add the current partial derivatives to the global jacobian matrix
        //             for (int eqIdx = 0; eqIdx < numEq; eqIdx++)
        //                 A[globalI][globalJ][eqIdx][pvIdx] += partialDeriv[eqIdx];
        //         }
        //     }
        // }

        // return the original residual
        return residual;
    }
private:
    template<class T = TypeTag>
    static typename std::enable_if<!GET_PROP_VALUE(T, EnableGridVolumeVariablesCache), VolumeVariables&>::type
    getVolVarAccess(GridVolumeVariables& gridVolVars, ElementVolumeVariables& elemVolVars, const SubControlVolume& scv)
    { return elemVolVars[scv]; }

    template<class T = TypeTag>
    static typename std::enable_if<GET_PROP_VALUE(T, EnableGridVolumeVariablesCache), VolumeVariables&>::type
    getVolVarAccess(GridVolumeVariables& gridVolVars, ElementVolumeVariables& elemVolVars, const SubControlVolume& scv)
    { return gridVolVars.volVars(scv); }
};


/*!
 * \ingroup Assembly
 * \ingroup CCDiscretization
 * \brief Cell-centered scheme local assembler using numeric differentiation and explicits time discretization
 */
template<class TypeTag>
class CCLocalAssembler<TypeTag,
                       DiffMethod::numeric,
                       /*implicit=*/false>
{
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using NumEqVector = typename GET_PROP_TYPE(TypeTag, NumEqVector);
    using ElementBoundaryTypes = typename GET_PROP_TYPE(TypeTag, ElementBoundaryTypes);
    using Element = typename GET_PROP_TYPE(TypeTag, GridView)::template Codim<0>::Entity;
    using SolutionVector = typename GET_PROP_TYPE(TypeTag, SolutionVector);
    using ElementSolutionVector = typename GET_PROP_TYPE(TypeTag, ElementSolutionVector);
    using ElementVolumeVariables = typename GET_PROP_TYPE(TypeTag, ElementVolumeVariables);
    using GridVolumeVariables = typename GET_PROP_TYPE(TypeTag, GridVolumeVariables);
    using VolumeVariables = typename GET_PROP_TYPE(TypeTag, VolumeVariables);
    using SubControlVolume = typename GET_PROP_TYPE(TypeTag, SubControlVolume);
    using JacobianMatrix = typename GET_PROP_TYPE(TypeTag, JacobianMatrix);

    enum { numEq = GET_PROP_VALUE(TypeTag, NumEq) };

public:

    /*!
     * \brief Computes the derivatives with respect to the given element and adds them
     *        to the global matrix. The element residual is written into the right hand side.
     */
    template<class Assembler>
    static void assemble(Assembler& assembler, JacobianMatrix& jac, SolutionVector& res,
                         const Element& element, const SolutionVector& curSol)
    {
        const auto globalI = assembler.fvGridGeometry().elementMapper().index(element);
        res[globalI] = assemble_(assembler, jac, element, curSol);
    }

    /*!
     * \brief Computes the derivatives with respect to the given element and adds them
     *        to the global matrix.
     */
    template<class Assembler>
    static void assemble(Assembler& assembler, JacobianMatrix& jac,
                         const Element& element, const SolutionVector& curSol)
    {
        assemble_(assembler, jac, element, curSol);
    }

    /*!
     * \brief Assemble the residual only
     */
    template<class Assembler>
    static void assemble(Assembler& assembler, SolutionVector& res,
                         const Element& element, const SolutionVector& curSol)
    {
        const auto globalI = assembler.fvGridGeometry().elementMapper().index(element);
        res[globalI] = assemble_(assembler, element, curSol);
    }

    /*!
     * \brief Computes the epsilon used for numeric differentiation
     *        for a given value of a primary variable.
     *
     * \param priVar The value of the primary variable
     */
    static Scalar numericEpsilon(const Scalar priVar)
    {
        // define the base epsilon as the geometric mean of 1 and the
        // resolution of the scalar type. E.g. for standard 64 bit
        // floating point values, the resolution is about 10^-16 and
        // the base epsilon is thus approximately 10^-8.
        /*
        static const Scalar baseEps
            = Dumux::geometricMean<Scalar>(std::numeric_limits<Scalar>::epsilon(), 1.0);
        */
        static const Scalar baseEps = 1e-10;
        assert(std::numeric_limits<Scalar>::epsilon()*1e4 < baseEps);
        // the epsilon value used for the numeric differentiation is
        // now scaled by the absolute value of the primary variable...
        return baseEps*(std::abs(priVar) + 1.0);
    }

private:

    /*!
     * \brief Computes the residual
     *
     * \return The element residual at the current solution.
     */
    template<class Assembler>
    static NumEqVector assemble_(Assembler& assembler,
                                 const Element& element, const SolutionVector& curSol)
    {
        // is the actual element a ghost element?
        const bool isGhost = (element.partitionType() == Dune::GhostEntity);
        if (isGhost) return NumEqVector(0.0);

        // get some references for convenience
        const auto& problem = assembler.problem();
        auto& localResidual = assembler.localResidual();
        auto& gridVariables = assembler.gridVariables();

        // using an explicit assembler doesn't make sense for stationary problems
        if (localResidual.isStationary())
            DUNE_THROW(Dune::InvalidStateException, "Using explicit jacobian assembler with stationary local residual");

        // prepare the local views
        auto fvGeometry = localView(assembler.fvGridGeometry());
        fvGeometry.bind(element);

        auto curElemVolVars = localView(gridVariables.curGridVolVars());
        curElemVolVars.bindElement(element, fvGeometry, curSol);

        auto prevElemVolVars = localView(gridVariables.prevGridVolVars());
        prevElemVolVars.bind(element, fvGeometry, localResidual.prevSol());

        auto elemFluxVarsCache = localView(gridVariables.gridFluxVarsCache());
        elemFluxVarsCache.bind(element, fvGeometry, prevElemVolVars);

        // compatibility with box method
        ElementBoundaryTypes elemBcTypes;

        // the actual element's previous time step residual
        auto residual = localResidual.eval(problem,
                                           element,
                                           fvGeometry,
                                           prevElemVolVars,
                                           elemBcTypes,
                                           elemFluxVarsCache)[0];

        auto storageResidual = localResidual.evalStorage(problem,
                                                         element,
                                                         fvGeometry,
                                                         prevElemVolVars,
                                                         curElemVolVars,
                                                         elemBcTypes,
                                                         elemFluxVarsCache)[0];

        residual += storageResidual;

        return residual;
    }

    /*!
     * \brief Computes the derivatives with respect to the given element and adds them
     *        to the global matrix.
     *
     * \return The element residual at the current solution.
     */
    template<class Assembler>
    static NumEqVector assemble_(Assembler& assembler, JacobianMatrix& A,
                                 const Element& element, const SolutionVector& curSol)
    {
        // get some references for convenience
        const auto& problem = assembler.problem();
        const auto& fvGridGeometry = assembler.fvGridGeometry();
        auto& localResidual = assembler.localResidual();
        auto& gridVariables = assembler.gridVariables();

        // using an explicit assembler doesn't make sense for stationary problems
        if (localResidual.isStationary())
            DUNE_THROW(Dune::InvalidStateException, "Using explicit jacobian assembler with stationary local residual");

        // prepare the local views
        auto fvGeometry = localView(assembler.fvGridGeometry());
        fvGeometry.bind(element);

        auto curElemVolVars = localView(gridVariables.curGridVolVars());
        curElemVolVars.bindElement(element, fvGeometry, curSol);

        auto prevElemVolVars = localView(gridVariables.prevGridVolVars());
        prevElemVolVars.bind(element, fvGeometry, localResidual.prevSol());

        auto elemFluxVarsCache = localView(gridVariables.gridFluxVarsCache());
        elemFluxVarsCache.bind(element, fvGeometry, prevElemVolVars);

        // the global dof of the actual element
        const auto globalI = fvGridGeometry.elementMapper().index(element);

        // check for boundaries on the element
        // TODO Do we need them for cell-centered models?
        ElementBoundaryTypes elemBcTypes;
        elemBcTypes.update(problem, element, fvGeometry);

        // is the actual element a ghost element?
        const bool isGhost = (element.partitionType() == Dune::GhostEntity);

        // the actual element's previous time step residual
        NumEqVector residual(0.0), storageResidual(0.0);
        if (!isGhost)
        {
            residual = localResidual.eval(problem,
                                          element,
                                          fvGeometry,
                                          prevElemVolVars,
                                          elemBcTypes,
                                          elemFluxVarsCache)[0];

            storageResidual = localResidual.evalStorage(problem,
                                                        element,
                                                        fvGeometry,
                                                        prevElemVolVars,
                                                        curElemVolVars,
                                                        elemBcTypes,
                                                        elemFluxVarsCache)[0];

            residual += storageResidual;
        }

        //////////////////////////////////////////////////////////////////////////////////////////////////
        // Calculate derivatives of all dofs in stencil with respect to the dofs in the element. In the //
        // neighboring elements all derivatives are zero. For the assembled element only the storage    //
        // derivatives are non-zero.                                                                    //
        //////////////////////////////////////////////////////////////////////////////////////////////////

        static const std::string group = GET_PROP_VALUE(TypeTag, ModelParameterGroup);
        static const int numericDifferenceMethod = getParamFromGroup<int>(group, "Implicit.NumericDifferenceMethod");

        // reference to the element's scv (needed later) and corresponding vol vars
        const auto& scv = fvGeometry.scv(globalI);
        auto& curVolVars = getVolVarAccess(gridVariables.curGridVolVars(), curElemVolVars, scv);

        // save a copy of the original privars and vol vars in order
        // to restore the original solution after deflection
        const auto origPriVars = curSol[globalI];
        const auto origVolVars = curVolVars;

        // element solution container to be deflected
        ElementSolutionVector elemSol({origPriVars});

        // derivatives in the neighbors with repect to the current elements
        for (int pvIdx = 0; pvIdx < numEq; pvIdx++)
        {
            // reset derivatives of element dof with respect to itself
            // as well as neighbor derivatives
            NumEqVector partialDeriv(0.0);

            if (isGhost)
                partialDeriv[pvIdx] = 1.0;

            Scalar eps = numericEpsilon(curVolVars.priVar(pvIdx));
            Scalar delta = 0;

            if (numericDifferenceMethod >= 0)
            {
                // we are not using backward differences, i.e. we need to
                // calculate f(x + \epsilon)

                // deflect primary variables
                elemSol[0][pvIdx] += eps;
                delta += eps;

                // update the volume variables and the flux var cache
                curVolVars.update(elemSol, problem, element, scv);

                // calculate the residual with the deflected primary variables
                if (!isGhost)
                {
                    partialDeriv = localResidual.evalStorage(problem,
                                                             element,
                                                             fvGeometry,
                                                             prevElemVolVars,
                                                             curElemVolVars,
                                                             elemBcTypes,
                                                             elemFluxVarsCache)[0];
                }
            }
            else
            {
                // we are using backward differences, i.e. we don't need
                // to calculate f(x + \epsilon) and we can recycle the
                // (already calculated) residual f(x)
                if (!isGhost)
                    partialDeriv = storageResidual;
            }

            if (numericDifferenceMethod <= 0)
            {
                // we are not using forward differences, i.e. we
                // need to calculate f(x - \epsilon)

                // deflect the primary variables
                elemSol[0][pvIdx] -= delta + eps;
                delta += eps;

                // update the volume variables and the flux var cache
                curVolVars.update(elemSol, problem, element, scv);

                // calculate the residual with the deflected primary variables and subtract it
                if (!isGhost)
                {
                   partialDeriv -= localResidual.evalStorage(problem,
                                                             element,
                                                             fvGeometry,
                                                             prevElemVolVars,
                                                             curElemVolVars,
                                                             elemBcTypes,
                                                             elemFluxVarsCache)[0];
                }
            }
            else
            {
                // we are using forward differences, i.e. we don't need to
                // calculate f(x - \epsilon) and we can recycle the
                // (already calculated) residual f(x)
                if (!isGhost)
                    partialDeriv -= storageResidual;
            }

            // divide difference in residuals by the magnitude of the
            // deflections between the two function evaluation
            if (!isGhost)
                partialDeriv /= delta;

            // restore the original state of the scv's volume variables
            curVolVars = origVolVars;

            // restore the current element solution
            elemSol[0][pvIdx] = origPriVars[pvIdx];

            // add the current partial derivatives to the global jacobian matrix
            for (int eqIdx = 0; eqIdx < numEq; eqIdx++)
            {
                // the diagonal entries
                A[globalI][globalI][eqIdx][pvIdx] += partialDeriv[eqIdx];
            }
        }

        // return the original residual
        return residual;
    }
private:
    template<class T = TypeTag>
    static typename std::enable_if<!GET_PROP_VALUE(T, EnableGridVolumeVariablesCache), VolumeVariables&>::type
    getVolVarAccess(GridVolumeVariables& gridVolVars, ElementVolumeVariables& elemVolVars, const SubControlVolume& scv)
    { return elemVolVars[scv]; }

    template<class T = TypeTag>
    static typename std::enable_if<GET_PROP_VALUE(T, EnableGridVolumeVariablesCache), VolumeVariables&>::type
    getVolVarAccess(GridVolumeVariables& gridVolVars, ElementVolumeVariables& elemVolVars, const SubControlVolume& scv)
    { return gridVolVars.volVars(scv); }
};

/*!
 * \ingroup Assembly
 * \ingroup CCDiscretization
 * \brief Cell-centered scheme local assembler using analytic (hand-coded) differentiation and implicit time discretization
 */
template<class TypeTag>
class CCLocalAssembler<TypeTag,
                       DiffMethod::analytic,
                       /*implicit=*/true>
{
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using NumEqVector = typename GET_PROP_TYPE(TypeTag, NumEqVector);
    using ElementBoundaryTypes = typename GET_PROP_TYPE(TypeTag, ElementBoundaryTypes);
    using Element = typename GET_PROP_TYPE(TypeTag, GridView)::template Codim<0>::Entity;
    using IndexType = typename GET_PROP_TYPE(TypeTag, GridView)::IndexSet::IndexType;
    using SolutionVector = typename GET_PROP_TYPE(TypeTag, SolutionVector);
    using JacobianMatrix = typename GET_PROP_TYPE(TypeTag, JacobianMatrix);
    using ElementSolutionVector = typename GET_PROP_TYPE(TypeTag, ElementSolutionVector);
    using ElementVolumeVariables = typename GET_PROP_TYPE(TypeTag, ElementVolumeVariables);
    using GridVolumeVariables = typename GET_PROP_TYPE(TypeTag, GridVolumeVariables);
    using VolumeVariables = typename GET_PROP_TYPE(TypeTag, VolumeVariables);
    using SubControlVolume = typename GET_PROP_TYPE(TypeTag, SubControlVolume);

    enum { numEq = GET_PROP_VALUE(TypeTag, NumEq) };

public:

    /*!
     * \brief Computes the derivatives with respect to the given element and adds them
     *        to the global matrix. The element residual is written into the right hand side.
     */
    template<class Assembler>
    static void assemble(Assembler& assembler, JacobianMatrix& jac, SolutionVector& res,
                         const Element& element, const SolutionVector& curSol)
    {
        const auto globalI = assembler.fvGridGeometry().elementMapper().index(element);
        res[globalI] = assemble_(assembler, jac, element, curSol);
    }

    /*!
     * \brief Computes the derivatives with respect to the given element and adds them
     *        to the global matrix.
     */
    template<class Assembler>
    static void assemble(Assembler& assembler, JacobianMatrix& jac,
                         const Element& element, const SolutionVector& curSol)
    {
        assemble_(assembler, jac, element, curSol);
    }

    /*!
     * \brief Assemble the residual only
     */
    template<class Assembler>
    static void assemble(Assembler& assembler, SolutionVector& res,
                         const Element& element, const SolutionVector& curSol)
    {
        const auto globalI = assembler.fvGridGeometry().elementMapper().index(element);
        res[globalI] = assemble_(assembler, element, curSol);
    }

private:

    /*!
     * \brief Computes the residual
     *
     * \return The element residual at the current solution.
     */
    template<class Assembler>
    static NumEqVector assemble_(Assembler& assembler,
                                 const Element& element, const SolutionVector& curSol)
    {
        // is the actual element a ghost element?
        const bool isGhost = (element.partitionType() == Dune::GhostEntity);
        if (isGhost) return NumEqVector(0.0);

        // get some references for convenience
        const auto& problem = assembler.problem();
        auto& localResidual = assembler.localResidual();
        auto& gridVariables = assembler.gridVariables();

        // prepare the local views
        auto fvGeometry = localView(assembler.fvGridGeometry());
        fvGeometry.bind(element);

        auto curElemVolVars = localView(gridVariables.curGridVolVars());
        curElemVolVars.bind(element, fvGeometry, curSol);

        auto elemFluxVarsCache = localView(gridVariables.gridFluxVarsCache());
        elemFluxVarsCache.bind(element, fvGeometry, curElemVolVars);

        auto prevElemVolVars = localView(gridVariables.prevGridVolVars());

        // check for boundaries on the element
        // TODO Do we need them for cell-centered models?
        ElementBoundaryTypes elemBcTypes;
        elemBcTypes.update(problem, element, fvGeometry);

        NumEqVector residual(0.0);
        if (localResidual.isStationary())
        {
            residual = localResidual.eval(problem,
                                          element,
                                          fvGeometry,
                                          curElemVolVars,
                                          elemBcTypes,
                                          elemFluxVarsCache)[0];
        }
        else
        {
            prevElemVolVars.bindElement(element, fvGeometry, localResidual.prevSol());
            residual = localResidual.eval(problem,
                                          element,
                                          fvGeometry,
                                          prevElemVolVars,
                                          curElemVolVars,
                                          elemBcTypes,
                                          elemFluxVarsCache)[0];
        }

        return residual;
    }

    /*!
     * \brief Computes the derivatives with respect to the given element and adds them
     *        to the global matrix.
     *
     * \return The element residual at the current solution.
     */
    template<class Assembler>
    static NumEqVector assemble_(Assembler& assembler, JacobianMatrix& A,
                                 const Element& element, const SolutionVector& curSol)
    {
        // get some references for convenience
        const auto& problem = assembler.problem();
        const auto& fvGridGeometry = assembler.fvGridGeometry();
        auto& localResidual = assembler.localResidual();
        auto& gridVariables = assembler.gridVariables();

        // prepare the local views
        auto fvGeometry = localView(assembler.fvGridGeometry());
        fvGeometry.bind(element);

        auto curElemVolVars = localView(gridVariables.curGridVolVars());
        curElemVolVars.bind(element, fvGeometry, curSol);

        auto elemFluxVarsCache = localView(gridVariables.gridFluxVarsCache());
        elemFluxVarsCache.bind(element, fvGeometry, curElemVolVars);

        auto prevElemVolVars = localView(gridVariables.prevGridVolVars());

        // the global dof of the actual element
        const auto globalI = fvGridGeometry.elementMapper().index(element);

        // check for boundaries on the element
        // TODO Do we need them for cell-centered models?
        ElementBoundaryTypes elemBcTypes;
        elemBcTypes.update(problem, element, fvGeometry);

        // is the actual element a ghost element?
        const bool isGhost = (element.partitionType() == Dune::GhostEntity);
        // is this a stationary simulation?
        const bool isStationary = localResidual.isStationary();

        // the actual element's current residual (will be returned by this function)
        NumEqVector residual(0.0);
        if (!isGhost)
        {
            if (isStationary)
            {
                residual = localResidual.eval(problem,
                                              element,
                                              fvGeometry,
                                              curElemVolVars,
                                              elemBcTypes,
                                              elemFluxVarsCache)[0];
            }
            else
            {
                prevElemVolVars.bindElement(element, fvGeometry, localResidual.prevSol());
                residual = localResidual.eval(problem,
                                              element,
                                              fvGeometry,
                                              prevElemVolVars,
                                              curElemVolVars,
                                              elemBcTypes,
                                              elemFluxVarsCache)[0];
            }
        }

        // get reference to the element's current vol vars
        const auto& scv = fvGeometry.scv(globalI);
        const auto& volVars = curElemVolVars[scv];

        // if the problem is instationary, add derivative of storage term
        if (!isStationary)
            localResidual.addStorageDerivatives(A[globalI][globalI],
                                                problem,
                                                element,
                                                fvGeometry,
                                                volVars,
                                                scv);

        // add source term derivatives
        localResidual.addSourceDerivatives(A[globalI][globalI],
                                           problem,
                                           element,
                                           fvGeometry,
                                           volVars,
                                           scv);

        // add flux derivatives for each scvf
        for (const auto& scvf : scvfs(fvGeometry))
        {
            if (!scvf.boundary())
            {
                localResidual.addFluxDerivatives(A[globalI],
                                                 problem,
                                                 element,
                                                 fvGeometry,
                                                 curElemVolVars,
                                                 elemFluxVarsCache,
                                                 scvf);
            }
            else
            {
                const auto& bcTypes = problem.boundaryTypes(element, scvf);

                // add Dirichlet boundary flux derivatives
                if (bcTypes.hasDirichlet() && !bcTypes.hasNeumann())
                {
                    localResidual.addCCDirichletFluxDerivatives(A[globalI],
                                                                problem,
                                                                element,
                                                                fvGeometry,
                                                                curElemVolVars,
                                                                elemFluxVarsCache,
                                                                scvf);
                }
                // add Robin ("solution dependent Neumann") boundary flux derivatives
                else if (bcTypes.hasNeumann() && !bcTypes.hasDirichlet())
                {
                    localResidual.addRobinFluxDerivatives(A[globalI],
                                                          problem,
                                                          element,
                                                          fvGeometry,
                                                          curElemVolVars,
                                                          elemFluxVarsCache,
                                                          scvf);
                }
                else
                    DUNE_THROW(Dune::NotImplemented, "Mixed boundary conditions. Use pure boundary conditions by converting Dirichlet BCs to Robin BCs");
            }
        }

        // TODO Do we really need this??????????
        // this->model_().updatePVWeights(fvGeometry);

        // TODO: Additional dof dependencies???

        // return element residual
        return residual;
    }
};

/*!
 * \ingroup Assembly
 * \ingroup CCDiscretization
 * \brief Cell-centered scheme local assembler using analytic (hand-coded) differentiation and explicit time discretization
 */
template<class TypeTag>
class CCLocalAssembler<TypeTag,
                       DiffMethod::analytic,
                       /*implicit=*/false>
{
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using NumEqVector = typename GET_PROP_TYPE(TypeTag, NumEqVector);
    using ElementBoundaryTypes = typename GET_PROP_TYPE(TypeTag, ElementBoundaryTypes);
    using Element = typename GET_PROP_TYPE(TypeTag, GridView)::template Codim<0>::Entity;
    using IndexType = typename GET_PROP_TYPE(TypeTag, GridView)::IndexSet::IndexType;
    using SolutionVector = typename GET_PROP_TYPE(TypeTag, SolutionVector);
    using JacobianMatrix = typename GET_PROP_TYPE(TypeTag, JacobianMatrix);
    using ElementSolutionVector = typename GET_PROP_TYPE(TypeTag, ElementSolutionVector);
    using ElementVolumeVariables = typename GET_PROP_TYPE(TypeTag, ElementVolumeVariables);
    using GridVolumeVariables = typename GET_PROP_TYPE(TypeTag, GridVolumeVariables);
    using VolumeVariables = typename GET_PROP_TYPE(TypeTag, VolumeVariables);
    using SubControlVolume = typename GET_PROP_TYPE(TypeTag, SubControlVolume);

    enum { numEq = GET_PROP_VALUE(TypeTag, NumEq) };

public:

    /*!
     * \brief Computes the derivatives with respect to the given element and adds them
     *        to the global matrix. The element residual is written into the right hand side.
     */
    template<class Assembler>
    static void assemble(Assembler& assembler, JacobianMatrix& jac, SolutionVector& res,
                         const Element& element, const SolutionVector& curSol)
    {
        const auto globalI = assembler.fvGridGeometry().elementMapper().index(element);
        res[globalI] = assemble_(assembler, jac, element, curSol);
    }

    /*!
     * \brief Computes the derivatives with respect to the given element and adds them
     *        to the global matrix.
     */
    template<class Assembler>
    static void assemble(Assembler& assembler, JacobianMatrix& jac,
                         const Element& element, const SolutionVector& curSol)
    {
        assemble_(assembler, jac, element, curSol);
    }

    /*!
     * \brief Assemble the residual only
     */
    template<class Assembler>
    static void assemble(Assembler& assembler, SolutionVector& res,
                         const Element& element, const SolutionVector& curSol)
    {
        const auto globalI = assembler.fvGridGeometry().elementMapper().index(element);
        res[globalI] = assemble_(assembler, element, curSol);
    }

private:

    /*!
     * \brief Computes the residual
     *
     * \return The element residual at the current solution.
     */
    template<class Assembler>
    static NumEqVector assemble_(Assembler& assembler,
                                 const Element& element, const SolutionVector& curSol)
    {
        // is the actual element a ghost element?
        const bool isGhost = (element.partitionType() == Dune::GhostEntity);
        if (isGhost) return NumEqVector(0.0);

        // get some references for convenience
        const auto& problem = assembler.problem();
        auto& localResidual = assembler.localResidual();
        auto& gridVariables = assembler.gridVariables();

        // using an explicit assembler doesn't make sense for stationary problems
        if (localResidual.isStationary())
            DUNE_THROW(Dune::InvalidStateException, "Using explicit jacobian assembler with stationary local residual");

        // prepare the local views
        auto fvGeometry = localView(assembler.fvGridGeometry());
        fvGeometry.bind(element);

        auto curElemVolVars = localView(gridVariables.curGridVolVars());
        curElemVolVars.bindElement(element, fvGeometry, curSol);

        auto prevElemVolVars = localView(gridVariables.prevGridVolVars());
        prevElemVolVars.bind(element, fvGeometry, localResidual.prevSol());

        auto elemFluxVarsCache = localView(gridVariables.gridFluxVarsCache());
        elemFluxVarsCache.bind(element, fvGeometry, prevElemVolVars);

        // check for boundaries on the element
        // TODO Do we need them for cell-centered models?
        ElementBoundaryTypes elemBcTypes;
        elemBcTypes.update(problem, element, fvGeometry);

        // the actual element's previous time step residual
        auto residual = localResidual.eval(problem,
                                           element,
                                           fvGeometry,
                                           prevElemVolVars,
                                           elemBcTypes,
                                           elemFluxVarsCache)[0];

        auto storageResidual = localResidual.evalStorage(problem,
                                                         element,
                                                         fvGeometry,
                                                         prevElemVolVars,
                                                         curElemVolVars,
                                                         elemBcTypes,
                                                         elemFluxVarsCache)[0];

        residual += storageResidual;

        return residual;
    }

    /*!
     * \brief Computes the derivatives with respect to the given element and adds them
     *        to the global matrix.
     *
     * \return The element residual at the current solution.
     */
    template<class Assembler>
    static NumEqVector assemble_(Assembler& assembler, JacobianMatrix& A,
                                 const Element& element, const SolutionVector& curSol)
    {
        // get some references for convenience
        const auto& problem = assembler.problem();
        const auto& fvGridGeometry = assembler.fvGridGeometry();
        auto& localResidual = assembler.localResidual();
        auto& gridVariables = assembler.gridVariables();

        // using an explicit assembler doesn't make sense for stationary problems
        if (localResidual.isStationary())
            DUNE_THROW(Dune::InvalidStateException, "Using explicit jacobian assembler with stationary local residual");

        // prepare the local views
        auto fvGeometry = localView(assembler.fvGridGeometry());
        fvGeometry.bind(element);

        auto curElemVolVars = localView(gridVariables.curGridVolVars());
        curElemVolVars.bindElement(element, fvGeometry, curSol);

        auto prevElemVolVars = localView(gridVariables.prevGridVolVars());
        prevElemVolVars.bind(element, fvGeometry, localResidual.prevSol());

        auto elemFluxVarsCache = localView(gridVariables.gridFluxVarsCache());
        elemFluxVarsCache.bind(element, fvGeometry, prevElemVolVars);

        // the global dof of the actual element
        const auto globalI = fvGridGeometry.elementMapper().index(element);

        // check for boundaries on the element
        // TODO Do we need them for cell-centered models?
        ElementBoundaryTypes elemBcTypes;
        elemBcTypes.update(problem, element, fvGeometry);

        // is the actual element a ghost element?
        const bool isGhost = (element.partitionType() == Dune::GhostEntity);

        // the actual element's previous time step residual
        NumEqVector residual(0.0), storageResidual(0.0);
        if (!isGhost)
        {
            residual = localResidual.eval(problem,
                                          element,
                                          fvGeometry,
                                          prevElemVolVars,
                                          elemBcTypes,
                                          elemFluxVarsCache)[0];

            storageResidual = localResidual.evalStorage(problem,
                                                        element,
                                                        fvGeometry,
                                                        prevElemVolVars,
                                                        curElemVolVars,
                                                        elemBcTypes,
                                                        elemFluxVarsCache)[0];

            residual += storageResidual;
        }

        // get reference to the element's current vol vars
        const auto& scv = fvGeometry.scv(globalI);
        const auto& volVars = curElemVolVars[scv];

        // add derivative of storage term
        localResidual.addStorageDerivatives(A[globalI][globalI],
                                            problem,
                                            element,
                                            fvGeometry,
                                            volVars,
                                            scv);

        // return the original residual
        return residual;
    }
};

} // end namespace Dumux

#endif
