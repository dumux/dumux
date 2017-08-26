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
#ifndef DUMUX_CC_EXPLICIT_LOCAL_ASSEMBLER_HH
#define DUMUX_CC_EXPLICIT_LOCAL_ASSEMBLER_HH

#include <dune/istl/matrixindexset.hh>
#include <dune/istl/bvector.hh>

#include <dumux/implicit/properties.hh>

namespace Dumux {

/*!
 * \ingroup ImplicitModel
 * \brief An assembler for the local contributions (per element) to the global
 *        linear system for fully implicit models and cell-centered discretization schemes.
 */
template<class TypeTag, DifferentiationMethods DM>
class CCExplicitLocalAssembler;


template<class TypeTag>
class CCExplicitLocalAssembler<TypeTag, DifferentiationMethods::numeric>
{
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using NumEqVector = typename GET_PROP_TYPE(TypeTag, NumEqVector);
    using ElementBoundaryTypes = typename GET_PROP_TYPE(TypeTag, ElementBoundaryTypes);
    using Element = typename GET_PROP_TYPE(TypeTag, GridView)::template Codim<0>::Entity;
    using SolutionVector = typename GET_PROP_TYPE(TypeTag, SolutionVector);
    using ElementSolutionVector = typename GET_PROP_TYPE(TypeTag, ElementSolutionVector);
    using ElementVolumeVariables = typename GET_PROP_TYPE(TypeTag, ElementVolumeVariables);
    using GridVolumeVariables = typename GET_PROP_TYPE(TypeTag, GlobalVolumeVariables);
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

        //////////////////////////////////////////////////////////////////////////////////////////////////
        // Calculate derivatives of all dofs in stencil with respect to the dofs in the element. In the //
        // neighboring elements all derivatives are zero. For the assembled element only the storage    //
        // derivatives are non-zero.                                                                    //
        //////////////////////////////////////////////////////////////////////////////////////////////////

        static const int numericDifferenceMethod = GET_PARAM_FROM_GROUP(TypeTag, int, Implicit, NumericDifferenceMethod);

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
                elemFluxVarsCache.update(element, fvGeometry, curElemVolVars);

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
    static typename std::enable_if<!GET_PROP_VALUE(T, EnableGlobalVolumeVariablesCache), VolumeVariables&>::type
    getVolVarAccess(GridVolumeVariables& gridVolVars, ElementVolumeVariables& elemVolVars, const SubControlVolume& scv)
    { return elemVolVars[scv]; }

    template<class T = TypeTag>
    static typename std::enable_if<GET_PROP_VALUE(T, EnableGlobalVolumeVariablesCache), VolumeVariables&>::type
    getVolVarAccess(GridVolumeVariables& gridVolVars, ElementVolumeVariables& elemVolVars, const SubControlVolume& scv)
    { return gridVolVars.volVars(scv); }
};

} // end namespace Dumux

#endif
