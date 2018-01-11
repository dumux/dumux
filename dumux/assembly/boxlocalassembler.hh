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
 * \brief An assembler for Jacobian and residual contribution per element (box method)
 */
#ifndef DUMUX_BOX_LOCAL_ASSEMBLER_HH
#define DUMUX_BOX_LOCAL_ASSEMBLER_HH

#include <dune/istl/matrixindexset.hh>
#include <dune/istl/bvector.hh>

#include <dumux/common/properties.hh>
#include <dumux/common/parameters.hh>
#include <dumux/assembly/diffmethod.hh>

#include "fvlocalassemblerbase.hh"

namespace Dumux {

/*!
 * \ingroup Assembly
 * \ingroup BoxDiscretization
 * \brief A base class for all local box assemblers
 * \tparam TypeTag the TypeTag
 * \tparam Assembler the assembler type
 * \tparam Assembler The acutal Implementation
 */
template<class TypeTag, class Assembler, class Implementation, bool implicit>
class BoxLocalAssemblerBase : public FVLocalAssemblerBase<TypeTag, Assembler, Implementation, implicit>
{
    using ParentType = FVLocalAssemblerBase<TypeTag, Assembler, Implementation, implicit>;
    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
    using JacobianMatrix = typename GET_PROP_TYPE(TypeTag, JacobianMatrix);
    using GridVariables = typename GET_PROP_TYPE(TypeTag, GridVariables);
    using SolutionVector = typename GET_PROP_TYPE(TypeTag, SolutionVector);
    using ElementVolumeVariables = typename GET_PROP_TYPE(TypeTag, ElementVolumeVariables);

    enum { numEq = GET_PROP_VALUE(TypeTag, NumEq) };

public:

    using ParentType::ParentType;

    void bindLocalViews()
    {
        ParentType::bindLocalViews();
        this->elemBcTypes().update(this->problem(), this->element(), this->fvGeometry());
    }


    /*!
     * \brief Computes the derivatives with respect to the given element and adds them
     *        to the global matrix. The element residual is written into the right hand side.
     */
    void assembleJacobianAndResidual(JacobianMatrix& jac, SolutionVector& res, GridVariables& gridVariables)
    {
        this->asImp_().bindLocalViews();

        const auto residual = this->asImp_().assembleJacobianAndResidualImpl(jac, gridVariables);

        for (const auto& scv : scvs(this->fvGeometry()))
            res[scv.dofIndex()] += residual[scv.indexInElement()];

        // enforce Dirichlet boundaries by overwriting partial derivatives with 1 or 0
        // and set the residual to (privar - dirichletvalue)
        // TODO: put this in separate method!
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

                            const auto& priVars = this->curElemVolVars()[scvI].priVars();
                            res[scvI.dofIndex()][eqIdx] = priVars[pvIdx] - dirichletValues[pvIdx];
                            for (const auto& scvJ : scvs(this->fvGeometry()))
                            {
                                jac[scvI.dofIndex()][scvJ.dofIndex()][eqIdx] = 0.0;
                                if (scvI.indexInElement() == scvJ.indexInElement())
                                    jac[scvI.dofIndex()][scvI.dofIndex()][eqIdx][pvIdx] = 1.0;
                            }
                        }
                    }
                }
            }
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
        const auto residual = this->evalLocalResidual();

        for (const auto& scv : scvs(this->fvGeometry()))
            res[scv.dofIndex()] += residual[scv.indexInElement()];
    }


    auto evalLocalFluxSourceResidual(const ElementVolumeVariables& elemVolVars) const
    {
        return this->evalLocalFluxSourceResidual(elemVolVars);
    }

    auto evalLocalStorageResidual() const
    {
        return this->evalLocalStorageResidual();
    }
};

/*!
 * \ingroup Assembly
 * \ingroup BoxDiscretization
 * \brief An assembler for Jacobian and residual contribution per element (box method)
 * \tparam TypeTag the TypeTag
 * \tparam DM the differentiation method to residual compute derivatives
 * \tparam implicit if to use an implicit or explicit time discretization
 */
template<class TypeTag, class Assembler, DiffMethod DM = DiffMethod::numeric, bool implicit = true>
class BoxLocalAssembler;

/*!
 * \ingroup Assembly
 * \ingroup BoxDiscretization
 * \brief Box local assembler using numeric differentiation and implicit time discretization
 */
template<class TypeTag, class Assembler>
class BoxLocalAssembler<TypeTag, Assembler, DiffMethod::numeric, /*implicit=*/true>
: public BoxLocalAssemblerBase<TypeTag, Assembler,
                              BoxLocalAssembler<TypeTag, Assembler, DiffMethod::numeric, true>, true>
{
    using ThisType = BoxLocalAssembler<TypeTag, Assembler, DiffMethod::numeric, true>;
    using ParentType = BoxLocalAssemblerBase<TypeTag, Assembler, ThisType, true>;
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using Element = typename GET_PROP_TYPE(TypeTag, GridView)::template Codim<0>::Entity;
    using SolutionVector = typename GET_PROP_TYPE(TypeTag, SolutionVector);
    using ElementSolutionVector = typename GET_PROP_TYPE(TypeTag, ElementSolutionVector);
    using ElementVolumeVariables = typename GET_PROP_TYPE(TypeTag, ElementVolumeVariables);
    using GridVolumeVariables = typename GET_PROP_TYPE(TypeTag, GridVolumeVariables);
    using GridVariables = typename GET_PROP_TYPE(TypeTag, GridVariables);
    using VolumeVariables = typename GET_PROP_TYPE(TypeTag, VolumeVariables);
    using JacobianMatrix = typename GET_PROP_TYPE(TypeTag, JacobianMatrix);

    enum { numEq = GET_PROP_VALUE(TypeTag, NumEq) };
    enum { dim = GET_PROP_TYPE(TypeTag, GridView)::dimension };

    static constexpr bool enableGridFluxVarsCache = GET_PROP_VALUE(TypeTag, EnableGridFluxVariablesCache);
    static constexpr int maxNeighbors = 4*(2*dim);

public:

    using ParentType::ParentType;

    /*!
     * \brief Computes the derivatives with respect to the given element and adds them
     *        to the global matrix.
     *
     * \return The element residual at the current solution.
     */
    auto assembleJacobianAndResidualImpl(JacobianMatrix& A, GridVariables& gridVariables)
    {
        // get some aliases for convenience
        const auto& element = this->element();
        const auto& fvGeometry = this->fvGeometry();
        const auto& curSol = this->curSol();
        auto&& curElemVolVars = this->curElemVolVars();
        auto&& elemFluxVarsCache = this->elemFluxVarsCache();

        // get the vecor of the acutal element residuals
        const auto origResiduals = this->evalLocalResidual();

        //////////////////////////////////////////////////////////////////////////////////////////////////
        //                                                                                              //
        // Calculate derivatives of all dofs in stencil with respect to the dofs in the element. In the //
        // neighboring elements we do so by computing the derivatives of the fluxes which depend on the //
        // actual element. In the actual element we evaluate the derivative of the entire residual.     //
        //                                                                                              //
        //////////////////////////////////////////////////////////////////////////////////////////////////

        // create the element solution
        ElementSolutionVector elemSol(element, curSol, fvGeometry);

        using ElementResidualVector = std::decay_t<decltype(origResiduals)>;
        ElementResidualVector partialDerivs(element.subEntities(dim));

        // calculation of the derivatives
        for (auto&& scv : scvs(fvGeometry))
        {
            // dof index and corresponding actual pri vars
            const auto dofIdx = scv.dofIndex();
            auto& curVolVars = this->getVolVarAccess(gridVariables.curGridVolVars(), curElemVolVars, scv);
            const VolumeVariables origVolVars(curVolVars);

            // calculate derivatives w.r.t to the privars at the dof at hand
            for (int pvIdx = 0; pvIdx < numEq; pvIdx++)
            {
                partialDerivs = 0.0;

                auto evalResiduals = [&](Scalar priVar)
                {
                    // auto partialDerivsTmp = partialDerivs;
                    // update the volume variables and the flux var cache
                    elemSol[scv.indexInElement()][pvIdx] = priVar;
                    curVolVars.update(elemSol, this->problem(), element, scv);
                    // if (enableGridFluxVarsCache)
                    //     gridVariables.gridFluxVarsCache().updateElement(element, fvGeometry, curElemVolVars);
                    // else
                    //     elemFluxVarsCache.update(element, fvGeometry, curElemVolVars);
                    //TODO: fluxvar caches!
                    //



                    return this->evalLocalResidual();
                };

                // derive the residuals numerically
                NumericDifferentiation::partialDerivative(evalResiduals, elemSol[scv.indexInElement()][pvIdx], partialDerivs, origResiduals);

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
                // TODO additional dof dependencies
            }
        }
        return origResiduals;
    }

}; // implicit BoxAssembler with numeric Jacobian

/*!
 * \ingroup Assembly
 * \ingroup BoxDiscretization
 * \brief Box local assembler using numeric differentiation and explicit time discretization
 */
template<class TypeTag, class Assembler>
class BoxLocalAssembler<TypeTag, Assembler, DiffMethod::numeric, /*implicit=*/false>
: public BoxLocalAssemblerBase<TypeTag, Assembler,
                              BoxLocalAssembler<TypeTag, Assembler, DiffMethod::numeric, false>, false>
{
    using ThisType = BoxLocalAssembler<TypeTag, Assembler, DiffMethod::numeric, false>;
    using ParentType = BoxLocalAssemblerBase<TypeTag, Assembler, ThisType, false>;
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using Element = typename GET_PROP_TYPE(TypeTag, GridView)::template Codim<0>::Entity;
    using SolutionVector = typename GET_PROP_TYPE(TypeTag, SolutionVector);
    using ElementSolutionVector = typename GET_PROP_TYPE(TypeTag, ElementSolutionVector);
    using ElementVolumeVariables = typename GET_PROP_TYPE(TypeTag, ElementVolumeVariables);
    using GridVolumeVariables = typename GET_PROP_TYPE(TypeTag, GridVolumeVariables);
    using GridVariables = typename GET_PROP_TYPE(TypeTag, GridVariables);
    using VolumeVariables = typename GET_PROP_TYPE(TypeTag, VolumeVariables);
    using JacobianMatrix = typename GET_PROP_TYPE(TypeTag, JacobianMatrix);

    enum { numEq = GET_PROP_VALUE(TypeTag, NumEq) };
    enum { dim = GET_PROP_TYPE(TypeTag, GridView)::dimension };

    static constexpr bool enableGridFluxVarsCache = GET_PROP_VALUE(TypeTag, EnableGridFluxVariablesCache);
    static constexpr int maxNeighbors = 4*(2*dim);

public:

    using ParentType::ParentType;

    /*!
     * \brief Computes the derivatives with respect to the given element and adds them
     *        to the global matrix.
     *
     * \return The element residual at the current solution.
     */
    auto assembleJacobianAndResidualImpl(JacobianMatrix& A, GridVariables& gridVariables)
    {
        if (this->assembler().isStationaryProblem())
            DUNE_THROW(Dune::InvalidStateException, "Using explicit jacobian assembler with stationary local residual");

        // get some aliases for convenience
        const auto& element = this->element();
        const auto& fvGeometry = this->fvGeometry();
        const auto& curSol = this->curSol();
        auto&& curElemVolVars = this->curElemVolVars();
        auto&& elemFluxVarsCache = this->elemFluxVarsCache();

        // get the vecor of the acutal element residuals
        const auto origResiduals = this->evalLocalResidual();

        //////////////////////////////////////////////////////////////////////////////////////////////////
        //                                                                                              //
        // Calculate derivatives of all dofs in stencil with respect to the dofs in the element. In the //
        // neighboring elements we do so by computing the derivatives of the fluxes which depend on the //
        // actual element. In the actual element we evaluate the derivative of the entire residual.     //
        //                                                                                              //
        //////////////////////////////////////////////////////////////////////////////////////////////////

        // create the element solution
        ElementSolutionVector elemSol(element, curSol, fvGeometry);

        using ElementResidualVector = std::decay_t<decltype(origResiduals)>;
        ElementResidualVector partialDerivs(element.subEntities(dim));

        // calculation of the derivatives
        for (auto&& scv : scvs(fvGeometry))
        {
            // dof index and corresponding actual pri vars
            const auto dofIdx = scv.dofIndex();
            auto& curVolVars = this->getVolVarAccess(gridVariables.curGridVolVars(), curElemVolVars, scv);
            const VolumeVariables origVolVars(curVolVars);

            // calculate derivatives w.r.t to the privars at the dof at hand
            for (int pvIdx = 0; pvIdx < numEq; pvIdx++)
            {
                partialDerivs = 0.0;

                auto evalStorage = [&](Scalar priVar)
                {
                    // auto partialDerivsTmp = partialDerivs;
                    // update the volume variables and the flux var cache
                    elemSol[scv.indexInElement()][pvIdx] = priVar;
                    curVolVars.update(elemSol, this->problem(), element, scv);
                    return this->evalLocalStorageResidual();
                };

                // derive the residuals numerically
                NumericDifferentiation::partialDerivative(evalStorage, elemSol[scv.indexInElement()][pvIdx], partialDerivs, origResiduals);

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
                // TODO additional dof dependencies
            }
        }
        return origResiduals;
    }
}; // explicit BoxAssembler with numeric Jacobian
//
//
// /*!
//  * \ingroup Assembly
//  * \ingroup BoxDiscretization
//  * \brief Box local assembler using numeric differentiation and explicit time discretization
//  */
// template<class TypeTag>
// class BoxLocalAssembler<TypeTag,
//                        DiffMethod::numeric,
//                        /*implicit=*/false>
// {
//     using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
//     using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
//     using ResidualVector = typename GET_PROP_TYPE(TypeTag, NumEqVector);
//     using ElementResidualVector = Dune::BlockVector<typename GET_PROP_TYPE(TypeTag, NumEqVector)>;
//     using ElementBoundaryTypes = typename GET_PROP_TYPE(TypeTag, ElementBoundaryTypes);
//     using Element = typename GET_PROP_TYPE(TypeTag, GridView)::template Codim<0>::Entity;
//     using SolutionVector = typename GET_PROP_TYPE(TypeTag, SolutionVector);
//     using ElementSolutionVector = typename GET_PROP_TYPE(TypeTag, ElementSolutionVector);
//     using ElementVolumeVariables = typename GET_PROP_TYPE(TypeTag, ElementVolumeVariables);
//     using GridVolumeVariables = typename GET_PROP_TYPE(TypeTag, GridVolumeVariables);
//     using VolumeVariables = typename GET_PROP_TYPE(TypeTag, VolumeVariables);
//     using SubControlVolume = typename GET_PROP_TYPE(TypeTag, SubControlVolume);
//     using JacobianMatrix = typename GET_PROP_TYPE(TypeTag, JacobianMatrix);
//
//     enum { numEq = GET_PROP_VALUE(TypeTag, NumEq) };
//
//     static constexpr int dim = GridView::dimension;
//
// public:
//
//     /*!
//      * \brief Computes the derivatives with respect to the given element and adds them
//      *        to the global matrix. The element residual is written into the right hand side.
//      */
//     template<class Assembler>
//     static void assemble(Assembler& assembler, JacobianMatrix& jac, SolutionVector& res,
//                          const Element& element, const SolutionVector& curSol)
//     {
//         assemble_(assembler, jac, res, element, curSol);
//
//     }
//
//     /*!
//      * \brief Computes the derivatives with respect to the given element and adds them
//      *        to the global matrix.
//      */
//     template<class Assembler>
//     static void assemble(Assembler& assembler, JacobianMatrix& jac,
//                          const Element& element, const SolutionVector& curSol)
//     {
//         auto dummyresidual = curSol;
//         assemble_(assembler, jac, dummyresidual, element, curSol);
//     }
//
//     /*!
//      * \brief Assemble the residual only
//      */
//     template<class Assembler>
//     static void assemble(Assembler& assembler, SolutionVector& res,
//                          const Element& element, const SolutionVector& curSol)
//     {
//         assemble_(assembler, res, element, curSol);
//     }
//
//     /*!
//      * \brief Computes the epsilon used for numeric differentiation
//      *        for a given value of a primary variable.
//      *
//      * \param priVar The value of the primary variable
//      */
//     static Scalar numericEpsilon(const Scalar priVar)
//     {
//         // define the base epsilon as the geometric mean of 1 and the
//         // resolution of the scalar type. E.g. for standard 64 bit
//         // floating point values, the resolution is about 10^-16 and
//         // the base epsilon is thus approximately 10^-8.
//         /*
//         static const Scalar baseEps
//             = Dumux::geometricMean<Scalar>(std::numeric_limits<Scalar>::epsilon(), 1.0);
//         */
//         static const Scalar baseEps = 1e-10;
//         assert(std::numeric_limits<Scalar>::epsilon()*1e4 < baseEps);
//         // the epsilon value used for the numeric differentiation is
//         // now scaled by the absolute value of the primary variable...
//         return baseEps*(std::abs(priVar) + 1.0);
//     }
//
// private:
//     /*!
//      * \brief Computes the residual
//      *
//      * \return The element residual at the current solution.
//      */
//     template<class Assembler>
//     static void assemble_(Assembler& assembler, SolutionVector& r,
//                           const Element& element, const SolutionVector& curSol)
//     {
//         // get some references for convenience
//         const auto& problem = assembler.problem();
//         auto& localResidual = assembler.localResidual();
//         auto& gridVariables = assembler.gridVariables();
//
//         // using an explicit assembler doesn't make sense for stationary problems
//         if (localResidual.isStationary())
//             DUNE_THROW(Dune::InvalidStateException, "Using explicit jacobian assembler with stationary local residual");
//
//         // prepare the local views
//         auto fvGeometry = localView(assembler.fvGridGeometry());
//         fvGeometry.bind(element);
//
//         auto curElemVolVars = localView(gridVariables.curGridVolVars());
//         curElemVolVars.bindElement(element, fvGeometry, curSol);
//
//         auto prevElemVolVars = localView(gridVariables.prevGridVolVars());
//         prevElemVolVars.bind(element, fvGeometry, localResidual.prevSol());
//
//         auto elemFluxVarsCache = localView(gridVariables.gridFluxVarsCache());
//         elemFluxVarsCache.bind(element, fvGeometry, prevElemVolVars);
//
//         // the element boundary types
//         ElementBoundaryTypes elemBcTypes;
//         elemBcTypes.update(problem, element, fvGeometry);
//
//         // the actual element's current residual
//         auto residual = localResidual.eval(problem,
//                                            element,
//                                            fvGeometry,
//                                            prevElemVolVars,
//                                            elemBcTypes,
//                                            elemFluxVarsCache);
//
//         auto storageResidual = localResidual.evalStorage(problem,
//                                                          element,
//                                                          fvGeometry,
//                                                          prevElemVolVars,
//                                                          curElemVolVars,
//                                                          elemBcTypes,
//                                                          elemFluxVarsCache);
//         residual += storageResidual;
//
//         for (const auto& scv : scvs(fvGeometry))
//             r[scv.dofIndex()] += residual[scv.indexInElement()];
//
//         // enforce Dirichlet boundaries by setting the residual to (privar - dirichletvalue)
//         if (elemBcTypes.hasDirichlet())
//         {
//             for (const auto& scvI : scvs(fvGeometry))
//             {
//                 const auto bcTypes = elemBcTypes[scvI.indexInElement()];
//                 if (bcTypes.hasDirichlet())
//                 {
//                     const auto dirichletValues = problem.dirichlet(element, scvI);
//
//                     // set the dirichlet conditions in residual and jacobian
//                     for (int eqIdx = 0; eqIdx < numEq; ++eqIdx)
//                     {
//                         if (bcTypes.isDirichlet(eqIdx))
//                         {
//                             const auto pvIdx = bcTypes.eqToDirichletIndex(eqIdx);
//                             assert(0 <= pvIdx && pvIdx < numEq);
//
//                             const auto& priVars = curElemVolVars[scvI].priVars();
//                             r[scvI.dofIndex()][eqIdx] = priVars[pvIdx] - dirichletValues[pvIdx];
//                         }
//                     }
//                 }
//             }
//         }
//     }
//
//     /*!
//      * \brief Computes the derivatives with respect to the given element and adds them
//      *        to the global matrix.
//      *
//      * \return The element residual at the current solution.
//      */
//     template<class Assembler>
//     static void assemble_(Assembler& assembler, JacobianMatrix& A, SolutionVector& r,
//                           const Element& element, const SolutionVector& curSol)
//     {
//         // get some references for convenience
//         const auto& problem = assembler.problem();
//         auto& localResidual = assembler.localResidual();
//         auto& gridVariables = assembler.gridVariables();
//
//         // using an explicit assembler doesn't make sense for stationary problems
//         if (localResidual.isStationary())
//             DUNE_THROW(Dune::InvalidStateException, "Using explicit jacobian assembler with stationary local residual");
//
//         // prepare the local views
//         auto fvGeometry = localView(assembler.fvGridGeometry());
//         fvGeometry.bind(element);
//
//         auto curElemVolVars = localView(gridVariables.curGridVolVars());
//         curElemVolVars.bindElement(element, fvGeometry, curSol);
//
//         auto prevElemVolVars = localView(gridVariables.prevGridVolVars());
//         prevElemVolVars.bind(element, fvGeometry, localResidual.prevSol());
//
//         auto elemFluxVarsCache = localView(gridVariables.gridFluxVarsCache());
//         elemFluxVarsCache.bind(element, fvGeometry, prevElemVolVars);
//
//         // the element boundary types
//         ElementBoundaryTypes elemBcTypes;
//         elemBcTypes.update(problem, element, fvGeometry);
//
//         // get the element solution
//         const auto numVert = element.subEntities(dim);
//         ElementSolutionVector elemSol(numVert);
//         for (const auto& scv : scvs(fvGeometry))
//             elemSol[scv.indexInElement()] = localResidual.prevSol()[scv.dofIndex()];
//
//         // the actual element's previous time step residual
//         auto residual = localResidual.eval(problem,
//                                            element,
//                                            fvGeometry,
//                                            prevElemVolVars,
//                                            elemBcTypes,
//                                            elemFluxVarsCache);
//
//         auto storageResidual = localResidual.evalStorage(problem,
//                                                          element,
//                                                          fvGeometry,
//                                                          prevElemVolVars,
//                                                          curElemVolVars,
//                                                          elemBcTypes,
//                                                          elemFluxVarsCache);
//
//         residual += storageResidual;
//
//         //////////////////////////////////////////////////////////////////////////////////////////////////
//         //                                                                                              //
//         // Calculate derivatives of all dofs in stencil with respect to the dofs in the element. In the //
//         // neighboring elements we do so by computing the derivatives of the fluxes which depend on the //
//         // actual element. In the actual element we evaluate the derivative of the entire residual.     //
//         //                                                                                              //
//         //////////////////////////////////////////////////////////////////////////////////////////////////
//
//         static const int numericDifferenceMethod = getParamFromGroup<int>(GET_PROP_VALUE(TypeTag, ModelParameterGroup), "Implicit.NumericDifferenceMethod");
//
//         // calculation of the derivatives
//         for (auto&& scv : scvs(fvGeometry))
//         {
//             // dof index and corresponding actual pri vars
//             const auto dofIdx = scv.dofIndex();
//             auto& volVars = getVolVarAccess(gridVariables.curGridVolVars(), curElemVolVars, scv);
//             VolumeVariables origVolVars(volVars);
//
//             // add precalculated residual for this scv into the global container
//             r[dofIdx] += residual[scv.indexInElement()];
//
//             // calculate derivatives w.r.t to the privars at the dof at hand
//             for (int pvIdx = 0; pvIdx < numEq; pvIdx++)
//             {
//                 ElementSolutionVector partialDeriv(element.subEntities(dim));
//                 Scalar eps = numericEpsilon(volVars.priVar(pvIdx));
//                 Scalar delta = 0;
//
//                 // calculate the residual with the forward deflected primary variables
//                 if (numericDifferenceMethod >= 0)
//                 {
//                     // we are not using backward differences, i.e. we need to
//                     // calculate f(x + \epsilon)
//
//                     // deflect primary variables
//                     elemSol[scv.indexInElement()][pvIdx] += eps;
//                     delta += eps;
//
//                     // update the volume variables connected to the dof
//                     volVars.update(elemSol, problem, element, scv);
//
//                     partialDeriv = localResidual.evalStorage(problem, element, fvGeometry,
//                                                              prevElemVolVars, curElemVolVars,
//                                                              elemBcTypes, elemFluxVarsCache);
//                 }
//                 else
//                 {
//                     // we are using backward differences, i.e. we don't need
//                     // to calculate f(x + \epsilon) and we can recycle the
//                     // (already calculated) residual f(x)
//                     partialDeriv = residual;
//                 }
//
//                 if (numericDifferenceMethod <= 0)
//                 {
//                     // we are not using forward differences, i.e. we
//                     // need to calculate f(x - \epsilon)
//
//                     // deflect the primary variables
//                     elemSol[scv.indexInElement()][pvIdx] -= delta + eps;
//                     delta += eps;
//
//                     // update the volume variables connected to the dof
//                     volVars.update(elemSol, problem, element, scv);
//
//                     partialDeriv -= localResidual.evalStorage(problem, element, fvGeometry,
//                                                               prevElemVolVars, curElemVolVars,
//                                                               elemBcTypes, elemFluxVarsCache);
//                 }
//                 else
//                 {
//                     // we are using forward differences, i.e. we don't need to
//                     // calculate f(x - \epsilon) and we can recycle the
//                     // (already calculated) residual f(x)
//                     partialDeriv -= residual;
//                 }
//
//                 // divide difference in residuals by the magnitude of the
//                 // deflections between the two function evaluation
//                 partialDeriv /= delta;
//
//                 // update the global stiffness matrix with the current partial derivatives
//                 for (int eqIdx = 0; eqIdx < numEq; eqIdx++)
//                 {
//                     // A[i][col][eqIdx][pvIdx] is the rate of change of
//                     // the residual of equation 'eqIdx' at dof 'i'
//                     // depending on the primary variable 'pvIdx' at dof
//                     // 'col'.
//                     A[dofIdx][dofIdx][eqIdx][pvIdx] += partialDeriv[scv.indexInElement()][eqIdx];
//                 }
//
//                 // restore the original state of the scv's volume variables
//                 volVars = origVolVars;
//
//                 // restore the original element solution
//                 elemSol[scv.indexInElement()][pvIdx] = curSol[scv.dofIndex()][pvIdx];
//             }
//
//             // TODO additional dof dependencies
//         }
//
//         // enforce Dirichlet boundaries by overwriting partial derivatives with 1 or 0
//         // and set the residual to (privar - dirichletvalue)
//         if (elemBcTypes.hasDirichlet())
//         {
//             for (const auto& scvI : scvs(fvGeometry))
//             {
//                 const auto bcTypes = elemBcTypes[scvI.indexInElement()];
//                 if (bcTypes.hasDirichlet())
//                 {
//                     const auto dirichletValues = problem.dirichlet(element, scvI);
//
//                     // set the dirichlet conditions in residual and jacobian
//                     for (int eqIdx = 0; eqIdx < numEq; ++eqIdx)
//                     {
//                         if (bcTypes.isDirichlet(eqIdx))
//                         {
//                             const auto pvIdx = bcTypes.eqToDirichletIndex(eqIdx);
//                             assert(0 <= pvIdx && pvIdx < numEq);
//
//                             const auto& priVars = curElemVolVars[scvI].priVars();
//                             r[scvI.dofIndex()][eqIdx] = priVars[pvIdx] - dirichletValues[pvIdx];
//                             A[scvI.dofIndex()][scvI.dofIndex()][eqIdx][pvIdx] = 1.0;
//                         }
//                     }
//                 }
//             }
//         }
//     }
// private:
//     template<class T = TypeTag>
//     static typename std::enable_if<!GET_PROP_VALUE(T, EnableGridVolumeVariablesCache), VolumeVariables&>::type
//     getVolVarAccess(GridVolumeVariables& gridVolVars, ElementVolumeVariables& elemVolVars, const SubControlVolume& scv)
//     { return elemVolVars[scv]; }
//
//     template<class T = TypeTag>
//     static typename std::enable_if<GET_PROP_VALUE(T, EnableGridVolumeVariablesCache), VolumeVariables&>::type
//     getVolVarAccess(GridVolumeVariables& gridVolVars, ElementVolumeVariables& elemVolVars, const SubControlVolume& scv)
//     { return gridVolVars.volVars(scv.elementIndex(), scv.indexInElement()); }
//
// }; // explicit BoxAssembler with numeric Jacobian

//
// /*!
//  * \ingroup Assembly
//  * \ingroup BoxDiscretization
//  * \brief Box local assembler using analytic (hard coded) derivatives and implicit time discretization
//  */
// template<class TypeTag>
// class BoxLocalAssembler<TypeTag,
//                        DiffMethod::analytic,
//                        /*implicit=*/true>
// {
//     using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
//     using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
//     using ResidualVector = typename GET_PROP_TYPE(TypeTag, NumEqVector);
//     using ElementResidualVector = Dune::BlockVector<typename GET_PROP_TYPE(TypeTag, NumEqVector)>;
//     using ElementBoundaryTypes = typename GET_PROP_TYPE(TypeTag, ElementBoundaryTypes);
//     using Element = typename GET_PROP_TYPE(TypeTag, GridView)::template Codim<0>::Entity;
//     using SolutionVector = typename GET_PROP_TYPE(TypeTag, SolutionVector);
//     using ElementSolutionVector = typename GET_PROP_TYPE(TypeTag, ElementSolutionVector);
//     using ElementVolumeVariables = typename GET_PROP_TYPE(TypeTag, ElementVolumeVariables);
//     using GridVolumeVariables = typename GET_PROP_TYPE(TypeTag, GridVolumeVariables);
//     using VolumeVariables = typename GET_PROP_TYPE(TypeTag, VolumeVariables);
//     using SubControlVolume = typename GET_PROP_TYPE(TypeTag, SubControlVolume);
//     using JacobianMatrix = typename GET_PROP_TYPE(TypeTag, JacobianMatrix);
//
//     enum { numEq = GET_PROP_VALUE(TypeTag, NumEq) };
//
//     static constexpr int dim = GridView::dimension;
//
// public:
//
//     /*!
//      * \brief Computes the derivatives with respect to the given element and adds them
//      *        to the global matrix. The element residual is written into the right hand side.
//      */
//     template<class Assembler>
//     static void assemble(Assembler& assembler, JacobianMatrix& jac, SolutionVector& res,
//                          const Element& element, const SolutionVector& curSol)
//     {
//         assemble_(assembler, jac, res, element, curSol);
//
//     }
//
//     /*!
//      * \brief Computes the derivatives with respect to the given element and adds them
//      *        to the global matrix.
//      */
//     template<class Assembler>
//     static void assemble(Assembler& assembler, JacobianMatrix& jac,
//                          const Element& element, const SolutionVector& curSol)
//     {
//         auto dummyresidual = curSol;
//         assemble_(assembler, jac, dummyresidual, element, curSol);
//     }
//
//     /*!
//      * \brief Assemble the residual only
//      */
//     template<class Assembler>
//     static void assemble(Assembler& assembler, SolutionVector& res,
//                          const Element& element, const SolutionVector& curSol)
//     {
//         assemble_(assembler, res, element, curSol);
//     }
//
// private:
//     /*!
//      * \brief Computes the residual
//      *
//      * \return The element residual at the current solution.
//      */
//     template<class Assembler>
//     static void assemble_(Assembler& assembler, SolutionVector& r,
//                           const Element& element, const SolutionVector& curSol)
//     {
//         // get some references for convenience
//         const auto& problem = assembler.problem();
//         auto& localResidual = assembler.localResidual();
//         auto& gridVariables = assembler.gridVariables();
//
//         // prepare the local views
//         auto fvGeometry = localView(assembler.fvGridGeometry());
//         fvGeometry.bind(element);
//
//         auto curElemVolVars = localView(gridVariables.curGridVolVars());
//         curElemVolVars.bind(element, fvGeometry, curSol);
//
//         auto elemFluxVarsCache = localView(gridVariables.gridFluxVarsCache());
//         elemFluxVarsCache.bind(element, fvGeometry, curElemVolVars);
//
//         const bool isStationary = localResidual.isStationary();
//         auto prevElemVolVars = localView(gridVariables.prevGridVolVars());
//         if (!isStationary)
//             prevElemVolVars.bindElement(element, fvGeometry, localResidual.prevSol());
//
//         // the element boundary types
//         ElementBoundaryTypes elemBcTypes;
//         elemBcTypes.update(problem, element, fvGeometry);
//
//         // the actual element's current residual
//         ElementResidualVector residual(0.0);
//         if (isStationary)
//         {
//             residual = localResidual.eval(problem,
//                                           element,
//                                           fvGeometry,
//                                           curElemVolVars,
//                                           elemBcTypes,
//                                           elemFluxVarsCache);
//         }
//         else
//         {
//             residual = localResidual.eval(problem,
//                                           element,
//                                           fvGeometry,
//                                           prevElemVolVars,
//                                           curElemVolVars,
//                                           elemBcTypes,
//                                           elemFluxVarsCache);
//         }
//
//         for (const auto& scv : scvs(fvGeometry))
//             r[scv.dofIndex()] += residual[scv.indexInElement()];
//
//         // enforce Dirichlet boundaries by setting the residual to (privar - dirichletvalue)
//         if (elemBcTypes.hasDirichlet())
//         {
//             for (const auto& scvI : scvs(fvGeometry))
//             {
//                 const auto bcTypes = elemBcTypes[scvI.indexInElement()];
//                 if (bcTypes.hasDirichlet())
//                 {
//                     const auto dirichletValues = problem.dirichlet(element, scvI);
//
//                     // set the dirichlet conditions in residual and jacobian
//                     for (int eqIdx = 0; eqIdx < numEq; ++eqIdx)
//                     {
//                         if (bcTypes.isDirichlet(eqIdx))
//                         {
//                             const auto pvIdx = bcTypes.eqToDirichletIndex(eqIdx);
//                             assert(0 <= pvIdx && pvIdx < numEq);
//
//                             const auto& priVars = curElemVolVars[scvI].priVars();
//                             r[scvI.dofIndex()][eqIdx] = priVars[pvIdx] - dirichletValues[pvIdx];
//                         }
//                     }
//                 }
//             }
//         }
//     }
//
//     /*!
//      * \brief Computes the derivatives with respect to the given element and adds them
//      *        to the global matrix.
//      *
//      * \return The element residual at the current solution.
//      */
//     template<class Assembler>
//     static void assemble_(Assembler& assembler, JacobianMatrix& A, SolutionVector& r,
//                           const Element& element, const SolutionVector& curSol)
//     {
//         // get some references for convenience
//         const auto& problem = assembler.problem();
//         const auto& fvGridGeometry = assembler.fvGridGeometry();
//         auto& localResidual = assembler.localResidual();
//         auto& gridVariables = assembler.gridVariables();
//
//         // prepare the local views
//         auto fvGeometry = localView(fvGridGeometry);
//         fvGeometry.bind(element);
//
//         auto curElemVolVars = localView(gridVariables.curGridVolVars());
//         curElemVolVars.bind(element, fvGeometry, curSol);
//
//         auto elemFluxVarsCache = localView(gridVariables.gridFluxVarsCache());
//         elemFluxVarsCache.bind(element, fvGeometry, curElemVolVars);
//
//         const bool isStationary = localResidual.isStationary();
//         auto prevElemVolVars = localView(gridVariables.prevGridVolVars());
//         if (!isStationary)
//             prevElemVolVars.bindElement(element, fvGeometry, localResidual.prevSol());
//
//         // check for boundaries on the element
//         ElementBoundaryTypes elemBcTypes;
//         elemBcTypes.update(problem, element, fvGeometry);
//
//         // the actual element's current residual
//         ElementResidualVector residual(0.0);
//         if (isStationary)
//         {
//             residual = localResidual.eval(problem,
//                                           element,
//                                           fvGeometry,
//                                           curElemVolVars,
//                                           elemBcTypes,
//                                           elemFluxVarsCache);
//         }
//         else
//         {
//             residual = localResidual.eval(problem,
//                                           element,
//                                           fvGeometry,
//                                           prevElemVolVars,
//                                           curElemVolVars,
//                                           elemBcTypes,
//                                           elemFluxVarsCache);
//         }
//
//         for (auto&& scv : scvs(fvGeometry))
//             r[scv.dofIndex()] += residual[scv.indexInElement()];
//
//         //////////////////////////////////////////////////////////////////////////////////////////////////
//         //                                                                                              //
//         // Calculate derivatives of all dofs in stencil with respect to the dofs in the element. In the //
//         // neighboring elements we do so by computing the derivatives of the fluxes which depend on the //
//         // actual element. In the actual element we evaluate the derivative of the entire residual.     //
//         //                                                                                              //
//         //////////////////////////////////////////////////////////////////////////////////////////////////
//
//         // calculation of the source and storage derivatives
//         for (const auto& scv : scvs(fvGeometry))
//         {
//             // dof index and corresponding actual pri vars
//             const auto dofIdx = scv.dofIndex();
//             const auto& volVars = curElemVolVars[scv];
//
//             // derivative of this scv residual w.r.t the d.o.f. of the same scv (because of mass lumping)
//             // only if the problem is instationary we add derivative of storage term
//             if (!isStationary)
//                 localResidual.addStorageDerivatives(A[dofIdx][dofIdx],
//                                                     problem,
//                                                     element,
//                                                     fvGeometry,
//                                                     volVars,
//                                                     scv);
//
//             // derivative of this scv residual w.r.t the d.o.f. of the same scv (because of mass lumping)
//             // add source term derivatives
//             localResidual.addSourceDerivatives(A[dofIdx][dofIdx],
//                                                problem,
//                                                element,
//                                                fvGeometry,
//                                                volVars,
//                                                scv);
//         }
//
//         // localJacobian[scvIdx][otherScvIdx][eqIdx][priVarIdx] of the fluxes
//         for (const auto& scvf : scvfs(fvGeometry))
//         {
//             if (!scvf.boundary())
//             {
//                 // add flux term derivatives
//                 localResidual.addFluxDerivatives(A,
//                                                  problem,
//                                                  element,
//                                                  fvGeometry,
//                                                  curElemVolVars,
//                                                  elemFluxVarsCache,
//                                                  scvf);
//             }
//
//             // the boundary gets special treatment to simplify
//             // for the user
//             else
//             {
//                 const auto& insideScv = fvGeometry.scv(scvf.insideScvIdx());
//                 if (elemBcTypes[insideScv.indexInElement()].hasNeumann())
//                 {
//                     // add flux term derivatives
//                     localResidual.addRobinFluxDerivatives(A[insideScv.dofIndex()],
//                                                           problem,
//                                                           element,
//                                                           fvGeometry,
//                                                           curElemVolVars,
//                                                           elemFluxVarsCache,
//                                                           scvf);
//                 }
//             }
//         }
//
//         // enforce Dirichlet boundaries by overwriting partial derivatives with 1 or 0
//         // and set the residual to (privar - dirichletvalue)
//         if (elemBcTypes.hasDirichlet())
//         {
//             for (const auto& scvI : scvs(fvGeometry))
//             {
//                 const auto bcTypes = elemBcTypes[scvI.indexInElement()];
//                 if (bcTypes.hasDirichlet())
//                 {
//                     const auto dirichletValues = problem.dirichlet(element, scvI);
//
//                     // set the dirichlet conditions in residual and jacobian
//                     for (int eqIdx = 0; eqIdx < numEq; ++eqIdx)
//                     {
//                         if (bcTypes.isDirichlet(eqIdx))
//                         {
//                             const auto pvIdx = bcTypes.eqToDirichletIndex(eqIdx);
//                             assert(0 <= pvIdx && pvIdx < numEq);
//
//                             const auto& priVars = curElemVolVars[scvI].priVars();
//                             r[scvI.dofIndex()][eqIdx] = priVars[pvIdx] - dirichletValues[pvIdx];
//                             for (const auto& scvJ : scvs(fvGeometry))
//                             {
//                                 A[scvI.dofIndex()][scvJ.dofIndex()][eqIdx] = 0.0;
//                                 if (scvI.indexInElement() == scvJ.indexInElement())
//                                     A[scvI.dofIndex()][scvI.dofIndex()][eqIdx][pvIdx] = 1.0;
//                             }
//                         }
//                     }
//                 }
//             }
//         }
//
//         // TODO additional dof dependencies
//     }
//
// }; // implicit BoxAssembler with analytic Jacobian
//
//
// /*!
//  * \ingroup Assembly
//  * \ingroup BoxDiscretization
//  * \brief Box local assembler using analytic (hard coded) derivatives and explicit time discretization
//  */
// template<class TypeTag>
// class BoxLocalAssembler<TypeTag,
//                        DiffMethod::analytic,
//                        /*implicit=*/false>
// {
//     using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
//     using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
//     using ResidualVector = typename GET_PROP_TYPE(TypeTag, NumEqVector);
//     using ElementResidualVector = Dune::BlockVector<typename GET_PROP_TYPE(TypeTag, NumEqVector)>;
//     using ElementBoundaryTypes = typename GET_PROP_TYPE(TypeTag, ElementBoundaryTypes);
//     using Element = typename GET_PROP_TYPE(TypeTag, GridView)::template Codim<0>::Entity;
//     using SolutionVector = typename GET_PROP_TYPE(TypeTag, SolutionVector);
//     using ElementSolutionVector = typename GET_PROP_TYPE(TypeTag, ElementSolutionVector);
//     using ElementVolumeVariables = typename GET_PROP_TYPE(TypeTag, ElementVolumeVariables);
//     using GridVolumeVariables = typename GET_PROP_TYPE(TypeTag, GridVolumeVariables);
//     using VolumeVariables = typename GET_PROP_TYPE(TypeTag, VolumeVariables);
//     using SubControlVolume = typename GET_PROP_TYPE(TypeTag, SubControlVolume);
//     using JacobianMatrix = typename GET_PROP_TYPE(TypeTag, JacobianMatrix);
//
//     enum { numEq = GET_PROP_VALUE(TypeTag, NumEq) };
//
//     static constexpr int dim = GridView::dimension;
//
// public:
//
//     /*!
//      * \brief Computes the derivatives with respect to the given element and adds them
//      *        to the global matrix. The element residual is written into the right hand side.
//      */
//     template<class Assembler>
//     static void assemble(Assembler& assembler, JacobianMatrix& jac, SolutionVector& res,
//                          const Element& element, const SolutionVector& curSol)
//     {
//         assemble_(assembler, jac, res, element, curSol);
//
//     }
//
//     /*!
//      * \brief Computes the derivatives with respect to the given element and adds them
//      *        to the global matrix.
//      */
//     template<class Assembler>
//     static void assemble(Assembler& assembler, JacobianMatrix& jac,
//                          const Element& element, const SolutionVector& curSol)
//     {
//         auto dummyresidual = curSol;
//         assemble_(assembler, jac, dummyresidual, element, curSol);
//     }
//
//     /*!
//      * \brief Assemble the residual only
//      */
//     template<class Assembler>
//     static void assemble(Assembler& assembler, SolutionVector& res,
//                          const Element& element, const SolutionVector& curSol)
//     {
//         assemble_(assembler, res, element, curSol);
//     }
//
// private:
//     /*!
//      * \brief Computes the residual
//      *
//      * \return The element residual at the current solution.
//      */
//     template<class Assembler>
//     static void assemble_(Assembler& assembler, SolutionVector& r,
//                           const Element& element, const SolutionVector& curSol)
//     {
//         // get some references for convenience
//         const auto& problem = assembler.problem();
//         auto& localResidual = assembler.localResidual();
//         auto& gridVariables = assembler.gridVariables();
//
//         // using an explicit assembler doesn't make sense for stationary problems
//         if (localResidual.isStationary())
//             DUNE_THROW(Dune::InvalidStateException, "Using explicit jacobian assembler with stationary local residual");
//
//         // prepare the local views
//         auto fvGeometry = localView(assembler.fvGridGeometry());
//         fvGeometry.bind(element);
//
//         auto curElemVolVars = localView(gridVariables.curGridVolVars());
//         curElemVolVars.bindElement(element, fvGeometry, curSol);
//
//         auto prevElemVolVars = localView(gridVariables.prevGridVolVars());
//         prevElemVolVars.bind(element, fvGeometry, localResidual.prevSol());
//
//         auto elemFluxVarsCache = localView(gridVariables.gridFluxVarsCache());
//         elemFluxVarsCache.bind(element, fvGeometry, prevElemVolVars);
//
//         // the element boundary types
//         ElementBoundaryTypes elemBcTypes;
//         elemBcTypes.update(problem, element, fvGeometry);
//
//         // the actual element's current residual
//         auto residual = localResidual.eval(problem,
//                                            element,
//                                            fvGeometry,
//                                            prevElemVolVars,
//                                            elemBcTypes,
//                                            elemFluxVarsCache);
//
//         auto storageResidual = localResidual.evalStorage(problem,
//                                                          element,
//                                                          fvGeometry,
//                                                          prevElemVolVars,
//                                                          curElemVolVars,
//                                                          elemBcTypes,
//                                                          elemFluxVarsCache);
//         residual += storageResidual;
//
//         for (const auto& scv : scvs(fvGeometry))
//             r[scv.dofIndex()] += residual[scv.indexInElement()];
//
//         // enforce Dirichlet boundaries by setting the residual to (privar - dirichletvalue)
//         if (elemBcTypes.hasDirichlet())
//         {
//             for (const auto& scvI : scvs(fvGeometry))
//             {
//                 const auto bcTypes = elemBcTypes[scvI.indexInElement()];
//                 if (bcTypes.hasDirichlet())
//                 {
//                     const auto dirichletValues = problem.dirichlet(element, scvI);
//
//                     // set the dirichlet conditions in residual and jacobian
//                     for (int eqIdx = 0; eqIdx < numEq; ++eqIdx)
//                     {
//                         if (bcTypes.isDirichlet(eqIdx))
//                         {
//                             const auto pvIdx = bcTypes.eqToDirichletIndex(eqIdx);
//                             assert(0 <= pvIdx && pvIdx < numEq);
//
//                             const auto& priVars = curElemVolVars[scvI].priVars();
//                             r[scvI.dofIndex()][eqIdx] = priVars[pvIdx] - dirichletValues[pvIdx];
//                         }
//                     }
//                 }
//             }
//         }
//     }
//
//     /*!
//      * \brief Computes the derivatives with respect to the given element and adds them
//      *        to the global matrix.
//      *
//      * \return The element residual at the current solution.
//      */
//     template<class Assembler>
//     static void assemble_(Assembler& assembler, JacobianMatrix& A, SolutionVector& r,
//                           const Element& element, const SolutionVector& curSol)
//     {
//         // get some references for convenience
//         const auto& problem = assembler.problem();
//         auto& localResidual = assembler.localResidual();
//         auto& gridVariables = assembler.gridVariables();
//
//         // using an explicit assembler doesn't make sense for stationary problems
//         if (localResidual.isStationary())
//             DUNE_THROW(Dune::InvalidStateException, "Using explicit jacobian assembler with stationary local residual");
//
//         // prepare the local views
//         auto fvGeometry = localView(assembler.fvGridGeometry());
//         fvGeometry.bind(element);
//
//         auto curElemVolVars = localView(gridVariables.curGridVolVars());
//         curElemVolVars.bindElement(element, fvGeometry, curSol);
//
//         auto prevElemVolVars = localView(gridVariables.prevGridVolVars());
//         prevElemVolVars.bind(element, fvGeometry, localResidual.prevSol());
//
//         auto elemFluxVarsCache = localView(gridVariables.gridFluxVarsCache());
//         elemFluxVarsCache.bind(element, fvGeometry, prevElemVolVars);
//
//         // the element boundary types
//         ElementBoundaryTypes elemBcTypes;
//         elemBcTypes.update(problem, element, fvGeometry);
//
//         // the actual element's current residual
//         auto residual = localResidual.eval(problem,
//                                            element,
//                                            fvGeometry,
//                                            prevElemVolVars,
//                                            elemBcTypes,
//                                            elemFluxVarsCache);
//
//         auto storageResidual = localResidual.evalStorage(problem,
//                                                          element,
//                                                          fvGeometry,
//                                                          prevElemVolVars,
//                                                          curElemVolVars,
//                                                          elemBcTypes,
//                                                          elemFluxVarsCache);
//         residual += storageResidual;
//
//         for (const auto& scv : scvs(fvGeometry))
//             r[scv.dofIndex()] += residual[scv.indexInElement()];
//
//         //////////////////////////////////////////////////////////////////////////////////////////////////
//         //                                                                                              //
//         // Calculate derivatives of all dofs in stencil with respect to the dofs in the element. In the //
//         // neighboring elements we do so by computing the derivatives of the fluxes which depend on the //
//         // actual element. In the actual element we evaluate the derivative of the entire residual.     //
//         //                                                                                              //
//         //////////////////////////////////////////////////////////////////////////////////////////////////
//
//         // calculation of the source and storage derivatives
//         for (const auto& scv : scvs(fvGeometry))
//         {
//             // dof index and corresponding actual pri vars
//             const auto dofIdx = scv.dofIndex();
//             const auto& volVars = curElemVolVars[scv];
//
//             // derivative of this scv residual w.r.t the d.o.f. of the same scv (because of mass lumping)
//             // only if the problem is instationary we add derivative of storage term
//             localResidual.addStorageDerivatives(A[dofIdx][dofIdx],
//                                                 problem,
//                                                 element,
//                                                 fvGeometry,
//                                                 volVars,
//                                                 scv);
//         }
//
//         // enforce Dirichlet boundaries by overwriting partial derivatives with 1 or 0
//         // and set the residual to (privar - dirichletvalue)
//         if (elemBcTypes.hasDirichlet())
//         {
//             for (const auto& scvI : scvs(fvGeometry))
//             {
//                 const auto bcTypes = elemBcTypes[scvI.indexInElement()];
//                 if (bcTypes.hasDirichlet())
//                 {
//                     const auto dirichletValues = problem.dirichlet(element, scvI);
//
//                     // set the dirichlet conditions in residual and jacobian
//                     for (int eqIdx = 0; eqIdx < numEq; ++eqIdx)
//                     {
//                         if (bcTypes.isDirichlet(eqIdx))
//                         {
//                             const auto pvIdx = bcTypes.eqToDirichletIndex(eqIdx);
//                             assert(0 <= pvIdx && pvIdx < numEq);
//
//                             const auto& priVars = curElemVolVars[scvI].priVars();
//                             r[scvI.dofIndex()][eqIdx] = priVars[pvIdx] - dirichletValues[pvIdx];
//                             A[scvI.dofIndex()][scvI.dofIndex()][eqIdx][pvIdx] = 1.0;
//                         }
//                     }
//                 }
//             }
//         }
//
//         // TODO additional dof dependencies
//     }
//
// }; // explicit BoxAssembler with analytic Jacobian

} // end namespace Dumux

#endif
