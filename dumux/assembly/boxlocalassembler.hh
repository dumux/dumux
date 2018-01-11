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

        this->asImp_().evalDirichletBoundaries(jac, res);
    }

    /*!
     * \brief Computes the derivatives with respect to the given element and adds them
     *        to the global matrix.
     */
    void assembleJacobian(JacobianMatrix& jac, GridVariables& gridVariables)
    {
        this->asImp_().bindLocalViews();
        this->asImp_().assembleJacobianAndResidualImpl(jac, gridVariables); // forward to the internal implementation
        // TODO: Dirichlet required here?
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

        // TODO: Dirichlet required here?
    }


    auto evalLocalFluxAndSourceResidual(const ElementVolumeVariables& elemVolVars) const
    {
        return this->evalLocalFluxAndSourceResidual(elemVolVars);
    }

    auto evalLocalStorageResidual() const
    {
        return this->evalLocalStorageResidual();
    }

    void evalDirichletBoundaries(JacobianMatrix& jac, SolutionVector& res)
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
    using ElementSolutionVector = typename GET_PROP_TYPE(TypeTag, ElementSolutionVector);
    using GridVariables = typename GET_PROP_TYPE(TypeTag, GridVariables);
    using VolumeVariables = typename GET_PROP_TYPE(TypeTag, VolumeVariables);
    using JacobianMatrix = typename GET_PROP_TYPE(TypeTag, JacobianMatrix);

    enum { numEq = GET_PROP_VALUE(TypeTag, NumEq) };
    enum { dim = GET_PROP_TYPE(TypeTag, GridView)::dimension };

    static constexpr bool enableGridFluxVarsCache = GET_PROP_VALUE(TypeTag, EnableGridFluxVariablesCache);

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
    using ElementSolutionVector = typename GET_PROP_TYPE(TypeTag, ElementSolutionVector);
    using GridVariables = typename GET_PROP_TYPE(TypeTag, GridVariables);
    using VolumeVariables = typename GET_PROP_TYPE(TypeTag, VolumeVariables);
    using JacobianMatrix = typename GET_PROP_TYPE(TypeTag, JacobianMatrix);

    enum { numEq = GET_PROP_VALUE(TypeTag, NumEq) };
    enum { dim = GET_PROP_TYPE(TypeTag, GridView)::dimension };

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

/*!
 * \ingroup Assembly
 * \ingroup BoxDiscretization
 * \brief Box local assembler using analytic differentiation and implicit time discretization
 */
template<class TypeTag, class Assembler>
class BoxLocalAssembler<TypeTag, Assembler, DiffMethod::analytic, /*implicit=*/true>
: public BoxLocalAssemblerBase<TypeTag, Assembler,
                              BoxLocalAssembler<TypeTag, Assembler, DiffMethod::analytic, true>, true>
{
    using ThisType = BoxLocalAssembler<TypeTag, Assembler, DiffMethod::analytic, true>;
    using ParentType = BoxLocalAssemblerBase<TypeTag, Assembler, ThisType, true>;
    using GridVariables = typename GET_PROP_TYPE(TypeTag, GridVariables);
    using JacobianMatrix = typename GET_PROP_TYPE(TypeTag, JacobianMatrix);

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
        const auto& problem = this->problem();
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

        // calculation of the source and storage derivatives
        for (const auto& scv : scvs(fvGeometry))
        {
            // dof index and corresponding actual pri vars
            const auto dofIdx = scv.dofIndex();
            const auto& volVars = curElemVolVars[scv];

            // derivative of this scv residual w.r.t the d.o.f. of the same scv (because of mass lumping)
            // only if the problem is instationary we add derivative of storage term
            if (!this->assembler().isStationaryProblem())
                this->localResidual().addStorageDerivatives(A[dofIdx][dofIdx],
                                                            problem,
                                                            element,
                                                            fvGeometry,
                                                            volVars,
                                                            scv);

            // derivative of this scv residual w.r.t the d.o.f. of the same scv (because of mass lumping)
            // add source term derivatives
            this->localResidual().addSourceDerivatives(A[dofIdx][dofIdx],
                                                       problem,
                                                       element,
                                                       fvGeometry,
                                                       volVars,
                                                       scv);
        }

        // localJacobian[scvIdx][otherScvIdx][eqIdx][priVarIdx] of the fluxes
        for (const auto& scvf : scvfs(fvGeometry))
        {
            if (!scvf.boundary())
            {
                // add flux term derivatives
                this->localResidual().addFluxDerivatives(A,
                                                         problem,
                                                         element,
                                                         fvGeometry,
                                                         curElemVolVars,
                                                         elemFluxVarsCache,
                                                         scvf);
            }

            // the boundary gets special treatment to simplify
            // for the user
            else
            {
                const auto& insideScv = fvGeometry.scv(scvf.insideScvIdx());
                if (this->elemBcTypes()[insideScv.indexInElement()].hasNeumann())
                {
                    // add flux term derivatives
                    this->localResidual().addRobinFluxDerivatives(A[insideScv.dofIndex()],
                                                                  problem,
                                                                  element,
                                                                  fvGeometry,
                                                                  curElemVolVars,
                                                                  elemFluxVarsCache,
                                                                  scvf);
                }
            }
        }

        return origResiduals;
    }

}; // implicit BoxAssembler with analytic Jacobian

/*!
 * \ingroup Assembly
 * \ingroup BoxDiscretization
 * \brief Box local assembler using analytic differentiation and explicit time discretization
 */
template<class TypeTag, class Assembler>
class BoxLocalAssembler<TypeTag, Assembler, DiffMethod::analytic, /*implicit=*/false>
: public BoxLocalAssemblerBase<TypeTag, Assembler,
                              BoxLocalAssembler<TypeTag, Assembler, DiffMethod::analytic, false>, false>
{
    using ThisType = BoxLocalAssembler<TypeTag, Assembler, DiffMethod::analytic, false>;
    using ParentType = BoxLocalAssemblerBase<TypeTag, Assembler, ThisType, false>;
    using GridVariables = typename GET_PROP_TYPE(TypeTag, GridVariables);
    using JacobianMatrix = typename GET_PROP_TYPE(TypeTag, JacobianMatrix);

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
        const auto& problem = this->problem();
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

        // calculation of the source and storage derivatives
        for (const auto& scv : scvs(fvGeometry))
        {
            // dof index and corresponding actual pri vars
            const auto dofIdx = scv.dofIndex();
            const auto& volVars = curElemVolVars[scv];

            // derivative of this scv residual w.r.t the d.o.f. of the same scv (because of mass lumping)
            // only if the problem is instationary we add derivative of storage term
            if (!this->assembler().isStationaryProblem())
                this->localResidual().addStorageDerivatives(A[dofIdx][dofIdx],
                                                            problem,
                                                            element,
                                                            fvGeometry,
                                                            volVars,
                                                            scv);
            else
                DUNE_THROW(Dune::InvalidStateException, "Using explicit jacobian assembler with stationary local residual");

        }

        return origResiduals;
    }

}; // explicit BoxAssembler with analytic Jacobian

} // end namespace Dumux

#endif
