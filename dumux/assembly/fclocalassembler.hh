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
 * \brief An assembler for Jacobian and residual contribution per element (face-centered methods)
 */
#ifndef DUMUX_FC_LOCAL_ASSEMBLER_HH
#define DUMUX_FC_LOCAL_ASSEMBLER_HH

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
#include <dumux/discretization/facecentered/staggered/elementsolution.hh>

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
class FaceCenteredLocalAssemblerBase : public FVLocalAssemblerBase<TypeTag, Assembler, Implementation, implicit>
{
    using ParentType = FVLocalAssemblerBase<TypeTag, Assembler, Implementation, implicit>;
    using GridView = typename GetPropType<TypeTag, Properties::GridGeometry>::GridView;
    using JacobianMatrix = GetPropType<TypeTag, Properties::JacobianMatrix>;
    using GridVariables = GetPropType<TypeTag, Properties::GridVariables>;
    using SolutionVector = GetPropType<TypeTag, Properties::SolutionVector>;
    using PrimaryVariables = GetPropType<TypeTag, Properties::PrimaryVariables>; // TODO
    using ElementVolumeVariables = typename GetPropType<TypeTag, Properties::GridVolumeVariables>::LocalView;
    using NumEqVector = GetPropType<TypeTag, Properties::NumEqVector>;

    static_assert(!Assembler::Problem::enableInternalDirichletConstraints(), "Internal Dirichlet constraints are currently not implemented for cc-methods!");

    static constexpr auto numEq = GetPropType<TypeTag, Properties::ModelTraits>::numEq();

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
    template <class PartialReassembler = DefaultPartialReassembler>
    void assembleJacobianAndResidual(JacobianMatrix& jac, SolutionVector& res, GridVariables& gridVariables,
                                     const PartialReassembler* partialReassembler)
    {
        this->asImp_().bindLocalViews();
        const auto eIdxGlobal = this->assembler().gridGeometry().elementMapper().index(this->element());
        if (partialReassembler
            && partialReassembler->elementColor(eIdxGlobal) == EntityColor::green)
        {
            const auto residual = this->asImp_().evalLocalResidual(); // forward to the internal implementation
            for (const auto& scv : scvs(this->fvGeometry()))
                res[scv.dofIndex()] += residual[scv.localDofIndex()];
        }
        else if (!this->elementIsGhost())
        {
            const auto residual = this->asImp_().assembleJacobianAndResidualImpl(jac, gridVariables, partialReassembler); // forward to the internal implementation
            for (const auto& scv : scvs(this->fvGeometry()))
                res[scv.dofIndex()] += residual[scv.localDofIndex()];
        }
        else
        {
            // using GridGeometry = typename GridVariables::GridGeometry;
            // using GridView = typename GridGeometry::GridView;
            // static constexpr auto dim = GridView::dimension;

            int numFacesLocal = this->element().subEntities(2);

            for (int i = 0; i < numFacesLocal; ++i)
            {
                auto face = this->element().template subEntity<2>(i);

                if (face.partitionType() == Dune::InteriorEntity ||
                    face.partitionType() == Dune::BorderEntity)
                {
                    // do not change the non-ghost vertices
                    continue;
                }

                // set main diagonal entries for the vertex
                int globalI = this->assembler().gridGeometry().gridView().indexSet().subIndex(face, i, 2);

                typedef typename JacobianMatrix::block_type BlockType;
                BlockType &J = jac[globalI][globalI];
                for (int j = 0; j < BlockType::rows; ++j)
                    J[j][j] = 1.0;

                // set residual for the face
                res[globalI] = 0;
            }
        }

        auto applyDirichlet = [&] (const auto& scvI,
                                   const auto& dirichletValues,
                                   const auto eqIdx,
                                   const auto pvIdx)
        {
            res[scvI.dofIndex()][eqIdx] = this->curElemVolVars()[scvI].priVars()[pvIdx] - dirichletValues[pvIdx];

            auto& row = jac[scvI.dofIndex()];
            for (auto col = row.begin(); col != row.end(); ++col)
                row[col.index()][eqIdx] = 0.0;

            jac[scvI.dofIndex()][scvI.dofIndex()][eqIdx][pvIdx] = 1.0;


            // if a periodic dof has Dirichlet values also apply the same Dirichlet values to the other dof TODO periodic
            // if (this->assembler().gridGeometry().dofOnPeriodicBoundary(scvI.dofIndex()))
            // {
            //     const auto periodicDof = this->assembler().gridGeometry().periodicallyMappedDof(scvI.dofIndex());
            //     res[periodicDof][eqIdx] = this->curElemVolVars()[scvI].priVars()[pvIdx] - dirichletValues[pvIdx];
            //     const auto end = jac[periodicDof].end();
            //     for (auto it = jac[periodicDof].begin(); it != end; ++it)
            //         (*it) = periodicDof != it.index() ? 0.0 : 1.0;
            // }
        };

        this->asImp_().enforceDirichletConstraints(applyDirichlet);
    }

    /*!
     * \brief Computes the derivatives with respect to the given element and adds them
     *        to the global matrix.
     */
    void assembleJacobian(JacobianMatrix& jac, GridVariables& gridVariables)
    {
        this->asImp_().bindLocalViews();
        this->asImp_().assembleJacobianAndResidualImpl(jac, gridVariables); // forward to the internal implementation

        auto applyDirichlet = [&] (const auto& scvI,
                                   const auto& dirichletValues,
                                   const auto eqIdx,
                                   const auto pvIdx)
        {
            auto& row = jac[scvI.dofIndex()];
            for (auto col = row.begin(); col != row.end(); ++col)
                row[col.index()][eqIdx] = 0.0;

            jac[scvI.dofIndex()][scvI.dofIndex()][eqIdx][pvIdx] = 1.0;
        };

        this->asImp_().enforceDirichletConstraints(applyDirichlet);
    }

    /*!
     * \brief Assemble the residual only
     */
    void assembleResidual(SolutionVector& res)
    {
        this->asImp_().bindLocalViews();
        const auto residual = this->evalLocalResidual();

        for (const auto& scv : scvs(this->fvGeometry()))
            res[scv.dofIndex()] += residual[scv.localDofIndex()];

        auto applyDirichlet = [&] (const auto& scvI,
                                   const auto& dirichletValues,
                                   const auto eqIdx,
                                   const auto pvIdx)
        {
            res[scvI.dofIndex()][eqIdx] = this->curElemVolVars()[scvI].priVars()[pvIdx] - dirichletValues[pvIdx];
        };

        this->asImp_().enforceDirichletConstraints(applyDirichlet);
    }

        //! Enforce Dirichlet constraints
    template<typename ApplyFunction>
    void enforceDirichletConstraints(const ApplyFunction& applyDirichlet)
    {
        // enforce Dirichlet boundary conditions
        this->asImp_().evalDirichletBoundaries(applyDirichlet);
        // take care of internal Dirichlet constraints (if enabled)
        this->asImp_().enforceInternalDirichletConstraints(applyDirichlet);
    }

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
            for (const auto& scvf : scvfs(this->fvGeometry()))
            {
                if (scvf.isFrontal() && scvf.boundary())
                {
                    const auto bcTypes = this->elemBcTypes()[scvf.localIndex()];
                    if (bcTypes.hasDirichlet())
                    {
                        const auto& scv = this->fvGeometry().scv(scvf.insideScvIdx());
                        const PrimaryVariables dirichletValues(this->problem().dirichlet(this->element(), scvf)[scv.directionIndex()]); // TODO

                        // set the Dirichlet conditions in residual and jacobian
                        for (int eqIdx = 0; eqIdx < numEq; ++eqIdx)
                        {
                            if (bcTypes.isDirichlet(eqIdx))
                            {
                                const auto pvIdx = bcTypes.eqToDirichletIndex(eqIdx);
                                assert(0 <= pvIdx && pvIdx < numEq);
                                applyDirichlet(scv, dirichletValues, eqIdx, pvIdx);
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
 * \ingroup CCDiscretization
 * \brief An assembler for Jacobian and residual contribution per element (Face-centered methods)
 * \tparam TypeTag The TypeTag
 * \tparam diffMethod The differentiation method to residual compute derivatives
 * \tparam implicit Specifies whether the time discretization is implicit or not not (i.e. explicit)
 */
template<class TypeTag, class Assembler, DiffMethod diffMethod = DiffMethod::numeric, bool implicit = true>
class FaceCenteredLocalAssembler;

/*!
 * \ingroup Assembly
 * \ingroup CCDiscretization
 * \brief Face-centered scheme local assembler using numeric differentiation and implicit time discretization
 */
template<class TypeTag, class Assembler>
class FaceCenteredLocalAssembler<TypeTag, Assembler, DiffMethod::numeric, /*implicit=*/true>
: public FaceCenteredLocalAssemblerBase<TypeTag, Assembler,
                                        FaceCenteredLocalAssembler<TypeTag, Assembler, DiffMethod::numeric, true>, true >
{
    using ThisType = FaceCenteredLocalAssembler<TypeTag, Assembler, DiffMethod::numeric, true>;
    using ParentType = FaceCenteredLocalAssemblerBase<TypeTag, Assembler, ThisType, true>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using NumEqVector = GetPropType<TypeTag, Properties::NumEqVector>;
    using Element = typename GetPropType<TypeTag, Properties::GridGeometry>::GridView::template Codim<0>::Entity;
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using FVElementGeometry = typename GridGeometry::LocalView;
    using SubControlVolume = typename FVElementGeometry::SubControlVolume;
    using SubControlVolumeFace = typename FVElementGeometry::SubControlVolumeFace;
    using GridVariables = GetPropType<TypeTag, Properties::GridVariables>;
    using JacobianMatrix = GetPropType<TypeTag, Properties::JacobianMatrix>;
    using LocalResidual = GetPropType<TypeTag, Properties::LocalResidual>;
    using ElementResidualVector = typename LocalResidual::ElementResidualVector;
    using VolumeVariables = GetPropType<TypeTag, Properties::VolumeVariables>;

    static constexpr auto numEq = GetPropType<TypeTag, Properties::ModelTraits>::numEq();
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
    template <class PartialReassembler = DefaultPartialReassembler>
    ElementResidualVector assembleJacobianAndResidualImpl(JacobianMatrix& A, GridVariables& gridVariables,
                                                          const PartialReassembler* partialReassembler = nullptr)
    {
        // get some aliases for convenience
        const auto& problem = this->problem();
        const auto& element = this->element();
        const auto& fvGeometry = this->fvGeometry();
        const auto& curSol = this->curSol();
        auto&& curElemVolVars = this->curElemVolVars();

        // get the vector of the actual element residuals
        const auto origResiduals = this->evalLocalResidual();

        //////////////////////////////////////////////////////////////////////////////////////////////////
        //                                                                                              //
        // Calculate derivatives of all dofs in stencil with respect to the dofs in the element. In the //
        // neighboring elements we do so by computing the derivatives of the fluxes which depend on the //
        // actual element. In the actual element we evaluate the derivative of the entire residual.     //
        //                                                                                              //
        //////////////////////////////////////////////////////////////////////////////////////////////////

        // create the element solution
        auto elemSol = elementSolution(element, curSol, fvGeometry.gridGeometry());

        // create the vector storing the partial derivatives
        ElementResidualVector partialDerivs(element.subEntities(2));

        // calculation of the derivatives
        for (auto&& scv : scvs(fvGeometry))
        {
            // NumEqVector residual(0.0);

            // dof index and corresponding actual pri vars
            const auto dofIdx = scv.dofIndex();
            const auto scvIdx = scv.index();
            auto& curVolVars = this->getVolVarAccess(gridVariables.curGridVolVars(), curElemVolVars, scv);
            const VolumeVariables origVolVars(curVolVars);

            auto evalSource = [&](ElementResidualVector& residual, const SubControlVolume& scv)
            {
                this->localResidual().evalSource(residual, problem, element, fvGeometry, curElemVolVars, scv);
            };

            auto evalStorage = [&](ElementResidualVector& residual, const SubControlVolume& scv)
            {
                this->localResidual().evalStorage(residual, problem, element, fvGeometry, this->prevElemVolVars(), curElemVolVars, scv);
            };

            auto evalFlux = [&](ElementResidualVector& residual, const SubControlVolumeFace& scvf)
            {
                this->localResidual().evalFlux(residual, problem, element, fvGeometry, curElemVolVars, this->elemBcTypes(), this->elemFluxVarsCache(), scvf);
            };

            // derivative w.r.t. own DOF
            for (int pvIdx = 0; pvIdx < numEq; pvIdx++)
            {
                partialDerivs = 0.0;

                auto evalResiduals = [&](Scalar priVar)
                {
                    // update the volume variables and compute element residual
                    elemSol[scv.localDofIndex()][pvIdx] = priVar;
                    curVolVars.update(elemSol, problem, element, scv);

                    ElementResidualVector residual(element.subEntities(2));
                    residual = 0;

                    evalSource(residual, scv);

                    if (!this->assembler().isStationaryProblem())
                        evalStorage(residual, scv);

                    for (const auto& scvf : scvfs(fvGeometry)) // TODO allow iteration over scvfs of scv
                    {
                        if (scvf.insideScvIdx() == scvIdx)
                            evalFlux(residual, scvf);
                    }

                    return residual;
                };

                // derive the residuals numerically
                static const NumericEpsilon<Scalar, numEq> eps_{this->problem().paramGroup()};
                static const int numDiffMethod = getParamFromGroup<int>(this->problem().paramGroup(), "Assembly.NumericDifferenceMethod");
                NumericDifferentiation::partialDerivative(evalResiduals, elemSol[scv.localDofIndex()][pvIdx], partialDerivs, origResiduals,
                                                          eps_(elemSol[scv.localDofIndex()][pvIdx], pvIdx), numDiffMethod);

                for (int eqIdx = 0; eqIdx < numEq; eqIdx++)
                {
                    // A[i][col][eqIdx][pvIdx] is the rate of change of
                    // the residual of equation 'eqIdx' at dof 'i'
                    // depending on the primary variable 'pvIdx' at dof
                    // 'col'.
                    A[dofIdx][dofIdx][eqIdx][pvIdx] += partialDerivs[scv.localDofIndex()][eqIdx];
                }

                // restore the original state of the scv's volume variables
                curVolVars = origVolVars;

                // restore the original element solution
                elemSol[scv.localDofIndex()][pvIdx] = curSol[scv.dofIndex()][pvIdx];
                // TODO additional dof dependencies

            }

            // derivative w.r.t. other DOFs
            const auto& otherScvIndices = fvGeometry.gridGeometry().connectivityMap()[scv.index()];
            for (const auto globalJ : otherScvIndices)
            {
                ElementResidualVector partialDerivsFluxOnly(element.subEntities(2));
                partialDerivsFluxOnly = 0;

                const auto scvJ = fvGeometry.scv(globalJ);
                auto& curOtherVolVars = this->getVolVarAccess(gridVariables.curGridVolVars(), curElemVolVars, scvJ);
                const VolumeVariables origOtherVolVars(curOtherVolVars);

                const auto& otherElement = fvGeometry.gridGeometry().element(scvJ.elementIndex());
                auto otherElemSol = elementSolution(otherElement, curSol, fvGeometry.gridGeometry()); // TODO allow selective creation of elemsol (for one scv)

                for (int pvIdx = 0; pvIdx < numEq; pvIdx++)
                {
                    partialDerivsFluxOnly = 0.0;

                    auto evalResiduals = [&](Scalar priVar)
                    {
                        // update the volume variables and compute element residual
                        otherElemSol[scvJ.localDofIndex()][pvIdx] = priVar;
                        curOtherVolVars.update(otherElemSol, problem, otherElement, scvJ);

                        ElementResidualVector residual(element.subEntities(2));
                        residual = 0;

                        for (const auto& scvf : scvfs(fvGeometry))
                        {
                            if (scvf.insideScvIdx() == scvIdx)
                            {
                                if (scvf.outsideScvIdx() == scvJ.index())
                                {
                                    evalFlux(residual, scvf);
                                    return residual;
                                }

                                if (scvf.isLateral())
                                {
                                    const auto& orthogonalScvf = fvGeometry.scvfWithCommonEntity(scvf);
                                    if (orthogonalScvf.insideScvIdx() == scvJ.index() || orthogonalScvf.outsideScvIdx() == scvJ.index())
                                    {
                                        evalFlux(residual, scvf);
                                        return residual;
                                    }
                                }
                            }
                        }

                        DUNE_THROW(Dune::InvalidStateException, "No scvf found");
                    };

                    // get original fluxes
                    const auto origFluxResidual = [&]()
                    {
                        ElementResidualVector result(element.subEntities(2));
                        result = 0;

                        for (const auto& scvf : scvfs(fvGeometry))
                        {
                            if (scvf.insideScvIdx() == scvIdx)
                            {
                                if (scvf.outsideScvIdx() == scvJ.index())
                                {
                                    evalFlux(result, scvf);
                                    return result;
                                }

                                if (scvf.isLateral())
                                {
                                    const auto& orthogonalScvf = fvGeometry.scvfWithCommonEntity(scvf);
                                    if (orthogonalScvf.insideScvIdx() == scvJ.index() || orthogonalScvf.outsideScvIdx() == scvJ.index())
                                    {
                                        evalFlux(result, scvf);
                                        return result;
                                    }
                                }
                            }
                        }
                        DUNE_THROW(Dune::InvalidStateException, "No scvf found");
                    }();

                    // derive the residuals numerically
                    static const NumericEpsilon<Scalar, numEq> eps_{this->problem().paramGroup()};
                    static const int numDiffMethod = getParamFromGroup<int>(this->problem().paramGroup(), "Assembly.NumericDifferenceMethod");
                    NumericDifferentiation::partialDerivative(evalResiduals, otherElemSol[scvJ.localDofIndex()][pvIdx], partialDerivsFluxOnly, origFluxResidual,
                                                            eps_(otherElemSol[scvJ.localDofIndex()][pvIdx], pvIdx), numDiffMethod);

                    for (int eqIdx = 0; eqIdx < numEq; eqIdx++)
                    {
                        // A[i][col][eqIdx][pvIdx] is the rate of change of
                        // the residual of equation 'eqIdx' at dof 'i'
                        // depending on the primary variable 'pvIdx' at dof
                        // 'col'.
                        A[dofIdx][scvJ.dofIndex()][eqIdx][pvIdx] += partialDerivsFluxOnly[scv.localDofIndex()][eqIdx];
                    }

                    // restore the original state of the scv's volume variables
                    curOtherVolVars = origOtherVolVars;

                    // restore the original element solution
                    otherElemSol[scvJ.localDofIndex()][pvIdx] = curSol[scvJ.dofIndex()][pvIdx];
                    // TODO additional dof dependencies
                }
            }
        }
        return origResiduals;
    }
};


/*!
 * \ingroup Assembly
 * \ingroup CCDiscretization
 * \brief Cell-centered scheme local assembler using numeric differentiation and explicit time discretization
 */
template<class TypeTag, class Assembler>
class FaceCenteredLocalAssembler<TypeTag, Assembler, DiffMethod::numeric, /*implicit=*/false>
: public FaceCenteredLocalAssemblerBase<TypeTag, Assembler,
            FaceCenteredLocalAssembler<TypeTag, Assembler, DiffMethod::numeric, false>, false>
{
    using ThisType = FaceCenteredLocalAssembler<TypeTag, Assembler, DiffMethod::numeric, false>;
    using ParentType = FaceCenteredLocalAssemblerBase<TypeTag, Assembler, ThisType, false>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using NumEqVector = GetPropType<TypeTag, Properties::NumEqVector>;
    using Element = typename GetPropType<TypeTag, Properties::GridGeometry>::GridView::template Codim<0>::Entity;
    using GridVariables = GetPropType<TypeTag, Properties::GridVariables>;
    using JacobianMatrix = GetPropType<TypeTag, Properties::JacobianMatrix>;
    using LocalResidual = GetPropType<TypeTag, Properties::LocalResidual>;
    using ElementResidualVector = typename LocalResidual::ElementResidualVector;
    using VolumeVariables = GetPropType<TypeTag, Properties::VolumeVariables>;

    enum { numEq = GetPropType<TypeTag, Properties::ModelTraits>::numEq() };

public:
    using ParentType::ParentType;

    /*!
     * \brief Computes the derivatives with respect to the given element and adds them
     *        to the global matrix.
     *
     * \return The element residual at the current solution.
     */
    template <class PartialReassembler = DefaultPartialReassembler>
    ElementResidualVector assembleJacobianAndResidualImpl(JacobianMatrix& A, GridVariables& gridVariables,
                                                          const PartialReassembler* partialReassembler = nullptr)
    {
        if (partialReassembler)
            DUNE_THROW(Dune::NotImplemented, "partial reassembly for explicit time discretization");

        // get some aliases for convenience
        const auto& element = this->element();
        const auto& fvGeometry = this->fvGeometry();
        const auto& curSol = this->curSol();
        auto&& curElemVolVars = this->curElemVolVars();

        // get the vecor of the acutal element residuals
        const auto origResiduals = this->evalLocalResidual();
        const auto origStorageResiduals = this->evalLocalStorageResidual();

        //////////////////////////////////////////////////////////////////////////////////////////////////
        //                                                                                              //
        // Calculate derivatives of all dofs in stencil with respect to the dofs in the element. In the //
        // neighboring elements we do so by computing the derivatives of the fluxes which depend on the //
        // actual element. In the actual element we evaluate the derivative of the entire residual.     //
        //                                                                                              //
        //////////////////////////////////////////////////////////////////////////////////////////////////

        // create the element solution
        auto elemSol = elementSolution(element, curSol, fvGeometry.gridGeometry());

        // create the vector storing the partial derivatives
        ElementResidualVector partialDerivs(element.subEntities(2));

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
                    elemSol[scv.localDofIndex()][pvIdx] = priVar;
                    curVolVars.update(elemSol, this->problem(), element, scv);
                    return this->evalLocalStorageResidual();
                };

                // derive the residuals numerically
                static const NumericEpsilon<Scalar, numEq> eps_{this->problem().paramGroup()};
                static const int numDiffMethod = getParamFromGroup<int>(this->problem().paramGroup(), "Assembly.NumericDifferenceMethod");
                NumericDifferentiation::partialDerivative(evalStorage, elemSol[scv.localDofIndex()][pvIdx], partialDerivs, origStorageResiduals,
                                                          eps_(elemSol[scv.localDofIndex()][pvIdx], pvIdx), numDiffMethod);

                // update the global stiffness matrix with the current partial derivatives
                for (int eqIdx = 0; eqIdx < numEq; eqIdx++)
                {
                    // A[i][col][eqIdx][pvIdx] is the rate of change of
                    // the residual of equation 'eqIdx' at dof 'i'
                    // depending on the primary variable 'pvIdx' at dof
                    // 'col'.
                    A[dofIdx][dofIdx][eqIdx][pvIdx] += partialDerivs[scv.localDofIndex()][eqIdx];
                }

                // restore the original state of the scv's volume variables
                curVolVars = origVolVars;

                // restore the original element solution
                elemSol[scv.localDofIndex()][pvIdx] = curSol[scv.dofIndex()][pvIdx];
                // TODO additional dof dependencies
            }
        }
        return origResiduals;
    }
};

/*!
 * \ingroup Assembly
 * \ingroup CCDiscretization
 * \brief Cell-centered scheme local assembler using analytic (hand-coded) differentiation and implicit time discretization
 */
template<class TypeTag, class Assembler>
class FaceCenteredLocalAssembler<TypeTag, Assembler, DiffMethod::analytic, /*implicit=*/true>
: public FaceCenteredLocalAssemblerBase<TypeTag, Assembler,
            FaceCenteredLocalAssembler<TypeTag, Assembler, DiffMethod::analytic, true>, true>
{
    using ThisType = FaceCenteredLocalAssembler<TypeTag, Assembler, DiffMethod::analytic, true>;
    using ParentType = FaceCenteredLocalAssemblerBase<TypeTag, Assembler, ThisType, true>;
    using NumEqVector = GetPropType<TypeTag, Properties::NumEqVector>;
    using JacobianMatrix = GetPropType<TypeTag, Properties::JacobianMatrix>;
    using GridVariables = GetPropType<TypeTag, Properties::GridVariables>;
    using LocalResidual = GetPropType<TypeTag, Properties::LocalResidual>;
    using ElementResidualVector = typename LocalResidual::ElementResidualVector;
    using VolumeVariables = GetPropType<TypeTag, Properties::VolumeVariables>;

    enum { numEq = GetPropType<TypeTag, Properties::ModelTraits>::numEq() };

public:
    using ParentType::ParentType;

    /*!
     * \brief Computes the derivatives with respect to the given element and adds them
     *        to the global matrix.
     *
     * \return The element residual at the current solution.
     */
    template <class PartialReassembler = DefaultPartialReassembler>
    ElementResidualVector assembleJacobianAndResidualImpl(JacobianMatrix& A, GridVariables& gridVariables,
                                                          const PartialReassembler* partialReassembler = nullptr)
    {
        if (partialReassembler)
            DUNE_THROW(Dune::NotImplemented, "partial reassembly for analytic differentiation");

        // get some aliases for convenience
        const auto& element = this->element();
        const auto& fvGeometry = this->fvGeometry();
        const auto& problem = this->problem();
        const auto& curElemVolVars = this->curElemVolVars();
        const auto& elemFluxVarsCache = this->elemFluxVarsCache();

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
            // TODO if e.g. porosity depends on all dofs in the element, we would have off-diagonal matrix entries!?
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
                if (this->elemBcTypes()[insideScv.localDofIndex()].hasNeumann())
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
};

/*!
 * \ingroup Assembly
 * \ingroup CCDiscretization
 * \brief Cell-centered scheme local assembler using analytic (hand-coded) differentiation and explicit time discretization
 */
template<class TypeTag, class Assembler>
class FaceCenteredLocalAssembler<TypeTag, Assembler, DiffMethod::analytic, /*implicit=*/false>
: public FaceCenteredLocalAssemblerBase<TypeTag, Assembler,
            FaceCenteredLocalAssembler<TypeTag, Assembler, DiffMethod::analytic, false>, false>
{
    using ThisType = FaceCenteredLocalAssembler<TypeTag, Assembler, DiffMethod::analytic, false>;
    using ParentType = FaceCenteredLocalAssemblerBase<TypeTag, Assembler, ThisType, false>;
    using NumEqVector = GetPropType<TypeTag, Properties::NumEqVector>;
    using JacobianMatrix = GetPropType<TypeTag, Properties::JacobianMatrix>;
    using GridVariables = GetPropType<TypeTag, Properties::GridVariables>;
    using LocalResidual = GetPropType<TypeTag, Properties::LocalResidual>;
    using ElementResidualVector = typename LocalResidual::ElementResidualVector;
    using VolumeVariables = GetPropType<TypeTag, Properties::VolumeVariables>;

    enum { numEq = GetPropType<TypeTag, Properties::ModelTraits>::numEq() };

public:
    using ParentType::ParentType;

    /*!
     * \brief Computes the derivatives with respect to the given element and adds them
     *        to the global matrix.
     *
     * \return The element residual at the current solution.
     */
    template <class PartialReassembler = DefaultPartialReassembler>
    ElementResidualVector assembleJacobianAndResidualImpl(JacobianMatrix& A, GridVariables& gridVariables,
                                                          const PartialReassembler* partialReassembler = nullptr)
    {
        if (partialReassembler)
            DUNE_THROW(Dune::NotImplemented, "partial reassembly for explicit time discretization");

        // get some aliases for convenience
        const auto& element = this->element();
        const auto& fvGeometry = this->fvGeometry();
        const auto& problem = this->problem();
        const auto& curElemVolVars = this->curElemVolVars();

        // get the vector of the actual element residuals
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
            this->localResidual().addStorageDerivatives(A[dofIdx][dofIdx],
                                                        problem,
                                                        element,
                                                        fvGeometry,
                                                        volVars,
                                                        scv);
        }

        return origResiduals;
    }
};

} // end namespace Dumux

#endif
