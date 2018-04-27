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
 * \ingroup StaggeredDiscretization
 * \ingroup Assembly
 * \brief An assembler for Jacobian and residual contribution per element (staggered FV method)
 */
#ifndef DUMUX_STAGGERED_LOCAL_ASSEMBLER_HH
#define DUMUX_STAGGERED_LOCAL_ASSEMBLER_HH

#include <dune/common/version.hh>
#include <dune/common/hybridutilities.hh>
#include <dune/common/rangeutilities.hh>

#include <dune/istl/matrixindexset.hh>
#include <dune/istl/bvector.hh>

#include <dumux/assembly/diffmethod.hh>
#include <dumux/assembly/numericepsilon.hh>
#include <dumux/common/numericdifferentiation.hh>
#include <dumux/common/properties.hh>
#include <dumux/discretization/staggered/facesolution.hh>
#include <dumux/discretization/staggered/elementsolution.hh>

namespace Dumux {

/*!
 * \ingroup Assembly
 * \ingroup StaggeredDiscretization
 * \brief An assembler for Jacobian and residual contribution per element (staggered FV method)
 * \tparam TypeTag the TypeTag
 * \tparam DM the differentiation method to residual compute derivatives
 * \tparam implicit if to use an implicit or explicit time discretization
 */
template<class TypeTag,
         DiffMethod diffMethod = DiffMethod::numeric,
         bool implicit = true>
class StaggeredLocalAssembler;

/*!
 * \ingroup Assembly
 * \brief Staggerd FV scheme local assembler using numeric differentiation and implicit time discretization
 */
template<class TypeTag>
class StaggeredLocalAssembler<TypeTag,
                       DiffMethod::numeric,
                       /*implicit=*/true>
{
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using NumEqVector = typename GET_PROP_TYPE(TypeTag, NumEqVector);
    using ElementBoundaryTypes = typename GET_PROP_TYPE(TypeTag, ElementBoundaryTypes);
    using FVGridGeometry = typename GET_PROP_TYPE(TypeTag, FVGridGeometry);
    using FVElementGeometry = typename FVGridGeometry::LocalView;
    using GridFaceVariables = typename GET_PROP_TYPE(TypeTag, GridFaceVariables);
    using ElementFluxVariablesCache = typename GET_PROP_TYPE(TypeTag, GridFluxVariablesCache)::LocalView;
    using Element = typename GET_PROP_TYPE(TypeTag, GridView)::template Codim<0>::Entity;
    using SolutionVector = typename GET_PROP_TYPE(TypeTag, SolutionVector);
    using ElementVolumeVariables = typename GET_PROP_TYPE(TypeTag, GridVolumeVariables)::LocalView;
    using GridVolumeVariables = typename GET_PROP_TYPE(TypeTag, GridVolumeVariables);
    using VolumeVariables = typename GET_PROP_TYPE(TypeTag, VolumeVariables);
    using SubControlVolume = typename FVElementGeometry::SubControlVolume;
    using JacobianMatrix = typename GET_PROP_TYPE(TypeTag, JacobianMatrix);

    using NumCellCenterEqVector =  typename GET_PROP_TYPE(TypeTag, CellCenterPrimaryVariables);
    using NumFaceEqVector = typename GET_PROP_TYPE(TypeTag, FacePrimaryVariables);

    using FaceSolutionVector = typename GET_PROP_TYPE(TypeTag, FaceSolutionVector);
    using FaceSolution = typename GET_PROP_TYPE(TypeTag, StaggeredFaceSolution);

    using Indices = typename GET_PROP_TYPE(TypeTag, ModelTraits)::Indices;
    using PrimaryVariables = typename GET_PROP_TYPE(TypeTag, PrimaryVariables);
    using FacePrimaryVariables = typename GET_PROP_TYPE(TypeTag, FacePrimaryVariables);
    using CellCenterPrimaryVariables = typename GET_PROP_TYPE(TypeTag, CellCenterPrimaryVariables);
    using FaceVariables = typename GET_PROP_TYPE(TypeTag, FaceVariables);
    using ElementFaceVariables = typename GET_PROP_TYPE(TypeTag, GridFaceVariables)::LocalView;
    using SubControlVolumeFace = typename FVElementGeometry::SubControlVolumeFace;

    static constexpr int numEq = GET_PROP_TYPE(TypeTag, ModelTraits)::numEq();
    static constexpr bool enableGridFluxVarsCache = GET_PROP_VALUE(TypeTag, EnableGridFluxVariablesCache);

public:

    /*!
     * \brief Assemble the residual only
     */
    template<class Assembler>
    static void assembleResidual(Assembler& assembler,
                                SolutionVector& res,
                                const Element& element,
                                const SolutionVector& curSol)
    {
        // get some references for convenience
        const auto& problem = assembler.problem();
        auto& localResidual = assembler.localResidual();
        auto& gridVariables = assembler.gridVariables();
        static constexpr auto cellCenterIdx = FVGridGeometry::cellCenterIdx();
        static constexpr auto faceIdx = FVGridGeometry::faceIdx();

        // prepare the local views
        auto fvGeometry = localView(assembler.fvGridGeometry());
        fvGeometry.bind(element);

        auto curElemVolVars = localView(gridVariables.curGridVolVars());
        curElemVolVars.bind(element, fvGeometry, curSol);

        auto elemFluxVarsCache = localView(gridVariables.gridFluxVarsCache());
        elemFluxVarsCache.bind(element, fvGeometry, curElemVolVars);

        auto curElemFaceVars = localView(gridVariables.curGridFaceVars());
        curElemFaceVars.bind(element, fvGeometry, curSol);

        const bool isStationary = localResidual.isStationary();
        auto prevElemVolVars = localView(gridVariables.prevGridVolVars());
        auto prevElemFaceVars = localView(gridVariables.prevGridFaceVars());
        if (!isStationary)
        {
            prevElemVolVars.bindElement(element, fvGeometry, localResidual.prevSol());
            prevElemFaceVars.bindElement(element, fvGeometry, localResidual.prevSol());
        }

        // for compatibility with box models
        ElementBoundaryTypes elemBcTypes;

        const auto cellCenterGlobalI = assembler.fvGridGeometry().elementMapper().index(element);
        res[cellCenterIdx][cellCenterGlobalI] = localResidual.evalCellCenter(problem,
                                                                             element,
                                                                             fvGeometry,
                                                                             prevElemVolVars,
                                                                             curElemVolVars,
                                                                             prevElemFaceVars,
                                                                             curElemFaceVars,
                                                                             elemBcTypes,
                                                                             elemFluxVarsCache);

        // treat the local residua of the face dofs:
        for(auto&& scvf : scvfs(fvGeometry))
        {
            res[faceIdx][scvf.dofIndex()] += localResidual.evalFace(problem,
                                                                    element,
                                                                    fvGeometry,
                                                                    scvf,
                                                                    prevElemVolVars,
                                                                    curElemVolVars,
                                                                    prevElemFaceVars,
                                                                    curElemFaceVars,
                                                                    elemBcTypes,
                                                                    elemFluxVarsCache);

        }
    }

    /*!
     * \brief Computes the derivatives with respect to the given element and adds them
     *        to the global matrix. The element residual is written into the right hand side.
     */
    template<class Assembler>
    static void assembleJacobianAndResidual(Assembler& assembler,
                                            JacobianMatrix& jac,
                                            SolutionVector& res,
                                            const Element& element,
                                            const SolutionVector& curSol)
    {
        // get some references for convenience
        const auto& problem = assembler.problem();
        auto& localResidual = assembler.localResidual();
        auto& gridVariables = assembler.gridVariables();
        static constexpr auto cellCenterIdx = FVGridGeometry::cellCenterIdx();
        static constexpr auto faceIdx = FVGridGeometry::faceIdx();

        // prepare the local views
        auto fvGeometry = localView(assembler.fvGridGeometry());
        fvGeometry.bind(element);

        auto curElemVolVars = localView(gridVariables.curGridVolVars());
        curElemVolVars.bind(element, fvGeometry, curSol);

        auto elemFluxVarsCache = localView(gridVariables.gridFluxVarsCache());
        elemFluxVarsCache.bind(element, fvGeometry, curElemVolVars);

        auto curElemFaceVars = localView(gridVariables.curGridFaceVars());
        curElemFaceVars.bind(element, fvGeometry, curSol);

        const bool isStationary = localResidual.isStationary();
        auto prevElemVolVars = localView(gridVariables.prevGridVolVars());
        auto prevElemFaceVars = localView(gridVariables.prevGridFaceVars());
        if (!isStationary)
        {
            prevElemVolVars.bindElement(element, fvGeometry, localResidual.prevSol());
            prevElemFaceVars.bindElement(element, fvGeometry, localResidual.prevSol());
        }

        // for compatibility with box models
        ElementBoundaryTypes elemBcTypes;

        const auto cellCenterGlobalI = assembler.fvGridGeometry().elementMapper().index(element);
        res[cellCenterIdx][cellCenterGlobalI] = localResidual.evalCellCenter(problem,
                                                                             element,
                                                                             fvGeometry,
                                                                             prevElemVolVars,
                                                                             curElemVolVars,
                                                                             prevElemFaceVars,
                                                                             curElemFaceVars,
                                                                             elemBcTypes,
                                                                             elemFluxVarsCache);

        // treat the local residua of the face dofs:
        // create a cache to reuse some results for the calculation of the derivatives
        FaceSolutionVector faceResidualCache;
        faceResidualCache.resize(fvGeometry.numScvf());
        faceResidualCache = 0.0;

        for(auto&& scvf : scvfs(fvGeometry))
        {
            faceResidualCache[scvf.localFaceIdx()] = localResidual.evalFace(problem,
                                                                            element,
                                                                            fvGeometry,
                                                                            scvf,
                                                                            prevElemVolVars,
                                                                            curElemVolVars,
                                                                            prevElemFaceVars,
                                                                            curElemFaceVars,
                                                                            elemBcTypes,
                                                                            elemFluxVarsCache);

            res[faceIdx][scvf.dofIndex()] += faceResidualCache[scvf.localFaceIdx()] ;
        }

        // calculate derivatives of all dofs in stencil with respect to the dofs in the element
        evalPartialDerivatives_(assembler,
                                element,
                                fvGeometry,
                                curSol,
                                prevElemVolVars,
                                curElemVolVars,
                                prevElemFaceVars,
                                curElemFaceVars,
                                elemFluxVarsCache,
                                elemBcTypes,
                                jac,
                                res[cellCenterIdx][cellCenterGlobalI],
                                faceResidualCache);

    }

protected:
    template<class Assembler>
    static void evalPartialDerivatives_(Assembler& assembler,
                                        const Element& element,
                                        const FVElementGeometry& fvGeometry,
                                        const SolutionVector& curSol,
                                        const ElementVolumeVariables& prevElemVolVars,
                                        ElementVolumeVariables& curElemVolVars,
                                        const ElementFaceVariables& prevElemFaceVars,
                                        ElementFaceVariables& curElemFaceVars,
                                        ElementFluxVariablesCache& elemFluxVarsCache,
                                        const ElementBoundaryTypes& elemBcTypes,
                                        JacobianMatrix& matrix,
                                        const NumCellCenterEqVector& ccResidual,
                                        const FaceSolutionVector& faceResidualCache)
    {
        // compute the derivatives of the cell center residuals with respect to cell center dofs
        dCCdCC_(assembler, element, fvGeometry, curSol, prevElemVolVars, curElemVolVars, prevElemFaceVars, curElemFaceVars, elemFluxVarsCache, elemBcTypes, matrix, ccResidual);

        // compute the derivatives of the cell center residuals with respect to face dofs
        dCCdFace_(assembler, element, fvGeometry, curSol, prevElemVolVars, curElemVolVars, prevElemFaceVars, curElemFaceVars, elemFluxVarsCache, elemBcTypes, matrix, ccResidual);

        // compute the derivatives of the face residuals with respect to cell center dofs
        dFacedCC_(assembler, element, fvGeometry, curSol, prevElemVolVars, curElemVolVars, prevElemFaceVars, curElemFaceVars, elemFluxVarsCache, elemBcTypes, matrix, faceResidualCache);

        // compute the derivatives of the face residuals with respect to face dofs
        dFacedFace_(assembler, element, fvGeometry, curSol, prevElemVolVars, curElemVolVars, prevElemFaceVars, curElemFaceVars, elemFluxVarsCache, elemBcTypes, matrix, faceResidualCache);
    }

    /*!
     * \brief Computes the derivatives of the cell center residuals with respect to cell center dofs
     */
    template<class Assembler>
    static void dCCdCC_(Assembler& assembler,
                 const Element& element,
                 const FVElementGeometry& fvGeometry,
                 const SolutionVector& curSol,
                 const ElementVolumeVariables& prevElemVolVars,
                 ElementVolumeVariables& curElemVolVars,
                 const ElementFaceVariables& prevElemFaceVars,
                 const ElementFaceVariables& curElemFaceVars,
                 ElementFluxVariablesCache& elemFluxVarsCache,
                 const ElementBoundaryTypes& elemBcTypes,
                 JacobianMatrix& matrix,
                 const NumCellCenterEqVector& ccResidual)
    {
        const auto& problem = assembler.problem();
        auto& localResidual = assembler.localResidual();
        auto& gridVariables = assembler.gridVariables();
        static constexpr auto cellCenterIdx = FVGridGeometry::cellCenterIdx();

        // build derivatives with for cell center dofs w.r.t. cell center dofs
        const auto cellCenterGlobalI = assembler.fvGridGeometry().elementMapper().index(element);

        const auto& connectivityMap = assembler.fvGridGeometry().connectivityMap();

        for(const auto& globalJ : connectivityMap(cellCenterIdx, cellCenterIdx, cellCenterGlobalI))
        {
            // get the volVars of the element with respect to which we are going to build the derivative
            auto&& scvJ = fvGeometry.scv(globalJ);
            const auto elementJ = fvGeometry.fvGridGeometry().element(globalJ);
            auto& curVolVars =  getVolVarAccess(gridVariables.curGridVolVars(), curElemVolVars, scvJ);
            VolumeVariables origVolVars(curVolVars);

            for(auto pvIdx : priVarIndices_(cellCenterIdx))
            {
                using PrimaryVariables = typename VolumeVariables::PrimaryVariables;
                const auto& cellCenterPriVars = curSol[cellCenterIdx][globalJ];
                PrimaryVariables priVars = makePriVarsFromCellCenterPriVars<PrimaryVariables>(cellCenterPriVars);

                constexpr auto offset = PrimaryVariables::dimension - CellCenterPrimaryVariables::dimension;
                NumCellCenterEqVector partialDeriv(0.0);

                auto evalResidual = [&](Scalar priVar)
                {
                      // update the volume variables and compute element residual
                      priVars[pvIdx + offset] = priVar;
                      auto elemSol = elementSolution<FVElementGeometry>(std::move(priVars));
                      curVolVars.update(elemSol, problem, elementJ, scvJ);
                      return localResidual.evalCellCenter(problem, element, fvGeometry,
                                                          prevElemVolVars, curElemVolVars,
                                                          prevElemFaceVars, curElemFaceVars,
                                                          elemBcTypes, elemFluxVarsCache);
                };

                // derive the residuals numerically
                static const int numDiffMethod = getParamFromGroup<int>(GET_PROP_VALUE(TypeTag, ModelParameterGroup), "Assembly.NumericDifferenceMethod");
                static const NumericEpsilon<Scalar, numEq> eps{GET_PROP_VALUE(TypeTag, ModelParameterGroup)};
                NumericDifferentiation::partialDerivative(evalResidual, priVars[pvIdx + offset], partialDeriv, ccResidual,
                                                          eps(priVars[pvIdx + offset], pvIdx + offset), numDiffMethod);

                // update the global jacobian matrix with the current partial derivatives
                updateGlobalJacobian_(matrix[cellCenterIdx][cellCenterIdx], cellCenterGlobalI, globalJ, pvIdx, partialDeriv);

                // restore the original volVars
                curVolVars = origVolVars;
            }
        }
    }

    /*!
     * \brief Computes the derivatives of the cell center residuals with respect to face dofs
     */
    template<class Assembler>
    static void dCCdFace_(Assembler& assembler,
                          const Element& element,
                          const FVElementGeometry& fvGeometry,
                          const SolutionVector& curSol,
                          const ElementVolumeVariables& prevElemVolVars,
                          const ElementVolumeVariables& curElemVolVars,
                          const ElementFaceVariables& prevElemFaceVars,
                          ElementFaceVariables& curElemFaceVars,
                          ElementFluxVariablesCache& elemFluxVarsCache,
                          const ElementBoundaryTypes& elemBcTypes,
                          JacobianMatrix& matrix,
                          const NumCellCenterEqVector& ccResidual)
    {
        const auto& problem = assembler.problem();
        auto& localResidual = assembler.localResidual();
        auto& gridVariables = assembler.gridVariables();
        static constexpr auto cellCenterIdx = FVGridGeometry::cellCenterIdx();
        static constexpr auto faceIdx = FVGridGeometry::faceIdx();

        // build derivatives with for cell center dofs w.r.t. cell center dofs
        const auto cellCenterGlobalI = assembler.fvGridGeometry().elementMapper().index(element);

        for(const auto& scvfJ : scvfs(fvGeometry))
        {
            const auto globalJ = scvfJ.dofIndex();

            // get the faceVars of the face with respect to which we are going to build the derivative
            auto& faceVars = getFaceVarAccess(gridVariables.curGridFaceVars(), curElemFaceVars, scvfJ);
            const auto origFaceVars(faceVars);

            for(auto pvIdx : priVarIndices_(faceIdx))
            {
                FacePrimaryVariables facePriVars(curSol[faceIdx][globalJ]);
                NumCellCenterEqVector partialDeriv(0.0);

                auto evalResidual = [&](Scalar priVar)
                {
                      // update the face variables and compute element residual
                      facePriVars[pvIdx] = priVar;
                      faceVars.updateOwnFaceOnly(facePriVars);
                      return localResidual.evalCellCenter(problem, element, fvGeometry,
                                                          prevElemVolVars, curElemVolVars,
                                                          prevElemFaceVars, curElemFaceVars,
                                                          elemBcTypes, elemFluxVarsCache);
                };

                // derive the residuals numerically
                static const int numDiffMethod = getParamFromGroup<int>(GET_PROP_VALUE(TypeTag, ModelParameterGroup), "Assembly.NumericDifferenceMethod");
                static const NumericEpsilon<Scalar, numEq> eps{GET_PROP_VALUE(TypeTag, ModelParameterGroup)};
                NumericDifferentiation::partialDerivative(evalResidual, facePriVars[pvIdx], partialDeriv, ccResidual,
                                                          eps(facePriVars, Indices::velocity(scvfJ.directionIndex())), numDiffMethod);

                // update the global jacobian matrix with the current partial derivatives
                updateGlobalJacobian_(matrix[cellCenterIdx][faceIdx], cellCenterGlobalI, globalJ, pvIdx, partialDeriv);

                // restore the original faceVars
                faceVars = origFaceVars;
            }
        }
    }

    /*!
     * \brief Computes the derivatives of the face residuals with respect to cell center dofs
     */
    template<class Assembler>
    static void dFacedCC_(Assembler& assembler,
                          const Element& element,
                          const FVElementGeometry& fvGeometry,
                          const SolutionVector& curSol,
                          const ElementVolumeVariables& prevElemVolVars,
                          ElementVolumeVariables& curElemVolVars,
                          const ElementFaceVariables& prevElemFaceVars,
                          const ElementFaceVariables& curElemFaceVars,
                          ElementFluxVariablesCache& elemFluxVarsCache,
                          const ElementBoundaryTypes& elemBcTypes,
                          JacobianMatrix& matrix,
                          const FaceSolutionVector& cachedResidual)
    {
        for(auto&& scvf : scvfs(fvGeometry))
        {
            const auto& problem = assembler.problem();
            auto& localResidual = assembler.localResidual();
            auto& gridVariables = assembler.gridVariables();
            const auto& connectivityMap = assembler.fvGridGeometry().connectivityMap();
            static constexpr auto cellCenterIdx = FVGridGeometry::cellCenterIdx();
            static constexpr auto faceIdx = FVGridGeometry::faceIdx();

            // set the actual dof index
            const auto faceGlobalI = scvf.dofIndex();

            // build derivatives with for face dofs w.r.t. cell center dofs
            for(const auto& globalJ : connectivityMap(faceIdx, cellCenterIdx, scvf.index()))
            {
                // get the volVars of the element with respect to which we are going to build the derivative
                auto&& scvJ = fvGeometry.scv(globalJ);
                const auto elementJ = fvGeometry.fvGridGeometry().element(globalJ);
                auto& curVolVars = getVolVarAccess(gridVariables.curGridVolVars(), curElemVolVars, scvJ);
                VolumeVariables origVolVars(curVolVars);

                for(auto pvIdx : priVarIndices_(cellCenterIdx))
                {
                    using PrimaryVariables = typename VolumeVariables::PrimaryVariables;
                    const auto& cellCenterPriVars = curSol[cellCenterIdx][globalJ];
                    PrimaryVariables priVars = makePriVarsFromCellCenterPriVars<PrimaryVariables>(cellCenterPriVars);

                    constexpr auto offset = PrimaryVariables::dimension - CellCenterPrimaryVariables::dimension;
                    FacePrimaryVariables partialDeriv(0.0);

                    auto evalResidual = [&](Scalar priVar)
                    {
                          // update the volume variables and compute element residual
                          priVars[pvIdx + offset] = priVar;
                          auto elemSol = elementSolution<FVElementGeometry>(std::move(priVars));
                          curVolVars.update(elemSol, problem, elementJ, scvJ);
                          return localResidual.evalFace(problem, element, fvGeometry, scvf,
                                                        prevElemVolVars, curElemVolVars,
                                                        prevElemFaceVars, curElemFaceVars,
                                                        elemBcTypes, elemFluxVarsCache);
                    };

                    // derive the residuals numerically
                    static const int numDiffMethod = getParamFromGroup<int>(GET_PROP_VALUE(TypeTag, ModelParameterGroup), "Assembly.NumericDifferenceMethod");
                    static const NumericEpsilon<Scalar, numEq> eps{GET_PROP_VALUE(TypeTag, ModelParameterGroup)};
                    NumericDifferentiation::partialDerivative(evalResidual, priVars[pvIdx + offset], partialDeriv, cachedResidual[scvf.localFaceIdx()],
                                                              eps(priVars[pvIdx + offset], pvIdx + offset), numDiffMethod);

                    // update the global jacobian matrix with the current partial derivatives
                    updateGlobalJacobian_(matrix[faceIdx][cellCenterIdx], faceGlobalI, globalJ, pvIdx, partialDeriv);

                    // restore the original volVars
                    curVolVars = origVolVars;
                }
            }
        }
    }

    /*!
     * \brief Computes the derivatives of the face residuals with respect to cell center dofs
     */
    template<class Assembler>
    static void dFacedFace_(Assembler& assembler,
                            const Element& element,
                            const FVElementGeometry& fvGeometry,
                            const SolutionVector& curSol,
                            const ElementVolumeVariables& prevElemVolVars,
                            const ElementVolumeVariables& curElemVolVars,
                            const ElementFaceVariables& prevElemFaceVars,
                            ElementFaceVariables& curElemFaceVars,
                            ElementFluxVariablesCache& elemFluxVarsCache,
                            const ElementBoundaryTypes& elemBcTypes,
                            JacobianMatrix& matrix,
                            const FaceSolutionVector& cachedResidual)
    {
        const auto& problem = assembler.problem();
        auto& localResidual = assembler.localResidual();
        const auto& connectivityMap = assembler.fvGridGeometry().connectivityMap();
        auto& gridVariables = assembler.gridVariables();
        static constexpr auto faceIdx = FVGridGeometry::faceIdx();

        for(auto&& scvf : scvfs(fvGeometry))
        {
            // set the actual dof index
            const auto faceGlobalI = scvf.dofIndex();

            // build derivatives with for face dofs w.r.t. cell center dofs
            for(const auto& globalJ : connectivityMap(faceIdx, faceIdx, scvf.index()))
            {
                // get the faceVars of the face with respect to which we are going to build the derivative
                auto& faceVars = getFaceVarAccess(gridVariables.curGridFaceVars(), curElemFaceVars, scvf);
                const auto origFaceVars(faceVars);

                for(auto pvIdx : priVarIndices_(faceIdx))
                {
                    auto faceSolution = FaceSolution(scvf, curSol[faceIdx], assembler.fvGridGeometry());
                    FacePrimaryVariables partialDeriv(0.0);

                    auto evalResidual = [&](Scalar priVar)
                    {
                          // update the volume variables and compute element residual
                          faceSolution[globalJ][pvIdx] = priVar;
                          faceVars.update(faceSolution, problem, element, fvGeometry, scvf);
                          return localResidual.evalFace(problem, element, fvGeometry, scvf,
                                                        prevElemVolVars, curElemVolVars,
                                                        prevElemFaceVars, curElemFaceVars,
                                                        elemBcTypes, elemFluxVarsCache);
                    };

                    // derive the residuals numerically
                    static const int numDiffMethod = getParamFromGroup<int>(GET_PROP_VALUE(TypeTag, ModelParameterGroup), "Assembly.NumericDifferenceMethod");
                    static const NumericEpsilon<Scalar, numEq> eps{GET_PROP_VALUE(TypeTag, ModelParameterGroup)};
                    NumericDifferentiation::partialDerivative(evalResidual, faceSolution[globalJ][pvIdx], partialDeriv, cachedResidual[scvf.localFaceIdx()],
                                                              eps(faceSolution[globalJ], Indices::velocity(scvf.directionIndex())), numDiffMethod);

                    // update the global jacobian matrix with the current partial derivatives
                    updateGlobalJacobian_(matrix[faceIdx][faceIdx], faceGlobalI, globalJ, pvIdx, partialDeriv);

                    // restore the original faceVars
                    faceVars = origFaceVars;
                }
            }
        }
    }

    /*!
     * \brief Updates the current global Jacobian matrix with the
     *        partial derivatives of all equations in regard to the
     *        primary variable 'pvIdx' at dof 'col'. Specialization for cc methods.
     */
    template<class SubMatrix, class CCOrFacePrimaryVariables>
    static void updateGlobalJacobian_(SubMatrix& matrix,
                          const int globalI,
                          const int globalJ,
                          const int pvIdx,
                          const CCOrFacePrimaryVariables &partialDeriv)
    {
        for (int eqIdx = 0; eqIdx < partialDeriv.size(); eqIdx++)
        {
            // A[i][col][eqIdx][pvIdx] is the rate of change of
            // the residual of equation 'eqIdx' at dof 'i'
            // depending on the primary variable 'pvIdx' at dof
            // 'col'.

            assert(pvIdx >= 0);
            assert(eqIdx < matrix[globalI][globalJ].size());
            assert(pvIdx < matrix[globalI][globalJ][eqIdx].size());
            matrix[globalI][globalJ][eqIdx][pvIdx] += partialDeriv[eqIdx];
        }
    }

    //! Helper function that returns an iterable range of primary variable indices.
    //! Specialization for cell center dofs.
    static auto priVarIndices_(typename FVGridGeometry::DofTypeIndices::CellCenterIdx)
    {
        constexpr auto numEqCellCenter = GET_PROP_VALUE(TypeTag, NumEqCellCenter);
        return Dune::range(0, numEqCellCenter);
    }

    //! Helper function that returns an iterable range of primary variable indices.
    //! Specialization for face dofs.
    static auto priVarIndices_(typename FVGridGeometry::DofTypeIndices::FaceIdx)
    {
        constexpr auto numEqFace = GET_PROP_VALUE(TypeTag, NumEqFace);
        return Dune::range(0, numEqFace);
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

    template<class T = TypeTag>
    static typename std::enable_if<!GET_PROP_VALUE(T, EnableGridFaceVariablesCache), FaceVariables&>::type
    getFaceVarAccess(GridFaceVariables& GridFaceVariables, ElementFaceVariables& elemFaceVars, const SubControlVolumeFace& scvf)
    { return elemFaceVars[scvf]; }

    template<class T = TypeTag>
    static typename std::enable_if<GET_PROP_VALUE(T, EnableGridFaceVariablesCache), FaceVariables&>::type
    getFaceVarAccess(GridFaceVariables& GridFaceVariables, ElementFaceVariables& elemFaceVars, const SubControlVolumeFace& scvf)
    { return GridFaceVariables.faceVars(scvf.index()); }
};

} // end namespace Dumux

#endif
