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
 * \ingroup MultiDomain
 * \brief Caculates the Jacobian of the local residual for fully-implicit models
 */
#ifndef DUMUX_SUBPROBLEM_LOCAL_JACOBIAN_FOR_STAGGERED_HH
#define DUMUX_SUBPROBLEM_LOCAL_JACOBIAN_FOR_STAGGERED_HH

#include <dune/common/indices.hh>
#include <dune/istl/matrix.hh>

#include <dumux/common/math.hh>
#include <dumux/common/valgrind.hh>

#include <dumux/multidomain/properties.hh>

namespace Dumux
{

namespace Properties
{
    // forward declaration of property tags
    NEW_PROP_TAG(StokesProblemTypeTag);
    NEW_PROP_TAG(DarcyProblemTypeTag);
}

template <class TypeTag, class SubProblemTypeTag>
class SubProblemLocalJacobianForStaggered;

template<class TypeTag>
using StokesLocalJacobianForStaggered = SubProblemLocalJacobianForStaggered<TypeTag, typename GET_PROP_TYPE(TypeTag, StokesProblemTypeTag)>;
template<class TypeTag>
using DarcyLocalJacobianForStaggered = SubProblemLocalJacobianForStaggered<TypeTag, typename GET_PROP_TYPE(TypeTag, DarcyProblemTypeTag)>;

/*!
 * \ingroup MultiDomain
 * \brief Calculates the Jacobian of the local residual for one subdomain domain
 *
 * The default behavior is to use numeric differentiation, i.e.
 * forward or backward differences (2nd order), or central
 * differences (3rd order). The method used is determined by the
 * "NumericDifferenceMethod" property:
 *
 * - if the value of this property is smaller than 0, backward
 *   differences are used, i.e.:
 *   \f[
 \frac{\partial f(x)}{\partial x} \approx \frac{f(x) - f(x - \epsilon)}{\epsilon}
 *   \f]
 *
 * - if the value of this property is 0, central
 *   differences are used, i.e.:
 *   \f[
 \frac{\partial f(x)}{\partial x} \approx \frac{f(x + \epsilon) - f(x - \epsilon)}{2 \epsilon}
 *   \f]
 *
 * - if the value of this property is larger than 0, forward
 *   differences are used, i.e.:
 *   \f[
 \frac{\partial f(x)}{\partial x} \approx \frac{f(x + \epsilon) - f(x)}{\epsilon}
 *   \f]
 *
 * Here, \f$ f \f$ is the residual function for all equations, \f$x\f$
 * is the value of a sub-control volume's primary variable at the
 * evaluation point and \f$\epsilon\f$ is a small value larger than 0.
 */
template<class TypeTag, class SubProblemTypeTag>
class SubProblemLocalJacobianForStaggered : public GET_PROP_TYPE(SubProblemTypeTag, LocalJacobian)
{
    using ParentType = typename GET_PROP_TYPE(SubProblemTypeTag, LocalJacobian);
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using GlobalProblem = typename GET_PROP_TYPE(TypeTag, Problem);
    using StokesProblemTypeTag = typename GET_PROP_TYPE(TypeTag, StokesProblemTypeTag);
    using DarcyProblemTypeTag = typename GET_PROP_TYPE(TypeTag, DarcyProblemTypeTag);

    // types of the sub problem
    using Problem = typename GET_PROP_TYPE(SubProblemTypeTag, Problem);
    using FVElementGeometry = typename GET_PROP_TYPE(SubProblemTypeTag, FVElementGeometry);
    using ElementSolutionVector = typename GET_PROP_TYPE(SubProblemTypeTag, ElementSolutionVector);
    using PrimaryVariables = typename GET_PROP_TYPE(SubProblemTypeTag, PrimaryVariables);
    using ElementVolumeVariables = typename GET_PROP_TYPE(SubProblemTypeTag, ElementVolumeVariables);
    using ElementFluxVariablesCache = typename GET_PROP_TYPE(SubProblemTypeTag, ElementFluxVariablesCache);
    using ElementBoundaryTypes = typename GET_PROP_TYPE(SubProblemTypeTag, ElementBoundaryTypes);
    using GridView = typename GET_PROP_TYPE(SubProblemTypeTag, GridView);
    using Element = typename GridView::template Codim<0>::Entity;
    using SolutionVector = typename GET_PROP_TYPE(SubProblemTypeTag, SolutionVector);

    using DofTypeIndices = typename GET_PROP(StokesProblemTypeTag, DofTypeIndices);
    typename DofTypeIndices::CellCenterIdx cellCenterIdx;
    typename DofTypeIndices::FaceIdx faceIdx;

    using SubProblemBlockIndices = typename GET_PROP(TypeTag, SubProblemBlockIndices);
    typename SubProblemBlockIndices::StokesIdx stokesIdx;
    typename SubProblemBlockIndices::DarcyIdx darcyIdx;

    typename GET_PROP(TypeTag, JacobianMatrix)::localDarcyIdx localDarcyIdx;

    using CellCenterPrimaryVariables = typename GET_PROP_TYPE(StokesProblemTypeTag, CellCenterPrimaryVariables);
    using GlobalFaceVars = typename GET_PROP_TYPE(StokesProblemTypeTag, GlobalFaceVars);
    using FaceSolutionVector = typename GET_PROP_TYPE(StokesProblemTypeTag, FaceSolutionVector);

public:

    // copying a local jacobian is not a good idea
    SubProblemLocalJacobianForStaggered(const SubProblemLocalJacobianForStaggered &) = delete;

    SubProblemLocalJacobianForStaggered()
    { numericDifferenceMethod_ = GET_PARAM_FROM_GROUP(TypeTag, int, Implicit, NumericDifferenceMethod); }

    /*!
     * \brief Initialize the local Jacobian object.
     *
     * At this point we can assume that everything has been allocated,
     * although some objects may not yet be completely initialized.
     *
     * \param problem The problem which we want to simulate.
     */
    void init(GlobalProblem &globalProblem)
    {
        globalProblemPtr_ = &globalProblem;
        ParentType::init(problem_());
    }

    /*!
     * \brief Assemble an element's local Jacobian matrix of the
     *        defect.
     *
     * \param element The DUNE Codim<0> entity which we look at.
     */
    template<class JacobianMatrix, class JacobianMatrixCoupling>
    void assemble(const Element &element,
                  JacobianMatrix& matrix,
                  JacobianMatrixCoupling& couplingMatrix,
                  SolutionVector& residual)
    {
        // prepare the volvars/fvGeometries in case caching is disabled
        auto fvGeometry = localView(this->model_().globalFvGeometry());
        fvGeometry.bind(element);

        auto curElemVolVars = localView(this->model_().curGlobalVolVars());
        curElemVolVars.bind(element, fvGeometry, this->model_().curSol());

        auto prevElemVolVars = localView(this->model_().prevGlobalVolVars());
        prevElemVolVars.bindElement(element, fvGeometry, this->model_().prevSol());

        auto elemFluxVarsCache = localView(this->model_().globalFluxVarsCache());
        elemFluxVarsCache.bind(element, fvGeometry, curElemVolVars);

        // check for boundaries on the element
        ElementBoundaryTypes elemBcTypes;
        elemBcTypes.update(this->problem_(), element, fvGeometry);

        assemble_(element, fvGeometry, prevElemVolVars, curElemVolVars, elemFluxVarsCache, elemBcTypes, matrix, couplingMatrix, residual);
    }

    protected:
    // for cell-centered models
    template<class JacobianMatrix, class JacobianMatrixCoupling, class T = SubProblemTypeTag,
    typename std::enable_if<((GET_PROP_VALUE(T, DiscretizationMethod) == DiscretizationMethods::CCTpfa)), bool>::type = 0>
    void assemble_(const Element& element,
                     const FVElementGeometry& fvGeometry,
                     ElementVolumeVariables& prevElemVolVars,
                     ElementVolumeVariables& curElemVolVars,
                     ElementFluxVariablesCache& elemFluxVarsCache,
                     const ElementBoundaryTypes& elemBcTypes,
                     JacobianMatrix& matrix,
                     JacobianMatrixCoupling& couplingMatrix,
                     SolutionVector& residual)
    {
        const bool isGhost = (element.partitionType() == Dune::GhostEntity);

        // set the actual dof index
        this->globalI_ = this->problem_().elementMapper().index(element);

        // calculate the local residual
        if (isGhost)
        {
            this->residual_ = 0.0;
            residual[this->globalI_] = 0.0;
        }
        else
        {
            // Evaluate the undeflected element local residual
            this->localResidual().eval(element,
                                       fvGeometry,
                                       prevElemVolVars,
                                       curElemVolVars,
                                       elemBcTypes,
                                       elemFluxVarsCache);
            this->residual_ = this->localResidual().residual();
            // store residual in global container as well
            residual[this->globalI_] = this->localResidual().residual(0);
        }

        this->model_().updatePVWeights(fvGeometry);

        // calculate derivatives of all dofs in stencil with respect to the dofs in the element
        this->evalPartialDerivatives_(element,
                                      fvGeometry,
                                      prevElemVolVars,
                                      curElemVolVars,
                                      elemFluxVarsCache,
                                      elemBcTypes,
                                      matrix,
                                      residual,
                                      isGhost);

//        // compute derivatives with respect to additional user defined DOF dependencies
//        const auto& additionalDofDepedencies = this->problem_().getAdditionalDofDependencies(this->globalI_);
//        if (!additionalDofDepedencies.empty() && !isGhost)
//        {
//            this->evalAdditionalDerivatives_(additionalDofDepedencies,
//                                             element,
//                                             fvGeometry,
//                                             curElemVolVars,
//                                             matrix,
//                                             residual);
//        }

        evalCCFVPartialDerivativeCoupling_(element,
                                       fvGeometry,
                                       curElemVolVars,
                                       elemFluxVarsCache,
                                       elemBcTypes,
                                       couplingMatrix,
                                       residual);
    }

    // for staggered grid model
    template<class JacobianMatrix, class JacobianMatrixCoupling, class T = SubProblemTypeTag,
    typename std::enable_if<((GET_PROP_VALUE(T, DiscretizationMethod) == DiscretizationMethods::Staggered)), bool>::type = 0>
    void assemble_(const Element& element,
                   const FVElementGeometry& fvGeometry,
                   ElementVolumeVariables& prevElemVolVars,
                   ElementVolumeVariables& curElemVolVars,
                   ElementFluxVariablesCache& elemFluxVarsCache,
                   const ElementBoundaryTypes& elemBcTypes,
                   JacobianMatrix& matrix,
                   JacobianMatrixCoupling& couplingMatrix,
                   SolutionVector& residual)
    {
        const bool isGhost = (element.partitionType() == Dune::GhostEntity);
        if(isGhost)
            DUNE_THROW(Dune::NotImplemented, "Support for ghost cells not implemented");


        // set the actual dof index
        this->ccGlobalI_ = this->problem_().elementMapper().index(element);

        auto curGlobalFaceVars = this->model_().curGlobalFaceVars();
        auto& prevGlobalFaceVars = this->model_().prevGlobalFaceVars();

        // Evaluate the undeflected element local residual
        this->localResidual().eval(element,
                                   fvGeometry,
                                   prevElemVolVars,
                                   curElemVolVars,
                                   prevGlobalFaceVars,
                                   curGlobalFaceVars,
                                   elemBcTypes,
                                   elemFluxVarsCache);

        // store the cell center residual in global container
        residual[cellCenterIdx][this->ccGlobalI_] = this->localResidual().ccResidual();

        // treat the local residua of the face dofs:
        // create a cache to reuse some results for the calculation of the derivatives
        FaceSolutionVector faceResidualCache;
        faceResidualCache.resize(fvGeometry.numScvf());
        faceResidualCache = 0.0;
        for(auto&& scvf : scvfs(fvGeometry))
        {
            residual[faceIdx][scvf.dofIndex()] += this->localResidual().faceResidual(scvf.localFaceIdx());
            faceResidualCache[scvf.localFaceIdx()] = this->localResidual().faceResidual(scvf.localFaceIdx());
        }

        this->model_().updatePVWeights(fvGeometry);

        // calculate derivatives of all dofs in stencil with respect to the dofs in the element
        this->evalPartialDerivatives_(element,
                                      fvGeometry,
                                      prevElemVolVars,
                                      curElemVolVars,
                                      prevGlobalFaceVars,
                                      curGlobalFaceVars,
                                      elemFluxVarsCache,
                                      elemBcTypes,
                                      matrix,
                                      residual[cellCenterIdx][this->ccGlobalI_],
                                      faceResidualCache);

        evalStaggeredPartialDerivativeCoupling_(element,
                                           fvGeometry,
                                           prevElemVolVars,
                                           curElemVolVars,
                                           prevGlobalFaceVars,
                                           curGlobalFaceVars,
                                           elemFluxVarsCache,
                                           elemBcTypes,
                                           couplingMatrix,
                                           residual[cellCenterIdx][this->ccGlobalI_],
                                           faceResidualCache);
    }


    /*!
     * \brief Returns a reference to the problem.
     */
    const GlobalProblem &globalProblem_() const
    { return *globalProblemPtr_; }

    GlobalProblem &globalProblem_()
    { return *globalProblemPtr_; }

    //! Return this subproblem
    template<class T = TypeTag>
    // static_assert(std::is_same<typename GET_PROP_TYPE(TypeTag, StokesProblemTypeTag), SubProblemTypeTag>::value, "class name: " + className<typename GET_PROP_TYPE(TypeTag, StokesProblemTypeTag)>());
    decltype(auto) problem_(typename std::enable_if<std::is_same<typename GET_PROP_TYPE(T, StokesProblemTypeTag), SubProblemTypeTag>::value, void>::type* x = nullptr)
    { return globalProblem_().stokesProblem(); }

    template<class T = TypeTag>
    decltype(auto) problem_(typename std::enable_if<std::is_same<typename GET_PROP_TYPE(T, DarcyProblemTypeTag), SubProblemTypeTag>::value, void>::type* x = nullptr)
    { return globalProblem_().darcyProblem(); }

    //! Return the other subproblem
    template<class T = TypeTag>
    decltype(auto) otherProblem_(typename std::enable_if<!std::is_same<typename GET_PROP_TYPE(T, StokesProblemTypeTag), SubProblemTypeTag>::value, void>::type* x = nullptr)
    { return globalProblem_().stokesProblem(); }

    template<class T = TypeTag>
    decltype(auto) otherProblem_(typename std::enable_if<!std::is_same<typename GET_PROP_TYPE(T, DarcyProblemTypeTag), SubProblemTypeTag>::value, void>::type* x = nullptr)
    { return globalProblem_().darcyProblem(); }

    // cell-centered
    template<class JacobianMatrixCoupling>
    void evalCCFVPartialDerivativeCoupling_(const Element& element,
                                        const FVElementGeometry& fvGeometry,
                                        ElementVolumeVariables& curElemVolVars,
                                        ElementFluxVariablesCache& elemFluxVarsCache,
                                        const ElementBoundaryTypes& elemBcTypes,
                                        JacobianMatrixCoupling& couplingMatrix,
                                        SolutionVector& residual)
    {
        const auto& couplingStencilCC = globalProblem_().couplingManager().couplingStencil(element, cellCenterIdx);
        const auto& couplingStencilFace = globalProblem_().couplingManager().couplingStencil(element, faceIdx);
        auto scv = fvGeometry.scv(0); // TODO 0 ?? idx

        // Darcy to cc
        for (auto globalJ : couplingStencilCC)
        {
            const auto originalResidual = globalProblem_().couplingManager().evalDarcyCouplingResidual(element,
                    fvGeometry,
                    curElemVolVars,
                    elemBcTypes,
                    elemFluxVarsCache);

            auto& otherPriVars = otherProblem_().model().curSol()[cellCenterIdx][globalJ];
            auto originalOtherPriVars = otherPriVars;

            // derivatives in the neighbors with respect to the current elements
            std::decay_t<decltype(originalResidual)> partialDeriv;
            for (int pvIdx = 0; pvIdx < otherPriVars.size(); pvIdx++)
            {
                const Scalar eps = this->numericEpsilon(otherPriVars[pvIdx]);
                Scalar delta = 0;

                if (numericDifferenceMethod_ >= 0)
                {
                    // we are not using backward differences, i.e. we need to
                    // calculate f(x + \epsilon)

                    // deflect primary variables
                    otherPriVars[pvIdx] += eps;
                    delta += eps;

                    // calculate the residual with the deflected primary variables
                    partialDeriv = globalProblem_().couplingManager().evalDarcyCouplingResidual(element,
                            fvGeometry,
                            curElemVolVars,
                            elemBcTypes,
                            elemFluxVarsCache);
                }
                else
                {
                    // we are using backward differences, i.e. we don't need
                    // to calculate f(x + \epsilon) and we can recycle the
                    // (already calculated) residual f(x)
                    partialDeriv = originalResidual;
                }


                if (numericDifferenceMethod_ <= 0)
                {
                    // we are not using forward differences, i.e. we
                    // need to calculate f(x - \epsilon)

                    // deflect the primary variables
                    otherPriVars[pvIdx] -= 2*eps;
                    delta += eps;

                    // calculate the residual with the deflected primary variables
                    partialDeriv -= globalProblem_().couplingManager().evalDarcyCouplingResidual(element,
                            fvGeometry,
                            curElemVolVars,
                            elemBcTypes,
                            elemFluxVarsCache);
                }
                else
                {
                    // we are using forward differences, i.e. we don't need to
                    // calculate f(x - \epsilon) and we can recycle the
                    // (already calculated) residual f(x)
                    partialDeriv -= originalResidual;

                }

                // divide difference in residuals by the magnitude of the
                // deflections between the two function evaluation
                partialDeriv /= delta;

                // restore the original state of the element solution vector
                otherPriVars = originalOtherPriVars;

                // update the global jacobian matrix (coupling block)
                for (auto&& scv : scvs(fvGeometry))
                    this->updateGlobalJacobian_(couplingMatrix[localDarcyIdx][cellCenterIdx], scv.dofIndex(), globalJ, pvIdx, partialDeriv[scv.indexInElement()]);
            }
        }

        // Darcy to face
        // TODO only horizontal faces on coupling interface? (compare matrix entries to Fetzer2017c for 2x2 cells)
        for (auto globalJ : couplingStencilFace)
        {
            const auto originalResidual = globalProblem_().couplingManager().evalDarcyCouplingResidual(element,
                    fvGeometry,
                    curElemVolVars,
                    elemBcTypes,
                    elemFluxVarsCache);

            auto& otherPriVars = otherProblem_().model().curSol()[faceIdx][globalJ];
            auto originalOtherPriVars = otherPriVars;

            // derivatives in the neighbors with respect to the current elements
            std::decay_t<decltype(originalResidual)> partialDeriv;
            for (int pvIdx = 0; pvIdx < otherPriVars.size(); pvIdx++)
            {
                const Scalar eps = this->numericEpsilon(otherPriVars[pvIdx]);
                Scalar delta = 0;

                if (numericDifferenceMethod_ >= 0)
                {
                    // we are not using backward differences, i.e. we need to
                    // calculate f(x + \epsilon)

                    // deflect primary variables
                    otherPriVars[pvIdx] += eps;
                    delta += eps;

                    // calculate the residual with the deflected primary variables
                    partialDeriv = globalProblem_().couplingManager().evalDarcyCouplingResidual(element,
                            fvGeometry,
                            curElemVolVars,
                            elemBcTypes,
                            elemFluxVarsCache);
                }
                else
                {
                    // we are using backward differences, i.e. we don't need
                    // to calculate f(x + \epsilon) and we can recycle the
                    // (already calculated) residual f(x)
                    partialDeriv = originalResidual;
                }


                if (numericDifferenceMethod_ <= 0)
                {
                    // we are not using forward differences, i.e. we
                    // need to calculate f(x - \epsilon)

                    // deflect the primary variables
                    otherPriVars[pvIdx] -= 2*eps;
                    delta += eps;

                    // calculate the residual with the deflected primary variables
                    partialDeriv -= globalProblem_().couplingManager().evalDarcyCouplingResidual(element,
                            fvGeometry,
                            curElemVolVars,
                            elemBcTypes,
                            elemFluxVarsCache);
                }
                else
                {
                    // we are using forward differences, i.e. we don't need to
                    // calculate f(x - \epsilon) and we can recycle the
                    // (already calculated) residual f(x)
                    partialDeriv -= originalResidual;

                }

                // divide difference in residuals by the magnitude of the
                // deflections between the two function evaluation
                partialDeriv /= delta;

                // restore the original state of the element solution vector
                otherPriVars = originalOtherPriVars;

                // update the global jacobian matrix (coupling block)
                for (auto&& scv : scvs(fvGeometry))
                    this->updateGlobalJacobian_(couplingMatrix[localDarcyIdx][faceIdx], scv.dofIndex(), globalJ, pvIdx, partialDeriv[scv.indexInElement()]);
            }
        }

    }

    // staggered grid
    template<class JacobianMatrixCoupling>
    void evalStaggeredPartialDerivativeCoupling_(const Element& element,
                                            const FVElementGeometry& fvGeometry,
                                            const ElementVolumeVariables& prevElemVolVars,
                                            ElementVolumeVariables& curElemVolVars,
                                            const GlobalFaceVars& prevGlobalFaceVars,
                                            GlobalFaceVars& curGlobalFaceVars,
                                            ElementFluxVariablesCache& elemFluxVarsCache,
                                            const ElementBoundaryTypes& elemBcTypes,
                                            JacobianMatrixCoupling& couplingMatrix,
                                            const CellCenterPrimaryVariables& ccResidual,
                                            const FaceSolutionVector& faceResidualCache)
    {
        const auto& couplingStencilCC = globalProblem_().couplingManager().couplingStencil(element);
        const auto ccGlobalI = globalProblem_().couplingManager().stokesProblem().elementMapper().index(element);

        for (auto globalJ : couplingStencilCC)
        {
            const auto originalResidual = globalProblem_().couplingManager().evalStokesCCCouplingResidual(element,
                    fvGeometry,
                    curElemVolVars,
                    prevGlobalFaceVars,
                    curGlobalFaceVars,
                    elemBcTypes,
                    elemFluxVarsCache);


            auto& otherPriVars = otherProblem_().model().curSol()[globalJ];
            auto originalOtherPriVars = otherPriVars;

            // derivatives in the neighbors with repect to the current elements
            std::decay_t<decltype(originalResidual)> partialDeriv;
            for (int pvIdx = 0; pvIdx < partialDeriv.size(); pvIdx++)
            {
                const Scalar eps = this->numericEpsilon(otherPriVars[pvIdx], cellCenterIdx, cellCenterIdx);
                Scalar delta = 0;

                if (numericDifferenceMethod_ >= 0)
                {
                    // we are not using backward differences, i.e. we need to
                    // calculate f(x + \epsilon)

                    // deflect primary variables
                    otherPriVars[pvIdx] += eps;
                    delta += eps;

                    // calculate the residual with the deflected primary variables
                    partialDeriv = globalProblem_().couplingManager().evalStokesCCCouplingResidual(element,
                            fvGeometry,
                            curElemVolVars,
                            prevGlobalFaceVars,
                            curGlobalFaceVars,
                            elemBcTypes,
                            elemFluxVarsCache);
                }
                else
                {
                    // we are using backward differences, i.e. we don't need
                    // to calculate f(x + \epsilon) and we can recycle the
                    // (already calculated) residual f(x)
                    partialDeriv = originalResidual;
                }


                if (numericDifferenceMethod_ <= 0)
                {
                    // we are not using forward differences, i.e. we
                    // need to calculate f(x - \epsilon)

                    // deflect the primary variables
                    otherPriVars[pvIdx] -= 2*eps;
                    delta += eps;

                    // calculate the residual with the deflected primary variables
                    partialDeriv -= globalProblem_().couplingManager().evalStokesCCCouplingResidual(element,
                            fvGeometry,
                            curElemVolVars,
                            prevGlobalFaceVars,
                            curGlobalFaceVars,
                            elemBcTypes,
                            elemFluxVarsCache);
                }
                else
                {
                    // we are using forward differences, i.e. we don't need to
                    // calculate f(x - \epsilon) and we can recycle the
                    // (already calculated) residual f(x)
                    partialDeriv -= originalResidual;
                }

                // divide difference in residuals by the magnitude of the
                // deflections between the two function evaluation
                partialDeriv /= delta;

                // restore the original state of the element solution vector
                otherPriVars = originalOtherPriVars;

                // update the global jacobian matrix (coupling block)
                this->updateGlobalJacobian_(couplingMatrix[cellCenterIdx][localDarcyIdx], ccGlobalI, globalJ, pvIdx, partialDeriv);

            }
        }

        // treat face derivatives
        for(auto&& scvf : scvfs(fvGeometry))
        {
            const auto originalResidual = globalProblem_().couplingManager().evalStokesFaceCouplingResidual(element,
                    fvGeometry,
                    scvf,
                    curElemVolVars,
                    prevGlobalFaceVars,
                    curGlobalFaceVars,
                    elemBcTypes,
                    elemFluxVarsCache);
            const auto& couplingStencilFace = globalProblem_().couplingManager().couplingStencil(scvf);
            const auto faceGlobalI = scvf.dofIndex();

            for(auto globalJ : couplingStencilFace)
            {
                auto& otherPriVars = otherProblem_().model().curSol()[globalJ];
                auto originalOtherPriVars = otherPriVars;

                // derivatives in the neighbors with repect to the current elements
                std::decay_t<decltype(originalResidual)> partialDeriv;
                for (int pvIdx = 0; pvIdx < partialDeriv.size(); pvIdx++)
                {
                    const Scalar eps = this->numericEpsilon(otherPriVars[pvIdx], cellCenterIdx, cellCenterIdx);
                    Scalar delta = 0;

                    if (numericDifferenceMethod_ >= 0)
                    {
                        // we are not using backward differences, i.e. we need to
                        // calculate f(x + \epsilon)

                        // deflect primary variables
                        otherPriVars[pvIdx] += eps;
                        delta += eps;

                        // calculate the residual with the deflected primary variables
                        partialDeriv = globalProblem_().couplingManager().evalStokesFaceCouplingResidual(element,
                                fvGeometry,
                                scvf,
                                curElemVolVars,
                                prevGlobalFaceVars,
                                curGlobalFaceVars,
                                elemBcTypes,
                                elemFluxVarsCache);
                    }
                    else
                    {
                        // we are using backward differences, i.e. we don't need
                        // to calculate f(x + \epsilon) and we can recycle the
                        // (already calculated) residual f(x)
                        partialDeriv = originalResidual;
                    }


                    if (numericDifferenceMethod_ <= 0)
                    {
                        // we are not using forward differences, i.e. we
                        // need to calculate f(x - \epsilon)

                        // deflect the primary variables
                        otherPriVars[pvIdx] -= 2*eps;
                        delta += eps;

                        // calculate the residual with the deflected primary variables
                        partialDeriv -= globalProblem_().couplingManager().evalStokesFaceCouplingResidual(element,
                                fvGeometry,
                                scvf,
                                curElemVolVars,
                                prevGlobalFaceVars,
                                curGlobalFaceVars,
                                elemBcTypes,
                                elemFluxVarsCache);
                    }
                    else
                    {
                        // we are using forward differences, i.e. we don't need to
                        // calculate f(x - \epsilon) and we can recycle the
                        // (already calculated) residual f(x)
                        partialDeriv -= originalResidual;
                    }

                    // divide difference in residuals by the magnitude of the
                    // deflections between the two function evaluation
                    partialDeriv /= delta;

                    // restore the original state of the element solution vector
                    otherPriVars = originalOtherPriVars;

                    // update the global jacobian matrix (coupling block)
                    this->updateGlobalJacobian_(couplingMatrix[faceIdx][localDarcyIdx], faceGlobalI, globalJ, pvIdx, partialDeriv);

                }
            }
        }


    }

    // The problem we would like to solve
    GlobalProblem *globalProblemPtr_;
    // The type of the numeric difference method (forward, center, backward)
    int numericDifferenceMethod_;
};

} // end namespace Dumux

#endif
