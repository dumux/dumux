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
 * \ingroup MixedDimension
 * \brief Caculates the Jacobian of the local residual for fully-implicit models
 */
#ifndef DUMUX_LOWDIM_LOCAL_JACOBIAN_HH
#define DUMUX_LOWDIM_LOCAL_JACOBIAN_HH

#include <dune/common/indices.hh>
#include <dune/istl/matrix.hh>

#include <dumux/common/math.hh>
#include <dumux/common/valgrind.hh>

#include <dumux/mixeddimension/properties.hh>

namespace Dumux
{
/*!
 * \ingroup MixedDimension
 * \brief Calculates the Jacobian of the local residual for the low dim domain
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
template<class TypeTag>
class LowDimLocalJacobian : public GET_PROP_TYPE(TypeTag, LowDimLocalJacobianBase)
{
    using ParentType = typename GET_PROP_TYPE(TypeTag, LowDimLocalJacobianBase);
    using Implementation = typename GET_PROP_TYPE(TypeTag, LowDimLocalJacobian);
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using GlobalProblem = typename GET_PROP_TYPE(TypeTag, Problem);
    using SubProblemBlockIndices = typename GET_PROP(TypeTag, SubProblemBlockIndices);

    // obtain the type tag of the lowdim sub problem
    using LowDimProblemTypeTag = typename GET_PROP_TYPE(TypeTag, LowDimProblemTypeTag);

    // types of the lowdim sub problem
    using Problem = typename GET_PROP_TYPE(LowDimProblemTypeTag, Problem);
    using FVElementGeometry = typename GET_PROP_TYPE(LowDimProblemTypeTag, FVElementGeometry);
    using ElementSolutionVector = typename GET_PROP_TYPE(LowDimProblemTypeTag, ElementSolutionVector);
    using PrimaryVariables = typename GET_PROP_TYPE(LowDimProblemTypeTag, PrimaryVariables);
    using VolumeVariables = typename GET_PROP_TYPE(LowDimProblemTypeTag, VolumeVariables);
    using ElementVolumeVariables = typename GET_PROP_TYPE(LowDimProblemTypeTag, ElementVolumeVariables);
    using ElementFluxVariablesCache = typename GET_PROP_TYPE(LowDimProblemTypeTag, ElementFluxVariablesCache);
    using ElementBoundaryTypes = typename GET_PROP_TYPE(LowDimProblemTypeTag, ElementBoundaryTypes);
    using LocalResidual = typename GET_PROP_TYPE(LowDimProblemTypeTag, LocalResidual);
    using GridView = typename GET_PROP_TYPE(LowDimProblemTypeTag, GridView);
    using IndexType = typename GridView::IndexSet::IndexType:
    using Element = typename GridView::template Codim<0>::Entity;

    enum { dim = GridView::dimension };
    enum { isBox = GET_PROP_VALUE(LowDimProblemTypeTag, ImplicitIsBox) };

    // matrix and solution types
    using SolutionVector = typename GET_PROP_TYPE(LowDimProblemTypeTag, SolutionVector);
    using JacobianMatrix = typename GET_PROP_TYPE(TypeTag, JacobianMatrix)::MatrixBlockLowDim;
    using JacobianMatrixCoupling = typename GET_PROP_TYPE(TypeTag, JacobianMatrix)::MatrixBlockLowDimCoupling;

public:

    // copying a local jacobian is not a good idea
    LowDimLocalJacobian(const LowDimLocalJacobian &) = delete;

    LowDimLocalJacobian()
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
        ParentType::init(globalProblem.lowDimProblem());
    }

    /*!
     * \brief Assemble an element's local Jacobian matrix of the
     *        defect.
     *
     * \param element The DUNE Codim<0> entity which we look at.
     */
    void assemble(const Element &element,
                  JacobianMatrix& matrix,
                  JacobianMatrixCoupling& couplingMatrix,
                  SolutionVector& residual)
    {
        const bool isGhost = (element.partitionType() == Dune::GhostEntity);

        // prepare the volvars/fvGeometries in case caching is disabled
        auto fvGeometry = localView(this->model_().globalFvGeometry());
        fvGeometry.bind(element);

        auto curElemVolVars = localView(this->model_().curGlobalVolVars());
        curElemVolVars.bind(element, fvGeometry, this->model_().curSol());

        auto prevElemVolVars = localView(this->model_().prevGlobalVolVars());
        prevElemVolVars.bindElement(element, fvGeometry, this->model_().prevSol());

        auto elemFluxVarsCache = localView(this->model_().globalFluxVarsCache());
        elemFluxVarsCache.bind(element, fvGeometry, curElemVolVars);

        // set the actual dof index
        globalI_ = this->problem_().elementMapper().index(element);

        // check for boundaries on the element
        ElementBoundaryTypes elemBcTypes;
        elemBcTypes.update(this->problem_(), element, fvGeometry);

        // The actual solution in the element
        auto curElemSol = this->model_().elementSolution(element, this->model_().curSol());

        // Evaluate the undeflected element local residual
        this->localResidual().eval(element,
                                   fvGeometry,
                                   prevElemVolVars,
                                   curElemVolVars,
                                   elemBcTypes,
                                   elemFluxVarsCache);
        this->residual_ = this->localResidual().residual();

        // calculate derivatives of all dofs in stencil with respect to the dofs in the element
        this->evalPartialDerivatives_(element,
                                      fvGeometry,
                                      prevElemVolVars,
                                      curElemVolVars,
                                      curElemSol,
                                      elemFluxVarsCache,
                                      elemBcTypes,
                                      matrix,
                                      residual,
                                      isGhost);

        // compute derivatives with respect to additional user defined DOF dependencies
        const auto& additionalDofDepedencies = this->problem_().getAdditionalDofDependencies(globalI_);
        if (!additionalDofDepedencies.empty() && !isGhost)
        {
            this->evalAdditionalDerivatives_(additionalDofDepedencies,
                                             element,
                                             fvGeometry,
                                             curElemVolVars,
                                             matrix,
                                             residual);
        }

        evalPartialDerivativeCoupling_(element,
                                       fvGeometry,
                                       curElemVolVars,
                                       elemFluxVarsCache,
                                       elemBcTypes,
                                       couplingMatrix,
                                       residual)
    }

protected:

    /*!
     * \brief Returns a reference to the problem.
     */
    const GlobalProblem &globalProblem_() const
    { return *globalProblemPtr_; }

    GlobalProblem &globalProblem_()
    { return *globalProblemPtr_; }

    void evalPartialDerivativeCoupling_(const Element& element,
                                        const FVElementGeometry& fvGeometry,
                                        ElementVolumeVariables& curElemVolVars,
                                        ElementFluxVariablesCache& elemFluxVarsCache,
                                        const ElementBoundaryTypes& elemBcTypes,
                                        JacobianMatrixCoupling& couplingMatrix,
                                        SolutionVector& residual)
    {
        const auto& couplingStencil = globalProblem_().couplingManager().couplingStencil(element);
        const auto numDofsCoupling = couplingStencil.size();

        for (auto globalJ : couplingStencil)
        {
            const auto lowDimElement = globalProblem_().lowDimProblem().model().globalFvGeometry().element(globalJ);
            const auto& lowDimSolution = globalProblem_().lowDimProblem().model().curSol();
            auto curElemSolLowDim = globalProblem_().lowDimProblem().model().elementSolution(lowDimElement, lowDimSolution);
            const auto lowDimResidual = globalProblem_().couplingManager().evalCouplingResidual(element,
                                                                                                fvGeometry,
                                                                                                curElemVolVars,
                                                                                                elemBcTypes,
                                                                                                elemFluxVarsCache,
                                                                                                lowDimElement,
                                                                                                curElemSolLowDim);
             // derivatives in the neighbors with repect to the current elements
            PrimaryVariables partialDeriv;
            for (int pvIdx = 0; pvIdx < numEq; pvIdx++)
            {
                const Scalar eps = numericEpsilon(curElemSolLowDim[0][pvIdx]);
                Scalar delta = 0;

                if (numericDifferenceMethod_ >= 0)
                {
                    // we are not using backward differences, i.e. we need to
                    // calculate f(x + \epsilon)

                    // deflect primary variables
                    curElemSolLowDim[0][pvIdx] += eps;
                    delta += eps;

                    // calculate the residual with the deflected primary variables
                    partialDeriv = globalProblem_().couplingManager().evalCouplingResidual(element,
                                                                                           fvGeometry,
                                                                                           curElemVolVars,
                                                                                           elemBcTypes,
                                                                                           elemFluxVarsCache,
                                                                                           lowDimElement,
                                                                                           curElemSolLowDim);
                }
                else
                {
                    // we are using backward differences, i.e. we don't need
                    // to calculate f(x + \epsilon) and we can recycle the
                    // (already calculated) residual f(x)
                    partialDeriv = lowDimResidual;
                }


                if (numericDifferenceMethod_ <= 0)
                {
                    // we are not using forward differences, i.e. we
                    // need to calculate f(x - \epsilon)

                    // deflect the primary variables
                    curElemSolLowDim[0][pvIdx] -= 2*eps;
                    delta += eps;

                    // calculate the residual with the deflected primary variables
                    partialDeriv -= globalProblem_().couplingManager().evalCouplingResidual(element,
                                                                                            fvGeometry,
                                                                                            curElemVolVars,
                                                                                            elemBcTypes,
                                                                                            elemFluxVarsCache,
                                                                                            lowDimElement,
                                                                                            curElemSolLowDim);
                }
                else
                {
                    // we are using forward differences, i.e. we don't need to
                    // calculate f(x - \epsilon) and we can recycle the
                    // (already calculated) residual f(x)
                    partialDeriv -= lowDimResidual;
                }

                // divide difference in residuals by the magnitude of the
                // deflections between the two function evaluation
                partialDeriv /= delta;

                // restore the original state of the element solution vector
                curElemSolLowDim[0][pvIdx] = lowDimSolution[globalJ][pvIdx];

                // update the global jacobian matrix (coupling block)
                this->updateGlobalJacobian_(couplingMatrix, globalI_, globalJ, pvIdx, partialDeriv);
            }
        }
    }

    // The global index of the current element
    IndexType globalI_;
    // The problem we would like to solve
    GlobalProblem *globalProblemPtr_;
    // The type of the numeric difference method (forward, center, backward)
    int numericDifferenceMethod_;
};

} // end namespace Dumux

#endif
