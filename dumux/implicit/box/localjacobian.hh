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
 * \brief Caculates the Jacobian of the local residual for fully-implicit box models
 */
#ifndef DUMUX_IMPLICIT_BOX_LOCAL_JACOBIAN_HH
#define DUMUX_IMPLICIT_BOX_LOCAL_JACOBIAN_HH

#include <dune/istl/io.hh>
#include <dune/istl/matrix.hh>

#include <dumux/common/math.hh>
#include <dumux/common/valgrind.hh>

#include <dumux/implicit/properties.hh>
#include <dumux/implicit/localjacobian.hh>

namespace Dumux
{
/*!
 * \ingroup ImplicitLocalJacobian
 * \brief Calculates the Jacobian of the local residual for fully-implicit models
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
class BoxLocalJacobian : public ImplicitLocalJacobian<TypeTag>
{
    using ParentType = ImplicitLocalJacobian<TypeTag>;
    using Implementation = typename GET_PROP_TYPE(TypeTag, LocalJacobian);
    using LocalResidual = typename GET_PROP_TYPE(TypeTag, LocalResidual);
    using Problem = typename GET_PROP_TYPE(TypeTag, Problem);
    using Model = typename GET_PROP_TYPE(TypeTag, Model);
    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
    using Element = typename GridView::template Codim<0>::Entity;
    using JacobianAssembler = typename GET_PROP_TYPE(TypeTag, JacobianAssembler);

    enum {
        numEq = GET_PROP_VALUE(TypeTag, NumEq),
        dim = GridView::dimension,
    };

    using JacobianMatrix = typename GET_PROP_TYPE(TypeTag, JacobianMatrix);
    using SolutionVector = typename GET_PROP_TYPE(TypeTag, SolutionVector);
    using FVElementGeometry = typename GET_PROP_TYPE(TypeTag, FVElementGeometry);
    using SubControlVolume = typename GET_PROP_TYPE(TypeTag, SubControlVolume);
    using VertexMapper = typename GET_PROP_TYPE(TypeTag, VertexMapper);
    using ElementSolutionVector = typename GET_PROP_TYPE(TypeTag, ElementSolutionVector);
    using PrimaryVariables = typename GET_PROP_TYPE(TypeTag, PrimaryVariables);
    using VolumeVariables = typename GET_PROP_TYPE(TypeTag, VolumeVariables);
    using ElementVolumeVariables = typename GET_PROP_TYPE(TypeTag, ElementVolumeVariables);
    using ElementBoundaryTypes = typename GET_PROP_TYPE(TypeTag, ElementBoundaryTypes);
    using ElementFluxVariablesCache = typename GET_PROP_TYPE(TypeTag, ElementFluxVariablesCache);
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);

public:
    BoxLocalJacobian()
    {
        numericDifferenceMethod_ = GET_PARAM_FROM_GROUP(TypeTag, int, Implicit, NumericDifferenceMethod);
    }

    /*!
     * \brief Assemble an element's local Jacobian matrix of the
     *        defect.
     *
     * \param element The DUNE Codim<0> entity which we look at.
     */
    void assemble(const Element& element,
                  JacobianMatrix& matrix,
                  SolutionVector& residual)
    {
        // prepare the volvars/fvGeometries for the case when caching is disabled
        auto fvGeometry = localView(this->problem_().model().globalFvGeometry());
        fvGeometry.bind(element);

        auto curElemVolVars = localView(this->model_().curGlobalVolVars());
        curElemVolVars.bind(element, fvGeometry, this->model_().curSol());

        auto prevElemVolVars = localView(this->model_().prevGlobalVolVars());
        prevElemVolVars.bindElement(element, fvGeometry, this->model_().prevSol());

        auto elemFluxVarsCache = localView(this->model_().globalFluxVarsCache());
        elemFluxVarsCache.bind(element, fvGeometry, curElemVolVars);

        // the element boundary types
        ElementBoundaryTypes elemBcTypes;
        elemBcTypes.update(this->problem_(), element, fvGeometry);

        // the element solution
        auto curElemSol = this->model_().elementSolution(element, this->model_().curSol());

        // calculate the actual element residual
        this->localResidual().eval(element, fvGeometry, prevElemVolVars, curElemVolVars, elemBcTypes, elemFluxVarsCache);
        this->residual_ = this->localResidual().residual();

        this->model_().updatePVWeights(fvGeometry);

        // calculation of the derivatives
        for (auto&& scv : scvs(fvGeometry))
        {
            // dof index and corresponding actual pri vars
            const auto dofIdx = scv.dofIndex();
            auto& curVolVars = getCurVolVars(curElemVolVars, scv);
            VolumeVariables origVolVars(curVolVars);

            // add precalculated residual for this scv into the global container
            residual[dofIdx] += this->residual_[scv.indexInElement()];

            // calculate derivatives w.r.t to the privars at the dof at hand
            for (int pvIdx = 0; pvIdx < numEq; pvIdx++)
            {
                evalPartialDerivative_(matrix,
                                       element,
                                       fvGeometry,
                                       prevElemVolVars,
                                       curElemVolVars,
                                       curElemSol,
                                       scv,
                                       elemBcTypes,
                                       elemFluxVarsCache,
                                       pvIdx);

                // restore the original state of the scv's volume variables
                curVolVars = origVolVars;

                // restore the original element solution
                curElemSol[scv.indexInElement()][pvIdx] = this->model_().curSol()[scv.dofIndex()][pvIdx];
            }

            // TODO: what if we have an extended source stencil????
        }
    }
protected:
    /*!
     * \brief Compute the partial derivatives to a primary variable at
     *        an degree of freedom.
     *
     * This method can be overwritten by the implementation if a
     * better scheme than numerical differentiation is available.
     *
     * The default implementation of this method uses numeric
     * differentiation, i.e. forward or backward differences (2nd
     * order), or central differences (3rd order). The method used is
     * determined by the "NumericDifferenceMethod" property:
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
     *
     * \param partialDeriv The vector storing the partial derivatives of all
     *              equations
     * \param storageDeriv the mass matrix contributions
     * \param col The block column index of the degree of freedom
     *            for which the partial derivative is calculated.
     *            Box: a sub-control volume index.
     *            Cell centered: a neighbor index.
     * \param pvIdx The index of the primary variable
     *              for which the partial derivative is calculated
     */
    void evalPartialDerivative_(JacobianMatrix& matrix,
                                const Element& element,
                                const FVElementGeometry& fvGeometry,
                                const ElementVolumeVariables& prevElemVolVars,
                                ElementVolumeVariables& curElemVolVars,
                                ElementSolutionVector& curElemSol,
                                const SubControlVolume& scv,
                                const ElementBoundaryTypes& elemBcTypes,
                                const ElementFluxVariablesCache& elemFluxVarsCache,
                                const int pvIdx)
    {
        const auto dofIdx = scv.dofIndex();
        auto& volVars = getCurVolVars(curElemVolVars, scv);

        ElementSolutionVector partialDeriv(element.subEntities(dim));
        Scalar eps = this->numericEpsilon(volVars.priVar(pvIdx));
        Scalar delta = 0;

        // calculate the residual with the forward deflected primary variables
        if (numericDifferenceMethod_ >= 0)
        {
            // we are not using backward differences, i.e. we need to
            // calculate f(x + \epsilon)

            // deflect primary variables
            curElemSol[scv.indexInElement()][pvIdx] += eps;
            delta += eps;

            // update the volume variables connected to the dof
            volVars.update(curElemSol, this->problem_(), element, scv);

            // calculate the deflected residual
            this->localResidual().eval(element, fvGeometry, prevElemVolVars, curElemVolVars, elemBcTypes, elemFluxVarsCache);

            // store the residual
            partialDeriv = this->localResidual().residual();
        }
        else
        {
            // we are using backward differences, i.e. we don't need
            // to calculate f(x + \epsilon) and we can recycle the
            // (already calculated) residual f(x)
            partialDeriv = this->residual_;
        }

        if (numericDifferenceMethod_ <= 0)
        {
            // we are not using forward differences, i.e. we
            // need to calculate f(x - \epsilon)

            // deflect the primary variables
            curElemSol[scv.indexInElement()][pvIdx] -= delta + eps;
            delta += eps;

            // update the volume variables connected to the dof
            volVars.update(curElemSol, this->problem_(), element, scv);

            // calculate the deflected residual
            this->localResidual().eval(element, fvGeometry, prevElemVolVars, curElemVolVars, elemBcTypes, elemFluxVarsCache);

            // subtract the residual from the derivative storage
            partialDeriv -= this->localResidual().residual();
        }
        else
        {
            // we are using forward differences, i.e. we don't need to
            // calculate f(x - \epsilon) and we can recycle the
            // (already calculated) residual f(x)
            partialDeriv -= this->residual_;
        }

        // divide difference in residuals by the magnitude of the
        // deflections between the two function evaluation
        partialDeriv /= delta;

        // update the global stiffness matrix with the current partial derivatives
        for (auto&& scvJ : scvs(fvGeometry))
            this->updateGlobalJacobian_(matrix, scvJ.dofIndex(), dofIdx, pvIdx, partialDeriv[scvJ.index()]);
    }

    //! If the global vol vars caching is enabled we have to modify the global volvar object
    template<class T = TypeTag>
    typename std::enable_if<GET_PROP_VALUE(T, EnableGlobalVolumeVariablesCache), VolumeVariables>::type&
    getCurVolVars(ElementVolumeVariables& elemVolVars, const SubControlVolume& scv)
    { return this->model_().nonConstCurGlobalVolVars().volVars(scv.elementIndex(), scv.indexInElement()); }

    //! When global volume variables caching is disabled, return the local volvar object
    template<class T = TypeTag>
    typename std::enable_if<!GET_PROP_VALUE(T, EnableGlobalVolumeVariablesCache), VolumeVariables>::type&
    getCurVolVars(ElementVolumeVariables& elemVolVars, const SubControlVolume& scv)
    { return elemVolVars[scv]; }

private:
    int numericDifferenceMethod_;
};

}

#endif
