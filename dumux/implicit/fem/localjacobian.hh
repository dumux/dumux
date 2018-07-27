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
 * \brief Caculates the Jacobian of the local residual for fully-implicit finite element models
 */
#ifndef DUMUX_IMPLICIT_FEM_LOCAL_JACOBIAN_HH
#define DUMUX_IMPLICIT_FEM_LOCAL_JACOBIAN_HH

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
class FemLocalJacobian : public ImplicitLocalJacobian<TypeTag>
{
    using ParentType = ImplicitLocalJacobian<TypeTag>;

    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using FEBasis = typename GET_PROP_TYPE(TypeTag, FeBasis);
    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
    using JacobianMatrix = typename GET_PROP_TYPE(TypeTag, JacobianMatrix);
    using SolutionVector = typename GET_PROP_TYPE(TypeTag, SolutionVector);
    using PrimaryVariables = typename GET_PROP_TYPE(TypeTag, PrimaryVariables);
    using ElementSolutionVector = typename GET_PROP_TYPE(TypeTag, ElementSolutionVector);

    static constexpr int dim = GridView::dimension;
    static constexpr int numEq = GET_PROP_VALUE(TypeTag, NumEq);

    using Element = typename GridView::template Codim<0>::Entity;
    using LocalView = typename FEBasis::LocalView;
    using LocalIndexSet = typename FEBasis::LocalIndexSet;

public:
    FemLocalJacobian()
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

//std::cout << "Wir sind in femLocalJacobian" << std::endl;

        // prepare the element solutions etc...
        const auto& curSol = this->model_().curSol();
        const auto& prevSol = this->model_().prevSol();

        //                printvector(std::cout, curSol, "LocJacCurSol", "");
        //                printvector(std::cout, prevSol, "LocJacPrevSol", "");


        // prepare the current and previous element solutions
        auto localView = this->model_().feBasis().localView();
        auto localIndexSet = this->model_().feBasis().localIndexSet();
        localView.bind(element);
        localIndexSet.bind(localView);


        auto numLocalDofs = localView.tree().finiteElement().localBasis().size();
        ElementSolutionVector curElemSol(numLocalDofs);
        ElementSolutionVector prevElemSol(numLocalDofs);
        for (unsigned int i = 0; i < numLocalDofs; i++)
        {
            auto dofIdxGlobal = localIndexSet.index(i);
            curElemSol[i] = curSol[dofIdxGlobal];
            prevElemSol[i] = prevSol[dofIdxGlobal];
        }

        //                printvector(std::cout, curElemSol, "LocJacCurElemSol", "");
        //                printvector(std::cout, prevElemSol, "LocJacPrevElemSol", "");

//std::cout << "locJac nach localResidual().eval()" << std::endl;

        // calculate the actual element residual
        this->localResidual().eval(element, localView, localIndexSet, curElemSol, prevElemSol);
        this->residual_ = this->localResidual().residual();

            //std::cout << "locJac nach localResidual().eval()" << std::endl;

        //printvector(std::cout, residual, "locResLocJacobian","");

        // TODO: WHAT IS THIS?
        // this->model_().updatePVWeights(fvGeometry);

//std::cout << "locJacMark: " << std::endl;

        // calculation of the derivatives
        for (unsigned int i = 0; i < numLocalDofs; ++i)
        {
            // store actual pri vars
            auto dofIdxGlobal = localIndexSet.index(i);

            // add precalculated residual for this scv into the global container
            residual[dofIdxGlobal] += this->residual(i);


                //            printvector(std::cout, residual, "LocJacRes", "");


            // calculate derivatives w.r.t to the privars at the dof at hand
            for (int pvIdx = 0; pvIdx < numEq; pvIdx++)
            {
//std::cout << "locJac - numLocalDof: " << i << "    numEq: " << pvIdx <<  std::endl;
                evalPartialDerivative_(matrix,
                                       element,
                                       localView,
                                       localIndexSet,
                                       prevElemSol,
                                       curElemSol,
                                       i,
                                       pvIdx);


                    //                printvector(std::cout, curElemSol, "curElemSol", "");

                // restore the original state of the primary variables
                curElemSol[i][pvIdx] = curSol[dofIdxGlobal][pvIdx];

            }

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
                                const LocalView& localView,
                                const LocalIndexSet& localIndexSet,
                                ElementSolutionVector& prevElemSol,
                                ElementSolutionVector& curElemSol,
                                const int localDofIdx,
                                const int pvIdx)
    {
        const auto numLocalDofs = localView.tree().size();
        const auto dofIdxGlobal = localIndexSet.index(localDofIdx);

        ElementSolutionVector partialDeriv(numLocalDofs);
        Scalar eps = this->numericEpsilon(curElemSol[localDofIdx][pvIdx]);
        Scalar delta = 0;

        //        std::cout << "LocJacEps: " << eps << std::endl;

        //printmatrix(std::cout, matrix, "LocJacVorAsssembleMatrix", "", 7, 0);


        // calculate the residual with the forward deflected primary variables
        if (numericDifferenceMethod_ >= 0)
        {
            // we are not using backward differences, i.e. we need to
            // calculate f(x + \epsilon)

//printvector(std::cout, curElemSol, "curElemSolVorEpsVorwaerts", "");


            // deflect primary variables
            curElemSol[localDofIdx][pvIdx] += eps;
            delta += eps;

                 //printvector(std::cout, curElemSol, "LocJacCurElemSolnachEpsVorwaerts", "");


            // calculate the deflected residual
            this->localResidual().eval(element, localView, localIndexSet, curElemSol, prevElemSol);

            // store the residual
            partialDeriv = this->localResidual().residual();

//printvector(std::cout, partialDeriv, "LocJacPartialDerivVorwaerts", "");
        }
        else
        {
            // we are using backward differences, i.e. we don't need
            // to calculate f(x + \epsilon) and we can recycle the
            // (already calculated) residual f(x)
            partialDeriv = this->residual_;

//printvector(std::cout, partialDeriv, "LocJacPartialDerivElse1", "");
        }

        if (numericDifferenceMethod_ <= 0)
        {
            // we are not using forward differences, i.e. we
            // need to calculate f(x - \epsilon)


                 //    printvector(std::cout, curElemSol, "curElemSolVorEpsRueck", "");


            // deflect the primary variables
            curElemSol[localDofIdx][pvIdx] -= delta + eps;
            delta += eps;

                //printvector(std::cout, curElemSol, "locJacCurElemSolnachEpsRueck", "");


            // calculate the deflected residual
            this->localResidual().eval(element, localView, localIndexSet, curElemSol, prevElemSol);

            // subtract the residual from the derivative storage
            partialDeriv -= this->localResidual().residual();

//printvector(std::cout, partialDeriv, "partialDerivRueck", "");
        }
        else
        {
            // we are using forward differences, i.e. we don't need to
            // calculate f(x - \epsilon) and we can recycle the
            // (already calculated) residual f(x)
            partialDeriv -= this->residual_;

//printvector(std::cout, partialDeriv, "partialDerivElse2", "");
        }

//printvector(std::cout, partialDeriv, "LocJacPartialDeriv", "");


        // divide difference in residuals by the magnitude of the
        // deflections between the two function evaluation
        partialDeriv /= delta;

//printvector(std::cout, this->residual_, "LocJacPartialthis->residual_", "");
//printvector(std::cout, partialDeriv, "LocJacPartialDeriv/delta", "");


        //printmatrix(std::cout, matrix, "matrix", "");
        //std::cout << "delta: " << delta << std::endl;
        //std::cout << "eps: " << eps << std::endl;

        //        std::cout << std::endl;

        // update the global stiffness matrix with the current partial derivatives
        for (unsigned int i = 0; i < numLocalDofs; ++i)
            this->updateGlobalJacobian_(matrix, localIndexSet.index(i), dofIdxGlobal, pvIdx, partialDeriv[i]);

//printmatrix(std::cout, matrix, "locJacMatrix", "", 6, 0);

    }

private:
    int numericDifferenceMethod_;
};

}

#endif
