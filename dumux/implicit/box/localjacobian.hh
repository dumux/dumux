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
 * \brief Caculates the Jacobian of the local residual for fully-implicit models
 */
#ifndef DUMUX_IMPLICIT_LOCAL_JACOBIAN_HH
#define DUMUX_IMPLICIT_LOCAL_JACOBIAN_HH

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
    typedef ImplicitLocalJacobian<TypeTag> ParentType;
    typedef typename GET_PROP_TYPE(TypeTag, LocalJacobian) Implementation;
    typedef typename GET_PROP_TYPE(TypeTag, LocalResidual) LocalResidual;
    typedef typename GET_PROP_TYPE(TypeTag, Problem) Problem;
    typedef typename GET_PROP_TYPE(TypeTag, Model) Model;
    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
    typedef typename GridView::template Codim<0>::Entity Element;
    typedef typename GET_PROP_TYPE(TypeTag, JacobianAssembler) JacobianAssembler;

    enum {
        numEq = GET_PROP_VALUE(TypeTag, NumEq),
        dim = GridView::dimension,
    };

    typedef typename GET_PROP_TYPE(TypeTag, FVElementGeometry) FVElementGeometry;
    typedef typename GET_PROP_TYPE(TypeTag, SubControlVolume) SubControlVolume;
    typedef typename GET_PROP_TYPE(TypeTag, VertexMapper) VertexMapper;
    typedef typename GET_PROP_TYPE(TypeTag, ElementSolutionVector) ElementSolutionVector;
    typedef typename GET_PROP_TYPE(TypeTag, PrimaryVariables) PrimaryVariables;

    typedef typename GET_PROP_TYPE(TypeTag, VolumeVariables) VolumeVariables;
    typedef typename GET_PROP_TYPE(TypeTag, ElementBoundaryTypes) ElementBoundaryTypes;

    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef Dune::FieldMatrix<Scalar, numEq, numEq> MatrixBlock;
    typedef Dune::Matrix<MatrixBlock> LocalBlockMatrix;

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
    void assemble(const Element& element, JacobianMatrix& matrix)
    {
        // set the current grid element and update the element's
        // finite volume geometry
        elemPtr_ = &element;
        reset_();

        bcTypes_.update(problem_(), element, fvElemGeom_());

        // calculate the local residual
        this->localResidual().eval(element, bcTypes_);
        this->residual_ = this->localResidual().residual();

        this->model_().updatePVWeights(fvElemGeom_());

        // get stencil informations
        const auto& elementStencil = this->model_().stencils(element).elementStencil();

        // calculate derivatives for the dofs inside the element
        ElementSolutionVector partialDeriv(numRows);
        for (auto&& scv : fvElemGeom_().scvs())
        {
            for (int pvIdx = 0; pvIdx < numEq; pvIdx++)
            {
                asImp_().evalPartialDerivative_(partialDeriv, element, scv, pvIdx);

                // update the local stiffness matrix with the current partial derivatives
                for (auto&& scvJ : fvElemGeom_().scvs())
                    this->updateGlobalJacobian_(matrix, scv.dofIndex(), scvJ.dofIndex(), pvIdx, partialDeriv[scvJ.indexInElement()]);
            }
        }

        // TODO: calculate derivatives in the case of an extended source stencil
        // const auto& extendedSourceStencil = model_().stencils(element).extendedSourceStencil();
        // for (auto&& globalJ : extendedSourceStencil)
        // {
        //     for (int pvIdx = 0; pvIdx < numEq; pvIdx++)
        //         {
        //             evalPartialDerivativeSource_(partialDeriv, globalJ, pvIdx, fluxVarIndices);

        //             // update the local stiffness matrix with the partial derivatives
        //             updateLocalJacobian_(j, pvIdx, partialDeriv);
        //         }
               // ++j;
        // }

        // for cellcentered methods, calculate the derivatives w.r.t cells in stencil
    }

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
     * \param col The block column index of the degree of freedom
     *            for which the partial derivative is calculated.
     *            Box: a sub-control volume index.
     *            Cell centered: a neighbor index.
     * \param pvIdx The index of the primary variable
     *              for which the partial derivative is calculated
     */
    void evalPartialDerivative_(ElementSolutionVector &partialDeriv,
                                const Element& element,
                                const SubControlVolume &scv,
                                const int pvIdx)
    {
        auto dofIdxGlobal = scv.dofIndex();

        PrimaryVariables priVars(model_().curSol()[dofIdxGlobal]);
        VolumeVariables origVolVars(model_().curVolVars(scv));

        Scalar eps = asImp_().numericEpsilon(scv, pvIdx);
        Scalar delta = 0;

        if (numericDifferenceMethod_ >= 0)
        {
            // we are not using backward differences, i.e. we need to
            // calculate f(x + \epsilon)

            // deflect primary variables
            priVars[pvIdx] += eps;
            delta += eps;

            // update the volume variables
            model_().curVolVars_(scv).update(priVars, problem_(), element, scv);

            // calculate the residual with the deflected primary variables
            localResidual().eval(element, bcTypes_);

            // store the residual and the storage term
            partialDeriv = localResidual().residual();
        }
        else
        {
            // we are using backward differences, i.e. we don't need
            // to calculate f(x + \epsilon) and we can recycle the
            // (already calculated) residual f(x)
            partialDeriv = residual_;
        }

        if (numericDifferenceMethod_ <= 0)
        {
            // we are not using forward differences, i.e. we
            // need to calculate f(x - \epsilon)

            // deflect the primary variables
            priVars[pvIdx] -= delta + eps;
            delta += eps;

            // update the volume variables
            model_().curVolVars_(scv).update(priVars, problem_(), element, scv);

            // calculate the residual with the deflected primary variables
            localResidual().eval(element, bcTypes_);

            // subtract the residual from the derivative storage
            partialDeriv -= localResidual().residual();
        }
        else
        {
            // we are using forward differences, i.e. we don't need to
            // calculate f(x - \epsilon) and we can recycle the
            // (already calculated) residual f(x)
            partialDeriv -= residual_;
        }

        // divide difference in residuals by the magnitude of the
        // deflections between the two function evaluation
        partialDeriv /= delta;

        // restore the original state of the scv's volume variables
        model_().curVolVars_(scv) = origVolVars;

#if HAVE_VALGRIND
        for (unsigned i = 0; i < partialDeriv.size(); ++i)
            Valgrind::CheckDefined(partialDeriv[i]);
#endif
    }

    const FVElementGeometry& fvElemGeom_() const
    {
        return model_().fvGeometries(eIdxGlobal_);
    }

    IndexType eIdxGlobal_;
    ElementBoundaryTypes bcTypes_;

    int numericDifferenceMethod_;
};
}

#endif
