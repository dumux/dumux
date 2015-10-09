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
 * \brief Calculates the partial derivatives of the local residual for the Jacobian of the
 *           one-phase two-component linear elasticity model.
 */
#ifndef DUMUX_EL1P2C_LOCAL_JACOBIAN_HH
#define DUMUX_EL1P2C_LOCAL_JACOBIAN_HH

#include <dune/common/version.hh>
#include <dumux/implicit/common/implicitlocaljacobian.hh>

namespace Dumux
{
/*!
 * \ingroup ElOnePTwoCBoxModel
 * \brief Calculates the partial derivatives of the local residual for the Jacobian
 *
 *  Except for the evalPartialDerivatives function all functions are taken from the
 *  base class ImplicitLocalJacobian
 */
template<class TypeTag>
class ElOnePTwoCLocalJacobian : public ImplicitLocalJacobian<TypeTag>
{
private:
    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
    enum {
        dim = GridView::dimension,
    };
    typedef typename GET_PROP_TYPE(TypeTag, FVElementGeometry) FVElementGeometry;
    typedef typename GET_PROP_TYPE(TypeTag, ElementSolutionVector) ElementSolutionVector;
    typedef typename GET_PROP_TYPE(TypeTag, PrimaryVariables) PrimaryVariables;
    typedef typename GET_PROP_TYPE(TypeTag, VolumeVariables) VolumeVariables;
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    enum { isBox = GET_PROP_VALUE(TypeTag, ImplicitIsBox) };

    // copying a local jacobian is not a good idea
    ElOnePTwoCLocalJacobian(const ElOnePTwoCLocalJacobian &);

public:
    ElOnePTwoCLocalJacobian()
    {}

    /*!
     * \brief Compute the partial derivatives to a primary variable at
     *        an degree of freedom.
     *
     * This method is overwritten here since this model requires a call of the model specific
     * elementvolumevariables which updates the effective porosities correctly.
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
    void evalPartialDerivative_(ElementSolutionVector &partialDeriv,
                                PrimaryVariables &storageDeriv,
                                const int col,
                                const int pvIdx)
    {
        int dofIdxGlobal;
        FVElementGeometry neighborFVGeom;
        auto neighbor = this->element_();
        if (isBox)
        {
#if DUNE_VERSION_NEWER(DUNE_COMMON, 2, 4)
            dofIdxGlobal = this->vertexMapper_().subIndex(this->element_(), col, dim);
#else
            dofIdxGlobal = this->vertexMapper_().map(this->element_(), col, dim);
#endif

        }
        else
        {
            neighbor = this->fvElemGeom_.neighbors[col];
            neighborFVGeom.updateInner(neighbor);
#if DUNE_VERSION_NEWER(DUNE_COMMON, 2, 4)
            dofIdxGlobal = this->problemPtr_->elementMapper().index(neighbor);
#else
            dofIdxGlobal = this->problemPtr_->elementMapper().map(neighbor);
#endif

        }

        PrimaryVariables priVars(this->model_().curSol()[dofIdxGlobal]);
        VolumeVariables origVolVars(this->curVolVars_[col]);

        this->curVolVars_[col].setEvalPoint(&origVolVars);
        Scalar eps = this->numericEpsilon(col, pvIdx);
        Scalar delta = 0;

        if (this->numericDifferenceMethod_ >= 0) {
            // we are not using backward differences, i.e. we need to
            // calculate f(x + \epsilon)

            // deflect primary variables
            priVars[pvIdx] += eps;
            delta += eps;

            // calculate the residual
            if (isBox){
                this->curVolVars_[col].update(priVars,
                        this->problem_(),
                        this->element_(),
                        this->fvElemGeom_,
                                        col,
                                        false);
                // update the effective porosities
                this->curVolVars_.updateEffPorosity(this->problem_(),
                        this->element_(),
                        this->fvElemGeom_);
            }
            else{
                this->curVolVars_[col].update(priVars,
                        this->problem_(),
                                        neighbor,
                                        neighborFVGeom,
                                        /*scvIdx=*/0,
                                        false);
                // update the effective porosities
                this->curVolVars_.updateEffPorosity(this->problem_(),
                        this->element_(),
                        this->fvElemGeom_);
            }

            this->localResidual().eval(this->element_(),
                    this->fvElemGeom_,
                    this->prevVolVars_,
                                 this->curVolVars_,
                                 this->bcTypes_);

            // store the residual and the storage term
            partialDeriv = this->localResidual().residual();
            if (isBox || col == 0)
                storageDeriv = this->localResidual().storageTerm()[col];
        }
        else {
            // we are using backward differences, i.e. we don't need
            // to calculate f(x + \epsilon) and we can recycle the
            // (already calculated) residual f(x)
            partialDeriv = this->residual_;
            storageDeriv = this->storageTerm_[col];
        }


        if (this->numericDifferenceMethod_ <= 0) {
            // we are not using forward differences, i.e. we don't
            // need to calculate f(x - \epsilon)

            // deflect the primary variables
            priVars[pvIdx] -= delta + eps;
            delta += eps;

            // calculate residual again
            if (isBox){
                this->curVolVars_[col].update(priVars,
                        this->problem_(),
                        this->element_(),
                        this->fvElemGeom_,
                                        col,
                                        false);
                // update the effective porosities
                this->curVolVars_.updateEffPorosity(this->problem_(),
                        this->element_(),
                        this->fvElemGeom_);
            }
            else{
                this->curVolVars_[col].update(priVars,
                        this->problem_(),
                                        neighbor,
                                        neighborFVGeom,
                                        /*scvIdx=*/0,
                                        false);
                // update the effective porosities
                this->curVolVars_.updateEffPorosity(this->problem_(),
                        this->element_(),
                        this->fvElemGeom_);
            }

            this->localResidual().eval(this->element_(),
                    this->fvElemGeom_,
                    this->prevVolVars_,
                                 this->curVolVars_,
                                 this->bcTypes_);
            partialDeriv -= this->localResidual().residual();
            if (isBox || col == 0)
                storageDeriv -= this->localResidual().storageTerm()[col];
        }
        else {
            // we are using forward differences, i.e. we don't need to
            // calculate f(x - \epsilon) and we can recycle the
            // (already calculated) residual f(x)
            partialDeriv -= this->residual_;
            if (isBox || col == 0)
                storageDeriv -= this->storageTerm_[col];
        }

        // divide difference in residuals by the magnitude of the
        // deflections between the two function evaluation
        partialDeriv /= delta;
        storageDeriv /= delta;

        // restore the original state of the element's volume variables
        this->curVolVars_[col] = origVolVars;
        // update the effective porosities
        this->curVolVars_.updateEffPorosity(this->problem_(),
                this->element_(),
                this->fvElemGeom_);

#if HAVE_VALGRIND
        for (unsigned i = 0; i < partialDeriv.size(); ++i)
            Valgrind::CheckDefined(partialDeriv[i]);
#endif
    }
};
}

#endif
