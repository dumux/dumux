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
#ifndef DUMUX_CC_LOCAL_JACOBIAN_HH
#define DUMUX_CC_LOCAL_JACOBIAN_HH

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
class CCLocalJacobian : public ImplicitLocalJacobian<TypeTag>
{
    typedef ImplicitLocalJacobian<TypeTag> ParentType;
    typedef typename GET_PROP_TYPE(TypeTag, LocalJacobian) Implementation;
    typedef typename GET_PROP_TYPE(TypeTag, LocalResidual) LocalResidual;
    typedef typename GET_PROP_TYPE(TypeTag, Problem) Problem;
    typedef typename GET_PROP_TYPE(TypeTag, Model) Model;
    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
    typedef typename GET_PROP_TYPE(TypeTag, JacobianMatrix) JacobianMatrix;

    typedef typename GET_PROP_TYPE(TypeTag, FVElementGeometry) FVElementGeometry;
    typedef typename GET_PROP_TYPE(TypeTag, SubControlVolume) SubControlVolume;
    typedef typename GET_PROP_TYPE(TypeTag, VertexMapper) VertexMapper;
    typedef typename GET_PROP_TYPE(TypeTag, PrimaryVariables) PrimaryVariables;
    typedef typename GET_PROP_TYPE(TypeTag, VolumeVariables) VolumeVariables;
    typedef typename GET_PROP_TYPE(TypeTag, ElementBoundaryTypes) ElementBoundaryTypes;

    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;

    typedef typename GridView::template Codim<0>::Entity Element;
    typedef typename GridView::IndexSet::IndexType IndexType;

    enum { numEq = GET_PROP_VALUE(TypeTag, NumEq) };

    typedef std::vector<std::vector<std::vector<IndexType>>> AssemblyMap;

public:

    CCLocalJacobian() : ParentType()
    {
        numericDifferenceMethod_ = GET_PARAM_FROM_GROUP(TypeTag, int, Implicit, NumericDifferenceMethod);
    }

    /*!
     * \brief Initialize the local Jacobian object.
     *
     * At this point we can assume that everything has been allocated,
     * although some objects may not yet be completely initialized.
     *
     * \param problem The problem which we want to simulate.
     */
    void init(Problem &problem)
    {
        ParentType::init(problem);

        assemblyMap_.resize(problem.gridView().size(0));
        for (const auto& element : elements(problem.gridView()))
        {
            auto globalI = problem.elementMapper().index(element);
            const auto& neighborStencil = this->model_().stencils(element).neighborStencil();

            assemblyMap_[globalI].reserve(neighborStencil.size());

            for (auto&& globalJ : neighborStencil)
            {
                // find the flux vars needed for the calculation of the flux into element
                std::vector<IndexType> fluxVarIndices;
                for (auto&& scvFaceJ : this->model_().fvGeometries(globalJ).scvfs())
                {
                    auto fluxVarsIdx = scvFaceJ.index();

                    // if globalI is in flux var stencil, add to list
                    const auto& fluxStencil = this->model_().fluxVars(fluxVarsIdx).stencil();
                    for (auto globalIdx : fluxStencil)
                        if (globalIdx == globalI)
                            fluxVarIndices.push_back(fluxVarsIdx);
                }

                assemblyMap_[globalI].emplace_back(std::move(fluxVarIndices));
            }
        }
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
        globalI_ = this->problem_().elementMapper().index(element);

        bcTypes_.update(this->problem_(), element, fvElemGeom_());

        // calculate the local residual
        this->localResidual().eval(element, bcTypes_);
        this->residual_ = this->localResidual().residual();

        this->model_().updatePVWeights(fvElemGeom_());

        // calculate derivatives of all dofs in stencil with respect to the dofs in the element
        evalPartialDerivatives_(element, matrix);

        // TODO: calculate derivatives in the case of an extended source stencil
        // const auto& extendedSourceStencil = model_().stencils(element).extendedSourceStencil();
        // for (auto&& globalJ : extendedSourceStencil)
        // {
        //     for (int pvIdx = 0; pvIdx < numEq; pvIdx++)
        //         {
        //             evalPartialDerivativeSource_(partialDeriv, globalJ, pvIdx, neighborToFluxVars[globalJ]);

        //             // update the local stiffness matrix with the partial derivatives
        //             updateLocalJacobian_(j, pvIdx, partialDeriv);
        //         }
               // ++j;
        // }
    }

    void evalPartialDerivatives_(const Element& element, JacobianMatrix& matrix)
    {
        // get stencil informations
        const auto& neighborStencil = this->model_().stencils(element).neighborStencil();
        const auto numNeighbors = neighborStencil.size();

        // derivatives in the element and the neighbors
        PrimaryVariables partialDeriv(0.0);
        Dune::BlockVector<PrimaryVariables> neighborDeriv(numNeighbors);

        // the localresidual class used for the flux calculations
        LocalResidual localRes;
        localRes.init(this->problem_());

        // container to store the neighboring elements
        std::vector<Element> neighborElements;
        neighborElements.reserve(numNeighbors);

        // get the elements and calculate the flux into the element in the undeflected state
        Dune::BlockVector<PrimaryVariables> origFlux(numNeighbors);
        origFlux = 0.0;

        unsigned int j = 0;
        for (auto globalJ : neighborStencil)
        {
            auto elementJ = this->problem_().model().fvGeometries().element(globalJ);
            neighborElements.push_back(elementJ);

            for (auto fluxVarIdx : assemblyMap_[globalI_][j])
                origFlux[j] += localRes.evalFlux_(elementJ, fluxVarIdx);

            ++j;
        }

        const auto& scv = this->model_().fvGeometries().subControlVolume(globalI_);
        VolumeVariables origVolVars(this->model_().curVolVars(scv));

        for (int pvIdx = 0; pvIdx < numEq; pvIdx++)
        {
            neighborDeriv = 0.0;
            PrimaryVariables priVars(this->model_().curSol()[globalI_]);

            Scalar eps = this->numericEpsilon(scv, pvIdx);
            Scalar delta = 0;

            if (numericDifferenceMethod_ >= 0)
            {
                // we are not using backward differences, i.e. we need to
                // calculate f(x + \epsilon)

                // deflect primary variables
                priVars[pvIdx] += eps;
                delta += eps;

                // update the volume variables
                this->model_().curVolVars_(scv).update(priVars, this->problem_(), element, scv);

                // calculate the residual with the deflected primary variables
                this->localResidual().eval(element, bcTypes_);

                // store the residual and the storage term
                partialDeriv = this->localResidual().residual(0);

                // calculate the fluxes into element with the deflected primary variables
                for (std::size_t k = 0; k < numNeighbors; ++k)
                {
                    for (auto fluxVarIdx : assemblyMap_[globalI_][k])
                        neighborDeriv[k] += localRes.evalFlux_(neighborElements[k], fluxVarIdx);
                }
            }
            else
            {
                // we are using backward differences, i.e. we don't need
                // to calculate f(x + \epsilon) and we can recycle the
                // (already calculated) residual f(x)
                partialDeriv = this->residual_[0];
                neighborDeriv = origFlux;
            }

            if (numericDifferenceMethod_ <= 0)
            {
                // we are not using forward differences, i.e. we
                // need to calculate f(x - \epsilon)

                // deflect the primary variables
                priVars[pvIdx] -= delta + eps;
                delta += eps;

                // update the volume variables
                this->model_().curVolVars_(scv).update(priVars, this->problem_(), element, scv);

                // calculate the residual with the deflected primary variables
                this->localResidual().eval(element, bcTypes_);

                // subtract the residual from the derivative storage
                partialDeriv -= this->localResidual().residual(0);

                // calculate the fluxes into element with the deflected primary variables
                for (std::size_t k = 0; k < numNeighbors; ++k)
                {
                    for (auto fluxVarIdx : assemblyMap_[globalI_][k])
                        neighborDeriv[k] -= localRes.evalFlux_(neighborElements[k], fluxVarIdx);
                }
            }
            else
            {
                // we are using forward differences, i.e. we don't need to
                // calculate f(x - \epsilon) and we can recycle the
                // (already calculated) residual f(x)
                partialDeriv -= this->residual_[0];
                neighborDeriv -= origFlux;
            }

            // divide difference in residuals by the magnitude of the
            // deflections between the two function evaluation
            partialDeriv /= delta;
            neighborDeriv /= delta;

            // restore the original state of the scv's volume variables
            this->model_().curVolVars_(scv) = origVolVars;

#if HAVE_VALGRIND
            for (unsigned i = 0; i < partialDeriv.size(); ++i)
                Valgrind::CheckDefined(partialDeriv[i]);
            for (int i = 0; i < neighborDeriv.size(); ++idx)
                for(unsigned k = 0; k < neighborDeriv[i].size(); ++k)
                    Valgrind::CheckDefined(neighborDeriv[i][k]);
#endif

            // update the global jacobian matrix with the current partial derivatives
            this->updateGlobalJacobian_(matrix, globalI_, globalI_, pvIdx, partialDeriv);

            j = 0;
            for (auto globalJ : neighborStencil)
                this->updateGlobalJacobian_(matrix, globalJ, globalI_, pvIdx, neighborDeriv[j++]);
        }
    }


    const FVElementGeometry& fvElemGeom_() const
    {
        return this->model_().fvGeometries(globalI_);
    }

    IndexType globalI_;
    ElementBoundaryTypes bcTypes_;

    int numericDifferenceMethod_;

    AssemblyMap assemblyMap_;
};
}

#endif
