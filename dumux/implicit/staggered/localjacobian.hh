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
#ifndef DUMUX_STAGGERED_LOCAL_JACOBIAN_HH
#define DUMUX_STAGGERED_LOCAL_JACOBIAN_HH

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
class StaggeredLocalJacobian : public ImplicitLocalJacobian<TypeTag>
{
    using ParentType = ImplicitLocalJacobian<TypeTag>;
    using Implementation = typename GET_PROP_TYPE(TypeTag, LocalJacobian);
    using LocalResidual = typename GET_PROP_TYPE(TypeTag, LocalResidual);
    using Problem = typename GET_PROP_TYPE(TypeTag, Problem);
    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
    using JacobianMatrix = typename GET_PROP_TYPE(TypeTag, JacobianMatrix);

    using FVElementGeometry = typename GET_PROP_TYPE(TypeTag, FVElementGeometry);
    using SubControlVolume = typename GET_PROP_TYPE(TypeTag, SubControlVolume);
    using PrimaryVariables = typename GET_PROP_TYPE(TypeTag, PrimaryVariables);
    using FluxVariables = typename GET_PROP_TYPE(TypeTag, FluxVariables);
    using VolumeVariables = typename GET_PROP_TYPE(TypeTag, VolumeVariables);
    using ElementVolumeVariables = typename GET_PROP_TYPE(TypeTag, ElementVolumeVariables);
    using ElementBoundaryTypes = typename GET_PROP_TYPE(TypeTag, ElementBoundaryTypes);
    using ElementFluxVariablesCache = typename GET_PROP_TYPE(TypeTag, ElementFluxVariablesCache);
    using SolutionVector = typename GET_PROP_TYPE(TypeTag, SolutionVector);
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using Element = typename GridView::template Codim<0>::Entity;
    using IndexType = typename GridView::IndexSet::IndexType;

    enum { numEq = GET_PROP_VALUE(TypeTag, NumEq) };

    using AssemblyMap = std::vector<std::vector<std::vector<IndexType>>>;

public:

    StaggeredLocalJacobian() : ParentType()
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
    }

    /*!
     * \brief Assemble an element's local Jacobian matrix of the
     *        defect.
     *
     * \param element The DUNE Codim<0> entity which we look at.
     */
    void assemble(const Element& element, JacobianMatrix& matrix, SolutionVector& residual)
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

        ElementBoundaryTypes elemBcTypes;
        elemBcTypes.update(this->problem_(), element, fvGeometry);

        // calculate the local residual
        if (isGhost)
        {
            this->residual_ = 0.0;
            residual[globalI_] = 0.0;
        }
        else
        {
            this->localResidual().eval(element, fvGeometry, prevElemVolVars, curElemVolVars, elemBcTypes, elemFluxVarsCache);
            this->residual_ = this->localResidual().residual();
            // store residual in global container as well
            residual[globalI_] = this->localResidual().residual(0);
        }

        this->model_().updatePVWeights(fvGeometry);

        // calculate derivatives of all dofs in stencil with respect to the dofs in the element
        evalPartialDerivatives_(element, fvGeometry, prevElemVolVars, curElemVolVars, elemFluxVarsCache, elemBcTypes, matrix, residual, isGhost);

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

    const AssemblyMap& assemblyMap() const
    { return assemblyMap_; }

private:
    void evalPartialDerivatives_(const Element& element,
                                 const FVElementGeometry& fvGeometry,
                                 const ElementVolumeVariables& prevElemVolVars,
                                 ElementVolumeVariables& curElemVolVars,
                                 ElementFluxVariablesCache& elemFluxVarsCache,
                                 const ElementBoundaryTypes& elemBcTypes,
                                 JacobianMatrix& matrix,
                                 SolutionVector& residual,
                                 const bool isGhost)
    {
        // get stencil informations
//         const auto& neighborStencil = this->model_().stencils(element).neighborStencil();
//         const auto numNeighbors = neighborStencil.size();

        // the localresidual class used for the flux calculations
//         LocalResidual localRes;
//         localRes.init(this->problem_());
//
//         // container to store the neighboring elements
//         std::vector<Element> neighborElements;
// //         neighborElements.reserve(numNeighbors);
//
//         // get the elements and calculate the flux into the element in the undeflected state
//         Dune::BlockVector<PrimaryVariables> origFlux(numNeighbors);
//         origFlux = 0.0;
//
//         unsigned int j = 0;
//         for (auto globalJ : neighborStencil)
//         {
//             auto elementJ = fvGeometry.globalFvGeometry().element(globalJ);
//             neighborElements.push_back(elementJ);
//
//             for (auto fluxVarIdx : assemblyMap_[globalI_][j])
//             {
//                 auto&& scvf = fvGeometry.scvf(fluxVarIdx);
//                 origFlux[j] += localRes.evalFlux_(elementJ, fvGeometry, curElemVolVars, scvf, elemFluxVarsCache[scvf]);
//             }
//
//             ++j;
//         }
//
//         auto&& scv = fvGeometry.scv(globalI_);
//         auto& curVolVars = getCurVolVars(curElemVolVars, scv);
//         // save a copy of the original vol vars
//         VolumeVariables origVolVars(curVolVars);
//
//         // derivatives in the neighbors with repect to the current elements
//         Dune::BlockVector<PrimaryVariables> neighborDeriv(numNeighbors);
//         for (int pvIdx = 0; pvIdx < numEq; pvIdx++)
//         {
//             // derivatives of element dof with respect to itself
//             PrimaryVariables partialDeriv(0.0);
//
//             if (isGhost)
//                 partialDeriv[pvIdx] = 1.0;
//
//             neighborDeriv = 0.0;
//             PrimaryVariables priVars(this->model_().curSol()[globalI_]);
//
//             Scalar eps = this->numericEpsilon(scv, curVolVars, pvIdx);
//             Scalar delta = 0;
//
//             if (numericDifferenceMethod_ >= 0)
//             {
//                 // we are not using backward differences, i.e. we need to
//                 // calculate f(x + \epsilon)
//
//                 // deflect primary variables
//                 priVars[pvIdx] += eps;
//                 delta += eps;
//
//                 // update the volume variables and bind the flux var cache again
//                 curVolVars.update(priVars, this->problem_(), element, scv);
//                 elemFluxVarsCache.bind(element, fvGeometry, curElemVolVars);
//
//                 if (!isGhost)
//                 {
//                     // calculate the residual with the deflected primary variables
//                     this->localResidual().eval(element, fvGeometry, prevElemVolVars, curElemVolVars, elemBcTypes, elemFluxVarsCache);
//
//                     // store the residual and the storage term
//                     partialDeriv = this->localResidual().residual(0);
//                 }
//
//                 // calculate the fluxes into element with the deflected primary variables
//                 for (std::size_t k = 0; k < numNeighbors; ++k)
//                 {
//                     for (auto fluxVarIdx : assemblyMap_[globalI_][k])
//                     {
//                         auto&& scvf = fvGeometry.scvf(fluxVarIdx);
//                         neighborDeriv[k] += localRes.evalFlux_(neighborElements[k], fvGeometry, curElemVolVars, scvf, elemFluxVarsCache[scvf]);
//                     }
//                 }
//             }
//             else
//             {
//                 // we are using backward differences, i.e. we don't need
//                 // to calculate f(x + \epsilon) and we can recycle the
//                 // (already calculated) residual f(x)
//                 if (!isGhost)
//                     partialDeriv = this->residual(0);
//                 neighborDeriv = origFlux;
//             }
//
//             if (numericDifferenceMethod_ <= 0)
//             {
//                 // we are not using forward differences, i.e. we
//                 // need to calculate f(x - \epsilon)
//
//                 // deflect the primary variables
//                 priVars[pvIdx] -= delta + eps;
//                 delta += eps;
//
//                 // update the volume variables and bind the flux var cache again
//                 curVolVars.update(priVars, this->problem_(), element, scv);
//                 elemFluxVarsCache.bind(element, fvGeometry, curElemVolVars);
//
//                 if (!isGhost)
//                 {
//                     // calculate the residual with the deflected primary variables
//                     this->localResidual().eval(element, fvGeometry, prevElemVolVars, curElemVolVars, elemBcTypes, elemFluxVarsCache);
//
//                     // subtract the residual from the derivative storage
//                     partialDeriv -= this->localResidual().residual(0);
//                 }
//
//                 // calculate the fluxes into element with the deflected primary variables
//                 for (std::size_t k = 0; k < numNeighbors; ++k)
//                 {
//                     for (auto fluxVarIdx : assemblyMap_[globalI_][k])
//                     {
//                         auto&& scvf = fvGeometry.scvf(fluxVarIdx);
//                         neighborDeriv[k] -= localRes.evalFlux_(neighborElements[k], fvGeometry, curElemVolVars, scvf, elemFluxVarsCache[scvf]);
//                     }
//                 }
//             }
//             else
//             {
//                 // we are using forward differences, i.e. we don't need to
//                 // calculate f(x - \epsilon) and we can recycle the
//                 // (already calculated) residual f(x)
//                 if (!isGhost)
//                     partialDeriv -= this->residual(0);
//                 neighborDeriv -= origFlux;
//             }
//
//             // divide difference in residuals by the magnitude of the
//             // deflections between the two function evaluation
//             if (!isGhost)
//                 partialDeriv /= delta;
//             neighborDeriv /= delta;
//
//             // restore the original state of the scv's volume variables
//             curVolVars = origVolVars;
//
//             // update the global jacobian matrix with the current partial derivatives
//             this->updateGlobalJacobian_(matrix, globalI_, globalI_, pvIdx, partialDeriv);
//
//             j = 0;
//             for (auto globalJ : neighborStencil)
//                 this->updateGlobalJacobian_(matrix, globalJ, globalI_, pvIdx, neighborDeriv[j++]);
//         }
    }

    //! If the global vol vars caching is enabled we have to modify the global volvar object
    template<class T = TypeTag>
    typename std::enable_if<GET_PROP_VALUE(T, EnableGlobalVolumeVariablesCache), VolumeVariables>::type&
    getCurVolVars(ElementVolumeVariables& elemVolVars, const SubControlVolume& scv)
    { return this->model_().nonConstCurGlobalVolVars().volVars(scv.index()); }

    //! When global volume variables caching is disabled, return the local volvar object
    template<class T = TypeTag>
    typename std::enable_if<!GET_PROP_VALUE(T, EnableGlobalVolumeVariablesCache), VolumeVariables>::type&
    getCurVolVars(ElementVolumeVariables& elemVolVars, const SubControlVolume& scv)
    { return elemVolVars[scv]; }

    IndexType globalI_;
    int numericDifferenceMethod_;
    AssemblyMap assemblyMap_;
};

}

#endif
