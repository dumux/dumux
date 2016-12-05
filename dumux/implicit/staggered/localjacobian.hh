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

    using DofTypeIndices = typename GET_PROP(TypeTag, DofTypeIndices);
    typename DofTypeIndices::CellCenterIdx cellCenterIdx;
    typename DofTypeIndices::FaceIdx faceIdx;

    using CellCenterPrimaryVariables = typename GET_PROP_TYPE(TypeTag, CellCenterPrimaryVariables);
    using FacePrimaryVariables = typename GET_PROP_TYPE(TypeTag, FacePrimaryVariables);
    using GlobalFaceVars = typename GET_PROP_TYPE(TypeTag, GlobalFaceVars);


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

        auto curGlobalFaceVars = this->model_().curGlobalFaceVars();
        auto& prevGlobalFaceVars = this->model_().prevGlobalFaceVars();

        // calculate the local residual
        if (isGhost)
        {
            DUNE_THROW(Dune::NotImplemented, "Support for ghost cells not implemented");
        }
        else
        {
            this->localResidual().eval(element, fvGeometry,
                                       prevElemVolVars, curElemVolVars,
                                       prevGlobalFaceVars, curGlobalFaceVars,
                                       elemBcTypes, elemFluxVarsCache);
            // store residual in global container as well
            residual[cellCenterIdx][globalI_] = this->localResidual().ccResidual();

            for(auto&& scvf : scvfs(fvGeometry))
            {
                residual[faceIdx][scvf.dofIndexSelf()] += this->localResidual().faceResidual(scvf.localFaceIdx());
            }
        }

        this->model_().updatePVWeights(fvGeometry);

        // calculate derivatives of all dofs in stencil with respect to the dofs in the element
        evalPartialDerivatives_(element, fvGeometry, prevElemVolVars, curElemVolVars, prevGlobalFaceVars, curGlobalFaceVars, elemFluxVarsCache, elemBcTypes, matrix, residual, isGhost);
    }

private:
    void evalPartialDerivatives_(const Element& element,
                                 const FVElementGeometry& fvGeometry,
                                 const ElementVolumeVariables& prevElemVolVars,
                                 ElementVolumeVariables& curElemVolVars,
                                 const GlobalFaceVars& prevGlobalFaceVars,
                                 GlobalFaceVars& curGlobalFaceVars,
                                 ElementFluxVariablesCache& elemFluxVarsCache,
                                 const ElementBoundaryTypes& elemBcTypes,
                                 JacobianMatrix& matrix,
                                 SolutionVector& residual,
                                 const bool isGhost)
    {
        // compute the derivatives of the cell center dofs with respect to cell center dofs
        dCCdCC_(element, fvGeometry, prevElemVolVars, curElemVolVars, prevGlobalFaceVars, curGlobalFaceVars, elemFluxVarsCache, elemBcTypes, matrix, residual, isGhost);


        printmatrix(std::cout, matrix[cellCenterIdx][cellCenterIdx], "A11 neu", "");
    }

     /*!
     * \brief Computes the derivatives of the cell center dofs with respect to cell center dofs
     */
    void dCCdCC_(const Element& element,
                 const FVElementGeometry& fvGeometry,
                 const ElementVolumeVariables& prevElemVolVars,
                 ElementVolumeVariables& curElemVolVars,
                 const GlobalFaceVars& prevGlobalFaceVars,
                 const GlobalFaceVars& curGlobalFaceVars,
                 ElementFluxVariablesCache& elemFluxVarsCache,
                 const ElementBoundaryTypes& elemBcTypes,
                 JacobianMatrix& matrix,
                 SolutionVector& residual,
                 const bool isGhost)
    {
        // set the actual dof index
        auto globalI = this->problem_().elementMapper().index(element);


        // build derivatives with for cell center dofs w.r.t. cell center dofs
        const auto& cellCenterToCellCenterStencil = this->model_().stencils(element).cellCenterToCellCenterStencil();

        for(const auto& globalJ : cellCenterToCellCenterStencil)
        {
            // get the volVars of the element with respect to which we are going to build the derivative
            auto&& scv = fvGeometry.scv(globalJ);
            const auto elementJ = fvGeometry.globalFvGeometry().element(globalJ);
            curElemVolVars.bind(elementJ, fvGeometry, this->model_().curSol());
            auto& curVolVars = getCurVolVars(curElemVolVars, scv);
            VolumeVariables origVolVars(curVolVars);

            CellCenterPrimaryVariables priVars(this->model_().curSol()[cellCenterIdx][globalJ]);


            for(int pvIdx = 0; pvIdx < priVars.size(); ++pvIdx)
            {
                const Scalar eps = 1e-4; // TODO: do properly
                priVars += eps;

                curVolVars.update(priVars, this->problem_(), elementJ, scv);

                this->localResidual().eval(element, fvGeometry,
                                        prevElemVolVars, curElemVolVars,
                                        prevGlobalFaceVars, curGlobalFaceVars,
                                        elemBcTypes, elemFluxVarsCache);

                auto partialDeriv = (this->localResidual().ccResidual() - residual[cellCenterIdx][globalI]);
                partialDeriv /= eps;

                // update the global jacobian matrix with the current partial derivatives
                this->updateGlobalJacobian_(matrix[cellCenterIdx][cellCenterIdx], globalI, globalJ, pvIdx, partialDeriv);

                // restore the original volVars
                curVolVars = origVolVars;
            }
        }
    }

    /*!
     * \brief Updates the current global Jacobian matrix with the
     *        partial derivatives of all equations in regard to the
     *        primary variable 'pvIdx' at dof 'col'. Specialization for cc methods.
     */
    template<class SubMatrix, class CCOrFacePrimaryVariables>
    void updateGlobalJacobian_(SubMatrix& matrix,
                          const int globalI,
                          const int globalJ,
                          const int pvIdx,
                          const CCOrFacePrimaryVariables &partialDeriv)
    {
        for (int eqIdx = 0; eqIdx < numEq; eqIdx++)
        {
            // A[i][col][eqIdx][pvIdx] is the rate of change of
            // the residual of equation 'eqIdx' at dof 'i'
            // depending on the primary variable 'pvIdx' at dof
            // 'col'.
            matrix[globalI][globalJ][eqIdx][pvIdx] = partialDeriv[eqIdx];
            Valgrind::CheckDefined(matrix[globalI][globalJ][eqIdx][pvIdx]);
        }
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
};

}

#endif
