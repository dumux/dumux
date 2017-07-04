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
#include "primaryvariables.hh"
#include "assemblymap.hh"

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
class StaggeredLocalJacobian
{
    using Implementation = typename GET_PROP_TYPE(TypeTag, LocalJacobian);
    using JacobianAssembler = typename GET_PROP_TYPE(TypeTag, JacobianAssembler);
    using LocalResidual = typename GET_PROP_TYPE(TypeTag, LocalResidual);
    using Problem = typename GET_PROP_TYPE(TypeTag, Problem);
    using Model = typename GET_PROP_TYPE(TypeTag, Model);
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
    using Indices = typename GET_PROP_TYPE(TypeTag, Indices);
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using Element = typename GridView::template Codim<0>::Entity;
    using IndexType = typename GridView::IndexSet::IndexType;

    enum {
        numEqCellCenter = GET_PROP_VALUE(TypeTag, NumEqCellCenter),
        numEqFace = GET_PROP_VALUE(TypeTag, NumEqFace),
    };

    using DofTypeIndices = typename GET_PROP(TypeTag, DofTypeIndices);
    typename DofTypeIndices::CellCenterIdx cellCenterIdx;
    typename DofTypeIndices::FaceIdx faceIdx;

    using CellCenterPrimaryVariables = typename GET_PROP_TYPE(TypeTag, CellCenterPrimaryVariables);
    using FacePrimaryVariables = typename GET_PROP_TYPE(TypeTag, FacePrimaryVariables);
    using GlobalFaceVars = typename GET_PROP_TYPE(TypeTag, GlobalFaceVars);

    using FaceSolutionVector = typename GET_PROP_TYPE(TypeTag, FaceSolutionVector);

    using PriVarIndices = typename Dumux::PriVarIndices<TypeTag>;

    using ElementSolutionVector = typename GET_PROP_TYPE(TypeTag, ElementSolutionVector);

    using AssemblyMap = Dumux::StaggeredAssemblyMap<TypeTag>;


public:

    StaggeredLocalJacobian(const StaggeredLocalJacobian &) = delete;

    StaggeredLocalJacobian()
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
        problemPtr_ = &problem;
        localResidual_.init(problem);
        assemblyMap_.init(problem);
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
        if(isGhost)
            DUNE_THROW(Dune::NotImplemented, "Support for ghost cells not implemented");

        // prepare the volvars/fvGeometries in case caching is disabled
        auto fvGeometry = localView(this->model_().globalFvGeometry());
        fvGeometry.bind(element);

        auto curElemVolVars = localView(this->model_().curGlobalVolVars());
        curElemVolVars.bind(element, fvGeometry, this->model_().curSol());

        auto prevElemVolVars = localView(this->model_().prevGlobalVolVars());
        prevElemVolVars.bindElement(element, fvGeometry, this->model_().prevSol());

        auto elemFluxVarsCache = localView(this->model_().globalFluxVarsCache());
        elemFluxVarsCache.bind(element, fvGeometry, curElemVolVars);

        // set the actual cell center dof index
        ccGlobalI_ = this->problem_().elementMapper().index(element);

        ElementBoundaryTypes elemBcTypes;
        elemBcTypes.update(this->problem_(), element, fvGeometry);

        auto& curGlobalFaceVars = this->model_().nonConstCurFaceVars();
        auto& prevGlobalFaceVars = this->model_().prevGlobalFaceVars();

        // calculate the local residual for all dofs of this element
        this->localResidual().eval(element, fvGeometry,
                                   prevElemVolVars, curElemVolVars,
                                   prevGlobalFaceVars, curGlobalFaceVars,
                                   elemBcTypes, elemFluxVarsCache);
        // store the cell center residual in global container
        residual[cellCenterIdx][ccGlobalI_] = this->localResidual().ccResidual();

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
        evalPartialDerivatives_(element,
                                fvGeometry,
                                prevElemVolVars,
                                curElemVolVars,
                                prevGlobalFaceVars,
                                curGlobalFaceVars,
                                elemFluxVarsCache,
                                elemBcTypes,
                                matrix,
                                residual[cellCenterIdx][ccGlobalI_],
                                faceResidualCache);
    }

     /*!
     * \brief Returns a reference to the object which calculates the
     *        local residual.
     */
    const LocalResidual &localResidual() const
    { return localResidual_; }

    /*!
     * \brief Returns a reference to the object which calculates the
     *        local residual.
     */
    LocalResidual &localResidual()
    { return localResidual_; }

    const AssemblyMap& assemblyMap() const
    { return assemblyMap_; }

protected:
    void evalPartialDerivatives_(const Element& element,
                                 const FVElementGeometry& fvGeometry,
                                 const ElementVolumeVariables& prevElemVolVars,
                                 ElementVolumeVariables& curElemVolVars,
                                 const GlobalFaceVars& prevGlobalFaceVars,
                                 GlobalFaceVars& curGlobalFaceVars,
                                 ElementFluxVariablesCache& elemFluxVarsCache,
                                 const ElementBoundaryTypes& elemBcTypes,
                                 JacobianMatrix& matrix,
                                 const CellCenterPrimaryVariables& ccResidual,
                                 const FaceSolutionVector& faceResidualCache)
    {
        // compute the derivatives of the cell center residuals with respect to cell center dofs
        dCCdCC_(element, fvGeometry, prevElemVolVars, curElemVolVars, prevGlobalFaceVars, curGlobalFaceVars, elemFluxVarsCache, elemBcTypes, matrix, ccResidual);

        // compute the derivatives of the cell center residuals with respect to face dofs
        dCCdFace_(element, fvGeometry, prevElemVolVars, curElemVolVars, prevGlobalFaceVars, curGlobalFaceVars, elemFluxVarsCache, elemBcTypes, matrix, ccResidual);

        // compute the derivatives of the face residuals with respect to cell center dofs
        dFacedCC_(element, fvGeometry, prevElemVolVars, curElemVolVars, prevGlobalFaceVars, curGlobalFaceVars, elemFluxVarsCache, elemBcTypes, matrix, faceResidualCache);

        // compute the derivatives of the face residuals with respect to face dofs
        dFacedFace_(element, fvGeometry, prevElemVolVars, curElemVolVars, prevGlobalFaceVars, curGlobalFaceVars, elemFluxVarsCache, elemBcTypes, matrix, faceResidualCache);
    }

     /*!
     * \brief Computes the derivatives of the cell center residuals with respect to cell center dofs
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
                 const CellCenterPrimaryVariables& ccResidual)
    {
        // build derivatives with for cell center dofs w.r.t. cell center dofs
        auto&& scvI = fvGeometry.scv(ccGlobalI_);

        for(const auto& globalJ : assemblyMap_(cellCenterIdx, cellCenterIdx, ccGlobalI_))
        {
            // get the volVars of the element with respect to which we are going to build the derivative
            auto&& scvJ = fvGeometry.scv(globalJ);
            const auto elementJ = fvGeometry.globalFvGeometry().element(globalJ);
            auto& curVolVars = getCurVolVars(curElemVolVars, scvJ);
            VolumeVariables origVolVars(curVolVars);

            for(auto pvIdx : PriVarIndices(cellCenterIdx))
            {
                PrimaryVariables priVars(CellCenterPrimaryVariables(this->model_().curSol()[cellCenterIdx][globalJ]),
                                         FacePrimaryVariables(0.0));

                const Scalar eps = numericEpsilon(priVars[pvIdx], cellCenterIdx, cellCenterIdx);
                priVars[pvIdx] += eps;
                ElementSolutionVector elemSol{std::move(priVars)};
                curVolVars.update(elemSol, this->problem_(), elementJ, scvJ);

                this->localResidual().evalCellCenter(element, fvGeometry, scvI,
                                        prevElemVolVars, curElemVolVars,
                                        prevGlobalFaceVars, curGlobalFaceVars,
                                        elemBcTypes, elemFluxVarsCache);

                auto partialDeriv = (this->localResidual().ccResidual() - ccResidual);
                partialDeriv /= eps;

                // update the global jacobian matrix with the current partial derivatives
                this->updateGlobalJacobian_(matrix[cellCenterIdx][cellCenterIdx], ccGlobalI_, globalJ, pvIdx, partialDeriv);

                // restore the original volVars
                curVolVars = origVolVars;
            }
        }
    }

     /*!
     * \brief Computes the derivatives of the cell center residuals with respect to face dofs
     */
    void dCCdFace_(const Element& element,
                   const FVElementGeometry& fvGeometry,
                   const ElementVolumeVariables& prevElemVolVars,
                   const ElementVolumeVariables& curElemVolVars,
                   const GlobalFaceVars& prevGlobalFaceVars,
                   GlobalFaceVars& curGlobalFaceVars,
                   ElementFluxVariablesCache& elemFluxVarsCache,
                   const ElementBoundaryTypes& elemBcTypes,
                   JacobianMatrix& matrix,
                   const CellCenterPrimaryVariables& ccResidual)
    {
        // build derivatives with for cell center dofs w.r.t. face dofs
        auto&& scvI = fvGeometry.scv(ccGlobalI_);

        for(const auto& globalJ : assemblyMap_(cellCenterIdx, faceIdx, ccGlobalI_))
        {
            // get the faceVars of the face with respect to which we are going to build the derivative
            auto origFaceVars = curGlobalFaceVars.faceVars(globalJ);
            auto& curFaceVars = curGlobalFaceVars.faceVars(globalJ);

            for(auto pvIdx : PriVarIndices(faceIdx))
            {
                PrimaryVariables priVars(CellCenterPrimaryVariables(0.0), FacePrimaryVariables(this->model_().curSol()[faceIdx][globalJ]));

                const Scalar eps = numericEpsilon(priVars[pvIdx], cellCenterIdx, faceIdx);
                priVars[pvIdx] += eps;
                curFaceVars.update(priVars[faceIdx]);

                this->localResidual().evalCellCenter(element, fvGeometry, scvI,
                                        prevElemVolVars, curElemVolVars,
                                        prevGlobalFaceVars, curGlobalFaceVars,
                                        elemBcTypes, elemFluxVarsCache);

                auto partialDeriv = (this->localResidual().ccResidual() - ccResidual);
                partialDeriv /= eps;

                // update the global jacobian matrix with the current partial derivatives
                this->updateGlobalJacobian_(matrix[cellCenterIdx][faceIdx], ccGlobalI_, globalJ, pvIdx - Indices::faceOffset, partialDeriv);

                // restore the original faceVars
                curFaceVars = origFaceVars;
            }
        }
    }

     /*!
     * \brief Computes the derivatives of the face residuals with respect to cell center dofs
     */
    void dFacedCC_(const Element& element,
                   const FVElementGeometry& fvGeometry,
                   const ElementVolumeVariables& prevElemVolVars,
                   ElementVolumeVariables& curElemVolVars,
                   const GlobalFaceVars& prevGlobalFaceVars,
                   const GlobalFaceVars& curGlobalFaceVars,
                   ElementFluxVariablesCache& elemFluxVarsCache,
                   const ElementBoundaryTypes& elemBcTypes,
                   JacobianMatrix& matrix,
                   const FaceSolutionVector& cachedResidual)
    {
        for(auto&& scvf : scvfs(fvGeometry))
        {
            // set the actual dof index
            const auto faceGlobalI = scvf.dofIndex();

            // build derivatives with for face dofs w.r.t. cell center dofs
            for(const auto& globalJ : assemblyMap_(faceIdx, cellCenterIdx, scvf.index()))
            {
                // get the volVars of the element with respect to which we are going to build the derivative
                auto&& scvJ = fvGeometry.scv(globalJ);
                const auto elementJ = fvGeometry.globalFvGeometry().element(globalJ);
                auto& curVolVars = getCurVolVars(curElemVolVars, scvJ);
                VolumeVariables origVolVars(curVolVars);

                for(auto pvIdx : PriVarIndices(cellCenterIdx))
                {
                    PrimaryVariables priVars(CellCenterPrimaryVariables(this->model_().curSol()[cellCenterIdx][globalJ]),
                                             FacePrimaryVariables(0.0));

                    const Scalar eps = numericEpsilon(priVars[pvIdx], faceIdx, cellCenterIdx);
                    priVars[pvIdx] += eps;
                    ElementSolutionVector elemSol{std::move(priVars)};
                    curVolVars.update(elemSol, this->problem_(), elementJ, scvJ);

                    this->localResidual().evalFace(element, fvGeometry, scvf,
                                            prevElemVolVars, curElemVolVars,
                                            prevGlobalFaceVars, curGlobalFaceVars,
                                            elemBcTypes, elemFluxVarsCache);

                    auto partialDeriv = (this->localResidual().faceResidual(scvf.localFaceIdx()) - cachedResidual[scvf.localFaceIdx()]);
                    partialDeriv /= eps;
                    // update the global jacobian matrix with the current partial derivatives
                    this->updateGlobalJacobian_(matrix[faceIdx][cellCenterIdx], faceGlobalI, globalJ, pvIdx, partialDeriv);

                    // restore the original volVars
                    curVolVars = origVolVars;
                }
            }
        }
    }

     /*!
     * \brief Computes the derivatives of the face residuals with respect to cell center dofs
     */
    void dFacedFace_(const Element& element,
                     const FVElementGeometry& fvGeometry,
                     const ElementVolumeVariables& prevElemVolVars,
                     const ElementVolumeVariables& curElemVolVars,
                     const GlobalFaceVars& prevGlobalFaceVars,
                     GlobalFaceVars& curGlobalFaceVars,
                     ElementFluxVariablesCache& elemFluxVarsCache,
                     const ElementBoundaryTypes& elemBcTypes,
                     JacobianMatrix& matrix,
                     const FaceSolutionVector& cachedResidual)
    {
        for(auto&& scvf : scvfs(fvGeometry))
        {
            // set the actual dof index
            const auto faceGlobalI = scvf.dofIndex();

            // build derivatives with for face dofs w.r.t. cell center dofs
            for(const auto& globalJ : assemblyMap_(faceIdx, faceIdx, scvf.index()))
            {
                // get the faceVars of the face with respect to which we are going to build the derivative
                auto origFaceVars = curGlobalFaceVars.faceVars(globalJ);
                auto& curFaceVars = curGlobalFaceVars.faceVars(globalJ);

                for(auto pvIdx : PriVarIndices(faceIdx))
                {
                    PrimaryVariables priVars(CellCenterPrimaryVariables(0.0), FacePrimaryVariables(this->model_().curSol()[faceIdx][globalJ]));

                    const Scalar eps = numericEpsilon(priVars[pvIdx], faceIdx, faceIdx);
                    priVars[pvIdx] += eps;
                    curFaceVars.update(priVars[faceIdx]);

                    this->localResidual().evalFace(element, fvGeometry, scvf,
                                            prevElemVolVars, curElemVolVars,
                                            prevGlobalFaceVars, curGlobalFaceVars,
                                            elemBcTypes, elemFluxVarsCache);

                    auto partialDeriv = (this->localResidual().faceResidual(scvf.localFaceIdx()) - cachedResidual[scvf.localFaceIdx()]);
                    partialDeriv /= eps;

                    // update the global jacobian matrix with the current partial derivatives
                    this->updateGlobalJacobian_(matrix[faceIdx][faceIdx], faceGlobalI, globalJ, pvIdx - Indices::faceOffset, partialDeriv);

                    // restore the original faceVars
                    curFaceVars = origFaceVars;
                }
            }
        }
    }

    /*!
     * \brief Returns a reference to the problem.
     */
    const Problem &problem_() const
    {
        Valgrind::CheckDefined(problemPtr_);
        return *problemPtr_;
    }

    /*!
     * \brief Returns a reference to the problem.
     */
    Problem &problem_()
    {
        Valgrind::CheckDefined(problemPtr_);
        return *problemPtr_;
    }

    /*!
     * \brief Returns a reference to the grid view.
     */
    const GridView &gridView_() const
    { return problem_().gridView(); }

    /*!
     * \brief Returns a reference to the model.
     */
    const Model &model_() const
    { return problem_().model(); }

    /*!
     * \brief Returns a reference to the model.
     */
    Model &model_()
    { return problem_().model(); }

    /*!
     * \brief Returns a reference to the jacobian assembler.
     */
    const JacobianAssembler &jacAsm_() const
    { return model_().jacobianAssembler(); }

     /*!
     * \brief Return the epsilon used to calculate the numeric derivative of a localResidual w.r.t a certain priVar
     *
     * \param priVar The the value of primary varible w.r.t which to derivative of the localResidual is calculated
     * \param idx1 Indicates whether the the derivative is build for a cellCenter or face localResidual
     * \param idx2 Indicates whether the the derivative is build w.r.t a priVar living on a cellCenter or face
     */
    Scalar numericEpsilon(const Scalar priVar, const int idx1, const int idx2) const
    {
        // define the base epsilon as the geometric mean of 1 and the
        // resolution of the scalar type. E.g. for standard 64 bit
        // floating point values, the resolution is about 10^-16 and
        // the base epsilon is thus approximately 10^-8.
        /*
        static const Scalar baseEps
            = Dumux::geometricMean<Scalar>(std::numeric_limits<Scalar>::epsilon(), 1.0);
        */

        static const Scalar baseEps = baseEps_[idx1][idx2];
        assert(std::numeric_limits<Scalar>::epsilon()*1e4 < baseEps);
        // the epsilon value used for the numeric differentiation is
        // now scaled by the absolute value of the primary variable...
        return baseEps*(std::abs(priVar) + 1.0);
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
            Valgrind::CheckDefined(matrix[globalI][globalJ][eqIdx][pvIdx]);
        }
    }

    //! If the global vol vars caching is enabled we have to modify the global volvar object
    template<class T = TypeTag>
    typename std::enable_if<GET_PROP_VALUE(T, EnableGlobalVolumeVariablesCache), VolumeVariables>::type&
    getCurVolVars(ElementVolumeVariables& elemVolVars, const SubControlVolume& scv)
    { return this->model_().nonConstCurGlobalVolVars().volVars(scv); }

    //! When global volume variables caching is disabled, return the local volvar object
    template<class T = TypeTag>
    typename std::enable_if<!GET_PROP_VALUE(T, EnableGlobalVolumeVariablesCache), VolumeVariables>::type&
    getCurVolVars(ElementVolumeVariables& elemVolVars, const SubControlVolume& scv)
    { return elemVolVars[scv]; }

    IndexType ccGlobalI_;
    int numericDifferenceMethod_;

    // The problem we would like to solve
    Problem *problemPtr_;

    LocalResidual localResidual_;

    AssemblyMap assemblyMap_;

    using BaseEpsilon = typename GET_PROP(TypeTag, BaseEpsilon);
    const std::array<std::array<Scalar, 2>, 2> baseEps_ = BaseEpsilon::getEps();
};

}

#endif
