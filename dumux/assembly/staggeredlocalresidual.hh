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
 * \ingroup StaggeredDiscretization
 * \ingroup Assembly
 * \brief Calculates the element-wise residual for the staggered FV scheme
 */
#ifndef DUMUX_STAGGERED_LOCAL_RESIDUAL_HH
#define DUMUX_STAGGERED_LOCAL_RESIDUAL_HH

#include <dumux/common/timeloop.hh>
#include <dumux/common/properties.hh>

#include "simpleassemblystructs.hh"

namespace Dumux {

/*!
 * \ingroup StaggeredDiscretization
 * \ingroup Assembly
 * \brief Calculates the element-wise residual for the staggered FV scheme
 */
template<class TypeTag>
class StaggeredLocalResidual
{
    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);

    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using Implementation = typename GET_PROP_TYPE(TypeTag, LocalResidual);
    using Problem = typename GET_PROP_TYPE(TypeTag, Problem);
    using Element = typename GridView::template Codim<0>::Entity;
    using ElementBoundaryTypes = typename GET_PROP_TYPE(TypeTag, ElementBoundaryTypes);
    using ElementVolumeVariables = typename GET_PROP_TYPE(TypeTag, GridVolumeVariables)::LocalView;
    using ElementFluxVariablesCache = typename GET_PROP_TYPE(TypeTag, GridFluxVariablesCache)::LocalView;
    using FVElementGeometry = typename GET_PROP_TYPE(TypeTag, FVGridGeometry)::LocalView;
    using SubControlVolume = typename FVElementGeometry::SubControlVolume;
    using SubControlVolumeFace = typename FVElementGeometry::SubControlVolumeFace;
    using CellCenterPrimaryVariables = typename GET_PROP_TYPE(TypeTag, CellCenterPrimaryVariables);
    using SolutionVector = typename GET_PROP_TYPE(TypeTag, SolutionVector);
    using Indices = typename GET_PROP_TYPE(TypeTag, ModelTraits)::Indices;

    using CellCenterResidual = typename GET_PROP_TYPE(TypeTag, CellCenterPrimaryVariables);
    using FaceResidual = typename GET_PROP_TYPE(TypeTag, FacePrimaryVariables);
    using ElementFaceVariables = typename GET_PROP_TYPE(TypeTag, GridFaceVariables)::LocalView;

    using TimeLoop = TimeLoopBase<Scalar>;

public:
    using CellCenterResidualValue = typename GET_PROP_TYPE(TypeTag, CellCenterPrimaryVariables);
    using FaceResidualValue = typename GET_PROP_TYPE(TypeTag, FacePrimaryVariables);
    using ElementResidualVector = CellCenterResidualValue;

    using SimpleMassBalanceSummands = typename GET_PROP_TYPE(TypeTag, SimpleMassBalanceSummands);
    using SimpleMomentumBalanceSummands = typename GET_PROP_TYPE(TypeTag, SimpleMomentumBalanceSummands);

    //! the constructor
    StaggeredLocalResidual(const Problem* problem,
                           const TimeLoop* timeLoop = nullptr)
    : problem_(problem)
    , timeLoop_(timeLoop)
    {}

    //! Convenience function to evaluate the flux and source terms for the cell center residual
    void evalFluxAndSourceForCellCenter(const Element& element,
                                                           const FVElementGeometry& fvGeometry,
                                                           const ElementVolumeVariables& elemVolVars,
                                                           const ElementFaceVariables& elemFaceVars,
                                                           const ElementBoundaryTypes& bcTypes,
                                                           const ElementFluxVariablesCache& elemFluxVarsCache,
                                                           SimpleMassBalanceSummands& simpleMassBalanceSummands) const
    {
        simpleMassBalanceSummands.setToZero(element, fvGeometry);

        // evaluate the source term
        for (auto&& scv : scvs(fvGeometry))
            asImp().evalSourceForCellCenter(this->problem(), element, fvGeometry, elemVolVars, elemFaceVars, scv, simpleMassBalanceSummands);

        // evaluate the flux term
        for (auto&& scvf : scvfs(fvGeometry))
            asImp().evalFluxForCellCenter(this->problem(), element, fvGeometry, elemVolVars, elemFaceVars, bcTypes, elemFluxVarsCache, scvf, simpleMassBalanceSummands);
    }

    //! Evaluate the flux terms for a cell center residual
    void evalFluxForCellCenter(const Problem& problem,
                               const Element& element,
                               const FVElementGeometry& fvGeometry,
                               const ElementVolumeVariables& elemVolVars,
                               const ElementFaceVariables& elemFaceVars,
                               const ElementBoundaryTypes& bcTypes,
                               const ElementFluxVariablesCache& elemFluxVarsCache,
                               const SubControlVolumeFace& scvf,
                               SimpleMassBalanceSummands& simpleMassBalanceSummands) const
    {
        if(!scvf.boundary())
            simpleMassBalanceSummands.coefficients[scvf.localFaceIdx()] += asImp_().computeFluxForCellCenter(problem, element, fvGeometry, elemVolVars, elemFaceVars, scvf, elemFluxVarsCache);
    }

    //! Evaluate the source terms for a cell center residual
    void evalSourceForCellCenter(const Problem& problem,
                                 const Element& element,
                                 const FVElementGeometry& fvGeometry,
                                 const ElementVolumeVariables& curElemVolVars,
                                 const ElementFaceVariables& curElemFaceVars,
                                 const SubControlVolume& scv,
                                 SimpleMassBalanceSummands& simpleMassBalanceSummands) const
    {
            const auto curExtrusionFactor = curElemVolVars[scv].extrusionFactor();

            // subtract the source term from the local rate
            auto source = asImp_().computeSourceForCellCenter(problem, element, fvGeometry, curElemVolVars, curElemFaceVars, scv);
            source *= scv.volume()*curExtrusionFactor;
            simpleMassBalanceSummands.RHS += source;
    }

    //! Evaluate the storage terms for a cell center residual
    void evalStorageForCellCenter(const Element &element,
                                                     const FVElementGeometry& fvGeometry,
                                                     const ElementVolumeVariables& prevElemVolVars,
                                                     const ElementVolumeVariables& curElemVolVars,
                                                     SimpleMassBalanceSummands& simpleMassBalanceSummands) const
    {
        assert(timeLoop_ && "no time loop set for storage term evaluation");
        CellCenterResidualValue storage(0.0);

        for (auto&& scv : scvs(fvGeometry))
            asImp().evalStorageForCellCenter(problem(), element, fvGeometry, prevElemVolVars, curElemVolVars, scv, simpleMassBalanceSummands);
    }

    //! Evaluate the storage terms for a cell center residual
    void evalStorageForCellCenter(const Problem& problem,
                                  const Element &element,
                                  const FVElementGeometry& fvGeometry,
                                  const ElementVolumeVariables& prevElemVolVars,
                                  const ElementVolumeVariables& curElemVolVars,
                                  const SubControlVolume& scv,
                                  SimpleMassBalanceSummands& simpleMassBalanceSummands) const
    {
        CellCenterResidualValue storage(0.0);
        const auto& curVolVars = curElemVolVars[scv];
        const auto& prevVolVars = prevElemVolVars[scv];

        // mass balance within the element. this is the
        // \f$\frac{m}{\partial t}\f$ term if using implicit
        // euler as time discretization.

        // We might need a more explicit way for
        // doing the time discretization...
        auto prevCCStorage = asImp_().computeStorageForCellCenter(problem, scv, prevVolVars);
        auto curCCStorage = asImp_().computeStorageForCellCenter(problem, scv, curVolVars);

        prevCCStorage *= prevVolVars.extrusionFactor();
        curCCStorage *= curVolVars.extrusionFactor();

        storage = std::move(curCCStorage);
        storage -= std::move(prevCCStorage);
        storage *= scv.volume();
        storage /= timeLoop_->timeStepSize();

        simpleMassBalanceSummands.RHS -= storage;
    }

    //! Evaluate the boundary conditions for a cell center residual
    void evalBoundaryForCellCenter(const Problem& problem,
                                   const Element& element,
                                   const FVElementGeometry& fvGeometry,
                                   const ElementVolumeVariables& elemVolVars,
                                   const ElementFaceVariables& elemFaceVars,
                                   const ElementBoundaryTypes& bcTypes,
                                   const ElementFluxVariablesCache& elemFluxVarsCache,
                                   SimpleMassBalanceSummands& simpleMassBalanceSummands) const
    {
        asImp_().evalBoundaryForCellCenter_(problem, element, fvGeometry, elemVolVars, elemFaceVars, bcTypes, elemFluxVarsCache, simpleMassBalanceSummands);
    }

    //! for compatibility with FVLocalAssemblerBase
    template<class... Args>
    CellCenterResidualValue evalFluxAndSource(Args&&... args) const
    {
        return CellCenterResidualValue(0.0);
    }

    //! for compatibility with FVLocalAssemblerBase
    template<class... Args>
    CellCenterResidualValue evalStorage(Args&&... args) const
    {
        return CellCenterResidualValue(0.0);
    }

    /*!
     * \name User interface
     * \note The following methods are usually expensive to evaluate
     *       They are useful for outputting residual information.
     */
    // \{

    //! Convenience function to evaluate the flux and source terms for the face residual
    void evalFluxAndSourceForFace(const Element& element,
                                               const FVElementGeometry& fvGeometry,
                                               const ElementVolumeVariables& elemVolVars,
                                               const ElementFaceVariables& elemFaceVars,
                                               const ElementBoundaryTypes& bcTypes,
                                               const ElementFluxVariablesCache& elemFluxVarsCache,
                                               const SubControlVolumeFace& scvf,
                                               SimpleMomentumBalanceSummands& simpleMomentumBalanceSummands) const
    {
        simpleMomentumBalanceSummands.setToZero(scvf);
        asImp().evalSourceForFace(this->problem(), element, fvGeometry, elemVolVars, elemFaceVars, scvf, simpleMomentumBalanceSummands);
        asImp().evalFluxForFace(this->problem(), element, fvGeometry, elemVolVars, elemFaceVars, bcTypes, elemFluxVarsCache, scvf, simpleMomentumBalanceSummands);
    }

    //! Evaluate the flux terms for a face residual
    void evalFluxForFace(const Problem& problem,
                         const Element& element,
                         const FVElementGeometry& fvGeometry,
                         const ElementVolumeVariables& elemVolVars,
                         const ElementFaceVariables& elemFaceVars,
                         const ElementBoundaryTypes& bcTypes,
                         const ElementFluxVariablesCache& elemFluxVarsCache,
                         const SubControlVolumeFace& scvf,
                         SimpleMomentumBalanceSummands& simpleMomentumBalanceSummands) const
    {
        if(!scvf.boundary())
            asImp_().computeFluxForFace(problem, element, scvf, fvGeometry, elemVolVars, elemFaceVars, elemFluxVarsCache, simpleMomentumBalanceSummands);
    }

    //! Evaluate the source terms for a face residual
    void evalSourceForFace(const Problem& problem,
                           const Element& element,
                           const FVElementGeometry& fvGeometry,
                           const ElementVolumeVariables& elemVolVars,
                           const ElementFaceVariables& elemFaceVars,
                           const SubControlVolumeFace& scvf,
                           SimpleMomentumBalanceSummands& simpleMomentumBalanceSummands) const
    {
        // the source term:
        auto source = asImp_().computeSourceForFace(problem, element, fvGeometry, scvf, elemVolVars, elemFaceVars);
        const auto& scv = fvGeometry.scv(scvf.insideScvIdx());
        const auto extrusionFactor = elemVolVars[scv].extrusionFactor();

        // multiply by 0.5 because we only consider half of a staggered control volume here
        source *= 0.5*scv.volume()*extrusionFactor;
        simpleMomentumBalanceSummands.RHS += source;
    }

    //! Evaluate the storage terms for a face residual
    FaceResidualValue evalStorageForFace(const Element& element,
                                         const FVElementGeometry& fvGeometry,
                                         const ElementVolumeVariables& prevElemVolVars,
                                         const ElementVolumeVariables& curElemVolVars,
                                         const ElementFaceVariables& prevElemFaceVars,
                                         const ElementFaceVariables& curElemFaceVars,
                                         const SubControlVolumeFace& scvf,
                                         SimpleMomentumBalanceSummands& simpleMomentumBalanceSummands) const
    {
        assert(timeLoop_ && "no time loop set for storage term evaluation");
        FaceResidualValue storage(0.0);
        asImp().evalStorageForFace(problem(), element, fvGeometry, prevElemVolVars, curElemVolVars, prevElemFaceVars, curElemFaceVars, scvf, simpleMomentumBalanceSummands);
        return storage;
    }

    //! Evaluate the storage terms for a face residual
    void evalStorageForFace(const Problem& problem,
                            const Element& element,
                            const FVElementGeometry& fvGeometry,
                            const ElementVolumeVariables& prevElemVolVars,
                            const ElementVolumeVariables& curElemVolVars,
                            const ElementFaceVariables& prevElemFaceVars,
                            const ElementFaceVariables& curElemFaceVars,
                            const SubControlVolumeFace& scvf,
                            SimpleMomentumBalanceSummands& simpleMomentumBalanceSummands) const
    {
        const auto& scv = fvGeometry.scv(scvf.insideScvIdx());
        auto prevFaceStorage = prevElemVolVars[scv].density() * prevElemFaceVars[scvf].velocitySelf();
        const auto extrusionFactor = curElemVolVars[scv].extrusionFactor();

        auto RHSContribution = std::move(prevFaceStorage);
        // multiply by 0.5 because we only consider half of a staggered control volume here
        RHSContribution *= (0.5*scv.volume()*extrusionFactor);
        RHSContribution /= timeLoop_->timeStepSize();

        simpleMomentumBalanceSummands.RHS += RHSContribution;

        //compute contribution to selfCoefficient
        if(scvf.boundary() && problem.boundaryTypes(element, scvf).isDirichlet(Indices::velocity(scvf.directionIndex()))){
            auto curFaceStorage = - curElemVolVars[scv].density() * curElemFaceVars[scvf].velocitySelf();

            auto secondRHSContribution = std::move(curFaceStorage);
            secondRHSContribution *= (0.5*scv.volume()*extrusionFactor);
            secondRHSContribution /= timeLoop_->timeStepSize();

            simpleMomentumBalanceSummands.RHS += secondRHSContribution;
        }
        else {
            auto curFaceStorage = curElemVolVars[scv].density();

            auto selfCoefficientContribution = std::move(curFaceStorage);
            selfCoefficientContribution *= (0.5*scv.volume()*extrusionFactor);
            selfCoefficientContribution /= timeLoop_->timeStepSize();

            simpleMomentumBalanceSummands.selfCoefficient += selfCoefficientContribution;
        }
    }

    //! Evaluate the boundary conditions for a face residual
    void evalBoundaryForFace(const Problem& problem,
                             const Element& element,
                             const FVElementGeometry& fvGeometry,
                             const ElementVolumeVariables& elemVolVars,
                             const ElementFaceVariables& elemFaceVars,
                             const ElementBoundaryTypes& bcTypes,
                             const ElementFluxVariablesCache& elemFluxVarsCache,
                             const SubControlVolumeFace& scvf,
                             SimpleMomentumBalanceSummands& simpleMomentumBalanceSummands) const
    {
        asImp_().evalBoundaryForFace_(problem, element, fvGeometry, scvf, elemVolVars, elemFaceVars, bcTypes, elemFluxVarsCache, simpleMomentumBalanceSummands);
    }

    //! If no solution has been set, we treat the problem as stationary.
    bool isStationary() const
    { return !timeLoop_; }

    //! the problem
    const Problem& problem() const
    { return *problem_; }

protected:

    Implementation &asImp_()
    { return *static_cast<Implementation*>(this); }

    const Implementation &asImp_() const
    { return *static_cast<const Implementation*>(this); }


    TimeLoop& timeLoop()
    { return *timeLoop_; }

    const TimeLoop& timeLoop() const
    { return *timeLoop_; }

    Implementation &asImp()
    { return *static_cast<Implementation*>(this); }

    const Implementation &asImp() const
    { return *static_cast<const Implementation*>(this); }

private:
    const Problem* problem_; //!< the problem we are assembling this residual for
    const TimeLoop* timeLoop_;

};

}

#endif   // DUMUX_CC_LOCAL_RESIDUAL_HH
