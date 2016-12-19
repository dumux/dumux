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
 * \brief Calculates the residual of models based on the box scheme element-wise.
 */
#ifndef DUMUX_STAGGERED_NAVIERSTOKES_LOCAL_RESIDUAL_HH
#define DUMUX_STAGGERED_NAVIERSTOKES_LOCAL_RESIDUAL_HH

#include <dune/istl/matrix.hh>

#include <dumux/common/valgrind.hh>
#include <dumux/implicit/staggered/localresidual.hh>

#include "properties.hh"

namespace Dumux
{
/*!
 * \ingroup CCModel
 * \ingroup StaggeredLocalResidual
 * \brief Element-wise calculation of the residual for models
 *        based on the fully implicit cell-centered scheme.
 *
 * \todo Please doc me more!
 */
template<class TypeTag>
class StaggeredNavierStokesResidual : public Dumux::StaggeredLocalResidual<TypeTag>
{
    using ParentType = StaggeredLocalResidual<TypeTag>;
    friend class StaggeredLocalResidual<TypeTag>;
    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);

    enum { numEq = GET_PROP_VALUE(TypeTag, NumEq) };

    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using Implementation = typename GET_PROP_TYPE(TypeTag, LocalResidual);
    using Problem = typename GET_PROP_TYPE(TypeTag, Problem);
    using Element = typename GridView::template Codim<0>::Entity;
    using PrimaryVariables = typename GET_PROP_TYPE(TypeTag, PrimaryVariables);
    using BoundaryTypes = typename GET_PROP_TYPE(TypeTag, BoundaryTypes);
    using ElementBoundaryTypes = typename GET_PROP_TYPE(TypeTag, ElementBoundaryTypes);
    using ElementVolumeVariables = typename GET_PROP_TYPE(TypeTag, ElementVolumeVariables);
    using ElementFluxVariablesCache = typename GET_PROP_TYPE(TypeTag, ElementFluxVariablesCache);
    using FluxVariablesCache = typename GET_PROP_TYPE(TypeTag, FluxVariablesCache);
    using FVElementGeometry = typename GET_PROP_TYPE(TypeTag, FVElementGeometry);
    using SubControlVolume = typename GET_PROP_TYPE(TypeTag, SubControlVolume);
    using SubControlVolumeFace = typename GET_PROP_TYPE(TypeTag, SubControlVolumeFace);
    using CellCenterSolutionVector = typename GET_PROP_TYPE(TypeTag, CellCenterSolutionVector);
    using FaceSolutionVector = typename GET_PROP_TYPE(TypeTag, FaceSolutionVector);
    using CellCenterPrimaryVariables = typename GET_PROP_TYPE(TypeTag, CellCenterPrimaryVariables);
    using FacePrimaryVariables = typename GET_PROP_TYPE(TypeTag, FacePrimaryVariables);
    using Indices = typename GET_PROP_TYPE(TypeTag, Indices);


    using DofTypeIndices = typename GET_PROP(TypeTag, DofTypeIndices);
    typename DofTypeIndices::CellCenterIdx cellCenterIdx;
    typename DofTypeIndices::FaceIdx faceIdx;

    enum {
         // grid and world dimension
        dim = GridView::dimension,
        dimWorld = GridView::dimensionworld,

        pressureIdx = Indices::pressureIdx,

        massBalanceIdx = Indices::massBalanceIdx,
        momentumBalanceIdx = Indices::momentumBalanceIdx
    };

    using VolumeVariables = typename GET_PROP_TYPE(TypeTag, VolumeVariables);
    using GlobalFaceVars = typename GET_PROP_TYPE(TypeTag, GlobalFaceVars);

    static constexpr bool navierStokes = GET_PROP_VALUE(TypeTag, EnableInertiaTerms);

public:
    // copying the local residual class is not a good idea
    StaggeredNavierStokesResidual(const StaggeredNavierStokesResidual &) = delete;

    StaggeredNavierStokesResidual() = default;


    CellCenterPrimaryVariables computeFluxForCellCenter(const Element &element,
                                  const FVElementGeometry& fvGeometry,
                                  const ElementVolumeVariables& elemVolVars,
                                  const GlobalFaceVars& globalFaceVars,
                                  const SubControlVolumeFace &scvf,
                                  const FluxVariablesCache& fluxVarsCache)
    {
        const auto& insideScv = fvGeometry.scv(scvf.insideScvIdx());
        const auto& insideVolVars = elemVolVars[insideScv];

        const Scalar velocity = globalFaceVars.faceVars(scvf.dofIndexSelf()).velocity();

        // if we are on an inflow/outflow boundary, use the volVars of the element itself
        const auto& outsideVolVars = scvf.boundary() ?  insideVolVars : elemVolVars[scvf.outsideScvIdx()];

        CellCenterPrimaryVariables flux(0.0);

        if(sign(scvf.outerNormalScalar()) == sign(velocity))
            flux[0] = insideVolVars.density() * velocity;
        else
            flux[0] = outsideVolVars.density() * velocity;
        return flux * scvf.area() * sign(scvf.outerNormalScalar());
    }

    CellCenterPrimaryVariables computeSourceForCellCenter(const Element &element,
                                     const FVElementGeometry& fvGeometry,
                                     const ElementVolumeVariables& elemVolVars,
                                     const GlobalFaceVars& globalFaceVars,
                                     const SubControlVolume &scv)
    {
        return CellCenterPrimaryVariables(0.0);
    }


     /*!
     * \brief Evaluate the rate of change of all conservation
     *        quantites (e.g. phase mass) within a sub-control
     *        volume of a finite volume element for the immiscible models.
     * \param scv The sub control volume
     * \param volVars The current or previous volVars
     * \note This function should not include the source and sink terms.
     * \note The volVars can be different to allow computing
     *       the implicit euler time derivative here
     */
    CellCenterPrimaryVariables computeStorageForCellCenter(const SubControlVolume& scv,
                                    const VolumeVariables& volVars)
    {
        CellCenterPrimaryVariables storage;
        storage[0] = volVars.density(0);
        return storage;
    }

     /*!
     * \brief Evaluate the rate of change of all conservation
     *        quantites (e.g. phase mass) within a sub-control
     *        volume of a finite volume element for the immiscible models.
     * \param scvf The sub control volume
     * \param volVars The current or previous volVars
     * \note This function should not include the source and sink terms.
     * \note The volVars can be different to allow computing
     *       the implicit euler time derivative here
     */
    FacePrimaryVariables computeStorageForFace(const SubControlVolumeFace& scvf,
                                               const VolumeVariables& volVars,
                                               const GlobalFaceVars& globalFaceVars)
    {
        FacePrimaryVariables storage(0.0);
        const Scalar velocity = globalFaceVars.faceVars(scvf.dofIndexSelf()).velocity();
        storage[0] = volVars.density(0) * velocity;
        return storage;
    }

    FacePrimaryVariables computeSourceForFace(const SubControlVolumeFace& scvf,
                                              const ElementVolumeVariables& elemVolVars,
                                              const GlobalFaceVars& globalFaceVars)
    {
        FacePrimaryVariables gravityTerm(0.0);
        const auto insideScvIdx = scvf.insideScvIdx();
        const auto& insideVolVars = elemVolVars[insideScvIdx];
        gravityTerm += this->problem().gravity()[scvf.directionIndex()] * insideVolVars.density();
        return gravityTerm;
    }

     /*!
     * \brief Returns the complete momentum flux for a face
     * \param scvf The sub control volume face
     * \param fvGeometry The finite-volume geometry
     * \param elemVolVars All volume variables for the element
     * \param globalFaceVars The face variables
     */
    FacePrimaryVariables computeFluxForFace(const SubControlVolumeFace& scvf,
                                            const FVElementGeometry& fvGeometry,
                                            const ElementVolumeVariables& elemVolVars,
                                            const GlobalFaceVars& globalFaceVars)
    {
        FacePrimaryVariables flux(0.0);
        flux += computeNormalMomentumFlux_(scvf, fvGeometry, elemVolVars, globalFaceVars);
        flux += computeTangetialMomentumFlux_(scvf, fvGeometry, elemVolVars, globalFaceVars);
        flux += computePressureTerm_(scvf, fvGeometry, elemVolVars, globalFaceVars);
        return flux;
    }

protected:

     /*!
     * \brief Evaluate boundary conditions
     */
    void evalBoundary_(const Element& element,
                       const FVElementGeometry& fvGeometry,
                       const ElementVolumeVariables& elemVolVars,
                       const GlobalFaceVars& faceVars,
                       const ElementBoundaryTypes& elemBcTypes,
                       const ElementFluxVariablesCache& elemFluxVarsCache)
    {
        evalBoundaryForCellCenter_(element, fvGeometry, elemVolVars, faceVars, elemBcTypes, elemFluxVarsCache);
        for (auto&& scvf : scvfs(fvGeometry))
        {
            evalBoundaryForFace_(element, fvGeometry, scvf, elemVolVars, faceVars, elemBcTypes, elemFluxVarsCache);
        }
    }

     /*!
     * \brief Evaluate boundary conditions for a cell center dof
     */
    void evalBoundaryForCellCenter_(const Element& element,
                                    const FVElementGeometry& fvGeometry,
                                    const ElementVolumeVariables& elemVolVars,
                                    const GlobalFaceVars& faceVars,
                                    const ElementBoundaryTypes& elemBcTypes,
                                    const ElementFluxVariablesCache& elemFluxVarsCache)
    {
        for (auto&& scvf : scvfs(fvGeometry))
        {
            if (scvf.boundary())
            {
                // For the mass-balance residual, do the same as if the face was not on a boundary.This might need to be changed sometime...
                this->ccResidual_ += computeFluxForCellCenter(element, fvGeometry, elemVolVars, faceVars, scvf, elemFluxVarsCache[scvf]);

                // handle the actual boundary conditions:
                const auto bcTypes = this->problem().boundaryTypes(element, scvf);

                // set a fixed pressure for cells adjacent to a wall
                if(bcTypes.isDirichlet(massBalanceIdx) && bcTypes.isDirichlet(momentumBalanceIdx))
                {
                    const auto& insideScv = fvGeometry.scv(scvf.insideScvIdx());
                    const auto& insideVolVars = elemVolVars[insideScv];
                    this->ccResidual_[pressureIdx] = insideVolVars.pressure() - this->problem().dirichletAtPos(scvf.center())[pressureIdx];
                }
            }
        }
    }

     /*!
     * \brief Evaluate boundary conditions for a face dof
     */
    void evalBoundaryForFace_(const Element& element,
                              const FVElementGeometry& fvGeometry,
                              const SubControlVolumeFace& scvf,
                              const ElementVolumeVariables& elemVolVars,
                              const GlobalFaceVars& faceVars,
                              const ElementBoundaryTypes& elemBcTypes,
                              const ElementFluxVariablesCache& elemFluxVarsCache)
    {
        if (scvf.boundary())
        {
            // handle the actual boundary conditions:
            const auto bcTypes = this->problem().boundaryTypes(element, scvf);

            // set a fixed value for the velocity
            if(bcTypes.isDirichlet(momentumBalanceIdx))
            {
                const Scalar velocity = faceVars.faceVars(scvf.dofIndexSelf()).velocity();
                const Scalar dirichletValue = this->problem().faceDirichletAtPos(scvf.center(), scvf.directionIndex());
                this->faceResiduals_[scvf.localFaceIdx()] = velocity - dirichletValue;
            }

            // outflow condition for the momentum balance equation
            if(bcTypes.isOutflow(momentumBalanceIdx))
            {
                if(bcTypes.isDirichlet(massBalanceIdx))
                    this->faceResiduals_[scvf.localFaceIdx()] += computeFluxForFace(scvf, fvGeometry, elemVolVars, faceVars);
                else
                    DUNE_THROW(Dune::InvalidStateException, "Face at " << scvf.center()  << " has an outflow BC for the momentum balance but no Dirichlet BC for the pressure!");
            }
        }
    }


private:
     /*!
     * \brief Returns the normal part of the momentum flux
     * \param scvf The sub control volume face
     * \param fvGeometry The finite-volume geometry
     * \param elemVolVars All volume variables for the element
     * \param globalFaceVars The face variables
     */
    FacePrimaryVariables computeNormalMomentumFlux_(const SubControlVolumeFace& scvf,
                                          const FVElementGeometry& fvGeometry,
                                          const ElementVolumeVariables& elemVolVars,
                                          const GlobalFaceVars& globalFaceVars)
    {
        const auto insideScvIdx = scvf.insideScvIdx();
        const auto& insideVolVars = elemVolVars[insideScvIdx];
        const Scalar velocitySelf = globalFaceVars.faceVars(scvf.dofIndexSelf()).velocity() ;
        const Scalar velocityOpposite = globalFaceVars.faceVars(scvf.dofIndexOpposite()).velocity();
        FacePrimaryVariables normalFlux(0.0);

        if(navierStokes)
        {
            // advective part
            const Scalar vAvg = (velocitySelf + velocityOpposite) * 0.5;
            const Scalar vUp = (sign(scvf.outerNormalScalar()) == sign(vAvg)) ? velocityOpposite : velocitySelf;
            normalFlux += vAvg * vUp * insideVolVars.density();
        }

        // diffusive part
        const Scalar deltaV = scvf.normalInPosCoordDir() ?
                              (velocitySelf - velocityOpposite) :
                              (velocityOpposite - velocitySelf);

        const Scalar deltaX = scvf.selfToOppositeDistance();
        normalFlux -= insideVolVars.viscosity() * 2.0 * deltaV/deltaX;

        // account for the orientation of the face
        const Scalar sgn = -1.0 * sign(scvf.outerNormalScalar());

        Scalar result = normalFlux * sgn * scvf.area();

        // treat outflow conditions
        if(navierStokes && scvf.boundary())
        {
            const auto& upVolVars = (sign(scvf.outerNormalScalar()) == sign(velocitySelf)) ?
                                    elemVolVars[insideScvIdx] : elemVolVars[scvf.outsideScvIdx()] ;

            result += velocitySelf * velocitySelf * upVolVars.density() * sign(scvf.outerNormalScalar()) * scvf.area() ;
        }
        return result;
    }

     /*!
     * \brief Returns the tangential part of the momentum flux
     * \param scvf The sub control volume face
     * \param fvGeometry The finite-volume geometry
     * \param elemVolVars All volume variables for the element
     * \param globalFaceVars The face variables
     */
    FacePrimaryVariables computeTangetialMomentumFlux_(const SubControlVolumeFace& scvf,
                                      const FVElementGeometry& fvGeometry,
                                      const ElementVolumeVariables& elemVolVars,
                                      const GlobalFaceVars& globalFaceVars)
    {
        FacePrimaryVariables tangentialFlux(0.0);

        // convenience function to get the velocity on a face
        auto velocity = [&globalFaceVars](const int dofIdx)
        {
            return globalFaceVars.faceVars(dofIdx).velocity();
        };

        // account for all sub-faces
        for(auto subFaceData : scvf.pairData())
        {
            const int eIdx = scvf.insideScvIdx();
            const auto& normalFace = fvGeometry.scvf(eIdx, subFaceData.localNormalFaceIdx);

            if(navierStokes)
                tangentialFlux += computeAdvectivePartOfTangentialMomentumFlux_(scvf, normalFace, subFaceData, elemVolVars, velocity);

            tangentialFlux += computeDiffusivePartOfTangentialMomentumFlux_(scvf, normalFace, subFaceData, elemVolVars, velocity);
        }
        return tangentialFlux;
    }


    template<class SubFaceData, class VelocityHelper>
    FacePrimaryVariables computeAdvectivePartOfTangentialMomentumFlux_(const SubControlVolumeFace& scvf,
                                                                       const SubControlVolumeFace& normalFace,
                                                                       const SubFaceData& subFaceData,
                                                                       const ElementVolumeVariables& elemVolVars,
                                                                       VelocityHelper velocity)
    {
        const Scalar transportingVelocity = velocity(subFaceData.normalPair.first);
        const auto insideScvIdx = normalFace.insideScvIdx();
        const auto outsideScvIdx = normalFace.outsideScvIdx();

        const bool innerElementIsUpstream = ( sign(normalFace.outerNormalScalar()) == sign(transportingVelocity) );

        const auto& upVolVars = innerElementIsUpstream ? elemVolVars[insideScvIdx] : elemVolVars[outsideScvIdx];

        Scalar transportedVelocity(0.0);

        if(innerElementIsUpstream)
            transportedVelocity = velocity(scvf.dofIndexSelf());
        else
        {
            const int outerDofIdx = subFaceData.outerParallelFaceDofIdx;
            if(outerDofIdx >= 0)
                transportedVelocity = velocity(outerDofIdx);
            else // this is the case when the outer parallal dof would lie outside the domain
            {
                const auto boundaryVelocity = this->problem().dirichletVelocityAtPos(subFaceData.virtualOuterParallelFaceDofPos);
                transportedVelocity = boundaryVelocity[scvf.directionIndex()];
            }
        }

        const Scalar momentum = upVolVars.density() * transportedVelocity;
        const int sgn = sign(normalFace.outerNormalScalar());

        return transportingVelocity * momentum * sgn * normalFace.area() * 0.5;
    }


    template<class SubFaceData, class VelocityHelper>
    FacePrimaryVariables computeDiffusivePartOfTangentialMomentumFlux_(const SubControlVolumeFace& scvf,
                                                                       const SubControlVolumeFace& normalFace,
                                                                       const SubFaceData& subFaceData,
                                                                       const ElementVolumeVariables& elemVolVars,
                                                                       VelocityHelper velocity)
    {
        FacePrimaryVariables tangentialDiffusiveFlux(0.0);

        const auto normalDirIdx = directionIndex(std::move(normalFace.unitOuterNormal()));
        const auto insideScvIdx = normalFace.insideScvIdx();
        const auto outsideScvIdx = normalFace.outsideScvIdx();

        const auto& insideVolVars = elemVolVars[insideScvIdx];
        const auto& outsideVolVars = elemVolVars[outsideScvIdx];

        // the averaged viscosity at the face normal to our face of interest (where we assemble the face residual)
        const Scalar muAvg = (insideVolVars.viscosity() + outsideVolVars.viscosity()) * 0.5;

        // the normal derivative
        const int innerNormalVelocityIdx = subFaceData.normalPair.first;
        const int outerNormalVelocityIdx = subFaceData.normalPair.second;

        const Scalar innerNormalVelocity = velocity(innerNormalVelocityIdx);

        const Scalar outerNormalVelocity = outerNormalVelocityIdx >= 0 ?
                                    velocity(outerNormalVelocityIdx) :
                                    this->problem().dirichletVelocityAtPos(subFaceData.virtualOuterNormalFaceDofPos)[normalDirIdx];

        const Scalar normalDeltaV = scvf.normalInPosCoordDir() ?
                                      (outerNormalVelocity - innerNormalVelocity) :
                                      (innerNormalVelocity - outerNormalVelocity);

        const Scalar normalDerivative = normalDeltaV / subFaceData.normalDistance;
        tangentialDiffusiveFlux -= muAvg * normalDerivative;

        // the parallel derivative
        const Scalar innerParallelVelocity = velocity(scvf.dofIndexSelf());

        const int outerParallelFaceDofIdx = subFaceData.outerParallelFaceDofIdx;
        const Scalar outerParallelVelocity = outerParallelFaceDofIdx >= 0 ?
                                             velocity(outerParallelFaceDofIdx) :
                                             this->problem().dirichletVelocityAtPos(subFaceData.virtualOuterParallelFaceDofPos)[scvf.directionIndex()];

        const Scalar parallelDeltaV = normalFace.normalInPosCoordDir() ?
                                     (outerParallelVelocity - innerParallelVelocity) :
                                     (innerParallelVelocity - outerParallelVelocity);

        const Scalar parallelDerivative = parallelDeltaV / subFaceData.parallelDistance;
        tangentialDiffusiveFlux -= muAvg * parallelDerivative;

        const Scalar sgn = sign(normalFace.outerNormalScalar());
        return tangentialDiffusiveFlux * sgn * normalFace.area() * 0.5;
    }


     /*!
     * \brief Returns the pressure term
     * \param scvf The sub control volume face
     * \param fvGeometry The finite-volume geometry
     * \param elemVolVars All volume variables for the element
     * \param globalFaceVars The face variables
     */
    FacePrimaryVariables computePressureTerm_(const SubControlVolumeFace& scvf,
                                      const FVElementGeometry& fvGeometry,
                                      const ElementVolumeVariables& elemVolVars,
                                      const GlobalFaceVars& globalFaceVars)
    {
        const auto insideScvIdx = scvf.insideScvIdx();
        const auto& insideVolVars = elemVolVars[insideScvIdx];

        Scalar result = insideVolVars.pressure() * scvf.area() * -1.0 * sign(scvf.outerNormalScalar());

        // treat outflow BCs
        if(scvf.boundary())
        {
            const Scalar pressure = this->problem().dirichletAtPos(scvf.center())[pressureIdx];
            result += pressure * scvf.area() * sign(scvf.outerNormalScalar());
        }

        return result;
    }

};

}

#endif   // DUMUX_CC_LOCAL_RESIDUAL_HH
