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
 * \brief Base class for the flux variables
 */
#ifndef DUMUX_FREELOW_IMPLICIT_FLUXVARIABLES_HH
#define DUMUX_FREELOW_IMPLICIT_FLUXVARIABLES_HH

#include <dumux/implicit/properties.hh>
#include <dumux/discretization/fluxvariablesbase.hh>

namespace Dumux
{

namespace Properties
{
// forward declaration
NEW_PROP_TAG(EnableComponentTransport);
NEW_PROP_TAG(EnableInertiaTerms);
}

// forward declaration
template<class TypeTag, bool enableComponentTransport>
class FreeFlowFluxVariablesImpl;

/*!
 * \ingroup ImplicitModel
 * \brief The flux variables class
 *        specializations are provided for combinations of physical processes
 * \note  Not all specializations are currently implemented
 */
template<class TypeTag>
using FreeFlowFluxVariables = FreeFlowFluxVariablesImpl<TypeTag, GET_PROP_VALUE(TypeTag, EnableComponentTransport)>;

/*!
 * \ingroup Discretization
 * \brief Base class for the flux variables
 *        Actual flux variables inherit from this class
 */
// specialization for immiscible, isothermal flow
template<class TypeTag>
class FreeFlowFluxVariablesImpl<TypeTag, false>
: public FluxVariablesBase<TypeTag>
{
    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
    using Problem = typename GET_PROP_TYPE(TypeTag, Problem);
    using Element = typename GridView::template Codim<0>::Entity;
    using FVElementGeometry = typename GET_PROP_TYPE(TypeTag, FVElementGeometry);
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using Indices = typename GET_PROP_TYPE(TypeTag, Indices);
    using ElementVolumeVariables = typename GET_PROP_TYPE(TypeTag, ElementVolumeVariables);
    using GlobalFaceVars = typename GET_PROP_TYPE(TypeTag, GlobalFaceVars);
    using SubControlVolumeFace = typename GET_PROP_TYPE(TypeTag, SubControlVolumeFace);
    using FluxVariablesCache = typename GET_PROP_TYPE(TypeTag, FluxVariablesCache);
    using CellCenterPrimaryVariables = typename GET_PROP_TYPE(TypeTag, CellCenterPrimaryVariables);
    using FacePrimaryVariables = typename GET_PROP_TYPE(TypeTag, FacePrimaryVariables);
    using IndexType = typename GridView::IndexSet::IndexType;
    using Stencil = std::vector<IndexType>;

    static constexpr bool navierStokes = GET_PROP_VALUE(TypeTag, EnableInertiaTerms);

    enum {
         // grid and world dimension
        dim = GridView::dimension,
        dimWorld = GridView::dimensionworld,

        pressureIdx = Indices::pressureIdx,
        velocityIdx = Indices::velocityIdx,

        massBalanceIdx = Indices::massBalanceIdx,
        momentumBalanceIdx = Indices::momentumBalanceIdx
    };

    using GlobalPosition = Dune::FieldVector<Scalar, dimWorld>;

    using DofTypeIndices = typename GET_PROP(TypeTag, DofTypeIndices);
    typename DofTypeIndices::CellCenterIdx cellCenterIdx;
    typename DofTypeIndices::FaceIdx faceIdx;

public:

    CellCenterPrimaryVariables computeFluxForCellCenter(const Problem& problem,
                                                        const Element &element,
                                                        const FVElementGeometry& fvGeometry,
                                                        const ElementVolumeVariables& elemVolVars,
                                                        const GlobalFaceVars& globalFaceVars,
                                                        const SubControlVolumeFace &scvf,
                                                        const FluxVariablesCache& fluxVarsCache)
    {
        const auto& insideScv = fvGeometry.scv(scvf.insideScvIdx());
        const auto& insideVolVars = elemVolVars[insideScv];

        // if we are on an inflow/outflow boundary, use the volVars of the element itself
        const auto& outsideVolVars = scvf.boundary() ?  insideVolVars : elemVolVars[scvf.outsideScvIdx()];

        CellCenterPrimaryVariables flux(0.0);
        const Scalar velocity = globalFaceVars.faceVars(scvf.dofIndex()).velocity();

        const bool insideIsUpstream = sign(scvf.outerNormalScalar()) == sign(velocity) ? true : false;
        const auto& upstreamVolVars = insideIsUpstream ? insideVolVars : outsideVolVars;
        const auto& downstreamVolVars = insideIsUpstream ? insideVolVars : outsideVolVars;

        const Scalar upWindWeight = GET_PROP_VALUE(TypeTag, ImplicitUpwindWeight);

        flux = (upWindWeight * upstreamVolVars.density(Indices::phaseIdx) +
               (1.0 - upWindWeight) * downstreamVolVars.density(Indices::phaseIdx)) * velocity;

        flux *= scvf.area() * sign(scvf.outerNormalScalar());

        return flux;
    }

    void computeCellCenterToCellCenterStencil(Stencil& stencil,
                                              const Problem& problem,
                                              const Element& element,
                                              const FVElementGeometry& fvGeometry,
                                              const SubControlVolumeFace& scvf)
    {
        // the first entry is always the cc dofIdx itself
        if(stencil.empty())
            stencil.push_back(scvf.insideScvIdx());
        if(!scvf.boundary())
            stencil.push_back(scvf.outsideScvIdx());
    }

    void computeCellCenterToFaceStencil(Stencil& stencil,
                                        const Problem& problem,
                                        const Element& element,
                                        const FVElementGeometry& fvGeometry,
                                        const SubControlVolumeFace& scvf)
    {
        stencil.push_back(scvf.dofIndex());
    }

    void computeFaceToCellCenterStencil(Stencil& stencil,
                                        const Problem& problem,
                                        const FVElementGeometry& fvGeometry,
                                        const SubControlVolumeFace& scvf)
    {
        const int eIdx = scvf.insideScvIdx();
        stencil.push_back(scvf.insideScvIdx());

        for(const auto& data : scvf.pairData())
        {
            auto& normalFace = fvGeometry.scvf(eIdx, data.localNormalFaceIdx);
            const auto outerParallelElementDofIdx = normalFace.outsideScvIdx();
            if(!normalFace.boundary())
                stencil.push_back(outerParallelElementDofIdx);
        }
    }

    void computeFaceToFaceStencil(Stencil& stencil,
                                  const Problem& problem,
                                  const FVElementGeometry& fvGeometry,
                                  const SubControlVolumeFace& scvf)
    {
        // the first entries are always the face dofIdx itself and the one of the opposing face
        if(stencil.empty())
        {
            stencil.push_back(scvf.dofIndex());
            stencil.push_back(scvf.dofIndexOpposingFace());
        }

        for(const auto& data : scvf.pairData())
        {
            stencil.push_back(data.normalPair.first);
            const auto outerParallelFaceDofIdx = data.outerParallelFaceDofIdx;
            if(outerParallelFaceDofIdx >= 0)
                stencil.push_back(outerParallelFaceDofIdx);
            if(!scvf.boundary())
                stencil.push_back(data.normalPair.second);
        }
    }

    /*!
    * \brief Returns the normal part of the momentum flux
    * \param scvf The sub control volume face
    * \param fvGeometry The finite-volume geometry
    * \param elemVolVars All volume variables for the element
    * \param globalFaceVars The face variables
    */
   FacePrimaryVariables computeNormalMomentumFlux(const Problem& problem,
                                                  const Element& element,
                                                  const SubControlVolumeFace& scvf,
                                                  const FVElementGeometry& fvGeometry,
                                                  const ElementVolumeVariables& elemVolVars,
                                                  const GlobalFaceVars& globalFaceVars)
   {
       const auto insideScvIdx = scvf.insideScvIdx();
       const auto& insideVolVars = elemVolVars[insideScvIdx];
       const Scalar velocitySelf = globalFaceVars.faceVars(scvf.dofIndex()).velocity() ;
       const Scalar velocityOpposite = globalFaceVars.faceVars(scvf.dofIndexOpposingFace()).velocity();
       FacePrimaryVariables normalFlux(0.0);

       if(navierStokes)
       {
           // advective part
           const Scalar vAvg = (velocitySelf + velocityOpposite) * 0.5;
           const Scalar vUp = (sign(scvf.outerNormalScalar()) == sign(vAvg)) ? velocityOpposite : velocitySelf;
           normalFlux += vAvg * vUp * insideVolVars.density(Indices::phaseIdx);
       }

       // diffusive part
       const Scalar deltaV = scvf.normalInPosCoordDir() ?
                             (velocitySelf - velocityOpposite) :
                             (velocityOpposite - velocitySelf);

       const Scalar deltaX = scvf.selfToOppositeDistance();
       normalFlux -= insideVolVars.viscosity(Indices::phaseIdx) * 2.0 * deltaV/deltaX;

       // account for the orientation of the face
       const Scalar sgn = -1.0 * sign(scvf.outerNormalScalar());

       Scalar result = normalFlux * sgn * scvf.area();

       // treat outflow conditions
       if(navierStokes && scvf.boundary())
       {
           const auto& upVolVars = (sign(scvf.outerNormalScalar()) == sign(velocitySelf)) ?
                                   elemVolVars[insideScvIdx] : elemVolVars[scvf.outsideScvIdx()] ;

           result += velocitySelf * velocitySelf * upVolVars.density(Indices::phaseIdx) * sign(scvf.outerNormalScalar()) * scvf.area() ;
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
  FacePrimaryVariables computeTangetialMomentumFlux(const Problem& problem,
                                                    const Element& element,
                                                    const SubControlVolumeFace& scvf,
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
              tangentialFlux += computeAdvectivePartOfTangentialMomentumFlux_(problem, element, scvf, normalFace, subFaceData, elemVolVars, velocity);

          tangentialFlux += computeDiffusivePartOfTangentialMomentumFlux_(problem, element, scvf, normalFace, subFaceData, elemVolVars, velocity);
      }
      return tangentialFlux;
  }

private:

  template<class SubFaceData, class VelocityHelper>
  FacePrimaryVariables computeAdvectivePartOfTangentialMomentumFlux_(const Problem& problem,
                                                                     const Element& element,
                                                                     const SubControlVolumeFace& scvf,
                                                                     const SubControlVolumeFace& normalFace,
                                                                     const SubFaceData& subFaceData,
                                                                     const ElementVolumeVariables& elemVolVars,
                                                                     VelocityHelper velocity)
  {
      const Scalar transportingVelocity = velocity(subFaceData.normalPair.first);
      const auto insideScvIdx = normalFace.insideScvIdx();
      const auto outsideScvIdx = normalFace.outsideScvIdx();

      // convenience function to create a ghost face, outside of the domain
    //   auto makeGhostFace = [insideScvIdx, outsideScvIdx] (const GlobalPosition& pos)
    //   {
    //       return SubControlVolumeFace(pos, std::vector<unsigned int>{insideScvIdx,insideScvIdx});
    //   };

      const bool innerElementIsUpstream = ( sign(normalFace.outerNormalScalar()) == sign(transportingVelocity) );

      const auto& upVolVars = innerElementIsUpstream ? elemVolVars[insideScvIdx] : elemVolVars[outsideScvIdx];

      Scalar transportedVelocity(0.0);

      if(innerElementIsUpstream)
          transportedVelocity = velocity(scvf.dofIndex());
      else
      {
          const int outerDofIdx = subFaceData.outerParallelFaceDofIdx;
          if(outerDofIdx >= 0)
              transportedVelocity = velocity(outerDofIdx);
          else // this is the case when the outer parallal dof would lie outside the domain TODO: discuss which one is better
            //   transportedVelocity = problem.dirichlet(makeGhostFace(subFaceData.virtualOuterParallelFaceDofPos))[faceIdx][scvf.directionIndex()];
              transportedVelocity = problem.dirichlet(element, scvf)[faceIdx][scvf.directionIndex()];
      }

      const Scalar momentum = upVolVars.density(Indices::phaseIdx) * transportedVelocity;
      const int sgn = sign(normalFace.outerNormalScalar());

      return transportingVelocity * momentum * sgn * normalFace.area() * 0.5;
  }

  template<class SubFaceData, class VelocityHelper>
  FacePrimaryVariables computeDiffusivePartOfTangentialMomentumFlux_(const Problem& problem,
                                                                     const Element& element,
                                                                     const SubControlVolumeFace& scvf,
                                                                     const SubControlVolumeFace& normalFace,
                                                                     const SubFaceData& subFaceData,
                                                                     const ElementVolumeVariables& elemVolVars,
                                                                     VelocityHelper velocity)
  {
      FacePrimaryVariables tangentialDiffusiveFlux(0.0);

      const auto normalDirIdx = normalFace.directionIndex();
      const auto insideScvIdx = normalFace.insideScvIdx();
      const auto outsideScvIdx = normalFace.outsideScvIdx();

      const auto& insideVolVars = elemVolVars[insideScvIdx];
      const auto& outsideVolVars = elemVolVars[outsideScvIdx];

      // convenience function to create a ghost face, outside of the domain
      auto makeGhostFace = [insideScvIdx, outsideScvIdx] (const GlobalPosition& pos)
      {
          return SubControlVolumeFace(pos, std::vector<unsigned int>{insideScvIdx,insideScvIdx});
      };

      // the averaged viscosity at the face normal to our face of interest (where we assemble the face residual)
      const Scalar muAvg = (insideVolVars.viscosity(Indices::phaseIdx) + outsideVolVars.viscosity(Indices::phaseIdx)) * 0.5;

      // the normal derivative
      const int innerNormalVelocityIdx = subFaceData.normalPair.first;
      const int outerNormalVelocityIdx = subFaceData.normalPair.second;

      const Scalar innerNormalVelocity = velocity(innerNormalVelocityIdx);

      const Scalar outerNormalVelocity = outerNormalVelocityIdx >= 0 ?
                                  velocity(outerNormalVelocityIdx) :
                                  problem.dirichlet(element, makeGhostFace(subFaceData.virtualOuterNormalFaceDofPos))[faceIdx][normalDirIdx];

      const Scalar normalDeltaV = scvf.normalInPosCoordDir() ?
                                    (outerNormalVelocity - innerNormalVelocity) :
                                    (innerNormalVelocity - outerNormalVelocity);

      const Scalar normalDerivative = normalDeltaV / subFaceData.normalDistance;
      tangentialDiffusiveFlux -= muAvg * normalDerivative;

      // the parallel derivative
      const Scalar innerParallelVelocity = velocity(scvf.dofIndex());

      const int outerParallelFaceDofIdx = subFaceData.outerParallelFaceDofIdx;
      const Scalar outerParallelVelocity = outerParallelFaceDofIdx >= 0 ?
                                           velocity(outerParallelFaceDofIdx) :
                                           problem.dirichlet(element, makeGhostFace(subFaceData.virtualOuterParallelFaceDofPos))[faceIdx][scvf.directionIndex()];

      const Scalar parallelDeltaV = normalFace.normalInPosCoordDir() ?
                                   (outerParallelVelocity - innerParallelVelocity) :
                                   (innerParallelVelocity - outerParallelVelocity);

      const Scalar parallelDerivative = parallelDeltaV / subFaceData.parallelDistance;
      tangentialDiffusiveFlux -= muAvg * parallelDerivative;

      const Scalar sgn = sign(normalFace.outerNormalScalar());
      return tangentialDiffusiveFlux * sgn * normalFace.area() * 0.5;
  }
};

} // end namespace

#endif
