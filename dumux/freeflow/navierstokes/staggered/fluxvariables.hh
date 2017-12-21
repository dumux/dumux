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
 * \ingroup NavierStokesModel
 * \copydoc Dumux::NavierStokesFluxVariablesImpl
 */
#ifndef DUMUX_NAVIERSTOKES_STAGGERED_FLUXVARIABLES_HH
#define DUMUX_NAVIERSTOKES_STAGGERED_FLUXVARIABLES_HH

#include <dumux/common/properties.hh>
#include <dumux/discretization/fluxvariablesbase.hh>
#include <dumux/discretization/methods.hh>

namespace Dumux
{

namespace Properties
{
// forward declaration
NEW_PROP_TAG(EnableInertiaTerms);
}

// forward declaration
template<class TypeTag, DiscretizationMethods Method>
class NavierStokesFluxVariablesImpl;


/*!
 * \ingroup NavierStokesModel
 * \brief The flux variables class for the Navier-Stokes model using the staggered grid discretization.
 */
template<class TypeTag>
class NavierStokesFluxVariablesImpl<TypeTag, DiscretizationMethods::Staggered>
: public FluxVariablesBase<TypeTag>
{
    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
    using Problem = typename GET_PROP_TYPE(TypeTag, Problem);
    using Element = typename GridView::template Codim<0>::Entity;
    using FVElementGeometry = typename GET_PROP_TYPE(TypeTag, FVElementGeometry);
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using Indices = typename GET_PROP_TYPE(TypeTag, Indices);
    using ElementVolumeVariables = typename GET_PROP_TYPE(TypeTag, ElementVolumeVariables);
    using SubControlVolumeFace = typename GET_PROP_TYPE(TypeTag, SubControlVolumeFace);
    using FluxVariablesCache = typename GET_PROP_TYPE(TypeTag, FluxVariablesCache);
    using CellCenterPrimaryVariables = typename GET_PROP_TYPE(TypeTag, CellCenterPrimaryVariables);
    using FacePrimaryVariables = typename GET_PROP_TYPE(TypeTag, FacePrimaryVariables);
    using IndexType = typename GridView::IndexSet::IndexType;
    using Stencil = std::vector<IndexType>;

    using ElementFaceVariables = typename GET_PROP_TYPE(TypeTag, ElementFaceVariables);
    using FaceVariables = typename GET_PROP_TYPE(TypeTag, FaceVariables);

    static constexpr bool navierStokes = GET_PROP_VALUE(TypeTag, EnableInertiaTerms);

    enum {
         // grid and world dimension
        dim = GridView::dimension,
        dimWorld = GridView::dimensionworld,

        massBalanceIdx = Indices::massBalanceIdx,
    };

    using GlobalPosition = Dune::FieldVector<Scalar, dimWorld>;

    using DofTypeIndices = typename GET_PROP(TypeTag, DofTypeIndices);
    typename DofTypeIndices::CellCenterIdx cellCenterIdx;
    typename DofTypeIndices::FaceIdx faceIdx;

public:

    /*!
    * \brief Returns the advective flux over a sub control volume face
    * \param elemVolVars All volume variables for the element
    * \param elemFaceVars The face variables
    * \param scvf The sub control volume face
    * \param upwindTerm The uwind term (i.e. the advectively transported quantity)
    * \param isOutflow Determines if we are on an outflow boundary
    */
    template<class UpwindTerm>
    static Scalar advectiveFluxForCellCenter(const ElementVolumeVariables& elemVolVars,
                                             const ElementFaceVariables& elemFaceVars,
                                             const SubControlVolumeFace &scvf,
                                             UpwindTerm upwindTerm,
                                             bool isOutflow = false)
    {
        const Scalar velocity = elemFaceVars[scvf].velocitySelf();
        const bool insideIsUpstream = scvf.directionSign() == sign(velocity);
        static const Scalar upWindWeight = getParamFromGroup<Scalar>(GET_PROP_VALUE(TypeTag, ModelParameterGroup), "Implicit.UpwindWeight");

        const auto& insideVolVars = elemVolVars[scvf.insideScvIdx()];
        const auto& outsideVolVars = isOutflow ? insideVolVars : elemVolVars[scvf.outsideScvIdx()];

        const auto& upstreamVolVars = insideIsUpstream ? insideVolVars : outsideVolVars;
        const auto& downstreamVolVars = insideIsUpstream ? outsideVolVars : insideVolVars;

        const Scalar flux = (upWindWeight * upwindTerm(upstreamVolVars) +
                            (1.0 - upWindWeight) * upwindTerm(downstreamVolVars))
                            * velocity * scvf.area() * scvf.directionSign();

        return flux;
    }

    /*!
    * \brief Computes the flux for the cell center residual.
    */
    CellCenterPrimaryVariables computeFluxForCellCenter(const Problem& problem,
                                                        const Element &element,
                                                        const FVElementGeometry& fvGeometry,
                                                        const ElementVolumeVariables& elemVolVars,
                                                        const ElementFaceVariables& elemFaceVars,
                                                        const SubControlVolumeFace &scvf,
                                                        const FluxVariablesCache& fluxVarsCache)
    {
        auto upwindTerm = [](const auto& volVars) { return volVars.density(); };

        const Scalar flux = advectiveFluxForCellCenter(elemVolVars, elemFaceVars, scvf, upwindTerm, false);

        CellCenterPrimaryVariables result(0.0);
        result[massBalanceIdx] = flux;

        return result;
    }

    void computeCellCenterToCellCenterStencil(Stencil& stencil,
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
                                        const Element& element,
                                        const FVElementGeometry& fvGeometry,
                                        const SubControlVolumeFace& scvf)
    {
        stencil.push_back(scvf.dofIndex());
    }

    void computeFaceToCellCenterStencil(Stencil& stencil,
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
    */
   FacePrimaryVariables computeNormalMomentumFlux(const Problem& problem,
                                                  const Element& element,
                                                  const SubControlVolumeFace& scvf,
                                                  const FVElementGeometry& fvGeometry,
                                                  const ElementVolumeVariables& elemVolVars,
                                                  const ElementFaceVariables& elementFaceVars)
   {
       const auto insideScvIdx = scvf.insideScvIdx();
       const auto& insideVolVars = elemVolVars[insideScvIdx];
       const Scalar velocitySelf = elementFaceVars[scvf].velocitySelf() ;
       const Scalar velocityOpposite = elementFaceVars[scvf].velocityOpposite();
       FacePrimaryVariables normalFlux(0.0);

       if(navierStokes)
       {
           // advective part
           const Scalar vAvg = (velocitySelf + velocityOpposite) * 0.5;
           const Scalar vUp = (scvf.directionSign() == sign(vAvg)) ? velocityOpposite : velocitySelf;
           normalFlux += vAvg * vUp * insideVolVars.density();
       }

       // diffusive part
       const Scalar deltaV = scvf.normalInPosCoordDir() ?
                             (velocitySelf - velocityOpposite) :
                             (velocityOpposite - velocitySelf);

       const Scalar deltaX = scvf.selfToOppositeDistance();
       normalFlux -= insideVolVars.viscosity() * 2.0 * deltaV/deltaX;

       // account for the orientation of the face
       const Scalar sgn = -1.0 * scvf.directionSign();

       Scalar result = normalFlux * sgn * scvf.area();

       // treat outflow conditions
       if(navierStokes && scvf.boundary())
       {
           const auto& upVolVars = (scvf.directionSign() == sign(velocitySelf)) ?
                                   elemVolVars[insideScvIdx] : elemVolVars[scvf.outsideScvIdx()] ;

           result += velocitySelf * velocitySelf * upVolVars.density() * scvf.directionSign() * scvf.area() ;
       }
       return result;
   }

   /*!
   * \brief Returns the tangential part of the momentum flux
   */
  FacePrimaryVariables computeTangetialMomentumFlux(const Problem& problem,
                                                    const Element& element,
                                                    const SubControlVolumeFace& scvf,
                                                    const FVElementGeometry& fvGeometry,
                                                    const ElementVolumeVariables& elemVolVars,
                                                    const ElementFaceVariables& elementFaceVars)
  {
      FacePrimaryVariables tangentialFlux(0.0);
      auto& faceVars = elementFaceVars[scvf];
      const int numSubFaces = scvf.pairData().size();

      // account for all sub-faces
      for(int localSubFaceIdx = 0; localSubFaceIdx < numSubFaces; ++localSubFaceIdx)
      {
          const auto eIdx = scvf.insideScvIdx();
          const auto& normalFace = fvGeometry.scvf(eIdx, scvf.pairData()[localSubFaceIdx].localNormalFaceIdx);

          // Check if we have a symmetry boundary condition. If yes, the tangental part of the momentum flux can be neglected.
          if(scvf.pairData()[localSubFaceIdx].outerParallelFaceDofIdx < 0)
          {
              // lambda to conveniently create a ghost face which is outside the domain, parallel to the scvf of interest
              auto makeGhostFace = [eIdx] (const GlobalPosition& pos)
              {
                  return SubControlVolumeFace(pos, std::vector<unsigned int>{eIdx,eIdx});
              };

              // use the ghost face to check if there is a symmetry boundary condition and skip any further steps if yes
              const auto bcTypes = problem.boundaryTypes(element, makeGhostFace(scvf.pairData()[localSubFaceIdx].virtualOuterParallelFaceDofPos));
              if(bcTypes.isSymmetry())
                continue;
          }

          // if there is no symmetry boundary condition, proceed to calculate the tangential momentum flux
          if(navierStokes)
              tangentialFlux += computeAdvectivePartOfTangentialMomentumFlux_(problem, element, scvf, normalFace, elemVolVars, faceVars, localSubFaceIdx);

          tangentialFlux += computeDiffusivePartOfTangentialMomentumFlux_(problem, element, scvf, normalFace, elemVolVars, faceVars, localSubFaceIdx);
      }
      return tangentialFlux;
  }

private:

  FacePrimaryVariables computeAdvectivePartOfTangentialMomentumFlux_(const Problem& problem,
                                                                     const Element& element,
                                                                     const SubControlVolumeFace& scvf,
                                                                     const SubControlVolumeFace& normalFace,
                                                                     const ElementVolumeVariables& elemVolVars,
                                                                     const FaceVariables& faceVars,
                                                                     const int localSubFaceIdx)
  {
      const Scalar transportingVelocity = faceVars.velocityNormalInside(localSubFaceIdx);
      const auto insideScvIdx = normalFace.insideScvIdx();
      const auto outsideScvIdx = normalFace.outsideScvIdx();

      const bool innerElementIsUpstream = ( normalFace.directionSign() == sign(transportingVelocity) );

      const auto& upVolVars = innerElementIsUpstream ? elemVolVars[insideScvIdx] : elemVolVars[outsideScvIdx];

      const Scalar transportedVelocity = innerElementIsUpstream ?
                                         faceVars.velocitySelf() :
                                         faceVars.velocityParallel(localSubFaceIdx);

      const Scalar momentum = upVolVars.density() * transportedVelocity;

      return transportingVelocity * momentum * normalFace.directionSign() * normalFace.area() * 0.5;
  }

  FacePrimaryVariables computeDiffusivePartOfTangentialMomentumFlux_(const Problem& problem,
                                                                     const Element& element,
                                                                     const SubControlVolumeFace& scvf,
                                                                     const SubControlVolumeFace& normalFace,
                                                                     const ElementVolumeVariables& elemVolVars,
                                                                     const FaceVariables& faceVars,
                                                                     const int localSubFaceIdx)
  {
      FacePrimaryVariables tangentialDiffusiveFlux(0.0);

      const auto insideScvIdx = normalFace.insideScvIdx();
      const auto outsideScvIdx = normalFace.outsideScvIdx();

      const auto& insideVolVars = elemVolVars[insideScvIdx];
      const auto& outsideVolVars = elemVolVars[outsideScvIdx];

      // the averaged viscosity at the face normal to our face of interest (where we assemble the face residual)
      const Scalar muAvg = (insideVolVars.viscosity() + outsideVolVars.viscosity()) * 0.5;

      // the normal derivative
      const Scalar innerNormalVelocity = faceVars.velocityNormalInside(localSubFaceIdx);
      const Scalar outerNormalVelocity = faceVars.velocityNormalOutside(localSubFaceIdx);

      const Scalar normalDeltaV = scvf.normalInPosCoordDir() ?
                                    (outerNormalVelocity - innerNormalVelocity) :
                                    (innerNormalVelocity - outerNormalVelocity);

      const Scalar normalDerivative = normalDeltaV / scvf.pairData(localSubFaceIdx).normalDistance;
      tangentialDiffusiveFlux -= muAvg * normalDerivative;

      // the parallel derivative
      const Scalar innerParallelVelocity = faceVars.velocitySelf();

      const Scalar outerParallelVelocity = faceVars.velocityParallel(localSubFaceIdx);

      const Scalar parallelDeltaV = normalFace.normalInPosCoordDir() ?
                                   (outerParallelVelocity - innerParallelVelocity) :
                                   (innerParallelVelocity - outerParallelVelocity);

      const Scalar parallelDerivative = parallelDeltaV / scvf.pairData(localSubFaceIdx).parallelDistance;
      tangentialDiffusiveFlux -= muAvg * parallelDerivative;

      return tangentialDiffusiveFlux * normalFace.directionSign() * normalFace.area() * 0.5;
  }
};

} // end namespace

#endif
