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

#include <dumux/common/math.hh>
#include <dumux/common/properties.hh>

#include <dumux/discretization/fluxvariablesbase.hh>
#include <dumux/discretization/methods.hh>

namespace Dumux {

// forward declaration
template<class TypeTag, DiscretizationMethod discMethod>
class NavierStokesFluxVariablesImpl;


/*!
 * \ingroup NavierStokesModel
 * \brief The flux variables class for the Navier-Stokes model using the staggered grid discretization.
 */
template<class TypeTag>
class NavierStokesFluxVariablesImpl<TypeTag, DiscretizationMethod::staggered>
: public FluxVariablesBase<TypeTag>
{
    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
    using Problem = typename GET_PROP_TYPE(TypeTag, Problem);
    using Element = typename GridView::template Codim<0>::Entity;
    using FVElementGeometry = typename GET_PROP_TYPE(TypeTag, FVGridGeometry)::LocalView;
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using Indices = typename GET_PROP_TYPE(TypeTag, Indices);
    using ElementVolumeVariables = typename GET_PROP_TYPE(TypeTag, GridVolumeVariables)::LocalView;
    using VolumeVariables = typename GET_PROP_TYPE(TypeTag, VolumeVariables);
    using SubControlVolumeFace = typename FVElementGeometry::SubControlVolumeFace;
    using FluxVariablesCache = typename GET_PROP_TYPE(TypeTag, FluxVariablesCache);
    using CellCenterPrimaryVariables = typename GET_PROP_TYPE(TypeTag, CellCenterPrimaryVariables);
    using FacePrimaryVariables = typename GET_PROP_TYPE(TypeTag, FacePrimaryVariables);

    using ElementFaceVariables = typename GET_PROP_TYPE(TypeTag, ElementFaceVariables);
    using FaceVariables = typename GET_PROP_TYPE(TypeTag, FaceVariables);

    static constexpr bool enableInertiaTerms = GET_PROP_VALUE(TypeTag, EnableInertiaTerms);
    static constexpr bool normalizePressure = GET_PROP_VALUE(TypeTag, NormalizePressure);

    using GlobalPosition = Dune::FieldVector<Scalar, GridView::dimensionworld>;

    using DofTypeIndices = typename GET_PROP(TypeTag, DofTypeIndices);
    typename DofTypeIndices::CellCenterIdx cellCenterIdx;
    typename DofTypeIndices::FaceIdx faceIdx;

public:

    /*!
    * \brief Returns the advective flux over a sub control volume face.
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
    * \brief Computes the flux for the cell center residual (mass balance).
    *
    * \verbatim
    *                    scvf
    *              ----------------
    *              |              # current scvf    # scvf over which fluxes are calculated
    *              |              #
    *              |      x       #~~~~> vel.Self   x dof position
    *              |              #
    *        scvf  |              #                 -- element
    *              ----------------
    *                   scvf
    * \endverbatim
    */
    CellCenterPrimaryVariables computeFluxForCellCenter(const Problem& problem,
                                                        const Element& element,
                                                        const FVElementGeometry& fvGeometry,
                                                        const ElementVolumeVariables& elemVolVars,
                                                        const ElementFaceVariables& elemFaceVars,
                                                        const SubControlVolumeFace& scvf,
                                                        const FluxVariablesCache& fluxVarsCache)
    {
        // The advectively transported quantity (i.e density for a single-phase system).
        auto upwindTerm = [](const auto& volVars) { return volVars.density(); };

        // Check if we are on an outflow boundary.
        const bool isOutflow = scvf.boundary() ?
                               problem.boundaryTypesAtPos(scvf.center()).isOutflow(Indices::totalMassBalanceIdx)
                             : false;

        // Call the generic flux function.
        const Scalar flux = advectiveFluxForCellCenter(elemVolVars, elemFaceVars, scvf, upwindTerm, isOutflow);

        CellCenterPrimaryVariables result(0.0);
        result[Indices::totalMassBalanceIdx] = flux;

        return result;
    }

    /*!
    * \brief Returns the momentum flux over all staggered faces.
    *
    */
    FacePrimaryVariables computeMomentumFlux(const Problem& problem,
                                             const Element& element,
                                             const SubControlVolumeFace& scvf,
                                             const FVElementGeometry& fvGeometry,
                                             const ElementVolumeVariables& elemVolVars,
                                             const ElementFaceVariables& elemFaceVars)
    {
        return computeFrontalMomentumFlux(problem, element, scvf, fvGeometry, elemVolVars, elemFaceVars) +
               computeNormalMomentumFlux(problem, element, scvf, fvGeometry, elemVolVars, elemFaceVars);
    }

    /*!
    * \brief Returns the frontal part of the momentum flux.
    *        This treats the flux over the staggered face at the center of an element,
    *        parallel to the current scvf where the velocity dof of interest lives.
    *
    * \verbatim
    *                    scvf
    *              ---------=======                 == and # staggered half-control-volume
    *              |       #      | current scvf
    *              |       #      |                 # staggered face over wich fluxes are calculated
    *   vel.Opp <~~|       #~~>   x~~~~> vel.Self
    *              |       #      |                 x dof position
    *        scvf  |       #      |
    *              --------========                 -- element
    *                   scvf
    * \endverbatim
    */
    FacePrimaryVariables computeFrontalMomentumFlux(const Problem& problem,
                                                    const Element& element,
                                                    const SubControlVolumeFace& scvf,
                                                    const FVElementGeometry& fvGeometry,
                                                    const ElementVolumeVariables& elemVolVars,
                                                    const ElementFaceVariables& elemFaceVars)
    {
        FacePrimaryVariables normalFlux(0.0);

        // The velocities of the dof at interest and the one of the opposite scvf.
        const Scalar velocitySelf = elemFaceVars[scvf].velocitySelf();
        const Scalar velocityOpposite = elemFaceVars[scvf].velocityOpposite();

        // The volume variables within the current element. We only require those (and none of neighboring elements)
        // because the fluxes are calculated over the staggered face at the center of the element.
        const auto& insideVolVars = elemVolVars[scvf.insideScvIdx()];

        // Advective flux.
        if(enableInertiaTerms)
        {
            // Get the average velocity at the center of the element (i.e. the location of the staggered face).
            const Scalar transportingVelocity = (velocitySelf + velocityOpposite) * 0.5;

            // Check if the the velocity of the dof at interest lies up- or downstream w.r.t. to the transporting velocity.
            const bool selfIsUpstream = scvf.directionSign() != sign(transportingVelocity);

            // Lamba function to evaluate the transported momentum, regarding an user-specified upwind weight.
            auto computeMomentum = [&insideVolVars](const Scalar upstreamVelocity, const Scalar downstreamVelocity)
            {
                static const Scalar upwindWeight = getParamFromGroup<Scalar>(GET_PROP_VALUE(TypeTag, ModelParameterGroup), "Implicit.UpwindWeight");
                return (upwindWeight * upstreamVelocity + (1.0 - upwindWeight) * downstreamVelocity) * insideVolVars.density();
            };


            // Get the momentum that is advectively transported and account for the flow direction.
            const Scalar momentum = selfIsUpstream ? computeMomentum(velocitySelf, velocityOpposite)
                                                   : computeMomentum(velocityOpposite, velocitySelf);

            // Account for the orientation of the staggered face's normal outer normal vector
            // (pointing in opposite direction of the scvf's one).
            normalFlux += transportingVelocity * momentum * -1.0 * scvf.directionSign();
        }

        // Diffusive flux.
        // The velocity gradient already accounts for the orientation
        // of the staggered face's outer normal vector.
        const Scalar gradV = (velocityOpposite - velocitySelf) / scvf.selfToOppositeDistance();
        normalFlux -= insideVolVars.effectiveViscosity() * 2.0 * gradV;

        // The pressure term.
        // If specified, the pressure can be normalized using the initial value on the scfv of interest.
        // Can potentially help to improve the condition number of the system matrix.
        const Scalar pressure = normalizePressure ?
                                insideVolVars.pressure() - problem.initialAtPos(scvf.center())[Indices::pressureIdx]
                              : insideVolVars.pressure();

        // Account for the orientation of the staggered face's normal outer normal vector
        // (pointing in opposite direction of the scvf's one).
        normalFlux += pressure * -1.0 * scvf.directionSign();

        // Treat outflow conditions.
        if(scvf.boundary())
        {
            if(problem.boundaryTypesAtPos(scvf.center()).isOutflow(Indices::momentumBalanceIdx))
            {
                // Treat the staggered half-volume adjacent to the boundary as if it was on the opposite side of the boundary.
                // The respective face's outer normal vector will point in the same direction as the scvf's one.
                normalFlux += outflowBoundaryFlux_(problem, element, scvf, elemVolVars, elemFaceVars);
            }
        }

        // Account for the staggered face's area. For rectangular elements, this equals the area of the scvf
        // our velocity dof of interest lives on.
        return normalFlux *  scvf.area();
   }

    /*!
    * \brief Returns the momentum flux over the staggered faces
    *        perpendicular to the scvf where the velocity dof of interest
    *        lives (coinciding with the element's scvfs).
    *
    * \verbatim
    *                scvf
    *              ---------#######                 || and # staggered half-control-volume
    *              |      ||      | current scvf
    *              |      ||      |                 # normal staggered faces over wich fluxes are calculated
    *              |      ||      x~~~~> vel.Self
    *              |      ||      |                 x dof position
    *        scvf  |      ||      |
    *              --------########                -- element
    *                 scvf
    * \endverbatim
    */
    FacePrimaryVariables computeNormalMomentumFlux(const Problem& problem,
                                                   const Element& element,
                                                   const SubControlVolumeFace& scvf,
                                                   const FVElementGeometry& fvGeometry,
                                                   const ElementVolumeVariables& elemVolVars,
                                                   const ElementFaceVariables& elemFaceVars)
    {
        FacePrimaryVariables normalFlux(0.0);
        auto& faceVars = elemFaceVars[scvf];
        const int numSubFaces = scvf.pairData().size();

        // Account for all sub-faces.
        for(int localSubFaceIdx = 0; localSubFaceIdx < numSubFaces; ++localSubFaceIdx)
        {
            const auto eIdx = scvf.insideScvIdx();
            const auto& normalFace = fvGeometry.scvf(eIdx, scvf.pairData()[localSubFaceIdx].localNormalFaceIdx);

            // Check if we have a symmetry boundary condition. If yes, the tangental part of the momentum flux can be neglected.
            if(!scvf.hasParallelNeighbor(localSubFaceIdx))
            {
                // Use the ghost face to check if there is a symmetry boundary condition and skip any further steps if yes.
                const auto bcTypes = problem.boundaryTypes(element, makeParallelGhostFace_(scvf, localSubFaceIdx));
                if(bcTypes.isSymmetry())
                    continue;
            }

            // If there is no symmetry boundary condition, proceed to calculate the tangential momentum flux.
            if(enableInertiaTerms)
                normalFlux += computeAdvectivePartOfNormalMomentumFlux_(problem, element, scvf, normalFace, elemVolVars, faceVars, localSubFaceIdx);

            normalFlux += computeDiffusivePartOfNormalMomentumFlux_(problem, element, scvf, normalFace, elemVolVars, faceVars, localSubFaceIdx);
        }
        return normalFlux;
    }

private:

    /*!
    * \brief Returns the advective momentum flux over the staggered face perpendicular to the scvf
    *        where the velocity dof of interest lives (coinciding with the element's scvfs).
    *
    * \verbatim
    *              ----------------
    *              |              |
    *              |    transp.   |
    *              |      vel.    |~~~~> vel.Parallel
    *              |       ^      |
    *              |       |      |
    *       scvf   ---------#######                 || and # staggered half-control-volume
    *              |      ||      | current scvf
    *              |      ||      |                 # normal staggered faces over wich fluxes are calculated
    *              |      ||      x~~~~> vel.Self
    *              |      ||      |                 x dof position
    *        scvf  |      ||      |
    *              ---------#######                -- elements
    *                 scvf
    * \endverbatim
    */
    FacePrimaryVariables computeAdvectivePartOfNormalMomentumFlux_(const Problem& problem,
                                                                   const Element& element,
                                                                   const SubControlVolumeFace& scvf,
                                                                   const SubControlVolumeFace& normalFace,
                                                                   const ElementVolumeVariables& elemVolVars,
                                                                   const FaceVariables& faceVars,
                                                                   const int localSubFaceIdx)
    {
        // Get the transporting velocity, located at the scvf perpendicular to the current scvf where the dof
        // of interest is located.
        const Scalar transportingVelocity = faceVars.velocityNormalInside(localSubFaceIdx);

        // Check whether the own or the neighboring element is upstream.
        const bool ownElementIsUpstream = ( normalFace.directionSign() == sign(transportingVelocity) );

        // Get the velocities at the current (own) scvf and at the parallel one at the neighboring scvf.
        const Scalar velocitySelf = faceVars.velocitySelf();

        // Lambda to conveniently get the outer parallel velocity for normal faces that are on the boundary
        // and therefore have no neighbor. Calls the problem to retrieve a fixed value set on the boundary.
        auto getParallelVelocityFromBoundary = [&]()
        {
            const auto ghostFace = makeParallelGhostFace_(scvf, localSubFaceIdx);
            return problem.dirichlet(element, ghostFace)[Indices::velocity(scvf.directionIndex())];
        };

        const Scalar velocityParallel = scvf.hasParallelNeighbor(localSubFaceIdx) ?
                                        faceVars.velocityParallel(localSubFaceIdx)
                                      : getParallelVelocityFromBoundary();

        // Get the volume variables of the own and the neighboring element
        const auto& insideVolVars = elemVolVars[normalFace.insideScvIdx()];
        const auto& outsideVolVars = elemVolVars[normalFace.outsideScvIdx()];

        // Lamba function to evaluate the transported momentum, regarding an user-specified upwind weight.
        auto computeMomentum = [](const VolumeVariables& upstreamVolVars,
                                  const VolumeVariables& downstreamVolVars,
                                  const Scalar upstreamVelocity,
                                  const Scalar downstreamVelocity)
        {
            static const Scalar upWindWeight = getParamFromGroup<Scalar>(GET_PROP_VALUE(TypeTag, ModelParameterGroup), "Implicit.UpwindWeight");
            const Scalar density = upWindWeight * upstreamVolVars.density() + (1.0 - upWindWeight) * downstreamVolVars.density();
            const Scalar transportedVelocity =  upWindWeight * upstreamVelocity + (1.0 - upWindWeight) * downstreamVelocity;
            return transportedVelocity * density;
        };

        // Get the momentum that is advectively transported and account for the flow direction.
        const Scalar momentum = ownElementIsUpstream ?
                                computeMomentum(insideVolVars, outsideVolVars, velocitySelf, velocityParallel)
                              : computeMomentum(outsideVolVars, insideVolVars, velocityParallel, velocitySelf);

        // Account for the orientation of the staggered normal face's outer normal vector
        // and its area (0.5 of the coinciding scfv).
        return transportingVelocity * momentum * normalFace.directionSign() * normalFace.area() * 0.5;
    }

    /*!
    * \brief Returns the diffusive momentum flux over the staggered face perpendicular to the scvf
    *        where the velocity dof of interest lives (coinciding with the element's scvfs).
    *
    * \verbatim
    *              ----------------
    *              |              |vel.
    *              |    in.norm.  |Parallel
    *              |       vel.   |~~~~>
    *              |       ^      |        ^ out.norm.vel.
    *              |       |      |        |
    *       scvf   ---------#######:::::::::       || and # staggered half-control-volume (own element)
    *              |      ||      | curr. ::
    *              |      ||      | scvf  ::       :: staggered half-control-volume (neighbor element)
    *              |      ||      x~~~~>  ::
    *              |      ||      | vel.  ::       # normal staggered faces over wich fluxes are calculated
    *        scvf  |      ||      | Self  ::
    *              ---------#######:::::::::       x dof position
    *                 scvf
    *                                              -- elements
    * \endverbatim
    */
    FacePrimaryVariables computeDiffusivePartOfNormalMomentumFlux_(const Problem& problem,
                                                                   const Element& element,
                                                                   const SubControlVolumeFace& scvf,
                                                                   const SubControlVolumeFace& normalFace,
                                                                   const ElementVolumeVariables& elemVolVars,
                                                                   const FaceVariables& faceVars,
                                                                   const int localSubFaceIdx)
    {
        FacePrimaryVariables normalDiffusiveFlux(0.0);

        // Get the volume variables of the own and the neighboring element. The neighboring
        // element is adjacent to the staggered face normal to the current scvf
        // where the dof of interest is located.
        const auto& insideVolVars = elemVolVars[normalFace.insideScvIdx()];
        const auto& outsideVolVars = elemVolVars[normalFace.outsideScvIdx()];

        // Get the averaged viscosity at the staggered face normal to the current scvf.
        const Scalar muAvg = (insideVolVars.effectiveViscosity() + outsideVolVars.effectiveViscosity()) * 0.5;

        // For the normal gradient, get the velocities perpendicular to the velocity at the current scvf.
        // The inner one is located at staggered face within the own element,
        // the outer one at the respective staggered face of the element on the other side of the
        // current scvf.
        const Scalar innerNormalVelocity = faceVars.velocityNormalInside(localSubFaceIdx);

        // Lambda to conveniently get the outer normal velocity for faces that are on the boundary
        // and therefore have no neighbor. Calls the problem to retrieve a fixed value set on the boundary.
        auto getNormalVelocityFromBoundary = [&]()
        {
            const auto ghostFace = makeNormalGhostFace_(scvf, localSubFaceIdx);
            return problem.dirichlet(element, ghostFace)[Indices::velocity(normalFace.directionIndex())];
        };

        const Scalar outerNormalVelocity = scvf.hasFrontalNeighbor(localSubFaceIdx) ?
                                           faceVars.velocityNormalOutside(localSubFaceIdx)
                                         : getNormalVelocityFromBoundary();


        // Calculate the velocity gradient in positive coordinate direction.
        const Scalar normalDeltaV = scvf.normalInPosCoordDir() ?
                                    (outerNormalVelocity - innerNormalVelocity)
                                  : (innerNormalVelocity - outerNormalVelocity);

        const Scalar normalGradient = normalDeltaV / scvf.pairData(localSubFaceIdx).normalDistance;

        // Account for the orientation of the staggered normal face's outer normal vector.
        normalDiffusiveFlux -= muAvg * normalGradient * normalFace.directionSign();

        // For the parallel derivative, get the velocities at the current (own) scvf
        // and at the parallel one at the neighboring scvf.
        const Scalar innerParallelVelocity = faceVars.velocitySelf();

        // Lambda to conveniently get the outer parallel velocity for normal faces that are on the boundary
        // and therefore have no neighbor. Calls the problem to retrieve a fixed value set on the boundary.
        auto getParallelVelocityFromBoundary = [&]()
        {
            const auto ghostFace = makeParallelGhostFace_(scvf, localSubFaceIdx);
            return problem.dirichlet(element, ghostFace)[Indices::velocity(scvf.directionIndex())];
        };

        const Scalar outerParallelVelocity = scvf.hasParallelNeighbor(localSubFaceIdx) ?
                                             faceVars.velocityParallel(localSubFaceIdx)
                                           : getParallelVelocityFromBoundary();

        // The velocity gradient already accounts for the orientation
        // of the staggered face's outer normal vector.
        const Scalar parallelGradient = (outerParallelVelocity - innerParallelVelocity)
                                        / scvf.pairData(localSubFaceIdx).parallelDistance;

        normalDiffusiveFlux -= muAvg * parallelGradient;

        // Account for the area of the staggered normal face (0.5 of the coinciding scfv).
        return normalDiffusiveFlux * normalFace.area() * 0.5;
    }

    /*!
    * \brief Returns the momentum flux over an outflow boundary face.
    *
    * \verbatim
    *                    scvf      //
    *              ---------=======//               == and # staggered half-control-volume
    *              |      ||      #// current scvf
    *              |      ||      #//               # staggered boundary face over wich fluxes are calculated
    *              |      ||      x~~~~> vel.Self
    *              |      ||      #//               x dof position
    *        scvf  |      ||      #//
    *              --------========//               -- element
    *                   scvf       //
    *                                              // boundary
    * \endverbatim
    */
    FacePrimaryVariables outflowBoundaryFlux_(const Problem& problem,
                                              const Element& element,
                                              const SubControlVolumeFace& scvf,
                                              const ElementVolumeVariables& elemVolVars,
                                              const ElementFaceVariables& elemFaceVars)
    {
        FacePrimaryVariables outflow(0.0);

        // Advective momentum outflow.
        if(enableInertiaTerms)
        {
            const Scalar velocitySelf = elemFaceVars[scvf].velocitySelf();
            const auto& insideVolVars = elemVolVars[scvf.insideScvIdx()];
            const auto& outsideVolVars = elemVolVars[scvf.outsideScvIdx()];
            const auto& upVolVars = (scvf.directionSign() == sign(velocitySelf)) ?
                                    insideVolVars : outsideVolVars;

            outflow += velocitySelf * velocitySelf * upVolVars.density();
        }

        // Apply a pressure at the boudary.
        const Scalar boundaryPressure = normalizePressure ?
                                        (problem.dirichlet(element, scvf)[Indices::pressureIdx] -
                                         problem.initialAtPos(scvf.center())[Indices::pressureIdx])
                                      : problem.dirichlet(element, scvf)[Indices::pressureIdx];
        outflow += boundaryPressure;

        // Account for the orientation of the face at the boundary,
        return outflow * scvf.directionSign();
    }

private:

    //! helper functiuob to conveniently create a ghost face which is outside the domain, parallel to the scvf of interest
    SubControlVolumeFace makeGhostFace_(const SubControlVolumeFace& ownScvf, const GlobalPosition& pos) const
    {
        return SubControlVolumeFace(pos, std::vector<unsigned int>{ownScvf.insideScvIdx(), ownScvf.outsideScvIdx()}, ownScvf.dofIndex());
    };

    //! helper functiuob to conveniently create a ghost face which is outside the domain, parallel to the scvf of interest
    SubControlVolumeFace makeNormalGhostFace_(const SubControlVolumeFace& ownScvf, const int localSubFaceIdx) const
    {
        return makeGhostFace_(ownScvf, ownScvf.pairData(localSubFaceIdx).virtualOuterNormalFaceDofPos);
    };

    //! helper functiuob to conveniently create a ghost face which is outside the domain, parallel to the scvf of interest
    SubControlVolumeFace makeParallelGhostFace_(const SubControlVolumeFace& ownScvf, const int localSubFaceIdx) const
    {
        return makeGhostFace_(ownScvf, ownScvf.pairData(localSubFaceIdx).virtualOuterParallelFaceDofPos);
    };

};

} // end namespace

#endif
