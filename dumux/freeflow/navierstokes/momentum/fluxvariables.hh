// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   See the file COPYING for full copying permissions.                      *
 *                                                                           *
 *   This program is free software: you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation, either version 3 of the License, or       *
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
#ifndef DUMUX_NAVIERSTOKES_MOMENTUM_FLUXVARIABLES_HH
#define DUMUX_NAVIERSTOKES_MOMENTUM_FLUXVARIABLES_HH

#include <array>
#include <optional>

#include <dumux/common/math.hh>
#include <dumux/common/exceptions.hh>
#include <dumux/common/parameters.hh>
#include <dumux/common/properties.hh>

#include <dumux/flux/fluxvariablesbase.hh>
#include <dumux/discretization/method.hh>


// #include "staggeredupwindhelper.hh"
#include "velocitygradients.hh"

namespace Dumux {

/*!
 * \ingroup NavierStokesModel
 * \brief The flux variables class for the Navier-Stokes model using the staggered grid discretization.
 */
template<class TypeTag>
class NavierStokesMomentumFluxVariables
: public FluxVariablesBase<GetPropType<TypeTag, Properties::Problem>,
                           typename GetPropType<TypeTag, Properties::GridGeometry>::LocalView,
                           typename GetPropType<TypeTag, Properties::GridVolumeVariables>::LocalView,
                           typename GetPropType<TypeTag, Properties::GridFluxVariablesCache>::LocalView>
{
    using GridVariables = GetPropType<TypeTag, Properties::GridVariables>;

    using GridVolumeVariables = typename GridVariables::GridVolumeVariables;
    using ElementVolumeVariables = typename GridVolumeVariables::LocalView;
    using VolumeVariables = typename GridVolumeVariables::VolumeVariables;

    using GridFluxVariablesCache = typename GridVariables::GridFluxVariablesCache;
    using FluxVariablesCache = typename GridFluxVariablesCache::FluxVariablesCache;

    using GridGeometry = typename GridVariables::GridGeometry;
    using FVElementGeometry = typename GridGeometry::LocalView;
    using SubControlVolume = typename FVElementGeometry::SubControlVolume;
    using SubControlVolumeFace = typename FVElementGeometry::SubControlVolumeFace;


    using Problem = GetPropType<TypeTag, Properties::Problem>;
    using GridView = typename GridGeometry::GridView;
    using Element = typename GridView::template Codim<0>::Entity;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using ModelTraits = GetPropType<TypeTag, Properties::ModelTraits>;
    using Indices = typename ModelTraits::Indices;
    using VelocityGradients = StaggeredVelocityGradients;

    static constexpr bool normalizePressure = getPropValue<TypeTag, Properties::NormalizePressure>();

    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;

    using NumEqVector = GetPropType<TypeTag, Properties::NumEqVector>;


public:

    /*!
     * \brief Returns the diffusive momentum flux due to viscous forces
     */
    NumEqVector diffusiveMomentumFlux() const
    {
        if (this->scvFace().isFrontal())
            return frontalDiffusiveMomentumFlux();
        else
            return lateralDiffusiveMomentumFlux();
    }

        /*!
     * \brief Returns the frontal part of the momentum flux.
     *        This treats the flux over the staggered face at the center of an element,
     *        parallel to the current scvf where the velocity dof of interest lives.
     *
     * \verbatim
     *                    scvf
     *              ---------=======                 == and # staggered half-control-volume
     *              |       #      | current scv
     *              |       #      |                 # staggered face over which fluxes are calculated
     *   vel.Opp <~~|       #~~>   x~~~~> vel.Self
     *              |       #      |                 x dof position
     *        scvf  |       #      |
     *              --------========                 -- element
     *                   scvf
     * \endverbatim
     */
    NumEqVector frontalDiffusiveMomentumFlux() const
    {
        const auto& scvf = this->scvFace();
        assert(scvf.isFrontal());

        NumEqVector result(0.0);
        const auto& fvGeometry = this->fvGeometry();
        const auto& elemVolVars = this->elemVolVars();

        if (scvf.boundary())
            return result;

        const Scalar velocityGrad_ii = VelocityGradients::velocityGradII(fvGeometry, scvf, elemVolVars) * scvf.directionSign();

        static const bool enableUnsymmetrizedVelocityGradient
        = getParamFromGroup<bool>(this->problem().paramGroup(), "FreeFlow.EnableUnsymmetrizedVelocityGradient", false);
        const Scalar factor = enableUnsymmetrizedVelocityGradient ? 1.0 : 2.0;

        const auto mu = this->problem().effectiveViscosity(this->element(), this->fvGeometry(), this->scvFace());
        result -= factor * mu * velocityGrad_ii * scvf.area() * elemVolVars[scvf.insideScvIdx()].extrusionFactor();

        return result;
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
     *              |      ||      | vel.  ::       # lateral staggered faces over which fluxes are calculated
     *        scvf  |      ||      | Self  ::
     *              ---------#######:::::::::       x dof position
     *                 scvf
     *                                              -- elements
     * \endverbatim
     */
    NumEqVector lateralDiffusiveMomentumFlux() const
    {
        const auto& scvf = this->scvFace();
        assert(scvf.isLateral());

        NumEqVector result(0.0);
        const auto& fvGeometry = this->fvGeometry();
        const auto& elemVolVars = this->elemVolVars();
        const auto& problem = this->problem();

        static const bool enableUnsymmetrizedVelocityGradient
            = getParamFromGroup<bool>(problem.paramGroup(), "FreeFlow.EnableUnsymmetrizedVelocityGradient", false);

        const auto mu = this->problem().effectiveViscosity(this->element(), this->fvGeometry(), this->scvFace());


        // Consider the shear stress caused by the gradient of the velocities normal to our face of interest.
        if (!enableUnsymmetrizedVelocityGradient)
        {
            // if (!scvf.boundary() ||
            //     currentScvfBoundaryTypes->isDirichlet(Indices::velocity(staggeredScvf.directionIndex())) ||
            //     currentScvfBoundaryTypes->isBeaversJoseph(Indices::velocity(staggeredScvf.directionIndex())))
            // {
                const Scalar velocityGrad_ji = VelocityGradients::velocityGradJI(fvGeometry, scvf, elemVolVars);
                // Account for the orientation of the staggered normal face's outer normal vector.
                result -= mu * velocityGrad_ji * scvf.directionSign();
            // }
        }

        // Consider the shear stress caused by the gradient of the velocities parallel to our face of interest.
        const Scalar velocityGrad_ij = VelocityGradients::velocityGradIJ(fvGeometry, scvf, elemVolVars);
        result -= mu * velocityGrad_ij * scvf.directionSign();

        // Account for the area of the staggered lateral face.
        return result * scvf.area() * elemVolVars[scvf.insideScvIdx()].extrusionFactor();
    }

            /*!
     * \brief Returns the frontal part of the momentum flux.
     *        This treats the flux over the staggered face at the center of an element,
     *        parallel to the current scvf where the velocity dof of interest lives.
     *
     * \verbatim
     *                    scvf
     *              ---------=======                 == and # staggered half-control-volume
     *              |       #      | current scv
     *              |       #      |                 # staggered face over which fluxes are calculated
     *   vel.Opp <~~|       #~~>   x~~~~> vel.Self
     *              |       #      |                 x dof position
     *        scvf  |       #      |
     *              --------========                 -- element
     *                   scvf
     * \endverbatim
     */
    NumEqVector pressureContribution() const
    {
        NumEqVector result(0.0);
        const auto& scvf = this->scvFace();
        if (scvf.isLateral() || scvf.boundary())
            return result;

        const auto pressure = this->problem().pressure(this->element(), this->fvGeometry(), this->scvFace());

        // Account for the orientation of the staggered face's normal outer normal vector
        // (pointing in opposite direction of the scvf's one).
        result += pressure * scvf.directionSign();

        // Account for the staggered face's area. For rectangular elements, this equals the area of the scvf
        // our velocity dof of interest lives on.
        return result * scvf.area() * this->elemVolVars()[scvf.insideScvIdx()].extrusionFactor();
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
     *              |       #      |                 # staggered face over which fluxes are calculated
     *   vel.Opp <~~|       #~~>   x~~~~> vel.Self
     *              |       #      |                 x dof position
     *        scvf  |       #      |
     *              --------========                 -- element
     *                   scvf
     * \endverbatim
     */
//     FacePrimaryVariables computeFrontalMomentumFlux(const Problem& problem,
//                                                     const Element& element,
//                                                     const FaceFVElementGeometry& faceFVGeometry,
//                                                     const StaggeredScvf& scvf,
//                                                     const ElementFaceVariables& elemFaceVars,
//                                                     const GridFluxVariablesCache& gridFluxVarsCache)
//     {
//         FacePrimaryVariables frontalFlux(0.0);

//         // The velocities of the dof at interest and the one of the opposite scvf.
//         const auto& scv = faceFVGeometry.scv(scvf.insideScvIdx());
//         const auto& faceVars = elemFaceVars[scv];
//         const auto& context = problem.couplingManager().faceCouplingContext(element, scv);

//         // Advective flux.
//         if (problem.enableInertiaTerms())
//         {
//             // const Scalar velocitySelf = faceVars.velocitySelf();
//             // const Scalar velocityOpposite = faceVars.velocityOpposite();
//             // // Get the average velocity at the center of the element (i.e. the location of the staggered face).
//             // const Scalar transportingVelocity = (velocitySelf + velocityOpposite) * 0.5;
//             // const bool selfIsUpstream = scvf.directionSign() != sign(transportingVelocity);
//             //
//             // StaggeredUpwindHelper<TypeTag, upwindSchemeOrder> upwindHelper(element, fvGeometry, scvf, elemFaceVars, context.elemVolVars, gridFluxVarsCache.staggeredUpwindMethods());
//             // frontalFlux += upwindHelper.computeUpwindFrontalMomentum(selfIsUpstream)
//             //                * transportingVelocity * -1.0 * staggeredScvf.directionSign();
//         }

//         // The volume variables within the current element. We only require those (and none of neighboring elements)
//         // because the fluxes are calculated over the staggered face at the center of the element.
//         const auto eIdx = faceFVGeometry.gridGeometry().elementMapper().index(element);
//         const auto& cellCenterVolVars = context.elemVolVars[eIdx];

//         // Diffusive flux.
//         const Scalar velocityGrad_ii = VelocityGradients::velocityGradII(faceFVGeometry, scvf, faceVars) * scvf.directionSign();

//         static const bool enableUnsymmetrizedVelocityGradient
//             = getParamFromGroup<bool>(problem.paramGroup(), "FreeFlow.EnableUnsymmetrizedVelocityGradient", false);
//         const Scalar factor = enableUnsymmetrizedVelocityGradient ? 1.0 : 2.0;
//         frontalFlux -= factor * cellCenterVolVars.effectiveViscosity() * velocityGrad_ii;

//         // The pressure term.
//         // If specified, the pressure can be normalized using the initial value on the scfv of interest.
//         // The scvf is used to normalize by the same value from the left and right side.
//         // Can potentially help to improve the condition number of the system matrix.
//         const Scalar pressure = normalizePressure ?
//                                 cellCenterVolVars.pressure() - problem.initial(scv)[Indices::pressureIdx] // TODO problem signatures
//                               : cellCenterVolVars.pressure();

//         // Account for the orientation of the staggered face's normal outer normal vector
//         // (pointing in opposite direction of the scvf's one).
//         frontalFlux += pressure * -1.0 * scvf.directionSign();

//         // Account for the staggered face's area. For rectangular elements, this equals the area of the scvf
//         // our velocity dof of interest lives on.
//         return frontalFlux * scvf.area() * cellCenterVolVars.extrusionFactor();
//    }

//     /*!
//      * \brief Returns the momentum flux over the staggered faces
//      *        perpendicular to the scvf where the velocity dof of interest
//      *         lives (coinciding with the element's scvfs).
//      *
//      * \verbatim
//      *                scvf
//      *              ---------#######                 || and # staggered half-control-volume
//      *              |      ||      | current scvf
//      *              |      ||      |                 # normal staggered sub faces over which fluxes are calculated
//      *              |      ||      x~~~~> vel.Self
//      *              |      ||      |                 x dof position
//      *        scvf  |      ||      |
//      *              --------########                -- element
//      *                 scvf
//      * \endverbatim
//      */
//     FacePrimaryVariables computeLateralMomentumFlux(const Problem& problem,
//                                                     const Element& element,
//                                                     const SubControlVolumeFace& scvf,
//                                                     const FVElementGeometry& fvGeometry,
//                                                     const FaceFVElementGeometry& faceFVGeometry,
//                                                     const StaggeredScv& staggeredScv,
//                                                     const ElementVolumeVariables& elemVolVars,
//                                                     const ElementFaceVariables& elemFaceVars,
//                                                     const GridFluxVariablesCache& gridFluxVarsCache)
//     {
//         FacePrimaryVariables lateralFlux(0.0);
//         const auto& faceVars = elemFaceVars[staggeredScv];

//         // If the current scvf is on a boundary, check if there is a Neumann BC for the stress in tangential direction.
//         // Create a boundaryTypes object (will be empty if not at a boundary).
//         std::optional<BoundaryTypes> currentScvfBoundaryTypes;
//         if (staggeredScv.boundary())
//         {
//             // TODO
//             for (const auto& staggeredScvf : scvfs(faceFVGeometry))
//             {
//                 if (staggeredScvf.boundary() && staggeredScvf.isFrontal())
//                     currentScvfBoundaryTypes.emplace(problem.boundaryTypes(element, staggeredScvf));
//             }
//         }

//         // Account for all sub faces.
//         int localSubFaceIdx = 0;
//         for (const auto& staggeredScvf : scvfs(faceFVGeometry))
//         {
//             if (!staggeredScvf.isLateral() || staggeredScvf.insideScvIdx() != staggeredScv.dofIndex())
//                 continue;

//             // Create a boundaryTypes object (will be empty if not at a boundary).
//             std::optional<BoundaryTypes> lateralFaceBoundaryTypes;

//             // Check if there is face/element parallel to our face of interest where the dof lives on. If there is no parallel neighbor,
//             // we are on a boundary where we have to check for boundary conditions.
//             if (staggeredScvf.boundary())
//             {
//                 // Retrieve the boundary types that correspond to the center of the lateral scvf. As a convention, we always query
//                 // the type of BCs at the center of the element's "actual" lateral scvf (not the face of the staggered control volume).
//                 // The value of the BC will be evaluated at the center of the staggered face.
//                 //     --------###T##V                 || frontal face of staggered half-control-volume
//                 //     |      ||      | current scvf    #  lateral staggered face of interest (may lie on a boundary)
//                 //     |      ||      |                 x  dof position
//                 //     |      ||      x~~~~> vel.Self   -- element boundaries
//                 //     |      ||      |                 T  position at which the type of boundary conditions will be evaluated
//                 //     |      ||      |                    (center of lateral scvf)
//                 //     ----------------                 V  position at which the value of the boundary conditions will be evaluated
//                 //                                         (center of the staggered lateral face)
//                 lateralFaceBoundaryTypes.emplace(problem.boundaryTypes(element, staggeredScvf));
//             }

//             // If the current scvf is on a bounary and if there is a Neumann or Beavers-Joseph-(Saffmann) BC for the stress in tangential direction,
//             // assign this value for the lateral momentum flux. No further calculations are required. We assume that all lateral faces
//             // have the same type of BC (Neumann or Beavers-Joseph-(Saffmann)), but we sample the value at their actual positions.
//             if (currentScvfBoundaryTypes)
//             {
//                 // Handle Neumann BCs.
//                 if (currentScvfBoundaryTypes->isNeumann(Indices::velocity(staggeredScvf.directionIndex())))
//                 {
//                     // TODO pass elemFluxVarsCache
//                     lateralFlux += problem.neumann(element, faceFVGeometry, elemFaceVars, staggeredScvf)[Indices::velocity(staggeredScvf.directionIndex())] * staggeredScvf.area() * staggeredScvf.directionSign();
//                     ++localSubFaceIdx;
//                     continue;
//                 }
//             }

//             // Check if the lateral face (perpendicular to our current scvf) lies on a boundary. If yes, boundary conditions might need to be treated
//             // and further calculations can be skipped.
//             if (staggeredScvf.boundary())
//             {
//                 // Check if we have a symmetry boundary condition. If yes, the tangential part of the momentum flux can be neglected
//                 // and we may skip any further calculations for the given sub face.
//                 if (lateralFaceBoundaryTypes->isSymmetry())
//                 {
//                     ++localSubFaceIdx;
//                     continue;
//                 }

//                 // Handle Neumann boundary conditions. No further calculations are then required for the given sub face.
//                 if (lateralFaceBoundaryTypes->isNeumann(Indices::velocity(staggeredScv.directionIndex())))
//                 {
//                     lateralFlux +=  problem.neumann(element, faceFVGeometry, elemFaceVars, staggeredScvf)[Indices::velocity(staggeredScv.directionIndex())] * staggeredScvf.area() * staggeredScvf.directionSign();
//                     ++localSubFaceIdx;
//                     continue;
//                 }

//                 // Handle wall-function fluxes (only required for RANS models)
//                 if (incorporateWallFunction_(lateralFlux, problem, element, fvGeometry, scvf, elemVolVars, elemFaceVars, localSubFaceIdx))
//                 {
//                     ++localSubFaceIdx;
//                     continue;
//                 }
//             }

//             // Check the consistency of the boundary conditions, exactly one of the following must be set
//             if (lateralFaceBoundaryTypes)
//             {
//                 std::bitset<3> admittableBcTypes;
//                 admittableBcTypes.set(0, lateralFaceBoundaryTypes->isDirichlet(Indices::pressureIdx));
//                 admittableBcTypes.set(1, lateralFaceBoundaryTypes->isDirichlet(Indices::velocity(staggeredScv.directionIndex())));
//                 admittableBcTypes.set(2, lateralFaceBoundaryTypes->isBeaversJoseph(Indices::velocity(staggeredScv.directionIndex())));
//                 if (admittableBcTypes.count() != 1)
//                 {
//                     DUNE_THROW(Dune::InvalidStateException, "Invalid boundary conditions for lateral scvf "
//                     "for the momentum equations at global position " << lateralStaggeredFaceCenter_(scvf, localSubFaceIdx)
//                     << ", current scvf global position " << scvf.center());
//                 }
//             }

//             // If none of the above boundary conditions apply for the given sub face, proceed to calculate the tangential momentum flux.
//             // TODO advective fluxes
//             // if (problem.enableInertiaTerms())
//             //     lateralFlux += computeAdvectivePartOfLateralMomentumFlux_(problem, fvGeometry, element,
//             //                                                               scvf, faceFVGeometry, staggeredScvf,
//             //                                                               elemVolVars, faceVars,
//             //                                                               gridFluxVarsCache,
//             //                                                               currentScvfBoundaryTypes, lateralFaceBoundaryTypes,
//             //                                                               localSubFaceIdx);

//             lateralFlux += computeDiffusivePartOfLateralMomentumFlux_(problem, fvGeometry, element,
//                                                                       scvf, faceFVGeometry, staggeredScvf,
//                                                                       elemVolVars, faceVars,
//                                                                       currentScvfBoundaryTypes, lateralFaceBoundaryTypes,
//                                                                       localSubFaceIdx);

//             ++localSubFaceIdx;
//         }
//         return lateralFlux;
//     }

//     /*!
//      * \brief Returns the momentum flux over an inflow or outflow boundary face.
//      *
//      * \verbatim
//      *                    scvf      //
//      *              ---------=======//               == and # staggered half-control-volume
//      *              |      ||      #// current scvf
//      *              |      ||      #//               # staggered boundary face over which fluxes are calculated
//      *              |      ||      x~~~~> vel.Self
//      *              |      ||      #//               x dof position
//      *        scvf  |      ||      #//
//      *              --------========//               -- element
//      *                   scvf       //
//      *                                              // boundary
//      * \endverbatim
//      */
//     FacePrimaryVariables inflowOutflowBoundaryFlux(const Problem& problem,
//                                                    const Element& element,
//                                                    const FaceFVElementGeometry& faceFVGeometry,
//                                                    const StaggeredScvf& scvf,
//                                                    const ElementFaceVariables& elemFaceVars) const
//     {
//         FacePrimaryVariables inOrOutflow(0.0);
//         const auto& scv = faceFVGeometry.scv(scvf.insideScvIdx());
//         const auto& context = problem.couplingManager().faceCouplingContext(element, scv);

//         const auto eIdx = problem.gridGeometry().elementMapper().index(element);
//         const auto& insideVolVars = context.elemVolVars[eIdx];

//         // Advective momentum flux.
//         if (problem.enableInertiaTerms())
//         {
//             // TODO
//             // const Scalar velocitySelf = elemFaceVars[scv].velocitySelf();
//             // const auto& outsideVolVars = elemVolVars[scvf.outsideScvIdx()];
//             // const auto& upVolVars = (scvf.directionSign() == sign(velocitySelf)) ?
//             //                         insideVolVars : outsideVolVars;
//             //
//             // inOrOutflow += velocitySelf * velocitySelf * upVolVars.density();
//         }

//         // Apply a pressure at the boundary.
//         const Scalar boundaryPressure = normalizePressure
//                                         ? (problem.dirichlet(element, scvf)[Indices::pressureIdx] -
//                                            problem.initial(scvf)[Indices::pressureIdx])
//                                         : problem.dirichlet(element, scvf)[Indices::pressureIdx];
//         inOrOutflow += boundaryPressure;

//         // Account for the orientation of the face at the boundary,
//         return inOrOutflow * scvf.directionSign() * scvf.area() * insideVolVars.extrusionFactor();
//     }

// private:

//     /*!
//      * \brief Returns the advective momentum flux over the staggered face perpendicular to the scvf
//      *        where the velocity dof of interest lives (coinciding with the element's scvfs).
//      *
//      * \verbatim
//      *              ----------------
//      *              |              |
//      *              |    transp.   |
//      *              |      vel.    |~~~~> vel.Parallel
//      *              |       ^      |
//      *              |       |      |
//      *       scvf   ---------#######                 || and # staggered half-control-volume
//      *              |      ||      | current scvf
//      *              |      ||      |                 # normal staggered faces over which fluxes are calculated
//      *              |      ||      x~~~~> vel.Self
//      *              |      ||      |                 x dof position
//      *        scvf  |      ||      |
//      *              ---------#######                -- elements
//      *                 scvf
//      * \endverbatim
//      */
//     template<class StaggeredFVElementGeometry>
//     FacePrimaryVariables computeAdvectivePartOfLateralMomentumFlux_(const Problem& problem,
//                                                                     const FVElementGeometry& fvGeometry,
//                                                                     const Element& element,
//                                                                     const SubControlVolumeFace& scvf,
//                                                                     const StaggeredFVElementGeometry staggeredFVGeometry,
//                                                                     const typename StaggeredFVElementGeometry::StaggeredSubControlVolumeFace& staggeredScvf,
//                                                                     const ElementVolumeVariables& elemVolVars,
//                                                                     const ElementFaceVariables& elemFaceVars,
//                                                                     const GridFluxVariablesCache& gridFluxVarsCache,
//                                                                     const std::optional<BoundaryTypes>& currentScvfBoundaryTypes,
//                                                                     const std::optional<BoundaryTypes>& lateralFaceBoundaryTypes,
//                                                                     const int localSubFaceIdx)
//     {
//         const auto eIdx = scvf.insideScvIdx();
//         const auto& lateralFace = fvGeometry.scvf(eIdx, scvf.pairData(localSubFaceIdx).localLateralFaceIdx);

//         // Get the transporting velocity, located at the scvf perpendicular to the current scvf where the dof
//         // of interest is located.
//         const Scalar transportingVelocity = [&]()
//         {
//             const auto& faceVars = elemFaceVars[scvf];
//             if (!scvf.boundary())
//                 return faceVars.velocityLateralInside(localSubFaceIdx);
//             else
//             {
//                 // Create a boundaryTypes object. Get the boundary conditions. We sample the type of BC at the center of the current scvf.
//                 const BoundaryTypes bcTypes = problem.boundaryTypes(element, scvf);

//                 if (bcTypes.isDirichlet(Indices::velocity(lateralFace.directionIndex())))
//                 {
//                     // Construct a temporary scvf which corresponds to the staggered sub face, featuring the location
//                     // the staggered faces's center.
//                     const auto& lateralBoundaryFacePos = lateralStaggeredFaceCenter_(scvf, localSubFaceIdx);
//                     return problem.dirichlet(element, scvf.makeBoundaryFace(lateralBoundaryFacePos))[Indices::velocity(lateralFace.directionIndex())];
//                 }
//                 else if (bcTypes.isBeaversJoseph(Indices::velocity(lateralFace.directionIndex())))
//                 {
//                     return VelocityGradients::beaversJosephVelocityAtCurrentScvf(problem, element, staggeredFVGeometry, staggeredScvf, faceVars,
//                                                                                  currentScvfBoundaryTypes, lateralFaceBoundaryTypes);
//                 }
//                 else
//                     return faceVars.velocityLateralInside(localSubFaceIdx);
//             }
//         }();

//         const bool selfIsUpstream = lateralFace.directionSign() == sign(transportingVelocity);
//         StaggeredUpwindHelper<TypeTag, upwindSchemeOrder> upwindHelper(element, fvGeometry, scvf, elemFaceVars, elemVolVars, gridFluxVarsCache.staggeredUpwindMethods());
//         return upwindHelper.computeUpwindLateralMomentum(selfIsUpstream, lateralFace, localSubFaceIdx, currentScvfBoundaryTypes, lateralFaceBoundaryTypes)
//                * transportingVelocity * lateralFace.directionSign() * lateralFace.area() * 0.5 * extrusionFactor_(elemVolVars, lateralFace);
//     }

//     /*!
//      * \brief Returns the diffusive momentum flux over the staggered face perpendicular to the scvf
//      *        where the velocity dof of interest lives (coinciding with the element's scvfs).
//      *
//      * \verbatim
//      *              ----------------
//      *              |              |vel.
//      *              |    in.norm.  |Parallel
//      *              |       vel.   |~~~~>
//      *              |       ^      |        ^ out.norm.vel.
//      *              |       |      |        |
//      *       scvf   ---------#######:::::::::       || and # staggered half-control-volume (own element)
//      *              |      ||      | curr. ::
//      *              |      ||      | scvf  ::       :: staggered half-control-volume (neighbor element)
//      *              |      ||      x~~~~>  ::
//      *              |      ||      | vel.  ::       # lateral staggered faces over which fluxes are calculated
//      *        scvf  |      ||      | Self  ::
//      *              ---------#######:::::::::       x dof position
//      *                 scvf
//      *                                              -- elements
//      * \endverbatim
//      */
//     template<class StaggeredFVElementGeometry>
//     FacePrimaryVariables computeDiffusivePartOfLateralMomentumFlux_(const Problem& problem,
//                                                                     const FVElementGeometry& fvGeometry,
//                                                                     const Element& element,
//                                                                     const SubControlVolumeFace& scvf,
//                                                                     const StaggeredFVElementGeometry staggeredFVGeometry,
//                                                                     const typename StaggeredFVElementGeometry::StaggeredSubControlVolumeFace& staggeredScvf,
//                                                                     const ElementVolumeVariables& elemVolVars,
//                                                                     const FaceVariables& faceVars,
//                                                                     const std::optional<BoundaryTypes>& currentScvfBoundaryTypes,
//                                                                     const std::optional<BoundaryTypes>& lateralFaceBoundaryTypes,
//                                                                     const int localSubFaceIdx)
//     {
//         assert(staggeredScvf.isLateral());

//         const auto eIdx = scvf.insideScvIdx();
//         const auto& lateralFace = fvGeometry.scvf(eIdx, scvf.pairData(localSubFaceIdx).localLateralFaceIdx);

//         FacePrimaryVariables lateralDiffusiveFlux(0.0);

//         static const bool enableUnsymmetrizedVelocityGradient
//             = getParamFromGroup<bool>(problem.paramGroup(), "FreeFlow.EnableUnsymmetrizedVelocityGradient", false);

//         // Get the volume variables of the own and the neighboring element. The neighboring
//         // element is adjacent to the staggered face normal to the current scvf
//         // where the dof of interest is located.
//         // TODO use coupling context to get interpolated mu
//         const auto& insideVolVars = elemVolVars[lateralFace.insideScvIdx()];
//         const auto& outsideVolVars = elemVolVars[lateralFace.outsideScvIdx()];

//         // Get the averaged viscosity at the staggered face normal to the current scvf.
//         const Scalar muAvg = lateralFace.boundary()
//                              ? insideVolVars.effectiveViscosity()
//                              : (insideVolVars.effectiveViscosity() + outsideVolVars.effectiveViscosity()) * 0.5;

//         // Consider the shear stress caused by the gradient of the velocities normal to our face of interest.
//         if (!enableUnsymmetrizedVelocityGradient)
//         {
//             if (!scvf.boundary() ||
//                 currentScvfBoundaryTypes->isDirichlet(Indices::velocity(staggeredScvf.directionIndex())) ||
//                 currentScvfBoundaryTypes->isBeaversJoseph(Indices::velocity(staggeredScvf.directionIndex())))
//             {
//                 const Scalar velocityGrad_ji = VelocityGradients::velocityGradJI(problem, element, staggeredFVGeometry, staggeredScvf, faceVars, currentScvfBoundaryTypes, lateralFaceBoundaryTypes);
//                 // Account for the orientation of the staggered normal face's outer normal vector.
//                 lateralDiffusiveFlux -= muAvg * velocityGrad_ji * staggeredScvf.directionSign();
//             }
//         }

//         // Consider the shear stress caused by the gradient of the velocities parallel to our face of interest.
//         const Scalar velocityGrad_ij = VelocityGradients::velocityGradIJ(problem, element, staggeredFVGeometry, staggeredScvf, faceVars, currentScvfBoundaryTypes, lateralFaceBoundaryTypes);
//         lateralDiffusiveFlux -= muAvg * velocityGrad_ij * staggeredScvf.directionSign();

//         // Account for the area of the staggered lateral face (0.5 of the coinciding scfv).
//         return lateralDiffusiveFlux * lateralFace.area() * 0.5 * extrusionFactor_(elemVolVars, lateralFace);
//     }

//     /*!
//      * \brief Get the location of the lateral staggered face's center.
//      *        Only needed for boundary conditions if the current scvf or the lateral one is on a bounary.
//      *
//      * \verbatim
//      *      --------#######o                 || frontal face of staggered half-control-volume
//      *      |      ||      | current scvf    #  lateral staggered face of interest (may lie on a boundary)
//      *      |      ||      |                 x  dof position
//      *      |      ||      x~~~~> vel.Self   -- element boundaries, current scvf may lie on a boundary
//      *      |      ||      |                 o  position at which the boundary conditions will be evaluated
//      *      |      ||      |                    (lateralStaggeredFaceCenter)
//      *      ----------------
//      * \endverbatim
//      */
//     const GlobalPosition& lateralStaggeredFaceCenter_(const SubControlVolumeFace& scvf, const int localSubFaceIdx) const
//     {
//         return scvf.pairData(localSubFaceIdx).lateralStaggeredFaceCenter;
//     };

//     //! helper function to get the averaged extrusion factor for a face
//     static Scalar extrusionFactor_(const ElementVolumeVariables& elemVolVars, const SubControlVolumeFace& scvf)
//     {
//         const auto& insideVolVars = elemVolVars[scvf.insideScvIdx()];
//         const auto& outsideVolVars = elemVolVars[scvf.outsideScvIdx()];
//         return harmonicMean(insideVolVars.extrusionFactor(), outsideVolVars.extrusionFactor());
//     }

//     //! do nothing if no turbulence model is used
//     template<class ...Args, bool turbulenceModel = ModelTraits::usesTurbulenceModel(), std::enable_if_t<!turbulenceModel, int> = 0>
//     bool incorporateWallFunction_(Args&&... args) const
//     { return false; }

//     //! if a turbulence model is used, ask the problem is a wall function shall be employed and get the flux accordingly
//     template<bool turbulenceModel = ModelTraits::usesTurbulenceModel(), std::enable_if_t<turbulenceModel, int> = 0>
//     bool incorporateWallFunction_(FacePrimaryVariables& lateralFlux,
//                                   const Problem& problem,
//                                   const Element& element,
//                                   const FVElementGeometry& fvGeometry,
//                                   const SubControlVolumeFace& scvf,
//                                   const ElementVolumeVariables& elemVolVars,
//                                   const ElementFaceVariables& elemFaceVars,
//                                   const std::size_t localSubFaceIdx) const
//     {
//         const auto eIdx = scvf.insideScvIdx();
//         const auto& lateralScvf = fvGeometry.scvf(eIdx, scvf.pairData(localSubFaceIdx).localLateralFaceIdx);

//         if (problem.useWallFunction(element, lateralScvf, Indices::velocity(scvf.directionIndex())))
//         {
//             const auto& lateralStaggeredFaceCenter = lateralStaggeredFaceCenter_(scvf, localSubFaceIdx);
//             const auto lateralBoundaryFace = lateralScvf.makeBoundaryFace(lateralStaggeredFaceCenter);
//             lateralFlux += problem.wallFunction(element, fvGeometry, elemVolVars, elemFaceVars, scvf, lateralBoundaryFace)[Indices::velocity(scvf.directionIndex())]
//                                                * extrusionFactor_(elemVolVars, lateralScvf) * lateralScvf.area() * 0.5;
//             return true;
//         }
//         else
//             return false;
//     }
};

} // end namespace Dumux

#endif // DUMUX_NAVIERSTOKES_STAGGERED_FLUXVARIABLES_HH
