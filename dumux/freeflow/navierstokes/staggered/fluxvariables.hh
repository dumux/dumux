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
#ifndef DUMUX_NAVIERSTOKES_STAGGERED_FLUXVARIABLES_HH
#define DUMUX_NAVIERSTOKES_STAGGERED_FLUXVARIABLES_HH

#include <array>

#include <dumux/common/math.hh>
#include <dumux/common/exceptions.hh>
#include <dumux/common/parameters.hh>
#include <dumux/common/properties.hh>

#include <dumux/flux/fluxvariablesbase.hh>
#include <dumux/discretization/method.hh>

#include "staggeredupwindfluxvariables.hh"

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
: public FluxVariablesBase<GetPropType<TypeTag, Properties::Problem>,
                           typename GetPropType<TypeTag, Properties::FVGridGeometry>::LocalView,
                           typename GetPropType<TypeTag, Properties::GridVolumeVariables>::LocalView,
                           typename GetPropType<TypeTag, Properties::GridFluxVariablesCache>::LocalView>
{
    using GridVariables = GetPropType<TypeTag, Properties::GridVariables>;

    using GridVolumeVariables = typename GridVariables::GridVolumeVariables;
    using ElementVolumeVariables = typename GridVolumeVariables::LocalView;
    using VolumeVariables = typename GridVolumeVariables::VolumeVariables;

    using GridFluxVariablesCache = typename GridVariables::GridFluxVariablesCache;
    using FluxVariablesCache = typename GridFluxVariablesCache::FluxVariablesCache;

    using GridFaceVariables = typename GridVariables::GridFaceVariables;
    using ElementFaceVariables = typename GridFaceVariables::LocalView;
    using FaceVariables = typename GridFaceVariables::FaceVariables;

    using Problem = GetPropType<TypeTag, Properties::Problem>;
    using FVGridGeometry = GetPropType<TypeTag, Properties::FVGridGeometry>;
    using FVElementGeometry = typename FVGridGeometry::LocalView;
    using GridView = typename FVGridGeometry::GridView;
    using Element = typename GridView::template Codim<0>::Entity;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using ModelTraits = GetPropType<TypeTag, Properties::ModelTraits>;
    using Indices = typename ModelTraits::Indices;
    using SubControlVolumeFace = typename FVElementGeometry::SubControlVolumeFace;
    using CellCenterPrimaryVariables = GetPropType<TypeTag, Properties::CellCenterPrimaryVariables>;
    using FacePrimaryVariables = GetPropType<TypeTag, Properties::FacePrimaryVariables>;
    using BoundaryTypes = GetPropType<TypeTag, Properties::BoundaryTypes>;

    static constexpr bool normalizePressure = getPropValue<TypeTag, Properties::NormalizePressure>();

    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;

    static constexpr auto cellCenterIdx = FVGridGeometry::cellCenterIdx();
    static constexpr auto faceIdx = FVGridGeometry::faceIdx();

    static constexpr int upwindSchemeOrder = getPropValue<TypeTag, Properties::UpwindSchemeOrder>();
    static constexpr bool useHigherOrder = upwindSchemeOrder > 1;

public:

    using HeatConductionType = GetPropType<TypeTag, Properties::HeatConductionType>;

    /*!
     * \brief Returns the advective flux over a sub control volume face.
     * \param problem The object specifying the problem which ought to be simulated
     * \param elemVolVars All volume variables for the element
     * \param elemFaceVars The face variables
     * \param scvf The sub control volume face
     * \param upwindTerm The uwind term (i.e. the advectively transported quantity)
     */
    template<class UpwindTerm>
    static Scalar advectiveFluxForCellCenter(const Problem& problem,
                                             const ElementVolumeVariables& elemVolVars,
                                             const ElementFaceVariables& elemFaceVars,
                                             const SubControlVolumeFace &scvf,
                                             UpwindTerm upwindTerm)
    {
        const Scalar velocity = elemFaceVars[scvf].velocitySelf();
        const bool insideIsUpstream = scvf.directionSign() == sign(velocity);
        static const Scalar upwindWeight = getParamFromGroup<Scalar>(problem.paramGroup(), "Flux.UpwindWeight");

        const auto& insideVolVars = elemVolVars[scvf.insideScvIdx()];
        const auto& outsideVolVars = elemVolVars[scvf.outsideScvIdx()];

        const auto& upstreamVolVars = insideIsUpstream ? insideVolVars : outsideVolVars;
        const auto& downstreamVolVars = insideIsUpstream ? outsideVolVars : insideVolVars;

        const Scalar flux = (upwindWeight * upwindTerm(upstreamVolVars) +
                            (1.0 - upwindWeight) * upwindTerm(downstreamVolVars))
                            * velocity * scvf.area() * scvf.directionSign();

        return flux * extrusionFactor_(elemVolVars, scvf);
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
    CellCenterPrimaryVariables computeMassFlux(const Problem& problem,
                                               const Element& element,
                                               const FVElementGeometry& fvGeometry,
                                               const ElementVolumeVariables& elemVolVars,
                                               const ElementFaceVariables& elemFaceVars,
                                               const SubControlVolumeFace& scvf,
                                               const FluxVariablesCache& fluxVarsCache)
    {
        // The advectively transported quantity (i.e density for a single-phase system).
        auto upwindTerm = [](const auto& volVars) { return volVars.density(); };

        // Call the generic flux function.
        CellCenterPrimaryVariables result(0.0);
        result[Indices::conti0EqIdx - ModelTraits::dim()] = advectiveFluxForCellCenter(problem, elemVolVars, elemFaceVars, scvf, upwindTerm);

        return result;
    }

    /*!
     * \brief Returns the momentum flux over all staggered faces.
     */
    FacePrimaryVariables computeMomentumFlux(const Problem& problem,
                                             const Element& element,
                                             const SubControlVolumeFace& scvf,
                                             const FVElementGeometry& fvGeometry,
                                             const ElementVolumeVariables& elemVolVars,
                                             const ElementFaceVariables& elemFaceVars,
                                             const GridFluxVariablesCache& gridFluxVarsCache)
    {
        return computeFrontalMomentumFlux(problem, element, scvf, fvGeometry, elemVolVars, elemFaceVars, gridFluxVarsCache) +
               computeLateralMomentumFlux(problem, element, scvf, fvGeometry, elemVolVars, elemFaceVars, gridFluxVarsCache);
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
    FacePrimaryVariables computeFrontalMomentumFlux(const Problem& problem,
                                                    const Element& element,
                                                    const SubControlVolumeFace& scvf,
                                                    const FVElementGeometry& fvGeometry,
                                                    const ElementVolumeVariables& elemVolVars,
                                                    const ElementFaceVariables& elemFaceVars,
                                                    const GridFluxVariablesCache& gridFluxVarsCache)
    {
        FacePrimaryVariables frontalFlux(0.0);

        // The velocities of the dof at interest and the one of the opposite scvf.
        const Scalar velocitySelf = elemFaceVars[scvf].velocitySelf();
        const Scalar velocityOpposite = elemFaceVars[scvf].velocityOpposite();

        // Advective flux.
        if (problem.enableInertiaTerms())
        {
            // Get the average velocity at the center of the element (i.e. the location of the staggered face).
            const Scalar transportingVelocity = (velocitySelf + velocityOpposite) * 0.5;

            frontalFlux += StaggeredUpwindFluxVariables<TypeTag, upwindSchemeOrder>::computeUpwindedFrontalMomentum(scvf, elemFaceVars, elemVolVars, gridFluxVarsCache, transportingVelocity)
                           * transportingVelocity * -1.0 * scvf.directionSign();
        }

        // The volume variables within the current element. We only require those (and none of neighboring elements)
        // because the fluxes are calculated over the staggered face at the center of the element.
        const auto& insideVolVars = elemVolVars[scvf.insideScvIdx()];

        // Diffusive flux.
        // The velocity gradient already accounts for the orientation
        // of the staggered face's outer normal vector.
        const Scalar gradV = (velocityOpposite - velocitySelf) / scvf.selfToOppositeDistance();

        static const bool enableUnsymmetrizedVelocityGradient
            = getParamFromGroup<bool>(problem.paramGroup(), "FreeFlow.EnableUnsymmetrizedVelocityGradient", false);
        const Scalar factor = enableUnsymmetrizedVelocityGradient ? 1.0 : 2.0;
        frontalFlux -= factor * insideVolVars.effectiveViscosity() * gradV;

        // The pressure term.
        // If specified, the pressure can be normalized using the initial value on the scfv of interest.
        // The scvf is used to normalize by the same value from the left and right side.
        // Can potentially help to improve the condition number of the system matrix.
        const Scalar pressure = normalizePressure ?
                                insideVolVars.pressure() - problem.initial(scvf)[Indices::pressureIdx]
                              : insideVolVars.pressure();

        // Account for the orientation of the staggered face's normal outer normal vector
        // (pointing in opposite direction of the scvf's one).
        frontalFlux += pressure * -1.0 * scvf.directionSign();

        // Handle inflow or outflow conditions.
        // Treat the staggered half-volume adjacent to the boundary as if it was on the opposite side of the boundary.
        // The respective face's outer normal vector will point in the same direction as the scvf's one.
        if(scvf.boundary() && problem.boundaryTypes(element, scvf).isDirichlet(Indices::pressureIdx))
            frontalFlux += inflowOutflowBoundaryFlux_(problem, element, scvf, elemVolVars, elemFaceVars);

        // Account for the staggered face's area. For rectangular elements, this equals the area of the scvf
        // our velocity dof of interest lives on.
        return frontalFlux * scvf.area() * insideVolVars.extrusionFactor();
   }

    /*!
     * \brief Returns the momentum flux over the staggered faces
     *        perpendicular to the scvf where the velocity dof of interest
     *         lives (coinciding with the element's scvfs).
     *
     * \verbatim
     *                scvf
     *              ---------#######                 || and # staggered half-control-volume
     *              |      ||      | current scvf
     *              |      ||      |                 # normal staggered sub faces over which fluxes are calculated
     *              |      ||      x~~~~> vel.Self
     *              |      ||      |                 x dof position
     *        scvf  |      ||      |
     *              --------########                -- element
     *                 scvf
     * \endverbatim
     */
    FacePrimaryVariables computeLateralMomentumFlux(const Problem& problem,
                                                    const Element& element,
                                                    const SubControlVolumeFace& scvf,
                                                    const FVElementGeometry& fvGeometry,
                                                    const ElementVolumeVariables& elemVolVars,
                                                    const ElementFaceVariables& elemFaceVars,
                                                    const GridFluxVariablesCache& gridFluxVarsCache)
    {
        FacePrimaryVariables normalFlux(0.0);
        const auto& faceVars = elemFaceVars[scvf];
        const int numSubFaces = scvf.pairData().size();

        // Account for all sub faces.
        for(int localSubFaceIdx = 0; localSubFaceIdx < numSubFaces; ++localSubFaceIdx)
        {
            const auto eIdx = scvf.insideScvIdx();
            // Get the face normal to the face the dof lives on. The staggered sub face conincides with half of this normal face.
            const auto& normalFace = fvGeometry.scvf(eIdx, scvf.pairData(localSubFaceIdx).localNormalFaceIdx);

            bool lateralFaceHasDirichletPressure = false; // check for Dirichlet boundary condition for the pressure
            bool lateralFaceHasBJS = false; // check for Beavers-Joseph-Saffman boundary condition

            // Check if there is face/element parallel to our face of interest where the dof lives on. If there is no parallel neighbor,
            // we are on a boundary where we have to check for boundary conditions.
            if(!scvf.hasParallelNeighbor(localSubFaceIdx,0))
            {
                // Construct a temporary scvf which corresponds to the staggered sub face, featuring the location
                // the sub faces's center.
                auto localSubFaceCenter = scvf.pairData(localSubFaceIdx).virtualFirstParallelFaceDofPos + normalFace.center();
                localSubFaceCenter *= 0.5;

                //    ________________
                //    --------####o###                 || frontal face of staggered half-control-volume
                //    |      ||      | current scvf    #  localSubFace of interest, lies on lateral boundary
                //    |      ||      |                 x  dof position
                //    |      ||      x~~~~> vel.Self   -- element boundaries
                //    |      ||      |                 __ domain boundary
                //    |      ||      |                 o  position at which the boundary conditions will be evaluated
                //    ----------------                    (localSubFaceCenter)

                const auto localSubFace = normalFace.makeBoundaryFace(localSubFaceCenter);

                // Retrieve the boundary types that correspond to the sub face.
                const auto bcTypes = problem.boundaryTypes(element, localSubFace);

                // Check if we have a symmetry boundary condition. If yes, the tangential part of the momentum flux can be neglected
                // and we may skip any further calculations for the given sub face.
                if(bcTypes.isSymmetry())
                    continue;

                // Handle Neumann boundary conditions. No further calculations are then required for the given sub face.
                if(bcTypes.isNeumann(Indices::velocity(scvf.directionIndex())))
                {
                    normalFlux += problem.neumann(element, fvGeometry, elemVolVars, elemFaceVars, localSubFace)[Indices::velocity(scvf.directionIndex())]
                                                  * elemVolVars[normalFace.insideScvIdx()].extrusionFactor() * normalFace.area() * 0.5;
                    continue;
                }

                // Handle wall-function fluxes (only required for RANS models)
                if (incorporateWallFunction_(normalFlux, problem, element, fvGeometry, scvf, normalFace, localSubFace, elemVolVars, elemFaceVars))
                    continue;

                // Check if we have a Beavers-Joseph-Saffman condition or a Dirichlet condition for the velocity or a Dirichlet condition for the pressure.
                // Then the parallel velocity at the boundary is calculated accordingly for the advective part and the diffusive part of the normal momentum flux.
                if (bcTypes.isDirichlet(Indices::pressureIdx))
                    lateralFaceHasDirichletPressure = true;
                else if (bcTypes.isBJS(Indices::velocity(scvf.directionIndex())))
                    lateralFaceHasBJS = true;
                else if (bcTypes.isDirichlet(Indices::velocity(scvf.directionIndex())) == false)
                    DUNE_THROW(Dune::InvalidStateException,  "Something went wrong with the boundary conditions "
                           "for the momentum equations at global position " << localSubFaceCenter);
            }

            // If there is no symmetry or Neumann boundary condition for the given sub face, proceed to calculate the tangential momentum flux.
            if (problem.enableInertiaTerms())
                normalFlux += computeAdvectivePartOfLateralMomentumFlux_(problem, fvGeometry, element,
                                                                         scvf, normalFace, elemVolVars, faceVars,
                                                                         gridFluxVarsCache, localSubFaceIdx,
                                                                         lateralFaceHasDirichletPressure, lateralFaceHasBJS);

            normalFlux += computeDiffusivePartOfLateralMomentumFlux_(problem, fvGeometry, element,
                                                                     scvf, normalFace, elemVolVars, faceVars,
                                                                     localSubFaceIdx,
                                                                     lateralFaceHasDirichletPressure, lateralFaceHasBJS);
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
     *              |      ||      |                 # normal staggered faces over which fluxes are calculated
     *              |      ||      x~~~~> vel.Self
     *              |      ||      |                 x dof position
     *        scvf  |      ||      |
     *              ---------#######                -- elements
     *                 scvf
     * \endverbatim
     */
    FacePrimaryVariables computeAdvectivePartOfLateralMomentumFlux_(const Problem& problem,
                                                                    const FVElementGeometry& fvGeometry,
                                                                    const Element& element,
                                                                    const SubControlVolumeFace& scvf,
                                                                    const SubControlVolumeFace& normalFace,
                                                                    const ElementVolumeVariables& elemVolVars,
                                                                    const FaceVariables& faceVars,
                                                                    const GridFluxVariablesCache& gridFluxVarsCache,
                                                                    const int localSubFaceIdx,
                                                                    const bool lateralFaceHasDirichletPressure,
                                                                    const bool lateralFaceHasBJS)
    {
        // Get the transporting velocity, located at the scvf perpendicular to the current scvf where the dof
        // of interest is located.
        const Scalar transportingVelocity = faceVars.velocityNormalInside(localSubFaceIdx);

        return StaggeredUpwindFluxVariables<TypeTag, upwindSchemeOrder>::computeUpwindedLateralMomentum(problem, fvGeometry, element, scvf, normalFace, elemVolVars, faceVars,
                                                                     gridFluxVarsCache, localSubFaceIdx, lateralFaceHasDirichletPressure,
                                                                     lateralFaceHasBJS)
               * transportingVelocity * normalFace.directionSign() * normalFace.area() * 0.5 * extrusionFactor_(elemVolVars, normalFace);
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
     *              |      ||      | vel.  ::       # normal staggered faces over which fluxes are calculated
     *        scvf  |      ||      | Self  ::
     *              ---------#######:::::::::       x dof position
     *                 scvf
     *                                              -- elements
     * \endverbatim
     */
    FacePrimaryVariables computeDiffusivePartOfLateralMomentumFlux_(const Problem& problem,
                                                                    const FVElementGeometry& fvGeometry,
                                                                    const Element& element,
                                                                    const SubControlVolumeFace& scvf,
                                                                    const SubControlVolumeFace& normalFace,
                                                                    const ElementVolumeVariables& elemVolVars,
                                                                    const FaceVariables& faceVars,
                                                                    const int localSubFaceIdx,
                                                                    const bool lateralFaceHasDirichletPressure,
                                                                    const bool lateralFaceHasBJS)
    {
        FacePrimaryVariables normalDiffusiveFlux(0.0);

        static const bool enableUnsymmetrizedVelocityGradient
            = getParamFromGroup<bool>(problem.paramGroup(), "FreeFlow.EnableUnsymmetrizedVelocityGradient", false);

        // Get the volume variables of the own and the neighboring element. The neighboring
        // element is adjacent to the staggered face normal to the current scvf
        // where the dof of interest is located.
        const auto& insideVolVars = elemVolVars[normalFace.insideScvIdx()];
        const auto& outsideVolVars = elemVolVars[normalFace.outsideScvIdx()];

        // Get the averaged viscosity at the staggered face normal to the current scvf.
        const Scalar muAvg = normalFace.boundary()
                             ? insideVolVars.effectiveViscosity()
                             : (insideVolVars.effectiveViscosity() + outsideVolVars.effectiveViscosity()) * 0.5;

        if (!enableUnsymmetrizedVelocityGradient)
        {
            // If we are at a boundary, a gradient of zero is implictly assumed for all velocities,
            // thus no further calculations are required.
            if (!scvf.boundary())
            {
                // For the normal gradient, get the velocities perpendicular to the velocity at the current scvf.
                // The inner one is located at staggered face within the own element,
                // the outer one at the respective staggered face of the element on the other side of the
                // current scvf.
                const Scalar innerNormalVelocity = faceVars.velocityNormalInside(localSubFaceIdx);
                const Scalar outerNormalVelocity = faceVars.velocityNormalOutside(localSubFaceIdx);

                // Calculate the velocity gradient in positive coordinate direction.
                const Scalar normalDeltaV = scvf.normalInPosCoordDir()
                                            ? (outerNormalVelocity - innerNormalVelocity)
                                            : (innerNormalVelocity - outerNormalVelocity);

                const Scalar normalGradient = normalDeltaV / scvf.pairData(localSubFaceIdx).normalDistance;

                // Account for the orientation of the staggered normal face's outer normal vector.
                normalDiffusiveFlux -= muAvg * normalGradient * normalFace.directionSign();
            }
        }

        // If we have a Dirichlet condition for the pressure we assume to have zero parallel gradient
        // so we can skip the computation.
        if (!lateralFaceHasDirichletPressure)
        {
            // For the parallel derivative, get the velocities at the current (own) scvf
            // and at the parallel one at the neighboring scvf.
            const Scalar innerParallelVelocity = faceVars.velocitySelf();

            const Scalar velocityFirstParallel = scvf.hasParallelNeighbor(localSubFaceIdx,0)
                                               ? faceVars.velocityParallel(localSubFaceIdx,0)
                                               : getParallelVelocityFromBoundary_(problem, scvf, normalFace,
                                                                                  innerParallelVelocity, localSubFaceIdx,
                                                                                  element, false, lateralFaceHasBJS);

            // The velocity gradient already accounts for the orientation
            // of the staggered face's outer normal vector.
            const Scalar parallelGradient = (velocityFirstParallel - innerParallelVelocity)
                                          / scvf.cellCenteredParallelDistance(localSubFaceIdx,0);

            normalDiffusiveFlux -= muAvg * parallelGradient;
        }

        // Account for the area of the staggered normal face (0.5 of the coinciding scfv).
        return normalDiffusiveFlux * normalFace.area() * 0.5 * extrusionFactor_(elemVolVars, normalFace);
    }

    /*!
     * \brief Returns the momentum flux over an inflow or outflow boundary face.
     *
     * \verbatim
     *                    scvf      //
     *              ---------=======//               == and # staggered half-control-volume
     *              |      ||      #// current scvf
     *              |      ||      #//               # staggered boundary face over which fluxes are calculated
     *              |      ||      x~~~~> vel.Self
     *              |      ||      #//               x dof position
     *        scvf  |      ||      #//
     *              --------========//               -- element
     *                   scvf       //
     *                                              // boundary
     * \endverbatim
     */
    FacePrimaryVariables inflowOutflowBoundaryFlux_(const Problem& problem,
                                                    const Element& element,
                                                    const SubControlVolumeFace& scvf,
                                                    const ElementVolumeVariables& elemVolVars,
                                                    const ElementFaceVariables& elemFaceVars)
    {
        FacePrimaryVariables inOrOutflow(0.0);

        // Advective momentum flux.
        if (problem.enableInertiaTerms())
        {
            const Scalar velocitySelf = elemFaceVars[scvf].velocitySelf();
            const auto& insideVolVars = elemVolVars[scvf.insideScvIdx()];
            const auto& outsideVolVars = elemVolVars[scvf.outsideScvIdx()];
            const auto& upVolVars = (scvf.directionSign() == sign(velocitySelf)) ?
                                    insideVolVars : outsideVolVars;

            inOrOutflow += velocitySelf * velocitySelf * upVolVars.density();
        }

        // Apply a pressure at the boundary.
        const Scalar boundaryPressure = normalizePressure
                                        ? (problem.dirichlet(element, scvf)[Indices::pressureIdx] -
                                           problem.initial(scvf)[Indices::pressureIdx])
                                        : problem.dirichlet(element, scvf)[Indices::pressureIdx];
        inOrOutflow += boundaryPressure;

        // Account for the orientation of the face at the boundary,
        return inOrOutflow * scvf.directionSign();
    }

private:

    //! helper function to conveniently create a ghost face which is outside the domain, parallel to the scvf of interest
    SubControlVolumeFace makeParallelGhostFace_(const SubControlVolumeFace& ownScvf, const int localSubFaceIdx) const
    {
        return ownScvf.makeBoundaryFace(ownScvf.pairData(localSubFaceIdx).virtualFirstParallelFaceDofPos);
    };

    //! helper function to get the averaged extrusion factor for a face
    static Scalar extrusionFactor_(const ElementVolumeVariables& elemVolVars, const SubControlVolumeFace& scvf)
    {
        const auto& insideVolVars = elemVolVars[scvf.insideScvIdx()];
        const auto& outsideVolVars = elemVolVars[scvf.outsideScvIdx()];
        return harmonicMean(insideVolVars.extrusionFactor(), outsideVolVars.extrusionFactor());
    }
    //! do nothing if no turbulence model is used
    template<class ...Args, bool turbulenceModel = ModelTraits::usesTurbulenceModel(), std::enable_if_t<!turbulenceModel, int> = 0>
    bool incorporateWallFunction_(Args&&... args) const
    { return false; };

    //! if a turbulence model is used, ask the problem is a wall function shall be employed and get the flux accordingly
    template<bool turbulenceModel = ModelTraits::usesTurbulenceModel(), std::enable_if_t<turbulenceModel, int> = 0>
    bool incorporateWallFunction_(FacePrimaryVariables& normalFlux,
                                  const Problem& problem,
                                  const Element& element,
                                  const FVElementGeometry& fvGeometry,
                                  const SubControlVolumeFace& scvf,
                                  const SubControlVolumeFace& normalFace,
                                  const SubControlVolumeFace& localSubFace,
                                  const ElementVolumeVariables& elemVolVars,
                                  const ElementFaceVariables& elemFaceVars) const
    {
        if (problem.useWallFunction(element, localSubFace, Indices::velocity(scvf.directionIndex())))
        {
            normalFlux += problem.wallFunction(element, fvGeometry, elemVolVars, elemFaceVars, scvf, localSubFace)[Indices::velocity(scvf.directionIndex())]
                                               * elemVolVars[normalFace.insideScvIdx()].extrusionFactor() * normalFace.area() * 0.5;
            return true;
        }
        else
            return false;
    }

    /*!
     * \brief Return the outer parallel velocity for normal faces that are on the boundary and therefore have no neighbor.
     *
     * Calls the problem to retrieve a fixed value set on the boundary.
     *
     * \param problem The problem
     * \param scvf The SubControlVolumeFace that is normal to the boundary
     * \param normalFace The face at the boundary
     * \param velocitySelf the velocity at scvf
     * \param localSubFaceIdx The local index of the face that is on the boundary
     * \param element The element that is on the boundary
     * \param lateralFaceHasDirichletPressure @c true if there is a dirichlet condition for the pressure on the boundary
     * \param lateralFaceHasBJS @c true if there is a BJS condition fot the velocity on the boundary
     */
    Scalar getParallelVelocityFromBoundary_(const Problem& problem,
                                            const SubControlVolumeFace& scvf,
                                            const SubControlVolumeFace& normalFace,
                                            const Scalar velocitySelf,
                                            const int localSubFaceIdx,
                                            const Element& element,
                                            const bool lateralFaceHasDirichletPressure,
                                            const bool lateralFaceHasBJS) const
    {
        // If there is a Dirichlet condition for the pressure we assume zero gradient for the velocity,
        // so the velocity at the boundary equal to that on the scvf.
        if (lateralFaceHasDirichletPressure)
            return velocitySelf;

        //    ________________
        //    --------#######o                 || frontal face of staggered half-control-volume
        //    |      ||      | current scvf    #  localSubFace of interest, lies on lateral boundary
        //    |      ||      |                 x  dof position
        //    |      ||      x~~~~> vel.Self   -- element boundaries
        //    |      ||      |                 __ domain boundary
        //    |      ||      |                 o  position at which the boundary conditions will be evaluated
        //    ----------------

        const auto ghostFace = makeParallelGhostFace_(scvf, localSubFaceIdx);
        if (lateralFaceHasBJS)
            return problem.bjsVelocity(element, scvf, normalFace, localSubFaceIdx, velocitySelf);
        return problem.dirichlet(element, ghostFace)[Indices::velocity(scvf.directionIndex())];
    }
};

} // end namespace Dumux

#endif // DUMUX_NAVIERSTOKES_STAGGERED_FLUXVARIABLES_HH
