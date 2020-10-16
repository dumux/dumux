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
#include <optional>
#include <type_traits>

#include <dumux/common/math.hh>
#include <dumux/common/exceptions.hh>
#include <dumux/common/parameters.hh>
#include <dumux/common/properties.hh>
#include <dumux/common/typetraits/problem.hh>

#include <dumux/flux/fluxvariablesbase.hh>
#include <dumux/discretization/method.hh>
#include <dumux/discretization/extrusion.hh>

#include "staggeredupwindhelper.hh"
#include "velocitygradients.hh"

namespace Dumux {

/*!
 * \ingroup NavierStokesModel
 * \brief The flux variables class for the Navier-Stokes model using the staggered grid discretization.
 */

// forward declaration
template<class TypeTag, DiscretizationMethod discMethod>
class NavierStokesFluxVariablesImpl;

template<class TypeTag>
class NavierStokesFluxVariablesImpl<TypeTag, DiscretizationMethod::staggered>
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

    using GridFaceVariables = typename GridVariables::GridFaceVariables;
    using ElementFaceVariables = typename GridFaceVariables::LocalView;
    using FaceVariables = typename GridFaceVariables::FaceVariables;

    using Problem = GetPropType<TypeTag, Properties::Problem>;
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using FVElementGeometry = typename GridGeometry::LocalView;
    using GridView = typename GridGeometry::GridView;
    using Element = typename GridView::template Codim<0>::Entity;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using ModelTraits = GetPropType<TypeTag, Properties::ModelTraits>;
    using Indices = typename ModelTraits::Indices;
    using SubControlVolumeFace = typename FVElementGeometry::SubControlVolumeFace;
    using FaceFrontalSubControlVolumeFace = typename GridGeometry::Traits::FaceFrontalSubControlVolumeFace;
    using FaceLateralSubControlVolumeFace = typename GridGeometry::Traits::FaceLateralSubControlVolumeFace;
    using Extrusion = Extrusion_t<GridGeometry>;
    using CellCenterPrimaryVariables = GetPropType<TypeTag, Properties::CellCenterPrimaryVariables>;
    using FacePrimaryVariables = GetPropType<TypeTag, Properties::FacePrimaryVariables>;
    using BoundaryTypes = typename ProblemTraits<Problem>::BoundaryTypes;
    using VelocityGradients = StaggeredVelocityGradients<Scalar, GridGeometry, BoundaryTypes, Indices>;

    static constexpr bool normalizePressure = getPropValue<TypeTag, Properties::NormalizePressure>();

    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;

    static constexpr auto cellCenterIdx = GridGeometry::cellCenterIdx();
    static constexpr auto faceIdx = GridGeometry::faceIdx();

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
                            * velocity * Extrusion::area(scvf) * scvf.directionSign();

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
        const auto& faceVars =  elemFaceVars[scvf];

        // Advective flux.
        if (problem.enableInertiaTerms())
        {
            // Get the average velocity at the center of the element (i.e. the location of the staggered face).
            const Scalar transportingVelocity = (velocitySelf + velocityOpposite) * 0.5;
            const bool selfIsUpstream = scvf.directionSign() != sign(transportingVelocity);

            StaggeredUpwindHelper<TypeTag, upwindSchemeOrder> upwindHelper(element, fvGeometry, scvf, elemFaceVars, elemVolVars, gridFluxVarsCache.staggeredUpwindMethods());
            frontalFlux += upwindHelper.computeUpwindFrontalMomentum(selfIsUpstream)
                           * transportingVelocity * -1.0 * scvf.directionSign();
        }

        // The volume variables within the current element. We only require those (and none of neighboring elements)
        // because the fluxes are calculated over the staggered face at the center of the element.
        const auto& insideVolVars = elemVolVars[scvf.insideScvIdx()];

        // Diffusive flux.
        const Scalar velocityGrad_ii = VelocityGradients::velocityGradII(scvf, faceVars) * scvf.directionSign();

        static const bool enableUnsymmetrizedVelocityGradient
            = getParamFromGroup<bool>(problem.paramGroup(), "FreeFlow.EnableUnsymmetrizedVelocityGradient", false);
        const Scalar factor = enableUnsymmetrizedVelocityGradient ? 1.0 : 2.0;
        frontalFlux -= factor * insideVolVars.effectiveViscosity() * velocityGrad_ii;

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

        // Account for the staggered face's area. For rectangular elements, this equals the area of the scvf
        // our velocity dof of interest lives on but with adjusted centroid
        const auto& scv = fvGeometry.scv(scvf.insideScvIdx());
        FaceFrontalSubControlVolumeFace frontalFace(scv.center(), scvf.area());
        return frontalFlux * Extrusion::area(frontalFace) * insideVolVars.extrusionFactor();
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
        FacePrimaryVariables lateralFlux(0.0);
        const auto& faceVars = elemFaceVars[scvf];
        const std::size_t numSubFaces = scvf.pairData().size();

        // If the current scvf is on a boundary, check if there is a Neumann BC for the stress in tangential direction.
        // Create a boundaryTypes object (will be empty if not at a boundary).
        std::optional<BoundaryTypes> currentScvfBoundaryTypes;
        if (scvf.boundary())
            currentScvfBoundaryTypes.emplace(problem.boundaryTypes(element, scvf));

        // Account for all sub faces.
        for (int localSubFaceIdx = 0; localSubFaceIdx < numSubFaces; ++localSubFaceIdx)
        {
            const auto eIdx = scvf.insideScvIdx();
            // Get the face normal to the face the dof lives on. The staggered sub face conincides with half of this lateral face.
            const auto& lateralFace = fvGeometry.scvf(eIdx, scvf.pairData(localSubFaceIdx).localLateralFaceIdx);

            // Create a boundaryTypes object (will be empty if not at a boundary).
            std::optional<BoundaryTypes> lateralFaceBoundaryTypes;

            // Check if there is face/element parallel to our face of interest where the dof lives on. If there is no parallel neighbor,
            // we are on a boundary where we have to check for boundary conditions.
            if (lateralFace.boundary())
            {
                // Retrieve the boundary types that correspond to the center of the lateral scvf. As a convention, we always query
                // the type of BCs at the center of the element's "actual" lateral scvf (not the face of the staggered control volume).
                // The value of the BC will be evaluated at the center of the staggered face.
                //     --------T######V                 || frontal face of staggered half-control-volume
                //     |      ||      | current scvf    #  lateral staggered face of interest (may lie on a boundary)
                //     |      ||      |                 x  dof position
                //     |      ||      x~~~~> vel.Self   -- element boundaries
                //     |      ||      |                 T  position at which the type of boundary conditions will be evaluated
                //     |      ||      |                    (center of lateral scvf)
                //     ----------------                 V  position at which the value of the boundary conditions will be evaluated
                //                                         (center of the staggered lateral face)
                lateralFaceBoundaryTypes.emplace(problem.boundaryTypes(element, lateralFace));
            }

            // If the current scvf is on a bounary and if there is a Neumann or Beavers-Joseph-(Saffmann) BC for the stress in tangential direction,
            // assign this value for the lateral momentum flux. No further calculations are required. We assume that all lateral faces
            // have the same type of BC (Neumann or Beavers-Joseph-(Saffmann)), but we sample the value at their actual positions.
            if (currentScvfBoundaryTypes)
            {
                // Handle Neumann BCs.
                if (currentScvfBoundaryTypes->isNeumann(Indices::velocity(lateralFace.directionIndex())))
                {
                    // Get the location of the lateral staggered face's center.
                    FaceLateralSubControlVolumeFace lateralScvf(lateralStaggeredSCVFCenter_(lateralFace, scvf, localSubFaceIdx), 0.5*lateralFace.area());
                    const auto& lateralStaggeredFaceCenter = lateralStaggeredFaceCenter_(scvf, localSubFaceIdx);
                    const auto lateralBoundaryFace = makeStaggeredBoundaryFace(scvf, lateralStaggeredFaceCenter);
                    lateralFlux += problem.neumann(element, fvGeometry, elemVolVars, elemFaceVars, lateralBoundaryFace)[Indices::velocity(lateralFace.directionIndex())]
                        * extrusionFactor_(elemVolVars, lateralFace) * Extrusion::area(lateralScvf) * lateralFace.directionSign();

                    continue;
                }

                // Treat the edge case when our current scvf has a BJ condition while the lateral face has a Neumann BC.
                // It is not clear how to evaluate the BJ condition here.
                // For symmetry reasons, our own scvf should then have the same Neumann flux as the lateral face.
                // TODO: We should clarify if this is the correct approach.
                if (currentScvfBoundaryTypes->isBeaversJoseph(Indices::velocity(lateralFace.directionIndex())) && lateralFaceBoundaryTypes &&
                    lateralFaceBoundaryTypes->isNeumann(Indices::velocity(scvf.directionIndex())))
                {
                    FaceLateralSubControlVolumeFace lateralScvf(lateralStaggeredSCVFCenter_(lateralFace, scvf, localSubFaceIdx), 0.5*lateralFace.area());
                    const auto& lateralStaggeredFaceCenter = lateralStaggeredFaceCenter_(scvf, localSubFaceIdx);
                    const auto lateralBoundaryFace = makeStaggeredBoundaryFace(lateralFace, lateralStaggeredFaceCenter);
                    lateralFlux += problem.neumann(element, fvGeometry, elemVolVars, elemFaceVars, lateralBoundaryFace)[Indices::velocity(scvf.directionIndex())]
                        * extrusionFactor_(elemVolVars, lateralFace) * Extrusion::area(lateralScvf) * lateralFace.directionSign();

                    continue;
                }
            }

            // Check if the lateral face (perpendicular to our current scvf) lies on a boundary. If yes, boundary conditions might need to be treated
            // and further calculations can be skipped.
            if (lateralFace.boundary())
            {
                // Check if we have a symmetry boundary condition. If yes, the tangential part of the momentum flux can be neglected
                // and we may skip any further calculations for the given sub face.
                if (lateralFaceBoundaryTypes->isSymmetry())
                    continue;

                // Handle Neumann boundary conditions. No further calculations are then required for the given sub face.
                if (lateralFaceBoundaryTypes->isNeumann(Indices::velocity(scvf.directionIndex())))
                {
                    // Construct a temporary scvf which corresponds to the staggered sub face, featuring the location
                    // the staggered faces's center.
                    FaceLateralSubControlVolumeFace lateralScvf(lateralStaggeredSCVFCenter_(lateralFace, scvf, localSubFaceIdx), 0.5*lateralFace.area());
                    const auto& lateralStaggeredFaceCenter = lateralStaggeredFaceCenter_(scvf, localSubFaceIdx);
                    const auto lateralBoundaryFace = makeStaggeredBoundaryFace(lateralFace, lateralStaggeredFaceCenter);
                    lateralFlux += problem.neumann(element, fvGeometry, elemVolVars, elemFaceVars, lateralBoundaryFace)[Indices::velocity(scvf.directionIndex())]
                                                  * elemVolVars[lateralFace.insideScvIdx()].extrusionFactor() * Extrusion::area(lateralScvf);

                    continue;
                }

                // Handle wall-function fluxes (only required for RANS models)
                if (incorporateWallFunction_(lateralFlux, problem, element, fvGeometry, scvf, elemVolVars, elemFaceVars, localSubFaceIdx))
                    continue;
            }

            // Check the consistency of the boundary conditions, exactly one of the following must be set
            if (lateralFaceBoundaryTypes)
            {
                std::bitset<3> admittableBcTypes;
                admittableBcTypes.set(0, lateralFaceBoundaryTypes->isDirichlet(Indices::pressureIdx));
                admittableBcTypes.set(1, lateralFaceBoundaryTypes->isDirichlet(Indices::velocity(scvf.directionIndex())));
                admittableBcTypes.set(2, lateralFaceBoundaryTypes->isBeaversJoseph(Indices::velocity(scvf.directionIndex())));
                if (admittableBcTypes.count() != 1)
                {
                    DUNE_THROW(Dune::InvalidStateException, "Invalid boundary conditions for lateral scvf "
                    "for the momentum equations at global position " << lateralStaggeredFaceCenter_(scvf, localSubFaceIdx)
                    << ", current scvf global position " << scvf.center());
                }
            }

            // If none of the above boundary conditions apply for the given sub face, proceed to calculate the tangential momentum flux.
            if (problem.enableInertiaTerms())
                lateralFlux += computeAdvectivePartOfLateralMomentumFlux_(problem, fvGeometry, element,
                                                                          scvf, elemVolVars, elemFaceVars,
                                                                          gridFluxVarsCache,
                                                                          currentScvfBoundaryTypes, lateralFaceBoundaryTypes,
                                                                          localSubFaceIdx);

            lateralFlux += computeDiffusivePartOfLateralMomentumFlux_(problem, fvGeometry, element,
                                                                      scvf, elemVolVars, faceVars,
                                                                      currentScvfBoundaryTypes, lateralFaceBoundaryTypes,
                                                                      localSubFaceIdx);
        }
        return lateralFlux;
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
    FacePrimaryVariables inflowOutflowBoundaryFlux(const Problem& problem,
                                                   const Element& element,
                                                   const SubControlVolumeFace& scvf,
                                                   const ElementVolumeVariables& elemVolVars,
                                                   const ElementFaceVariables& elemFaceVars) const
    {
        FacePrimaryVariables inOrOutflow(0.0);
        const auto& insideVolVars = elemVolVars[scvf.insideScvIdx()];

        // Advective momentum flux.
        if (problem.enableInertiaTerms())
        {
            const Scalar velocitySelf = elemFaceVars[scvf].velocitySelf();
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
        return inOrOutflow * scvf.directionSign() * Extrusion::area(scvf) * insideVolVars.extrusionFactor();
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
                                                                    const ElementVolumeVariables& elemVolVars,
                                                                    const ElementFaceVariables& elemFaceVars,
                                                                    const GridFluxVariablesCache& gridFluxVarsCache,
                                                                    const std::optional<BoundaryTypes>& currentScvfBoundaryTypes,
                                                                    const std::optional<BoundaryTypes>& lateralFaceBoundaryTypes,
                                                                    const int localSubFaceIdx)
    {
        const auto eIdx = scvf.insideScvIdx();
        const auto& lateralFace = fvGeometry.scvf(eIdx, scvf.pairData(localSubFaceIdx).localLateralFaceIdx);

        // Get the transporting velocity, located at the scvf perpendicular to the current scvf where the dof
        // of interest is located.
        const Scalar transportingVelocity = [&]()
        {
            const auto& faceVars = elemFaceVars[scvf];
            if (!scvf.boundary())
                return faceVars.velocityLateralInside(localSubFaceIdx);
            else
            {
                // Create a boundaryTypes object. Get the boundary conditions. We sample the type of BC at the center of the current scvf.
                const auto bcTypes = problem.boundaryTypes(element, scvf);

                if (bcTypes.isDirichlet(Indices::velocity(lateralFace.directionIndex())))
                {
                    // Construct a temporary scvf which corresponds to the staggered sub face, featuring the location
                    // the staggered faces's center.
                    const auto& lateralBoundaryFacePos = lateralStaggeredFaceCenter_(scvf, localSubFaceIdx);
                    const auto lateralBoundaryFace = makeStaggeredBoundaryFace(scvf, lateralBoundaryFacePos);
                    return problem.dirichlet(element, lateralBoundaryFace)[Indices::velocity(lateralFace.directionIndex())];
                }
                else if (bcTypes.isBeaversJoseph(Indices::velocity(lateralFace.directionIndex())))
                {
                    return VelocityGradients::beaversJosephVelocityAtCurrentScvf(problem, element, fvGeometry, scvf,  faceVars,
                                                                                 currentScvfBoundaryTypes, lateralFaceBoundaryTypes, localSubFaceIdx);
                }
                else
                    return faceVars.velocityLateralInside(localSubFaceIdx);
            }
        }();

        const bool selfIsUpstream = lateralFace.directionSign() == sign(transportingVelocity);
        StaggeredUpwindHelper<TypeTag, upwindSchemeOrder> upwindHelper(element, fvGeometry, scvf, elemFaceVars, elemVolVars, gridFluxVarsCache.staggeredUpwindMethods());
        FaceLateralSubControlVolumeFace lateralScvf(lateralStaggeredSCVFCenter_(lateralFace, scvf, localSubFaceIdx), 0.5*lateralFace.area());
        return upwindHelper.computeUpwindLateralMomentum(selfIsUpstream, lateralFace, localSubFaceIdx, currentScvfBoundaryTypes, lateralFaceBoundaryTypes)
               * transportingVelocity * lateralFace.directionSign() * Extrusion::area(lateralScvf) * extrusionFactor_(elemVolVars, lateralFace);
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
    FacePrimaryVariables computeDiffusivePartOfLateralMomentumFlux_(const Problem& problem,
                                                                    const FVElementGeometry& fvGeometry,
                                                                    const Element& element,
                                                                    const SubControlVolumeFace& scvf,
                                                                    const ElementVolumeVariables& elemVolVars,
                                                                    const FaceVariables& faceVars,
                                                                    const std::optional<BoundaryTypes>& currentScvfBoundaryTypes,
                                                                    const std::optional<BoundaryTypes>& lateralFaceBoundaryTypes,
                                                                    const int localSubFaceIdx)
    {
        const auto eIdx = scvf.insideScvIdx();
        const auto& lateralFace = fvGeometry.scvf(eIdx, scvf.pairData(localSubFaceIdx).localLateralFaceIdx);

        FacePrimaryVariables lateralDiffusiveFlux(0.0);

        static const bool enableUnsymmetrizedVelocityGradient
            = getParamFromGroup<bool>(problem.paramGroup(), "FreeFlow.EnableUnsymmetrizedVelocityGradient", false);

        // Get the volume variables of the own and the neighboring element. The neighboring
        // element is adjacent to the staggered face normal to the current scvf
        // where the dof of interest is located.
        const auto& insideVolVars = elemVolVars[lateralFace.insideScvIdx()];
        const auto& outsideVolVars = elemVolVars[lateralFace.outsideScvIdx()];

        // Get the averaged viscosity at the staggered face normal to the current scvf.
        const Scalar muAvg = lateralFace.boundary()
                             ? insideVolVars.effectiveViscosity()
                             : (insideVolVars.effectiveViscosity() + outsideVolVars.effectiveViscosity()) * 0.5;

        // Consider the shear stress caused by the gradient of the velocities normal to our face of interest.
        if (!enableUnsymmetrizedVelocityGradient)
        {
            if (!scvf.boundary() ||
                currentScvfBoundaryTypes->isDirichlet(Indices::velocity(lateralFace.directionIndex())) ||
                currentScvfBoundaryTypes->isBeaversJoseph(Indices::velocity(lateralFace.directionIndex())))
            {
                const Scalar velocityGrad_ji = VelocityGradients::velocityGradJI(problem, element, fvGeometry, scvf, faceVars, currentScvfBoundaryTypes, lateralFaceBoundaryTypes, localSubFaceIdx);
                // Account for the orientation of the staggered normal face's outer normal vector.
                lateralDiffusiveFlux -= muAvg * velocityGrad_ji * lateralFace.directionSign();
            }
        }

        // Consider the shear stress caused by the gradient of the velocities parallel to our face of interest.
        // If we have a Dirichlet condition for the pressure at the lateral face we assume to have a zero velocityGrad_ij velocity gradient
        // so we can skip the computation.
        if (!lateralFace.boundary() || !lateralFaceBoundaryTypes->isDirichlet(Indices::pressureIdx))
        {
            const Scalar velocityGrad_ij = VelocityGradients::velocityGradIJ(problem, element, fvGeometry, scvf, faceVars, currentScvfBoundaryTypes, lateralFaceBoundaryTypes, localSubFaceIdx);
            lateralDiffusiveFlux -= muAvg * velocityGrad_ij * lateralFace.directionSign();
        }

        // Account for the area of the staggered lateral face (0.5 of the coinciding scfv).
        FaceLateralSubControlVolumeFace lateralScvf(lateralStaggeredSCVFCenter_(lateralFace, scvf, localSubFaceIdx), 0.5*lateralFace.area());
        return lateralDiffusiveFlux * Extrusion::area(lateralScvf) * extrusionFactor_(elemVolVars, lateralFace);
    }

    /*!
     * \brief Get the location of the lateral staggered face's center.
     *        Only needed for boundary conditions if the current scvf or the lateral one is on a bounary.
     *
     * \verbatim
     *      --------#######o                 || frontal face of staggered half-control-volume
     *      |      ||      | current scvf    #  lateral staggered face of interest (may lie on a boundary)
     *      |      ||      |                 x  dof position
     *      |      ||      x~~~~> vel.Self   -- element boundaries, current scvf may lie on a boundary
     *      |      ||      |                 o  position at which the boundary conditions will be evaluated
     *      |      ||      |                    (lateralStaggeredFaceCenter)
     *      ----------------
     * \endverbatim
     */
    const GlobalPosition& lateralStaggeredFaceCenter_(const SubControlVolumeFace& scvf, const int localSubFaceIdx) const
    {
        return scvf.pairData(localSubFaceIdx).lateralStaggeredFaceCenter;
    };

    /*!
     * \brief Get the location of the lateral staggered sub control volume face center.
     *        Only needed for boundary conditions if the current scvf or the lateral one is on a bounary.
     *
     * \verbatim
     *      --------###o####                 || frontal face of staggered half-control-volume
     *      |      ||      | current scvf    #  lateral staggered face of interest (may lie on a boundary)
     *      |      ||      |                 x  dof position
     *      |      ||      x~~~~> vel.Self   -- element boundaries, current scvf may lie on a boundary
     *      |      ||      |                 o  center of the lateral staggered scvf
     *      |      ||      |
     *      ----------------
     * \endverbatim
     */
    GlobalPosition lateralStaggeredSCVFCenter_(const SubControlVolumeFace& lateralFace,
                                               const SubControlVolumeFace& currentFace,
                                               const int localSubFaceIdx) const
    {
        auto pos = lateralStaggeredFaceCenter_(currentFace, localSubFaceIdx) + lateralFace.center();
        pos *= 0.5;
        return pos;
    }

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
    { return false; }

    //! if a turbulence model is used, ask the problem is a wall function shall be employed and get the flux accordingly
    template<bool turbulenceModel = ModelTraits::usesTurbulenceModel(), std::enable_if_t<turbulenceModel, int> = 0>
    bool incorporateWallFunction_(FacePrimaryVariables& lateralFlux,
                                  const Problem& problem,
                                  const Element& element,
                                  const FVElementGeometry& fvGeometry,
                                  const SubControlVolumeFace& scvf,
                                  const ElementVolumeVariables& elemVolVars,
                                  const ElementFaceVariables& elemFaceVars,
                                  const std::size_t localSubFaceIdx) const
    {
        const auto eIdx = scvf.insideScvIdx();
        const auto& lateralFace = fvGeometry.scvf(eIdx, scvf.pairData(localSubFaceIdx).localLateralFaceIdx);

        if (problem.useWallFunction(element, lateralFace, Indices::velocity(scvf.directionIndex())))
        {
            FaceLateralSubControlVolumeFace lateralScvf(lateralStaggeredSCVFCenter_(lateralFace, scvf, localSubFaceIdx), 0.5*lateralFace.area());
            const auto lateralBoundaryFace = makeStaggeredBoundaryFace(lateralFace, lateralStaggeredFaceCenter_(scvf, localSubFaceIdx));
            lateralFlux += problem.wallFunction(element, fvGeometry, elemVolVars, elemFaceVars, scvf, lateralBoundaryFace)[Indices::velocity(scvf.directionIndex())]
                                               * extrusionFactor_(elemVolVars, lateralFace) * Extrusion::area(lateralScvf);
            return true;
        }
        else
            return false;
    }
};

} // end namespace Dumux

#endif
