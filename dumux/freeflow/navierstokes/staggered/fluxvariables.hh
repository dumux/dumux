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
#include <dumux/common/parameters.hh>
#include <dumux/common/properties.hh>

#include <dumux/discretization/fluxvariablesbase.hh>
#include <dumux/discretization/methods.hh>

#include <dumux/assembly/simpleassemblystructs.hh>

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
: public FluxVariablesBase<typename GET_PROP_TYPE(TypeTag, Problem),
                           typename GET_PROP_TYPE(TypeTag, FVGridGeometry)::LocalView,
                           typename GET_PROP_TYPE(TypeTag, GridVolumeVariables)::LocalView,
                           typename GET_PROP_TYPE(TypeTag, GridFluxVariablesCache)::LocalView>
{
    using GridVariables = typename GET_PROP_TYPE(TypeTag, GridVariables);

    using GridVolumeVariables = typename GridVariables::GridVolumeVariables;
    using ElementVolumeVariables = typename GridVolumeVariables::LocalView;
    using VolumeVariables = typename GridVolumeVariables::VolumeVariables;

    using GridFluxVariablesCache = typename GridVariables::GridFluxVariablesCache;
    using FluxVariablesCache = typename GridFluxVariablesCache::FluxVariablesCache;

    using GridFaceVariables = typename GridVariables::GridFaceVariables;
    using ElementFaceVariables = typename GridFaceVariables::LocalView;
    using FaceVariables = typename GridFaceVariables::FaceVariables;

    using Problem = typename GET_PROP_TYPE(TypeTag, Problem);
    using FVGridGeometry = typename GET_PROP_TYPE(TypeTag, FVGridGeometry);
    using FVElementGeometry = typename FVGridGeometry::LocalView;
    using GridView = typename FVGridGeometry::GridView;
    using Element = typename GridView::template Codim<0>::Entity;
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using ModelTraits = typename GET_PROP_TYPE(TypeTag, ModelTraits);
    using Indices = typename ModelTraits::Indices;
    using SubControlVolumeFace = typename FVElementGeometry::SubControlVolumeFace;
    using CellCenterPrimaryVariables = typename GET_PROP_TYPE(TypeTag, CellCenterPrimaryVariables);
    using FacePrimaryVariables = typename GET_PROP_TYPE(TypeTag, FacePrimaryVariables);
    using BoundaryTypes = typename GET_PROP_TYPE(TypeTag, BoundaryTypes);

    static constexpr bool enableInertiaTerms = GET_PROP_VALUE(TypeTag, EnableInertiaTerms);
    static constexpr bool normalizePressure = GET_PROP_VALUE(TypeTag, NormalizePressure);

    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;

    static constexpr auto cellCenterIdx = FVGridGeometry::cellCenterIdx();
    static constexpr auto faceIdx = FVGridGeometry::faceIdx();

    using SimpleMomentumBalanceSummands = typename GET_PROP_TYPE(TypeTag, SimpleMomentumBalanceSummands);

    enum {
        pressureIdx = Indices::pressureIdx
    };

public:

    using HeatConductionType = typename GET_PROP_TYPE(TypeTag, HeatConductionType);

    /*!
    * \brief Returns the advective flux over a sub control volume face.
    * \param elemVolVars All volume variables for the element
    * \param elemFaceVars The face variables
    * \param scvf The sub control volume face
    * \param upwindTerm The uwind term (i.e. the advectively transported quantity)
    */
    template<class UpwindTerm>
    static Scalar advectiveFluxForCellCenter(const Problem& problem,
                                             const Element& element,
                                             const ElementVolumeVariables& elemVolVars,
                                             const ElementFaceVariables& elemFaceVars,
                                             const SubControlVolumeFace &scvf,
                                             UpwindTerm upwindTerm)
    {
        const Scalar velocity = elemFaceVars[scvf].velocitySelf();
        const bool insideIsUpstream = scvf.directionSign() == sign(velocity);
        static const Scalar upWindWeight = getParamFromGroup<Scalar>(problem.paramGroup(), "Implicit.UpwindWeight");

        const auto& insideVolVars = elemVolVars[scvf.insideScvIdx()];
        const auto& outsideVolVars = elemVolVars[scvf.outsideScvIdx()];

        const auto& upstreamVolVars = insideIsUpstream ? insideVolVars : outsideVolVars;
        const auto& downstreamVolVars = insideIsUpstream ? outsideVolVars : insideVolVars;

        Scalar flux = (upWindWeight * upwindTerm(upstreamVolVars) +
                            (1.0 - upWindWeight) * upwindTerm(downstreamVolVars))
                            * scvf.area() * scvf.directionSign();

        if (scvf.boundary()){
            const auto bcTypes = problem.boundaryTypes(element, scvf);
            if (bcTypes.isDirichlet(Indices::velocity(scvf.directionIndex()))){
                flux *= velocity;
            }
        }

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
        result[Indices::conti0EqIdx - ModelTraits::dim()] = advectiveFluxForCellCenter(problem, element, elemVolVars, elemFaceVars, scvf, upwindTerm);

        return result;
    }

    /*!
    * \brief Returns the momentum flux over all staggered faces.
    */
    void computeMomentumFlux(const Problem& problem,
                             const Element& element,
                             const SubControlVolumeFace& scvf,
                             const FVElementGeometry& fvGeometry,
                             const ElementVolumeVariables& elemVolVars,
                             const ElementFaceVariables& elemFaceVars,
                             SimpleMomentumBalanceSummands& simpleMomentumBalanceSummands)
    {
        computeFrontalMomentumFlux(problem, element, scvf, fvGeometry, elemVolVars, elemFaceVars, simpleMomentumBalanceSummands);
        computeLateralMomentumFlux(problem, element, scvf, fvGeometry, elemVolVars, elemFaceVars, simpleMomentumBalanceSummands);
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
    void computeFrontalMomentumFlux(const Problem& problem,
                                    const Element& element,
                                    const SubControlVolumeFace& scvf,
                                    const FVElementGeometry& fvGeometry,
                                    const ElementVolumeVariables& elemVolVars,
                                    const ElementFaceVariables& elemFaceVars,
                                    SimpleMomentumBalanceSummands& simpleMomentumBalanceSummands)
    {
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

            //upwinding
            static const Scalar upwindWeight = getParamFromGroup<Scalar>(problem.paramGroup(), "Implicit.UpwindWeight");

            // Check if the the velocity of the dof at interest lies up- or downstream w.r.t. to the transporting velocity.
            const bool selfIsUpstream = scvf.directionSign() != sign(transportingVelocity);

            Scalar selfUpwindFactor;
            Scalar oppositeUpwindFactor;

            if (selfIsUpstream){
                selfUpwindFactor = upwindWeight;
                oppositeUpwindFactor = (1.0 - upwindWeight);
            }
            else {
                selfUpwindFactor = (1.0 - upwindWeight);
                oppositeUpwindFactor = upwindWeight;
            }

            // self
            // * scvf.directionSign() accounts for the orientation of the staggered face's normal outer normal vector
            // (pointing in opposite direction of the scvf's one).
            // "* scvf.area()" accounts for the staggered face's area. For rectangular elements, this equals the area // of the scvf our velocity dof of interest lives on.
            if (scvf.boundary() && problem.boundaryTypes(element, scvf).isDirichlet(Indices::velocity(scvf.directionIndex()))){
                simpleMomentumBalanceSummands.RHS -= transportingVelocity * insideVolVars.density() * -1.0 * scvf.directionSign() * scvf.area() * insideVolVars.extrusionFactor() * selfUpwindFactor * velocitySelf;
            }
            else {
                simpleMomentumBalanceSummands.selfCoefficient += transportingVelocity * insideVolVars.density() * -1.0 * scvf.directionSign() * scvf.area() * insideVolVars.extrusionFactor() *  selfUpwindFactor;
            }

            // opposite
            const auto eIdx = scvf.insideScvIdx();
            const auto opposingFace = fvGeometry.scvf(eIdx, scvf.localIdxOpposingFace());

            if (opposingFace.boundary() && problem.boundaryTypes(element, opposingFace).isDirichlet(Indices::velocity(scvf.directionIndex())))
            {
                simpleMomentumBalanceSummands.RHS -= transportingVelocity * insideVolVars.density() * -1.0 * scvf.directionSign() * scvf.area() * insideVolVars.extrusionFactor() * oppositeUpwindFactor * velocityOpposite;
            }
            else {
                simpleMomentumBalanceSummands.oppositeCoefficient += transportingVelocity *  insideVolVars.density() * -1.0 * scvf.directionSign() * scvf.area() * insideVolVars.extrusionFactor() * oppositeUpwindFactor;
            }
        }

        // Diffusive flux.
        // The velocity gradient already accounts for the orientation
        // of the staggered face's outer normal vector.
        static const bool enableUnsymmetrizedVelocityGradient
            = getParamFromGroup<bool>(problem.paramGroup(), "FreeFlow.EnableUnsymmetrizedVelocityGradient", false);
        const Scalar factor = enableUnsymmetrizedVelocityGradient ? 1.0 : 2.0;

        if(scvf.boundary() && problem.boundaryTypes(element, scvf).isDirichlet(Indices::velocity(scvf.directionIndex()))){
            simpleMomentumBalanceSummands.RHS -= factor * insideVolVars.effectiveViscosity() * scvf.area() * insideVolVars.extrusionFactor() * velocitySelf / scvf.selfToOppositeDistance();
        }
        else {
            simpleMomentumBalanceSummands.selfCoefficient += factor * insideVolVars.effectiveViscosity() * scvf.area() * insideVolVars.extrusionFactor() / scvf.selfToOppositeDistance();
        }

        const auto eIdx = scvf.insideScvIdx();
        const auto opposingFace = fvGeometry.scvf(eIdx, scvf.localIdxOpposingFace());

        if(opposingFace.boundary() && problem.boundaryTypes(element, opposingFace).isDirichlet(Indices::velocity(scvf.directionIndex()))){
            simpleMomentumBalanceSummands.RHS += factor * insideVolVars.effectiveViscosity() * scvf.area() * insideVolVars.extrusionFactor() * velocityOpposite / scvf.selfToOppositeDistance();
        }
        else {
            simpleMomentumBalanceSummands.oppositeCoefficient -= factor * insideVolVars.effectiveViscosity() * scvf.area() * insideVolVars.extrusionFactor() / scvf.selfToOppositeDistance();
        }

        // The pressure term.
        // Account for the orientation of the staggered face's normal outer normal vector
        // (pointing in opposite direction of the scvf's one).
        const auto& fixedPressureScvsIndexSet = problem.fixedPressureScvsIndexSet();

        const auto& it = std::find(fixedPressureScvsIndexSet.begin(), fixedPressureScvsIndexSet.end(), scvf.insideScvIdx());
        if(it != fixedPressureScvsIndexSet.end() /*scvf.insideScvIdx() is in fixedPressureScvsIndexSet*/)
        {
            const auto& insideScv = fvGeometry.scv(scvf.insideScvIdx());
            simpleMomentumBalanceSummands.RHS -= -1.0 * scvf.directionSign() * scvf.area() * insideVolVars.extrusionFactor() * problem.dirichletAtPos(insideScv.center())[pressureIdx];
        }
        else {
            simpleMomentumBalanceSummands.pressureCoefficient += -1.0 * scvf.directionSign() * scvf.area() * insideVolVars.extrusionFactor();
        }

        // If specified, the pressure can be normalized using the initial value on the scfv of interest.
        // The scvf is used to normalize by the same value from the left and right side.
        // Can potentially help to improve the condition number of the system matrix.
        if (normalizePressure){
            simpleMomentumBalanceSummands.RHS -= problem.initial(scvf)[Indices::pressureIdx] * scvf.directionSign() * scvf.area() * insideVolVars.extrusionFactor();
        }

        // Handle inflow or outflow conditions.
        // Treat the staggered half-volume adjacent to the boundary as if it was on the opposite side of the boundary.
        // The respective face's outer normal vector will point in the same direction as the scvf's one.
        if(scvf.boundary() && problem.boundaryTypes(element, scvf).isDirichlet(Indices::pressureIdx))
            inflowOutflowBoundaryFlux_(problem, element, scvf, elemVolVars, elemFaceVars, simpleMomentumBalanceSummands);
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
    *              |      ||      |                 # normal staggered sub faces over which fluxes are calculated
    *              |      ||      x~~~~> vel.Self
    *              |      ||      |                 x dof position
    *        scvf  |      ||      |
    *              --------########                -- element
    *                 scvf
    * \endverbatim
    */
    void computeLateralMomentumFlux(const Problem& problem,
                                    const Element& element,
                                    const SubControlVolumeFace& scvf,
                                    const FVElementGeometry& fvGeometry,
                                    const ElementVolumeVariables& elemVolVars,
                                    const ElementFaceVariables& elemFaceVars,
                                    SimpleMomentumBalanceSummands& simpleMomentumBalanceSummands)
    {
        auto& faceVars = elemFaceVars[scvf];
        const int numSubFaces = scvf.pairData().size();

        // Account for all sub faces.
        for(int localSubFaceIdx = 0; localSubFaceIdx < numSubFaces; ++localSubFaceIdx)
        {
            const auto eIdx = scvf.insideScvIdx();
            // Get the face normal to the face the dof lives on. The staggered sub face conincides with half of this normal face.
            const auto& normalFace = fvGeometry.scvf(eIdx, scvf.pairData()[localSubFaceIdx].localNormalFaceIdx);

            // Construct a temporary scvf which corresponds to the staggered sub face, featuring the location
            // the sub faces's center.
            auto localSubFaceCenter = scvf.pairData(localSubFaceIdx).virtualOuterParallelFaceDofPos - normalFace.center();
            localSubFaceCenter *= 0.5;
            localSubFaceCenter += normalFace.center();
            const auto localSubFace = makeGhostFace_(normalFace, localSubFaceCenter);

            // Retrieve the boundary types that correspond to the sub face.
            const auto bcTypes = problem.boundaryTypes(element, localSubFace);

            // Check if there is face/element parallel to our face of interest where the dof lives on. If there is no parallel neighbor,
            // we are on a boundary where we have to check for boundary conditions.
            if(!scvf.hasParallelNeighbor(localSubFaceIdx))
            {
                // Check if we have a symmetry boundary condition. If yes, the tangential part of the momentum flux can be neglected
                // and we may skip any further calculations for the given sub face.
                if(bcTypes.isSymmetry())
                    continue;

                // Handle Neumann boundary conditions. No further calculations are then required for the given sub face.
                if(bcTypes.isNeumann(Indices::velocity(scvf.directionIndex())))
                {
                    simpleMomentumBalanceSummands.RHS -= problem.neumann(element, fvGeometry, elemVolVars, elemFaceVars, localSubFace)[Indices::velocity(scvf.directionIndex())]
                                                  * elemVolVars[normalFace.insideScvIdx()].extrusionFactor() * normalFace.area() * 0.5;
                    continue;
                }

                // Handle wall-function fluxes (only required for RANS models)
                if(problem.useWallFunction(element, localSubFace, Indices::velocity(scvf.directionIndex())))
                {
                    simpleMomentumBalanceSummands.RHS -= problem.wallFunction(element, fvGeometry, elemVolVars, elemFaceVars, scvf, localSubFace)[Indices::velocity(scvf.directionIndex())]
                                                       * elemVolVars[normalFace.insideScvIdx()].extrusionFactor() * normalFace.area() * 0.5;
                    continue;
                }
            }

            // If there is no symmetry or Neumann boundary condition for the given sub face, proceed to calculate the tangential momentum flux.
            if(enableInertiaTerms)
                computeAdvectivePartOfLateralMomentumFlux_(problem, element, scvf, normalFace, elemVolVars, faceVars, localSubFaceIdx, bcTypes, simpleMomentumBalanceSummands, fvGeometry);

            computeDiffusivePartOfLateralMomentumFlux_(problem, element, scvf, fvGeometry, normalFace, elemVolVars, faceVars, localSubFaceIdx, bcTypes, simpleMomentumBalanceSummands);
        }
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
    void computeAdvectivePartOfLateralMomentumFlux_(const Problem& problem,
                                                    const Element& element,
                                                    const SubControlVolumeFace& scvf,
                                                    const SubControlVolumeFace& normalFace,
                                                    const ElementVolumeVariables& elemVolVars,
                                                    const FaceVariables& faceVars,
                                                    const int localSubFaceIdx,
                                                    const BoundaryTypes& bcTypes,
                                                    SimpleMomentumBalanceSummands& simpleMomentumBalanceSummands,
                                                    const FVElementGeometry& fvGeometry)
    {
        // Get the transporting velocity, located at the scvf perpendicular to the current scvf where the dof
        // of interest is located.
        const Scalar transportingVelocity = faceVars.velocityNormalInside(localSubFaceIdx);

        // Get the volume variables of the own and the neighboring element
        const auto& insideVolVars = elemVolVars[normalFace.insideScvIdx()];
        const auto& outsideVolVars = elemVolVars[normalFace.outsideScvIdx()];

        //upwinding
        static const Scalar upwindWeight = getParamFromGroup<Scalar>(problem.paramGroup(), "Implicit.UpwindWeight");

        // Check whether the own or the neighboring element is upstream.
        const bool ownElementIsUpstream = ( normalFace.directionSign() == sign(transportingVelocity) );

        Scalar selfUpwindFactor;
        Scalar parallelUpwindFactor;

        if (ownElementIsUpstream) {
            selfUpwindFactor = upwindWeight;
            parallelUpwindFactor = 1.0 - upwindWeight;
        }
        else {
            selfUpwindFactor = 1.0 - upwindWeight;
            parallelUpwindFactor = upwindWeight;
        }

        const Scalar density = selfUpwindFactor * insideVolVars.density() + parallelUpwindFactor * outsideVolVars.density();

        //self
        // Account for the orientation of the staggered normal face's outer normal vector
        // and its area (0.5 of the coinciding scfv).
        if(scvf.boundary() && problem.boundaryTypes(element, scvf).isDirichlet(Indices::velocity(scvf.directionIndex()))){
            simpleMomentumBalanceSummands.RHS -= transportingVelocity * density * normalFace.directionSign() * normalFace.area()  * 0.5 * extrusionFactor_(elemVolVars, normalFace) * faceVars.velocitySelf() * selfUpwindFactor;
        }
        else {
            simpleMomentumBalanceSummands.selfCoefficient += transportingVelocity * density * normalFace.directionSign() * normalFace.area() * 0.5 * extrusionFactor_(elemVolVars, normalFace) * selfUpwindFactor;
        }

        //parallel
        if(scvf.hasParallelNeighbor(localSubFaceIdx)){
            const auto eIdx = normalFace.outsideScvIdx();
            const auto parallelFace = fvGeometry.scvf(eIdx, scvf.localFaceIdx());
            if (parallelFace.boundary() && problem.boundaryTypes(element, parallelFace).isDirichlet(Indices::velocity(scvf.directionIndex()))) {
                simpleMomentumBalanceSummands.RHS -= transportingVelocity * density * normalFace.directionSign() * normalFace.area() * 0.5 * extrusionFactor_(elemVolVars, normalFace) * faceVars.velocityParallel(localSubFaceIdx) * parallelUpwindFactor;
            }
            else {
                simpleMomentumBalanceSummands.parallelCoefficients[localSubFaceIdx] += transportingVelocity * density * normalFace.directionSign() * normalFace.area() * 0.5 * extrusionFactor_(elemVolVars, normalFace) * parallelUpwindFactor;
            }
        }
        //treat boundary
        else{
            // Lambda to conveniently get the outer parallel velocity for normal faces that are on the boundary
            // and therefore have no neighbor. Calls the problem to retrieve a fixed value set on the boundary.
            auto getParallelVelocityFromBoundary = [&]()
            {
                const auto ghostFace = makeParallelGhostFace_(scvf, localSubFaceIdx);

                // Check if we have a Beavers-Joseph-Saffman condition.
                // If yes, the parallel velocity at the boundary is calculated accordingly.
                if (bcTypes.isBJS(Indices::velocity(scvf.directionIndex())))
                    return problem.bjsVelocity(scvf, normalFace, localSubFaceIdx, faceVars.velocitySelf());
                return problem.dirichlet(element, ghostFace)[Indices::velocity(scvf.directionIndex())];
            };

            const Scalar velocityParallel = getParallelVelocityFromBoundary();

            simpleMomentumBalanceSummands.RHS -= transportingVelocity * density * velocityParallel * normalFace.directionSign() * normalFace.area() * 0.5 * extrusionFactor_(elemVolVars, normalFace) * parallelUpwindFactor;
        }
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
    void computeDiffusivePartOfLateralMomentumFlux_(const Problem& problem,
                                                    const Element& element,
                                                    const SubControlVolumeFace& scvf,
                                                    const FVElementGeometry& fvGeometry,
                                                    const SubControlVolumeFace& normalFace,
                                                    const ElementVolumeVariables& elemVolVars,
                                                    const FaceVariables& faceVars,
                                                    const int localSubFaceIdx,
                                                    const BoundaryTypes& bcTypes,
                                                    SimpleMomentumBalanceSummands& simpleMomentumBalanceSummands)
    {
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
                //normal gradient
                // Account for the orientation of the staggered normal face's outer normal vector.
                // Account for the area of the staggered normal face (0.5 of the coinciding scfv)
                if (normalFace.boundary() && problem.boundaryTypes(element, normalFace).isDirichlet(Indices::velocity(scvf.directionIndex()))){
                    const auto innerNormalVelocity =  faceVars.velocityNormalInside(localSubFaceIdx);
                    simpleMomentumBalanceSummands.RHS -= scvf.directionSign() * muAvg * normalFace.directionSign() *  extrusionFactor_(elemVolVars, normalFace) * normalFace.area() * 0.5 * innerNormalVelocity / scvf.pairData(localSubFaceIdx).normalDistance;
                }
                else{
                    simpleMomentumBalanceSummands.innerNormalCoefficients[localSubFaceIdx] += scvf.directionSign() * muAvg * normalFace.directionSign() * extrusionFactor_(elemVolVars, normalFace) * normalFace.area() * 0.5 / scvf.pairData(localSubFaceIdx).normalDistance;
                }

                const auto outerNormalFace = fvGeometry.scvf(scvf.outsideScvIdx(), normalFace.localFaceIdx());

                if(outerNormalFace.boundary() && problem.boundaryTypes(element, outerNormalFace).isDirichlet(Indices::velocity(scvf.directionIndex()))){
                    const auto outerNormalVelocity =  faceVars.velocityNormalOutside(localSubFaceIdx);
                    simpleMomentumBalanceSummands.RHS += scvf.directionSign() * muAvg * normalFace.directionSign() * extrusionFactor_(elemVolVars, normalFace) * normalFace.area() * 0.5 * outerNormalVelocity / scvf.pairData(localSubFaceIdx).normalDistance;
                }
                else{
                    simpleMomentumBalanceSummands.outerNormalCoefficients[localSubFaceIdx] -= scvf.directionSign() * muAvg * normalFace.directionSign() * extrusionFactor_(elemVolVars, normalFace) * normalFace.area() * 0.5 / scvf.pairData(localSubFaceIdx).normalDistance;
                }
            }
        }

        //parallel gradient
        if (!(normalFace.boundary() && problem.boundaryTypes(element, normalFace).isOutflow(Indices::velocity(scvf.directionIndex())))){
            const auto innerParallelVelocity = faceVars.velocitySelf();
            // Account for the area of the staggered normal face (0.5 of the coinciding scfv)
            if (scvf.boundary() && problem.boundaryTypes(element, scvf).isDirichlet(Indices::velocity(scvf.directionIndex()))){
                simpleMomentumBalanceSummands.RHS -= muAvg * extrusionFactor_(elemVolVars, normalFace) * normalFace.area() * 0.5 * innerParallelVelocity/ scvf.pairData(localSubFaceIdx).parallelDistance;
            }
            else{
                simpleMomentumBalanceSummands.selfCoefficient += muAvg * extrusionFactor_(elemVolVars, normalFace) * normalFace.area() * 0.5/ scvf.pairData(localSubFaceIdx).parallelDistance;
            }

            if(scvf.hasParallelNeighbor(localSubFaceIdx))
            {
                const auto parallelFace = fvGeometry.scvf(normalFace.outsideScvIdx(), scvf.localFaceIdx());
                if (parallelFace.boundary() && problem.boundaryTypes(element, parallelFace).isDirichlet(Indices::velocity(scvf.directionIndex()))){
                    const auto outerParallelVelocity = faceVars.velocityParallel(localSubFaceIdx);
                    simpleMomentumBalanceSummands.RHS += muAvg * extrusionFactor_(elemVolVars, normalFace) * normalFace.area() * 0.5 * outerParallelVelocity/ scvf.pairData(localSubFaceIdx).parallelDistance;
                }
                else{
                    simpleMomentumBalanceSummands.parallelCoefficients[localSubFaceIdx] -= muAvg * extrusionFactor_(elemVolVars, normalFace) * normalFace.area() * 0.5/ scvf.pairData(localSubFaceIdx).parallelDistance;
                }
            }
            else
            {
                // Lambda to conveniently get the outer parallel velocity for normal faces that are on the boundary
                // and therefore have no neighbor. Calls the problem to retrieve a fixed value set on the boundary.
                auto getParallelVelocityFromBoundary = [&]()
                {
                    const auto ghostFace = makeParallelGhostFace_(scvf, localSubFaceIdx);

                    // Check if we have a Beavers-Joseph-Saffman condition.
                    // If yes, the parallel velocity at the boundary is calculated accordingly.
                    if(bcTypes.isBJS(Indices::velocity(scvf.directionIndex())))
                        return problem.bjsVelocity(scvf, normalFace, localSubFaceIdx, innerParallelVelocity);

                    else if(bcTypes.isDirichlet(Indices::velocity(scvf.directionIndex())))
                        return problem.dirichlet(element, ghostFace)[Indices::velocity(scvf.directionIndex())];

                    //bcTypes.isDirichlet(Indices::pressureIdx) is assumed to be true now
                    else
                        return innerParallelVelocity;
                };

                const Scalar outerParallelVelocity = getParallelVelocityFromBoundary();

                // Account for the area of the staggered normal face (0.5 of the coinciding scfv)
                simpleMomentumBalanceSummands.RHS += muAvg * outerParallelVelocity * extrusionFactor_(elemVolVars, normalFace) * normalFace.area() * 0.5/ scvf.pairData(localSubFaceIdx).parallelDistance;
            }
        }
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
    void inflowOutflowBoundaryFlux_(const Problem& problem,
                                    const Element& element,
                                    const SubControlVolumeFace& scvf,
                                    const ElementVolumeVariables& elemVolVars,
                                    const ElementFaceVariables& elemFaceVars,
                                    SimpleMomentumBalanceSummands& simpleMomentumBalanceSummands)
    {
        FacePrimaryVariables inOrOutflow(0.0);

        const auto& insideVolVars = elemVolVars[scvf.insideScvIdx()];

        // Advective momentum flux.
        if(enableInertiaTerms)
        {
            const Scalar velocitySelf = elemFaceVars[scvf].velocitySelf();
            const auto& outsideVolVars = elemVolVars[scvf.outsideScvIdx()];
            const auto& upVolVars = (scvf.directionSign() == sign(velocitySelf)) ?
                                    insideVolVars : outsideVolVars;

            // "* scvf.directionSign()" is to account for the orientation of the face at the boundary
            simpleMomentumBalanceSummands.selfCoefficient += velocitySelf * upVolVars.density() * scvf.directionSign() * scvf.area() * insideVolVars.extrusionFactor();
        }

        // Apply a pressure at the boundary.
        const Scalar boundaryPressure = normalizePressure
                                        ? (problem.dirichlet(element, scvf)[Indices::pressureIdx] -
                                           problem.initial(scvf)[Indices::pressureIdx])
                                        : problem.dirichlet(element, scvf)[Indices::pressureIdx];

        // "* scvf.directionSign()" is to account for the orientation of the face at the boundary
        simpleMomentumBalanceSummands.RHS -= boundaryPressure * scvf.directionSign() * scvf.area() * insideVolVars.extrusionFactor();

        // contribution from grad v^T -> outflow straight out of pipe, assume linear interpolation of velocity,
        simpleMomentumBalanceSummands.oppositeCoefficient += insideVolVars.viscosity() * scvf.area() * insideVolVars.extrusionFactor() * scvf.directionSign() / scvf.selfToOppositeDistance();

        simpleMomentumBalanceSummands.selfCoefficient -= insideVolVars.viscosity() * scvf.area() * insideVolVars.extrusionFactor() * scvf.directionSign() / scvf.selfToOppositeDistance();
    }

private:

    //! helper function to conveniently create a ghost face used to retrieve boundary values from the problem
    SubControlVolumeFace makeGhostFace_(const SubControlVolumeFace& ownScvf, const GlobalPosition& pos) const
    {
        return SubControlVolumeFace(pos, std::vector<unsigned int>{ownScvf.insideScvIdx(), ownScvf.outsideScvIdx()}, ownScvf.directionIndex(), ownScvf.dofIndex(), ownScvf.index());
    };

    //! helper function to conveniently create a ghost face which is outside the domain, parallel to the scvf of interest
    SubControlVolumeFace makeParallelGhostFace_(const SubControlVolumeFace& ownScvf, const int localSubFaceIdx) const
    {
        return makeGhostFace_(ownScvf, ownScvf.pairData(localSubFaceIdx).virtualOuterParallelFaceDofPos);
    };

    //! helper function to get the averaged extrusion factor for a face
    static Scalar extrusionFactor_(const ElementVolumeVariables& elemVolVars, const SubControlVolumeFace& scvf)
    {
        const auto& insideVolVars = elemVolVars[scvf.insideScvIdx()];
        const auto& outsideVolVars = elemVolVars[scvf.outsideScvIdx()];
        return harmonicMean(insideVolVars.extrusionFactor(), outsideVolVars.extrusionFactor());
    }
};

} // end namespace

#endif
