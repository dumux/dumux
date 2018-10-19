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

#include <array>

#include <dumux/common/math.hh>
#include <dumux/common/exceptions.hh>
#include <dumux/common/parameters.hh>
#include <dumux/common/properties.hh>

#include <dumux/freeflow/higherorderapproximation.hh>

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

    static constexpr bool enableInertiaTerms = GET_PROP_VALUE(TypeTag, EnableInertiaTerms);
    static constexpr bool normalizePressure = GET_PROP_VALUE(TypeTag, NormalizePressure);

    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;

    static constexpr auto cellCenterIdx = FVGridGeometry::cellCenterIdx();
    static constexpr auto faceIdx = FVGridGeometry::faceIdx();

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

        const Scalar flux = (upWindWeight * upwindTerm(upstreamVolVars) +
                            (1.0 - upWindWeight) * upwindTerm(downstreamVolVars))
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
                                             const ElementFaceVariables& elemFaceVars)
    {
        return computeFrontalMomentumFlux(problem, element, scvf, fvGeometry, elemVolVars, elemFaceVars) +
               computeLateralMomentumFlux(problem, element, scvf, fvGeometry, elemVolVars, elemFaceVars);
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
                                                    const ElementFaceVariables& elemFaceVars)
    {
        FacePrimaryVariables frontalFlux(0.0);

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

            Scalar momentum = 0.0;

            // Variables that will store the velocities of interest:
            // velocities[0]: downstream velocity
            // velocities[1]: upsteram velocity
            // velocities[2]: upstream-upstream velocity
            std::array<Scalar, 3> velocities{0.0, 0.0, 0.0};

            // Variables that will store the distances between the dofs of interest:
            // distances[0]: upstream to downstream distance
            // distances[1]: upstream-upstream to upstream distance
            // distances[2]: downstream staggered cell size
            std::array<Scalar, 3> distances{0.0, 0.0, 0.0};

            const auto& highOrder = problem.higherOrderApproximation();

            // If a Tvd approac has been specified and I am not too near to the boundary I can use a second order
            // approximation for the velocity. In this frontal flux I use for the density always the value that I have on the scvf.
            if (highOrder.tvdApproach() != TvdApproach::none && canFrontalSecondOrder_(scvf, selfIsUpstream, velocities, distances, elemFaceVars[scvf]))
            {
                switch (highOrder.tvdApproach())
                {
                    case TvdApproach::uniform :
                    {
                        momentum = highOrder.tvd(velocities[0], velocities[1], velocities[2], insideVolVars.density());
                        break;
                    }
                    case TvdApproach::li :
                    {
                        momentum = highOrder.tvd(velocities[0], velocities[1], velocities[2], distances[0], distances[1], selfIsUpstream, insideVolVars.density());
                        break;
                    }
                    case TvdApproach::hou :
                    {
                        momentum = highOrder.tvd(velocities[0], velocities[1], velocities[2], distances[0], distances[1], distances[2], insideVolVars.density());
                        break;
                    }
                    default:
                    {
                        DUNE_THROW(ParameterException, "\nTvd approach " << static_cast<int>(highOrder.tvdApproach()) << " is not implemented.\n" <<
                                    static_cast<int>(TvdApproach::uniform) << ": Uniform Tvd\n" <<
                                    static_cast<int>(TvdApproach::li) << ": Li's approach\n" <<
                                    static_cast<int>(TvdApproach::hou) << ": Hou's approach");
                        break;
                    }
                }
            }
            else
                momentum = highOrder.upwind(velocities[0], velocities[1], insideVolVars.density());

            // Account for the orientation of the staggered face's normal outer normal vector
            // (pointing in opposite direction of the scvf's one).
            frontalFlux += transportingVelocity * momentum * -1.0 * scvf.directionSign();
        }

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
    FacePrimaryVariables computeLateralMomentumFlux(const Problem& problem,
                                                    const Element& element,
                                                    const SubControlVolumeFace& scvf,
                                                    const FVElementGeometry& fvGeometry,
                                                    const ElementVolumeVariables& elemVolVars,
                                                    const ElementFaceVariables& elemFaceVars)
    {
        FacePrimaryVariables normalFlux(0.0);
        auto& faceVars = elemFaceVars[scvf];
        const int numSubFaces = scvf.pairData().size();

        // Account for all sub faces.
        for(int localSubFaceIdx = 0; localSubFaceIdx < numSubFaces; ++localSubFaceIdx)
        {
            const auto eIdx = scvf.insideScvIdx();
            // Get the face normal to the face the dof lives on. The staggered sub face conincides with half of this normal face.
            const auto& normalFace = fvGeometry.scvf(eIdx, scvf.pairData(localSubFaceIdx).localNormalFaceIdx);

            bool isDirichletPressure = false; // check for Dirichlet boundary condition for the pressure
            bool isBJS = false; // check for Beavers-Joseph-Saffman boundary condition

            // Check if there is face/element parallel to our face of interest where the dof lives on. If there is no parallel neighbor,
            // we are on a boundary where we have to check for boundary conditions.
            if(!scvf.hasParallelNeighbor(localSubFaceIdx,0))
            {
                // Construct a temporary scvf which corresponds to the staggered sub face, featuring the location
                // the sub faces's center.
                auto localSubFaceCenter = scvf.pairData(localSubFaceIdx).virtualFirstParallelFaceDofPos - normalFace.center();
                localSubFaceCenter *= 0.5;
                localSubFaceCenter += normalFace.center();
                const auto localSubFace = makeGhostFace_(normalFace, localSubFaceCenter);

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
                if(problem.useWallFunction(element, localSubFace, Indices::velocity(scvf.directionIndex())))
                {
                    normalFlux += problem.wallFunction(element, fvGeometry, elemVolVars, elemFaceVars, scvf, localSubFace)[Indices::velocity(scvf.directionIndex())]
                                                       * elemVolVars[normalFace.insideScvIdx()].extrusionFactor() * normalFace.area() * 0.5;
                    continue;
                }

                // Check if we have a Beavers-Joseph-Saffman condition or a Dirichlet condition for the velocity or a Dirichlet condition for the pressure.
                // Then the parallel velocity at the boundary is calculated accordingly for the advective part and the diffusive part of the normal momentum flux.
                if (bcTypes.isDirichlet(Indices::pressureIdx))
                    isDirichletPressure = true;
                else if (bcTypes.isBJS(Indices::velocity(scvf.directionIndex())))
                    isBJS = true;
                else if (bcTypes.isDirichlet(Indices::velocity(scvf.directionIndex())) == false)
                    DUNE_THROW(Dune::InvalidStateException,  "Something went wrong with the boundary conditions "
                           "for the momentum equations at global position " << localSubFaceCenter);
            }

            // If there is no symmetry or Neumann boundary condition for the given sub face, proceed to calculate the tangential momentum flux.
            if(enableInertiaTerms)
                normalFlux += computeAdvectivePartOfLateralMomentumFlux_(problem, element, scvf, normalFace, elemVolVars, faceVars, localSubFaceIdx, isDirichletPressure, isBJS);

            normalFlux += computeDiffusivePartOfLateralMomentumFlux_(problem, element, scvf, normalFace, elemVolVars, faceVars, localSubFaceIdx, isDirichletPressure, isBJS);
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
                                                                    const Element& element,
                                                                    const SubControlVolumeFace& scvf,
                                                                    const SubControlVolumeFace& normalFace,
                                                                    const ElementVolumeVariables& elemVolVars,
                                                                    const FaceVariables& faceVars,
                                                                    const int localSubFaceIdx,
                                                                    const bool isDirichletPressure,
                                                                    const bool isBJS)
    {
        // Get the transporting velocity, located at the scvf perpendicular to the current scvf where the dof
        // of interest is located.
        const Scalar transportingVelocity = faceVars.velocityNormalInside(localSubFaceIdx);

        // Check whether the own or the neighboring element is upstream.
        const bool selfIsUpstream = ( normalFace.directionSign() == sign(transportingVelocity) );

        // Get the volume variables of the own and the neighboring element
        const auto& insideVolVars = elemVolVars[normalFace.insideScvIdx()];
        const auto& outsideVolVars = elemVolVars[normalFace.outsideScvIdx()];

        Scalar momentum = 0.0;

        // Variables that will store the velocities of interest:
        // velocities[0]: downstream velocity
        // velocities[1]: upsteram velocity
        // velocities[2]: upstream-upstream velocity
        std::array<Scalar, 3> velocities{0.0, 0.0, 0.0};

        // Variables that will store the distances between the dofs of interest:
        // distances[0]: upstream to downstream distance
        // distances[1]: upstream-upstream to upstream distance
        // distances[2]: downstream staggered cell size
        std::array<Scalar, 3> distances{0.0, 0.0, 0.0};

        const auto& highOrder = problem.higherOrderApproximation();

        // If a Tvd approach has been specified and I am not too near to the boundary I can use a second order approximation.
        if (highOrder.tvdApproach() != TvdApproach::none && canLateralSecondOrder_(scvf, selfIsUpstream, localSubFaceIdx, velocities, distances, problem, element, faceVars, isDirichletPressure, isBJS))
        {
            switch (highOrder.tvdApproach())
            {
                case TvdApproach::uniform :
                {
                    momentum = highOrder.tvd(velocities[0], velocities[1], velocities[2], selfIsUpstream ? insideVolVars.density() : outsideVolVars.density());
                    break;
                }
                case TvdApproach::li :
                {
                    momentum = highOrder.tvd(velocities[0], velocities[1], velocities[2], distances[0], distances[1], selfIsUpstream, selfIsUpstream ? insideVolVars.density() : outsideVolVars.density());
                    break;
                }
                case TvdApproach::hou :
                {
                    momentum = highOrder.tvd(velocities[0], velocities[1], velocities[2], distances[0], distances[1], distances[2], selfIsUpstream ? insideVolVars.density() : outsideVolVars.density());
                    break;
                }
                default:
                {
                    DUNE_THROW(ParameterException, "\nTvd approach " << static_cast<int>(highOrder.tvdApproach()) << " is not implemented.\n" <<
                                static_cast<int>(TvdApproach::uniform) << ": Uniform Tvd\n" <<
                                static_cast<int>(TvdApproach::li) << ": Li's approach\n" <<
                                static_cast<int>(TvdApproach::hou) << ": Hou's approach");
                    break;
                }
            }
        }
        else
            momentum = highOrder.upwind(velocities[0], velocities[1], selfIsUpstream ? insideVolVars.density() : outsideVolVars.density());

        // Account for the orientation of the staggered normal face's outer normal vector
        // and its area (0.5 of the coinciding scfv).
        return transportingVelocity * momentum * normalFace.directionSign() * normalFace.area() * 0.5 * extrusionFactor_(elemVolVars, normalFace);
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
                                                                    const Element& element,
                                                                    const SubControlVolumeFace& scvf,
                                                                    const SubControlVolumeFace& normalFace,
                                                                    const ElementVolumeVariables& elemVolVars,
                                                                    const FaceVariables& faceVars,
                                                                    const int localSubFaceIdx,
                                                                    const bool isDirichletPressure,
                                                                    const bool isBJS)
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
        if (!isDirichletPressure)
        {
            // For the parallel derivative, get the velocities at the current (own) scvf
            // and at the parallel one at the neighboring scvf.
            const Scalar innerParallelVelocity = faceVars.velocitySelf();

            // Lambda to conveniently get the outer parallel velocity for normal faces that are on the boundary
            // and therefore have no neighbor. Calls the problem to retrieve a fixed value set on the boundary.
            auto getParallelVelocityFromBoundary = [&]()
            {
                const auto ghostFace = makeParallelGhostFace_(scvf, localSubFaceIdx);
                if (isBJS)
                    return problem.bjsVelocity(scvf, normalFace, localSubFaceIdx, innerParallelVelocity);
                return problem.dirichlet(element, ghostFace)[Indices::velocity(scvf.directionIndex())];
            };

            const Scalar velocityFirstParallel = scvf.hasParallelNeighbor(localSubFaceIdx,0)
                                               ? faceVars.velocityParallel(localSubFaceIdx,0)
                                               : getParallelVelocityFromBoundary();

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
        if(enableInertiaTerms)
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

    //! helper function to conveniently create a ghost face used to retrieve boundary values from the problem
    SubControlVolumeFace makeGhostFace_(const SubControlVolumeFace& ownScvf, const GlobalPosition& pos) const
    {
        return SubControlVolumeFace(pos, std::vector<unsigned int>{ownScvf.insideScvIdx(), ownScvf.outsideScvIdx()}, ownScvf.directionIndex(), ownScvf.dofIndex(), ownScvf.index());
    };

    //! helper function to conveniently create a ghost face which is outside the domain, parallel to the scvf of interest
    SubControlVolumeFace makeParallelGhostFace_(const SubControlVolumeFace& ownScvf, const int localSubFaceIdx) const
    {
        return makeGhostFace_(ownScvf, ownScvf.pairData(localSubFaceIdx).virtualFirstParallelFaceDofPos);
    };

    //! helper function to get the averaged extrusion factor for a face
    static Scalar extrusionFactor_(const ElementVolumeVariables& elemVolVars, const SubControlVolumeFace& scvf)
    {
        const auto& insideVolVars = elemVolVars[scvf.insideScvIdx()];
        const auto& outsideVolVars = elemVolVars[scvf.outsideScvIdx()];
        return harmonicMean(insideVolVars.extrusionFactor(), outsideVolVars.extrusionFactor());
    }

    /*!
     * \brief Check if a second order approximation for the frontal part of the advective term can be used
     *
     * This helper function checks if the scvf of interest is not too near to the
     * boundary so that a dof upstream with respect to the upstream dof is available.
     *
     * \param ownScvf The SubControlVolumeFace we are considering
     * \param selfIsUpstream @c true if the velocity ownScvf is upstream wrt the transporting velocity
     * \param velocities Variable that will store the velocities of interest
     * \param distances Variable that will store the distances of interest
     * \param faceVars The face variables related to ownScvf
     */
    bool canFrontalSecondOrder_(const SubControlVolumeFace& ownScvf,
                                const bool selfIsUpstream,
                                std::array<Scalar, 3>& velocities,
                                std::array<Scalar, 3>& distances,
                                const FaceVariables& faceVars) const
    {
        // Depending on selfIsUpstream I can assign the downstream and the upstream velocities,
        // then I have to check if I have a forward or a backward neighbot to retrieve
        // an "upstream-upstream velocity" and be able to use a second order scheme.
        if (selfIsUpstream)
        {
            velocities[0] = faceVars.velocityOpposite();
            velocities[1] = faceVars.velocitySelf();

            if (ownScvf.hasForwardNeighbor(0))
            {
                velocities[2] = faceVars.velocityForward(0);
                distances[0] = ownScvf.selfToOppositeDistance();
                distances[1] = ownScvf.axisData().inAxisForwardDistances[0];
                distances[2] = 0.5 * (ownScvf.axisData().selfToOppositeDistance + ownScvf.axisData().inAxisForwardDistances[0]);
                return true;
            }
            else
                return false;
        }
        else
        {
            velocities[0] = faceVars.velocitySelf();
            velocities[1] = faceVars.velocityOpposite();

            if (ownScvf.hasBackwardNeighbor(0))
            {
                velocities[2] = faceVars.velocityBackward(0);
                distances[0] = ownScvf.selfToOppositeDistance();
                distances[1] = ownScvf.axisData().inAxisBackwardDistances[0];
                // pairData(0) or pairData(1) have the same normalDistance
                distances[2] = ownScvf.pairData(0).normalDistance;
                return true;
            }
            else
                return false;
        }
    }

    /*!
     * \brief Check if a second order approximation for the lateral part of the advective term can be used
     *
     * This helper function checks if the scvf of interest is not too near to the
     * boundary so that a dof upstream with respect to the upstream dof is available.
     *
     * \param ownScvf The SubControlVolumeFace we are considering
     * \param selfIsUpstream @c true if the velocity ownScvf is upstream wrt the transporting velocity
     * \param localSubFaceIdx The local subface index
     * \param velocities Variable that will store the velocities of interest
     * \param distances Variable that will store the distances of interest
     */
    bool canLateralSecondOrder_(const SubControlVolumeFace& ownScvf,
                                const bool selfIsUpstream,
                                const int localSubFaceIdx,
                                std::array<Scalar, 3>& velocities,
                                std::array<Scalar, 3>& distances,
                                const Problem& problem,
                                const Element& element,
                                const FaceVariables& faceVars,
                                const bool isDirichletPressure,
                                const bool isBJS) const
    {
        const SubControlVolumeFace& normalFace = problem.fvGridGeometry().scvf(ownScvf.insideScvIdx(), ownScvf.pairData(localSubFaceIdx).localNormalFaceIdx);

        // The local index of the faces that is opposite to localSubFaceIdx
        const int oppositeSubFaceIdx = localSubFaceIdx % 2 ? localSubFaceIdx - 1 : localSubFaceIdx + 1;

        // Lambda to conveniently get the outer parallel velocity for normal faces that are on the boundary
        // and therefore have no neighbor. Calls the problem to retrieve a fixed value set on the boundary.
        auto getParallelVelocityFromBoundary = [&]()
        {
            // If there is a Dirichlet condition for the pressure we assume zero gradient for the velocity,
            // so the velocity at the boundary equal to that on the scvf.
            if (isDirichletPressure)
                return faceVars.velocitySelf();

            const auto ghostFace = makeParallelGhostFace_(ownScvf, localSubFaceIdx);
            if (isBJS)
                return problem.bjsVelocity(ownScvf, normalFace, localSubFaceIdx, faceVars.velocitySelf());
            return problem.dirichlet(element, ghostFace)[Indices::velocity(ownScvf.directionIndex())];
        };

        if (selfIsUpstream)
        {
            // I can assign the upstream velocity. The downstream velocity can be assigned or retrieved
            // from the boundary if there is no parallel neighbor.
            velocities[1] = faceVars.velocitySelf();

            if(ownScvf.hasParallelNeighbor(localSubFaceIdx, 0))
            {
                velocities[0] = faceVars.velocityParallel(localSubFaceIdx, 0);
                distances[2] = ownScvf.pairData(localSubFaceIdx).parallelDistances[1];
            }
            else
            {
                velocities[0] = getParallelVelocityFromBoundary();
                distances[2] = ownScvf.area() / 2.0;
            }

            // The "upstream-upsteram" velocity is retrieved from the other parallel neighbor
            // or from the boundary.
            if (ownScvf.hasParallelNeighbor(oppositeSubFaceIdx, 0))
                velocities[2] = faceVars.velocityParallel(oppositeSubFaceIdx, 0);
            else
                velocities[2] = getParallelVelocityFromOtherBoundary_(problem, ownScvf, oppositeSubFaceIdx, element, velocities[1]);

            distances[0] = ownScvf.cellCenteredParallelDistance(localSubFaceIdx, 0);
            distances[1] = ownScvf.cellCenteredParallelDistance(oppositeSubFaceIdx, 0);

            return true;
        }
        else
        {
            // The self velocity is downstream, then if there is no parallel neighbor I can not use
            // a second order approximation beacuse I have only two velocities.
            velocities[0] = faceVars.velocitySelf();

            if (!ownScvf.hasParallelNeighbor(localSubFaceIdx, 0))
            {
                velocities[1] = getParallelVelocityFromBoundary();
                return false;
            }

            velocities[1] = faceVars.velocityParallel(localSubFaceIdx, 0);

            // If there is another parallel neighbor I can assign the "upstream-upstream"
            // velocity, otherwise I retrieve it from the boundary.
            if (ownScvf.hasParallelNeighbor(localSubFaceIdx, 1))
                velocities[2] = faceVars.velocityParallel(localSubFaceIdx, 1);
            else
            {
                const Element& elementParallel = problem.fvGridGeometry().element(problem.fvGridGeometry().scv(normalFace.outsideScvIdx()));
                const SubControlVolumeFace& firstParallelScvf = problem.fvGridGeometry().scvf(normalFace.outsideScvIdx(), ownScvf.localFaceIdx());
                velocities[2] = getParallelVelocityFromOtherBoundary_(problem, firstParallelScvf, localSubFaceIdx, elementParallel, velocities[1]);
            }

            distances[0] = ownScvf.cellCenteredParallelDistance(localSubFaceIdx, 0);
            distances[1] = ownScvf.cellCenteredParallelDistance(localSubFaceIdx, 1);
            distances[2] = ownScvf.area();

            return true;
        }
    }

    /*!
     * \brief Return a velocity value from a boundary for which the boundary conditions have to be checked.
     *
     * \param problem The problem
     * \param scvf The SubControlVolumeFace that is normal to the boundary
     * \param localIdx The local index of the face that is on the boundary
     * \param boundaryElement The element that is on the boundary
     * \param parallelVelocity The velocity over scvf
     */
    Scalar getParallelVelocityFromOtherBoundary_(const Problem& problem,
                                                 const SubControlVolumeFace& scvf,
                                                 const int localIdx,
                                                 const Element& boundaryElement,
                                                 const Scalar parallelVelocity) const
    {
        // A ghost subface at the boundary is created, featuring the location of the sub face's center
        const SubControlVolumeFace& boundaryNormalFace = problem.fvGridGeometry().scvf(scvf.insideScvIdx(), scvf.pairData(localIdx).localNormalFaceIdx);
        GlobalPosition boundarySubFaceCenter = scvf.pairData(localIdx).virtualFirstParallelFaceDofPos + boundaryNormalFace.center();
        boundarySubFaceCenter *= 0.5;
        const SubControlVolumeFace boundarySubFace = makeGhostFace_(boundaryNormalFace, boundarySubFaceCenter);

        // The boundary condition is checked, in case of symmetry or Dirichlet for the pressure
        // a gradient of zero is assumed in the direction normal to the bounadry, while if there is
        // Dirichlet of BJS for the velocity the related values are exploited.
        const auto bcTypes = problem.boundaryTypes(boundaryElement, boundarySubFace);

        if (bcTypes.isDirichlet(Indices::velocity(scvf.directionIndex())))
        {
            const SubControlVolumeFace ghostFace = makeParallelGhostFace_(scvf, localIdx);
            return problem.dirichlet(boundaryElement, ghostFace)[Indices::velocity(scvf.directionIndex())];
        }
        else if (bcTypes.isSymmetry() || bcTypes.isDirichlet(Indices::pressureIdx))
            return parallelVelocity;
        else if (bcTypes.isBJS(Indices::velocity(scvf.directionIndex())))
        {
            const SubControlVolumeFace ghostFace = makeParallelGhostFace_(scvf, localIdx);
            return problem.bjsVelocity(scvf, boundaryNormalFace, localIdx, parallelVelocity);
        }
        else
        {
            // Neumann conditions are not well implemented
            DUNE_THROW(Dune::InvalidStateException, "Something went wrong with the boundary conditions for the momentum equations at global position " << boundarySubFaceCenter);
        }
    }
};

} // end namespace Dumux

#endif // DUMUX_NAVIERSTOKES_STAGGERED_FLUXVARIABLES_HH
