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
 * \copydoc Dumux::StaggeredVelocityGradients
 */
#ifndef DUMUX_NAVIERSTOKES_STAGGERED_VELOCITYGRADIENTS_HH
#define DUMUX_NAVIERSTOKES_STAGGERED_VELOCITYGRADIENTS_HH

#include <optional>
#include <dumux/common/exceptions.hh>
#include <dumux/common/parameters.hh>

namespace Dumux {

/*!
 * \ingroup NavierStokesModel
 * \brief Helper class for calculating the velocity gradients for the Navier-Stokes model using the staggered grid discretization.
 */
template<class Scalar, class GridGeometry, class BoundaryTypes, class Indices>
class StaggeredVelocityGradients
{
    using FVElementGeometry = typename GridGeometry::LocalView;
    using GridView = typename GridGeometry::GridView;
    using Element = typename GridView::template Codim<0>::Entity;
    using SubControlVolumeFace = typename FVElementGeometry::SubControlVolumeFace;
    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;

public:

    /*!
     * \brief Returns the in-axis velocity gradient.
     *
     * \verbatim
     *              ---------=======                 == and # staggered half-control-volume
     *              |       #      | current scvf
     *              |       #      |                 # staggered face over which fluxes are calculated
     *   vel.Opp <~~|       O~~>   x~~~~> vel.Self
     *              |       #      |                 x dof position
     *        scvf  |       #      |
     *              --------========                 -- element
     *
     *                                               O position at which gradient is evaluated
     * \endverbatim
     */
    template<class FaceVariables>
    static Scalar velocityGradII(const SubControlVolumeFace& scvf,
                                 const FaceVariables& faceVars)
    {
        // The velocities of the dof at interest and the one of the opposite scvf.
        const Scalar velocitySelf = faceVars.velocitySelf();
        const Scalar velocityOpposite = faceVars.velocityOpposite();

        return ((velocityOpposite - velocitySelf) / scvf.selfToOppositeDistance()) * scvf.directionSign();
    }

    /*!
     * \brief Returns the velocity gradient perpendicular to the orientation of our current scvf.
     *
     * \verbatim
     *              ----------------
     *              |              |vel.
     *              |              |Parallel
     *              |              |~~~~>       ------->
     *              |              |             ------>     * gradient
     *              |              |              ----->
     *       scvf   ---------######O:::::::::      ---->     || and # staggered half-control-volume (own element)
     *              |      ||      | curr. ::       --->
     *              |      ||      | scvf  ::        -->     :: staggered half-control-volume (neighbor element)
     *              |      ||      x~~~~>  ::         ->
     *              |      ||      | vel.  ::                # lateral staggered faces over which fluxes are calculated
     *        scvf  |      ||      | Self  ::
     *              ---------#######:::::::::                x dof position
     *                 scvf
     *                                                       -- elements
     *
     *                                                       O position at which gradient is evaluated
     * \endverbatim
     *
     * ------------
     * |     xxxx s
     * |     xxxx a
     * |     xxxx s
     * -----------O-----------
     * |     yyyy s zzzz     |
     * |     yyyy b zzzz     |
     * |     yyyy s zzzz     |
     * -----------------------
     *
     * In a corner geometry (scvf is sas or sbs), we calculate the velocity gradient at O, by
     * (velocity(a)-velocity(b))/distance(a,b) for the half-control volumes x and y, but by
     * (velocity(O)-velocity(b))/distance(O,b) for z. This does not harm flux continuity (x and y use the same
     * formulation). We do this different formulation for y (approximate gradient by central differncing) and
     * z (approximate gradient by forward/backward differencing), because it is the natural way of implementing
     * it and it is not clear which gradient is the better approximation in this case anyway.
     * Particularly, for the frequent case of no-slip, no-flow boundaries, the velocity would be zero at O and a
     * and thus, the gradient within the flow domain might be better approximated by velocity(b)/distanc(O,b)
     * than by velocity(b)/distance(a,b).
     */
    template<class Problem, class FaceVariables>
    static Scalar velocityGradIJ(const Problem& problem,
                                 const Element& element,
                                 const FVElementGeometry& fvGeometry,
                                 const SubControlVolumeFace& scvf,
                                 const FaceVariables& faceVars,
                                 const std::optional<BoundaryTypes>& currentScvfBoundaryTypes,
                                 const std::optional<BoundaryTypes>& lateralFaceBoundaryTypes,
                                 const std::size_t localSubFaceIdx)
    {
        const auto eIdx = scvf.insideScvIdx();
        const auto& lateralScvf = fvGeometry.scvf(eIdx, scvf.pairData(localSubFaceIdx).localLateralFaceIdx);

        // For the velocityGrad_ij derivative, get the velocities at the current (own) scvf
        // and at the parallel one at the neighboring scvf.
        const Scalar innerParallelVelocity = faceVars.velocitySelf();

        const auto outerParallelVelocity = [&]()
        {
            if (!lateralScvf.boundary())
                return faceVars.velocityParallel(localSubFaceIdx, 0);
            else if (lateralFaceBoundaryTypes->isDirichlet(Indices::velocity(scvf.directionIndex())))
            {
                // Sample the value of the Dirichlet BC at the center of the staggered lateral face.
                const auto& lateralBoundaryFacePos = lateralStaggeredFaceCenter_(scvf, localSubFaceIdx);
                const auto lateralBoundaryFace = makeStaggeredBoundaryFace(lateralScvf, lateralBoundaryFacePos);
                return problem.dirichlet(element, lateralBoundaryFace)[Indices::velocity(scvf.directionIndex())];
            }
            else if (lateralFaceBoundaryTypes->isBeaversJoseph(Indices::velocity(scvf.directionIndex())))
            {
                return beaversJosephVelocityAtLateralScvf(problem, element, fvGeometry, scvf,  faceVars,
                                                          currentScvfBoundaryTypes, lateralFaceBoundaryTypes, localSubFaceIdx);
            }
            else
                DUNE_THROW(Dune::InvalidStateException, "Invalid lateral boundary type at " << lateralScvf.center());
        }();

        // The velocity gradient already accounts for the orientation
        // of the staggered face's outer normal vector. This also correctly accounts for the reduced
        // distance used in the gradient if the lateral scvf lies on a boundary.
        return (outerParallelVelocity - innerParallelVelocity)
               / scvf.parallelDofsDistance(localSubFaceIdx, 0) * lateralScvf.directionSign();
    }

    /*!
     * \brief Returns the velocity gradient in line with our current scvf.
     *
     * \verbatim
     *                      ^       gradient
     *                      |  ^
     *                      |  |  ^
     *                      |  |  |  ^
     *                      |  |  |  |  ^
     *                      |  |  |  |  |  ^
     *                      |  |  |  |  |  |
     *
     *              ----------------
     *              |              |
     *              |    in.norm.  |
     *              |       vel.   |
     *              |       ^      |        ^ out.norm.vel.
     *              |       |      |        |
     *       scvf   ---------######O:::::::::       || and # staggered half-control-volume (own element)
     *              |      ||      | curr. ::
     *              |      ||      | scvf  ::       :: staggered half-control-volume (neighbor element)
     *              |      ||      x~~~~>  ::
     *              |      ||      | vel.  ::       # lateral staggered faces over which fluxes are calculated
     *        scvf  |      ||      | Self  ::
     *              ---------#######:::::::::       x dof position
     *                 scvf
     *                                              -- elements
     *
     *                                              O position at which gradient is evaluated
     * \endverbatim
     */
    template<class Problem, class FaceVariables>
    static Scalar velocityGradJI(const Problem& problem,
                                 const Element& element,
                                 const FVElementGeometry& fvGeometry,
                                 const SubControlVolumeFace& scvf,
                                 const FaceVariables& faceVars,
                                 const std::optional<BoundaryTypes>& currentScvfBoundaryTypes,
                                 const std::optional<BoundaryTypes>& lateralFaceBoundaryTypes,
                                 const std::size_t localSubFaceIdx)
    {
        const auto eIdx = scvf.insideScvIdx();
        const auto& lateralScvf = fvGeometry.scvf(eIdx, scvf.pairData(localSubFaceIdx).localLateralFaceIdx);

        // Assume a zero velocity gradient for pressure boundary conditions.
        if (currentScvfBoundaryTypes && currentScvfBoundaryTypes->isDirichlet(Indices::pressureIdx))
            return 0.0;

        // For the velocityGrad_ji gradient, get the velocities perpendicular to the velocity at the current scvf.
        // The inner one is located at staggered face within the own element,
        // the outer one at the respective staggered face of the element on the other side of the
        // current scvf.
        const Scalar innerLateralVelocity = faceVars.velocityLateralInside(localSubFaceIdx);
        const Scalar outerLateralVelocity = [&]()
        {
            if (!scvf.boundary())
                return faceVars.velocityLateralOutside(localSubFaceIdx);
            else if (currentScvfBoundaryTypes->isDirichlet(Indices::velocity(lateralScvf.directionIndex())))
            {
                // Sample the value of the Dirichlet BC at the center of the lateral face intersecting with the boundary.
                const auto& lateralBoundaryFacePos = lateralStaggeredFaceCenter_(scvf, localSubFaceIdx);
                const auto lateralBoundaryFace = makeStaggeredBoundaryFace(scvf, lateralBoundaryFacePos);
                return problem.dirichlet(element, lateralBoundaryFace)[Indices::velocity(lateralScvf.directionIndex())];
            }
            else if (currentScvfBoundaryTypes->isBeaversJoseph(Indices::velocity(lateralScvf.directionIndex())))
            {
                return beaversJosephVelocityAtCurrentScvf(problem, element, fvGeometry, scvf,  faceVars,
                                                          currentScvfBoundaryTypes, lateralFaceBoundaryTypes, localSubFaceIdx);
            }
            else
                DUNE_THROW(Dune::InvalidStateException, "Invalid lateral boundary types at " << lateralScvf.center());
        }();

        // Calculate the velocity gradient in positive coordinate direction.
        const Scalar lateralDeltaV = scvf.normalInPosCoordDir()
                                    ? (outerLateralVelocity - innerLateralVelocity)
                                    : (innerLateralVelocity - outerLateralVelocity);

        return lateralDeltaV / scvf.pairData(localSubFaceIdx).lateralDistance;
    }

    /*!
     * \brief Returns the Beavers-Jospeh slip velocity for a scvf which lies on the boundary itself.
     *
     * \verbatim
     *                  in.norm.  B-J slip
     *                     vel.   vel.
     *                     ^       ^
     *                     |       |
     *       scvf   ---------######|*               * boundary
     *              |      ||      |* curr.
     *              |      ||      |* scvf          || and # staggered half-control-volume (own element)
     *              |      ||      x~~~~>
     *              |      ||      |* vel.          # lateral staggered faces
     *        scvf  |      ||      |* Self
     *              ---------#######*                x dof position
     *                 scvf
     *                                              -- element
     * \endverbatim
     *
     */
    template<class Problem, class FaceVariables>
    static Scalar beaversJosephVelocityAtCurrentScvf(const Problem& problem,
                                                     const Element& element,
                                                     const FVElementGeometry& fvGeometry,
                                                     const SubControlVolumeFace& scvf,
                                                     const FaceVariables& faceVars,
                                                     const std::optional<BoundaryTypes>& currentScvfBoundaryTypes,
                                                     const std::optional<BoundaryTypes>& lateralFaceBoundaryTypes,
                                                     const std::size_t localSubFaceIdx)
    {
        const auto eIdx = scvf.insideScvIdx();
        const auto& lateralScvf = fvGeometry.scvf(eIdx, scvf.pairData(localSubFaceIdx).localLateralFaceIdx);
        const Scalar innerLateralVelocity = faceVars.velocityLateralInside(localSubFaceIdx);

        const auto tangentialVelocityGradient = [&]()
        {
            // If the current scvf is on a boundary and if a Dirichlet BC for the pressure or a BJ condition for
            // the slip velocity is set there, assume a tangential velocity gradient of zero along the lateral face
            // (towards the current scvf).
            static const bool unsymmetrizedGradientForBJ = getParamFromGroup<bool>(problem.paramGroup(),
                                                           "FreeFlow.EnableUnsymmetrizedVelocityGradientForBeaversJoseph", false);

            if (unsymmetrizedGradientForBJ)
                return 0.0;

            if (lateralScvf.boundary())
            {
                if (lateralFaceBoundaryTypes->isDirichlet(Indices::pressureIdx) ||
                    lateralFaceBoundaryTypes->isBeaversJoseph(Indices::velocity(scvf.directionIndex())))
                    return 0.0;
            }

            return velocityGradIJ(problem, element, fvGeometry, scvf, faceVars, currentScvfBoundaryTypes, lateralFaceBoundaryTypes, localSubFaceIdx);
        }();

        return problem.beaversJosephVelocity(element,
                                             fvGeometry.scv(scvf.insideScvIdx()),
                                             lateralScvf,
                                             scvf, /*on boundary*/
                                             innerLateralVelocity,
                                             tangentialVelocityGradient);
    }

    /*!
     * \brief Returns the Beavers-Jospeh slip velocity for a lateral scvf which lies on the boundary.
     *
     * \verbatim
     *                             B-J slip                  * boundary
     *              ************** vel. *****
     *       scvf   ---------##### ~~~~> ::::                || and # staggered half-control-volume (own element)
     *              |      ||      | curr. ::
     *              |      ||      | scvf  ::                :: staggered half-control-volume (neighbor element)
     *              |      ||      x~~~~>  ::
     *              |      ||      | vel.  ::                # lateral staggered faces
     *        scvf  |      ||      | Self  ::
     *              ---------#######:::::::::                x dof position
     *                 scvf
     *                                                       -- elements
     * \endverbatim
     */
    template<class Problem, class FaceVariables>
    static Scalar beaversJosephVelocityAtLateralScvf(const Problem& problem,
                                                     const Element& element,
                                                     const FVElementGeometry& fvGeometry,
                                                     const SubControlVolumeFace& scvf,
                                                     const FaceVariables& faceVars,
                                                     const std::optional<BoundaryTypes>& currentScvfBoundaryTypes,
                                                     const std::optional<BoundaryTypes>& lateralFaceBoundaryTypes,
                                                     const std::size_t localSubFaceIdx)
    {
        const auto eIdx = scvf.insideScvIdx();
        const auto& lateralScvf = fvGeometry.scvf(eIdx, scvf.pairData(localSubFaceIdx).localLateralFaceIdx);
        const Scalar innerParallelVelocity = faceVars.velocitySelf();

        const auto tangentialVelocityGradient = [&]()
        {
            // If the current scvf is on a boundary and if a Dirichlet BC for the pressure or a BJ condition for
            // the slip velocity is set there, assume a tangential velocity gradient of zero along the lateral face
            // (towards the current scvf).
            static const bool unsymmetrizedGradientForBJ = getParamFromGroup<bool>(problem.paramGroup(),
                                                           "FreeFlow.EnableUnsymmetrizedVelocityGradientForBeaversJoseph", false);

            if (unsymmetrizedGradientForBJ)
                return 0.0;

            if (scvf.boundary())
            {
                if (currentScvfBoundaryTypes->isDirichlet(Indices::pressureIdx) ||
                    currentScvfBoundaryTypes->isBeaversJoseph(Indices::velocity(lateralScvf.directionIndex())))
                    return 0.0;
            }

            return velocityGradJI(problem, element, fvGeometry, scvf, faceVars, currentScvfBoundaryTypes, lateralFaceBoundaryTypes, localSubFaceIdx);
        }();

        return problem.beaversJosephVelocity(element,
                                             fvGeometry.scv(scvf.insideScvIdx()),
                                             scvf,
                                             lateralScvf, /*on boundary*/
                                             innerParallelVelocity,
                                             tangentialVelocityGradient);
    }

private:

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
    static const GlobalPosition& lateralStaggeredFaceCenter_(const SubControlVolumeFace& scvf, const int localSubFaceIdx)
    {
        return scvf.pairData(localSubFaceIdx).lateralStaggeredFaceCenter;
    };
};

} // end namespace Dumux

#endif
