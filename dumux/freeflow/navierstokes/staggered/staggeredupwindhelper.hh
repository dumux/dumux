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
 * \copydoc Dumux::StaggeredUpwindHelper
 */
#ifndef DUMUX_NAVIERSTOKES_STAGGERED_UPWINDVARIABLES_HH
#define DUMUX_NAVIERSTOKES_STAGGERED_UPWINDVARIABLES_HH

#include <array>
#include <optional>

#include <dumux/common/math.hh>
#include <dumux/common/exceptions.hh>
#include <dumux/common/parameters.hh>
#include <dumux/common/properties.hh>
#include <dumux/common/typetraits/problem.hh>

#include <dumux/discretization/method.hh>
#include <dumux/freeflow/staggeredupwindmethods.hh>
#include "velocitygradients.hh"

namespace Dumux {

/*!
 * \ingroup NavierStokesModel
 * \brief The upwinding variables class for the Navier-Stokes model using the staggered grid discretization.
 */
template<class TypeTag, int upwindSchemeOrder>
class StaggeredUpwindHelper
{
    using GridVariables = GetPropType<TypeTag, Properties::GridVariables>;

    using GridVolumeVariables = typename GridVariables::GridVolumeVariables;
    using ElementVolumeVariables = typename GridVolumeVariables::LocalView;
    using VolumeVariables = typename GridVolumeVariables::VolumeVariables;

    using GridFluxVariablesCache = typename GridVariables::GridFluxVariablesCache;
    using UpwindScheme = typename GridFluxVariablesCache::UpwindScheme;
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
    using SubControlVolume = typename FVElementGeometry::SubControlVolume;
    using SubControlVolumeFace = typename FVElementGeometry::SubControlVolumeFace;
    using CellCenterPrimaryVariables = GetPropType<TypeTag, Properties::CellCenterPrimaryVariables>;
    using FacePrimaryVariables = GetPropType<TypeTag, Properties::FacePrimaryVariables>;
    using BoundaryTypes = typename ProblemTraits<Problem>::BoundaryTypes;
    using VelocityGradients = StaggeredVelocityGradients<Scalar, GridGeometry, BoundaryTypes, Indices>;

    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;
    static constexpr bool useHigherOrder = upwindSchemeOrder > 1;

    static_assert(upwindSchemeOrder <= 2, "Not implemented: Order higher than 2!");
    static_assert(upwindSchemeOrder <= 1 || GridFluxVariablesCache::cachingEnabled,
        "Higher order upwind method requires caching enabled for the GridFluxVariablesCache!");
    static_assert(upwindSchemeOrder <= 1 || GridGeometry::cachingEnabled,
        "Higher order upwind method requires caching enabled for the GridGeometry!");
    static_assert(upwindSchemeOrder <= 1 || GridFaceVariables::cachingEnabled,
        "Higher order upwind method requires caching enabled for the GridFaceVariables!");
    static_assert(upwindSchemeOrder <= 1 || GridVolumeVariables::cachingEnabled,
        "Higher order upwind method requires caching enabled for the GridGeometry!");

public:
    StaggeredUpwindHelper(const Element& element,
                          const FVElementGeometry& fvGeometry,
                          const SubControlVolumeFace& scvf,
                          const ElementFaceVariables& elemFaceVars,
                          const ElementVolumeVariables& elemVolVars,
                          const UpwindScheme& upwindScheme)
    : element_(element)
    , fvGeometry_(fvGeometry)
    , scvf_(scvf)
    , elemFaceVars_(elemFaceVars)
    , faceVars_(elemFaceVars[scvf])
    , elemVolVars_(elemVolVars)
    , upwindScheme_(upwindScheme)
    {}

    /*!
     * \brief Returns the momentum in the frontal directon.
     *
     *        Checks if the model has higher order methods enabled and if the scvf in
     *        question is far enough from the boundary such that higher order methods can be employed.
     *        Then the corresponding set of momenta are collected and the prescribed
     *        upwinding method is used to calculate the momentum.
     */
    FacePrimaryVariables computeUpwindFrontalMomentum(const bool selfIsUpstream) const
    {
        const auto density = elemVolVars_[scvf_.insideScvIdx()].density();

        // for higher order schemes do higher order upwind reconstruction
        if constexpr (useHigherOrder)
        {
            // only second order is implemented so far
            if (canDoFrontalSecondOrder_(selfIsUpstream))
            {
                const auto distances = getFrontalDistances_(selfIsUpstream);
                const auto upwindMomenta = getFrontalSecondOrderUpwindMomenta_(density, selfIsUpstream);
                return upwindScheme_.tvd(upwindMomenta, distances, selfIsUpstream, upwindScheme_.tvdApproach());
            }
        }

        // otherwise apply first order upwind scheme
        const auto upwindMomenta = getFrontalFirstOrderUpwindMomenta_(density, selfIsUpstream);
        return upwindScheme_.upwind(upwindMomenta[0], upwindMomenta[1]);
    }

    /*!
     * \brief Returns the momentum in the lateral directon.
     *
     *        Evaluates which face is upstream.
     *        Checks if the model has higher order methods enabled and if the scvf in
     *        question is far enough from the boundary such that higher order methods can be employed.
     *        Then the corresponding set of momenta are collected and the prescribed
     *        upwinding method is used to calculate the momentum.
     */
    FacePrimaryVariables computeUpwindLateralMomentum(const bool selfIsUpstream,
                                                      const SubControlVolumeFace& lateralFace,
                                                      const int localSubFaceIdx,
                                                      const std::optional<BoundaryTypes>& currentScvfBoundaryTypes,
                                                      const std::optional<BoundaryTypes>& lateralFaceBoundaryTypes) const
    {
        // Check whether the own or the neighboring element is upstream.
        // Get the transporting velocity, located at the scvf perpendicular to the current scvf where the dof
        // of interest is located.
        const Scalar insideDensity = elemVolVars_[lateralFace.insideScvIdx()].density();
        const Scalar outsideDensity = elemVolVars_[lateralFace.outsideScvIdx()].density();

        // for higher order schemes do higher order upwind reconstruction
        if constexpr (useHigherOrder)
        {
            // only second order is implemented so far
            if (canDoLateralSecondOrder_(selfIsUpstream, localSubFaceIdx))
            {
                const auto upwindMomenta = getLateralSecondOrderUpwindMomenta_(insideDensity, outsideDensity, selfIsUpstream, localSubFaceIdx, currentScvfBoundaryTypes, lateralFaceBoundaryTypes);
                const auto distances = getLateralDistances_(localSubFaceIdx, selfIsUpstream);
                return upwindScheme_.tvd(upwindMomenta, distances, selfIsUpstream, upwindScheme_.tvdApproach());
            }
        }

        // otherwise apply first order upwind scheme
        const auto upwindMomenta = getLateralFirstOrderUpwindMomenta_(insideDensity, outsideDensity, selfIsUpstream, localSubFaceIdx, currentScvfBoundaryTypes, lateralFaceBoundaryTypes);
        return upwindScheme_.upwind(upwindMomenta[0], upwindMomenta[1]);
    }

private:
    /*!
     * \brief Returns whether or not the face in question is far enough from the wall to handle higher order methods.
     *
     *        Evaluates which face is upstream.
     *        If the face is upstream, and the scvf has a forward neighbor, higher order methods are possible.
     *        If the face is downstream, and the scvf has a backwards neighbor, higher order methods are possible.
     *        Otherwise, higher order methods are not possible.
     */
    bool canDoFrontalSecondOrder_(bool selfIsUpstream) const
    {
        // Depending on selfIsUpstream I have to check if I have a forward or a backward neighbor to retrieve
        return selfIsUpstream ? scvf_.hasForwardNeighbor(0) : scvf_.hasBackwardNeighbor(0);
    }

    /*!
     * \brief Returns an array of the three momenta needed for second order upwinding methods.
     */
    std::array<Scalar, 3> getFrontalSecondOrderUpwindMomenta_(const Scalar density, bool selfIsUpstream) const
    {
        static_assert(useHigherOrder, "Should only be reached if higher order methods are enabled");
        const Scalar momentumSelf = faceVars_.velocitySelf() * density;
        const Scalar momentumOpposite = faceVars_.velocityOpposite() * density;
        if (selfIsUpstream)
            return {momentumOpposite, momentumSelf, faceVars_.velocityForward(0)*density};
        else
            return {momentumSelf, momentumOpposite, faceVars_.velocityBackward(0)*density};
    }

    /*!
     * \brief Returns an array of the two momenta needed for first order upwinding method
     */
    std::array<Scalar, 2> getFrontalFirstOrderUpwindMomenta_(const Scalar density, bool selfIsUpstream) const
    {
        const Scalar momentumSelf = faceVars_.velocitySelf() * density;
        const Scalar momentumOpposite = faceVars_.velocityOpposite() * density;
        if (selfIsUpstream)
            return {momentumOpposite, momentumSelf};
        else
            return {momentumSelf, momentumOpposite};
    }

    /*!
     * \brief Returns an array of distances needed for non-uniform higher order upwind schemes
     * Depending on selfIsUpstream the downstream and the (up)upstream distances are saved.
     * distances {upstream to downstream distance, up-upstream to upstream distance, downstream staggered cell size}
     */
    std::array<Scalar, 3> getFrontalDistances_(const bool selfIsUpstream) const
    {
        static_assert(useHigherOrder, "Should only be reached if higher order methods are enabled");
        if (selfIsUpstream)
        {
            std::array<Scalar, 3> distances;
            distances[0] = scvf_.selfToOppositeDistance();
            distances[1] = scvf_.axisData().inAxisForwardDistances[0];
            distances[2] = 0.5 * (scvf_.axisData().selfToOppositeDistance + scvf_.axisData().inAxisBackwardDistances[0]);
            return distances;
        }
        else
        {
            std::array<Scalar, 3> distances;
            distances[0] = scvf_.selfToOppositeDistance();
            distances[1] = scvf_.axisData().inAxisBackwardDistances[0];
            distances[2] = 0.5 * (scvf_.axisData().selfToOppositeDistance + scvf_.axisData().inAxisForwardDistances[0]);
            return distances;
        }
    }

    /*!
     * \brief Check if a second order approximation for the lateral part of the advective term can be used
     *
     * This helper function checks if the scvf of interest is not too near to the
     * boundary so that a dof upstream with respect to the upstream dof is available.
     */
    bool canDoLateralSecondOrder_(const bool selfIsUpstream, const int localSubFaceIdx) const
    {
        // If the self velocity is upstream, the downstream velocity can be assigned or retrieved
        // from the boundary, even if there is no parallel neighbor.
        // If the self velocity is downstream and there is no parallel neighbor I cannot use a second order approximation.
        return selfIsUpstream || scvf_.hasParallelNeighbor(localSubFaceIdx, 0);
    }

    /*!
     * \brief Returns an array of momenta needed for higher order or calls a function to return an array for basic upwinding methods.
     * TODO: In order to get a second order momentum upwind scheme for compressible flow the densities have to be evaluated
     * at the same integration points / positions as the velocities. The currently implementation just takes the closest upwind density
     * to compute the momentum as a crude approximation.
     *
     * ------------
     * |     xxxx o
     * |     xxxx o
     * |     xxxx o
     * -----------*bbbbbbbbbbb
     * |     yyyy o zzzz     |
     * |     yyyy o zzzz     |
     * |     yyyy o zzzz     |
     * -----------------------
     * If scvf_ is touching a corner, at which there is a Dirichlet condition for face b (half-control
     * volumes x, y and z), the transported velocity at * is given. No upwinding or
     * higher-order approximation for the velocity at * are required. This also means the transported
     * velocity at * is the same for the half-control volumes y and z.
     *
     * ------------
     * |          |
     * |          |
     * |          |
     * -----------------------
     * |     xxxx o wwww     |
     * |     xxxx o wwww     |
     * |     xxxx o wwww     |
     * ------+++++*~~~~~------
     * |     yyyy o zzzz     |
     * |     yyyy o zzzz     |
     * |     yyyy o zzzz     |
     * -----------------------
     * If the flux over face + is calculated and a corner occurs a bit further away (here upper right),
     * no special treatment of the corner geometry is provided. This means, the transported velocity at
     * star (*) is obtained with an upwind scheme. In particularly, this also means
     * that while x and y see the same velocity * for the flux over face x (continuity OK), z sees
     * another velocity * for the flux over face ~ (continuity still OK, as only z and w have to use the
     * same velocity at *).
     */
    std::array<Scalar, 3> getLateralSecondOrderUpwindMomenta_(const Scalar insideDensity,
                                                              const Scalar outsideDensity,
                                                              const bool selfIsUpstream,
                                                              const int localSubFaceIdx,
                                                              const std::optional<BoundaryTypes>& currentScvfBoundaryTypes,
                                                              const std::optional<BoundaryTypes>& lateralFaceBoundaryTypes) const
    {
        static_assert(useHigherOrder, "Should only be reached if higher order methods are enabled");

        // If the lateral face lies on a boundary, we assume that the parallel velocity on the boundary is actually known,
        // thus we always use this value for the computation of the transported momentum.
        if (!scvf_.hasParallelNeighbor(localSubFaceIdx, 0) || scvf_.hasHalfParallelNeighbor(localSubFaceIdx) || scvf_.hasCornerParallelNeighbor(localSubFaceIdx))
        {
            if ((scvf_.hasHalfParallelNeighbor(localSubFaceIdx) || scvf_.hasCornerParallelNeighbor(localSubFaceIdx)) && dirichletParallelNeighbor_(localSubFaceIdx))
            {
                const Scalar boundaryVelocity = getParallelVelocityFromCorner_(localSubFaceIdx);
                const Scalar boundaryMomentum = boundaryVelocity*insideDensity;

                return {boundaryMomentum, boundaryMomentum, boundaryMomentum};
            }
            else if (!scvf_.hasParallelNeighbor(localSubFaceIdx, 0))
            {
                const Scalar boundaryVelocity = getParallelVelocityFromBoundary_(element_, scvf_, faceVars_, currentScvfBoundaryTypes, lateralFaceBoundaryTypes, localSubFaceIdx);
                const Scalar boundaryMomentum = boundaryVelocity*insideDensity;

                return {boundaryMomentum, boundaryMomentum, boundaryMomentum};
            }
        }

        if (selfIsUpstream)
        {
            std::array<Scalar, 3> momenta;
            momenta[1] = faceVars_.velocitySelf()*insideDensity;
            momenta[0] = faceVars_.velocityParallel(localSubFaceIdx, 0)*insideDensity;

            // The local index of the faces that is opposite to localSubFaceIdx
            const int oppositeSubFaceIdx = localSubFaceIdx % 2 ? localSubFaceIdx - 1 : localSubFaceIdx + 1;

            // The "upstream-upstream" velocity is retrieved from the other parallel neighbor or from the boundary.
            if (scvf_.hasParallelNeighbor(oppositeSubFaceIdx, 0))
                momenta[2] = faceVars_.velocityParallel(oppositeSubFaceIdx, 0)*insideDensity;
            else
                momenta[2] = getParallelVelocityFromOppositeBoundary_(element_, scvf_, faceVars_, currentScvfBoundaryTypes, oppositeSubFaceIdx)*insideDensity;
            return momenta;
        }
        else
        {
            std::array<Scalar, 3> momenta;
            momenta[0] = faceVars_.velocitySelf()*outsideDensity;
            momenta[1] = faceVars_.velocityParallel(localSubFaceIdx, 0)*outsideDensity;

            // If there is another parallel neighbor I can assign the "upstream-upstream" velocity, otherwise I retrieve it from the boundary.
            if (scvf_.hasParallelNeighbor(localSubFaceIdx, 1))
                momenta[2] = faceVars_.velocityParallel(localSubFaceIdx, 1)*outsideDensity;
            else
            {
                const auto& lateralFace = fvGeometry_.scvf(scvf_.insideScvIdx(), scvf_.pairData(localSubFaceIdx).localLateralFaceIdx);
                const auto& elementParallel = fvGeometry_.gridGeometry().element(lateralFace.outsideScvIdx());
                const auto& firstParallelScvf = fvGeometry_.scvf(lateralFace.outsideScvIdx(), scvf_.localFaceIdx());
                const auto& problem = elemVolVars_.gridVolVars().problem();
                const auto& boundaryTypes = problem.boundaryTypes(elementParallel, firstParallelScvf);
                momenta[2] = getParallelVelocityFromOppositeBoundary_(elementParallel, firstParallelScvf, elemFaceVars_[firstParallelScvf], boundaryTypes, localSubFaceIdx)*outsideDensity;
            }
            return momenta;
        }
    }

    /*!
     * \brief Returns an array of momenta needed for basic upwinding methods.
     */
    std::array<Scalar, 2> getLateralFirstOrderUpwindMomenta_(const Scalar insideDensity,
                                                             const Scalar outsideDensity,
                                                             const bool selfIsUpstream,
                                                             const int localSubFaceIdx,
                                                             const std::optional<BoundaryTypes>& currentScvfBoundaryTypes,
                                                             const std::optional<BoundaryTypes>& lateralFaceBoundaryTypes) const
    {
        // If the lateral face lies on a boundary, we assume that the parallel velocity on the boundary is actually known,
        // thus we always use this value for the computation of the transported momentum.
        if ((scvf_.hasHalfParallelNeighbor(localSubFaceIdx) || scvf_.hasCornerParallelNeighbor(localSubFaceIdx)) && dirichletParallelNeighbor_(localSubFaceIdx))
        {
            const Scalar boundaryVelocity =  getParallelVelocityFromCorner_(localSubFaceIdx);
            const Scalar boundaryMomentum = boundaryVelocity*outsideDensity;
            return {boundaryMomentum, boundaryMomentum};
        }
        else if (!scvf_.hasParallelNeighbor(localSubFaceIdx, 0))
        {
            const Scalar boundaryVelocity = getParallelVelocityFromBoundary_(element_, scvf_, faceVars_, currentScvfBoundaryTypes, lateralFaceBoundaryTypes, localSubFaceIdx);
            const Scalar boundaryMomentum = boundaryVelocity*insideDensity;
            return {boundaryMomentum, boundaryMomentum};
        }

        const Scalar momentumParallel = faceVars_.velocityParallel(localSubFaceIdx, 0)*outsideDensity;
        const Scalar momentumSelf = faceVars_.velocitySelf()*insideDensity;
        if (selfIsUpstream)
            return {momentumParallel, momentumSelf};
        else
            return {momentumSelf, momentumParallel};
    }

    /*!
     * \brief Returns an array of distances needed for non-uniform higher order upwind schemes
     * computes lateral distances {upstream to downstream distance, up-upstream to upstream distance, downstream staggered cell size}
     */
    std::array<Scalar, 3> getLateralDistances_(const int localSubFaceIdx, const bool selfIsUpstream) const
    {
        static_assert(useHigherOrder, "Should only be reached if higher order methods are enabled");
        if (selfIsUpstream)
        {
            // The local index of the faces that is opposite to localSubFaceIdx
            const int oppositeSubFaceIdx = localSubFaceIdx % 2 ? localSubFaceIdx - 1 : localSubFaceIdx + 1;

            std::array<Scalar, 3> distances;
            distances[0] = scvf_.parallelDofsDistance(localSubFaceIdx, 0);
            distances[1] = scvf_.parallelDofsDistance(oppositeSubFaceIdx, 0);
            if (scvf_.hasParallelNeighbor(localSubFaceIdx, 0))
                distances[2] = scvf_.pairData(localSubFaceIdx).parallelCellWidths[0];
            else
                distances[2] = scvf_.area() / 2.0;
            return distances;
        }
        else
        {
            std::array<Scalar, 3> distances;
            distances[0] = scvf_.parallelDofsDistance(localSubFaceIdx, 0);
            distances[1] = scvf_.parallelDofsDistance(localSubFaceIdx, 1);
            distances[2] = scvf_.pairData(localSubFaceIdx).parallelCellWidths[0];
            return distances;
        }
    }

    /*!
     * \brief Return the outer parallel velocity for normal faces that are on the boundary and therefore have no neighbor.
     * Calls the problem to retrieve a fixed value set on the boundary.
     */
    Scalar getParallelVelocityFromBoundary_(const Element& element,
                                            const SubControlVolumeFace& scvf,
                                            const FaceVariables& faceVars,
                                            const std::optional<BoundaryTypes>& currentScvfBoundaryTypes,
                                            const std::optional<BoundaryTypes>& lateralFaceBoundaryTypes,
                                            const int localSubFaceIdx) const
    {
        // If there is a Dirichlet condition for the pressure we assume zero gradient for the velocity,
        // so the velocity at the boundary equal to that on the scvf.
        const bool useZeroGradient = lateralFaceBoundaryTypes && (lateralFaceBoundaryTypes->isSymmetry() || lateralFaceBoundaryTypes->isDirichlet(Indices::pressureIdx));
        if (useZeroGradient)
            return faceVars.velocitySelf();

        const bool lateralFaceHasBJS = lateralFaceBoundaryTypes && lateralFaceBoundaryTypes->isBeaversJoseph(Indices::velocity(scvf.directionIndex()));
        if (lateralFaceHasBJS)
            return VelocityGradients::beaversJosephVelocityAtLateralScvf(elemVolVars_.gridVolVars().problem(), element, fvGeometry_, scvf, faceVars,
                                                                         currentScvfBoundaryTypes, lateralFaceBoundaryTypes, localSubFaceIdx);

        const bool lateralFaceHasDirichletVelocity = lateralFaceBoundaryTypes && lateralFaceBoundaryTypes->isDirichlet(Indices::velocity(scvf.directionIndex()));
        if (lateralFaceHasDirichletVelocity)
        {
            //     ________________
            //     ---------------o                 || frontal face of staggered half-control-volume
            //     |      ||      # current scvf    #  ghostFace of interest, lies on boundary
            //     |      ||      #                 x  dof position
            //     |      ||      x~~~~> vel.Self   -- element boundaries
            //     |      ||      #                 __ domain boundary
            //     |      ||      #                 o  position at which the boundary conditions will be evaluated
            //     ---------------#

            const auto& lateralFace = fvGeometry_.scvf(scvf.insideScvIdx(), scvf.pairData(localSubFaceIdx).localLateralFaceIdx);
            const auto ghostFace = makeStaggeredBoundaryFace(lateralFace, scvf.pairData(localSubFaceIdx).lateralStaggeredFaceCenter);
            const auto& problem = elemVolVars_.gridVolVars().problem();
            return problem.dirichlet(element, ghostFace)[Indices::velocity(scvf.directionIndex())];
        }

        // Neumann conditions are not well implemented
        DUNE_THROW(Dune::InvalidStateException, "Something went wrong with the boundary conditions for the momentum equations at global position "
                    << fvGeometry_.scvf(scvf.insideScvIdx(), scvf.pairData(localSubFaceIdx).localLateralFaceIdx).center());
    }

    /*!
     * \brief Return a velocity value from a boundary for which the boundary conditions have to be checked.
     */
    Scalar getParallelVelocityFromOppositeBoundary_(const Element& element,
                                                    const SubControlVolumeFace& scvf,
                                                    const FaceVariables& faceVars,
                                                    const std::optional<BoundaryTypes>& currentScvfBoundaryTypes,
                                                    const int localOppositeSubFaceIdx) const
    {
        // A ghost subface at the boundary is created, featuring the location of the sub face's center
        const auto& lateralOppositeScvf = fvGeometry_.scvf(scvf.insideScvIdx(), scvf.pairData(localOppositeSubFaceIdx).localLateralFaceIdx);
        GlobalPosition center = scvf.pairData(localOppositeSubFaceIdx).lateralStaggeredFaceCenter + lateralOppositeScvf.center();
        center *= 0.5;

        //          lateral face               #  lateralFace currently being assembled
        //    --------########                 || frontal face of staggered half-control-volume
        //    |      ||      | current scvf    %  lateralOppositeBoundaryFace of interest, lies on boundary
        //    |      ||      |                 x  dof position
        //    |      ||      x~~~~> vel.Self   -- element boundaries
        //    |      ||      |                 __ domain boundary
        //    |      ||      |                 o  position at which the boundary conditions will be evaluated
        //    --------%%%c%%%o
        //    ________________                 c  center of opposite boundary face


        // Get the boundary types of the lateral opposite boundary face
        const auto& problem = elemVolVars_.gridVolVars().problem();
        const auto lateralOppositeBoundaryFace = makeStaggeredBoundaryFace(lateralOppositeScvf, center);
        const auto lateralOppositeFaceBoundaryTypes = problem.boundaryTypes(element, lateralOppositeBoundaryFace);
        return getParallelVelocityFromBoundary_(element, scvf, faceVars, currentScvfBoundaryTypes, lateralOppositeFaceBoundaryTypes, localOppositeSubFaceIdx);
    }

    /*!
     * \brief Returns the boundary subcontrolvolumeface in a corner geometry.
     *
     * ------------
     * |     xxxx o
     * |     xxxx o
     * |     xxxx o
     * ------------bbbbbbbbbbb
     * |     yyyy o          |
     * |     yyyy o          |
     * |     yyyy o          |
     * -----------------------
     *
     * This function will be entered in such a corner geometry (there is no cell in the upper right, --- and |
     * stand for the grid cells). The scvf_ will be one of the two ones denoted by o (upper one
     * hasCornerParallelNeighbor, lower one hasHalfParallelNeighbor). x and y are the two possible corresponding
     * half-control volumes. In both cases, the returned boundaryScvf is the one marked by b. It needs to be the
     * same boundaryScvf returned for the sake of flux continuity.
     */
    const SubControlVolumeFace& boundaryScvf_(const int localSubFaceIdx) const
    {
        if (scvf_.hasHalfParallelNeighbor(localSubFaceIdx))
        {
            return fvGeometry_.scvf(scvf_.outsideScvIdx(), scvf_.pairData(localSubFaceIdx).localLateralFaceIdx);
        }
        else if (scvf_.hasCornerParallelNeighbor(localSubFaceIdx))
        {
           /*
            *      ------------
            *      | xxxxxxxx o
            *      | xxxxxxxx o
            *      | xxxxxxxx o
            *      lllllllllll-bbbbbbbbbbb
            *      | yyyyyyyy p          |
            *      | yyyyyyyy p          |
            *      | yyyyyyyy p          |
            *      -----------------------
            *
            * o: scvf_, l: lateralFace, p: parallelFace, b: returned scvf, x: scvf_ inside scv, y: lateralFace
            * outside scv
            */
            const SubControlVolumeFace& lateralFace = fvGeometry_.scvf(scvf_.insideScvIdx(), scvf_.pairData(localSubFaceIdx).localLateralFaceIdx);
            const SubControlVolumeFace& parallelFace = fvGeometry_.scvf(lateralFace.outsideScvIdx(), scvf_.localFaceIdx());

            const auto& localLateralIdx = scvf_.pairData(localSubFaceIdx).localLateralFaceIdx;
            const auto& localLateralOppositeIdx =  (localLateralIdx % 2) ? (localLateralIdx - 1) : (localLateralIdx + 1);

            return fvGeometry_.scvf(parallelFace.outsideScvIdx(), localLateralOppositeIdx);
        }
        else
        {
            DUNE_THROW(Dune::InvalidStateException, "The function boundaryScvf_ should only be called when hasHalfParallelNeighbor or hasCornerParallelNeighbor.");
        }
    }

   /*!
    * \brief Gets the boundary element in a corner geometry.
    *
    * ------------
    * |     xxxx o
    * |     xxxx o
    * |     xxxx o
    * -----------------------
    * |     yyyy o bbbbbbbb |
    * |     yyyy o bbbbbbbb |
    * |     yyyy o bbbbbbbb |
    * -----------------------
    *
    * This function will be entered in such a corner geometry (there is no cell in the upper right, --- and |
    * stand for the grid cells). The scvf_ will be one of the two ones denoted by o (upper one
    * hasCornerParallelNeighbor, lower one hasHalfParallelNeighbor). x and y are the two possible corresponding
    * half-control volumes. In both cases, the returned boundaryElement is the one marked by b.  It needs to be
    * the same boundaryScvf returned for the sake of flux continuity.
    */
    Element boundaryElement_(const int localSubFaceIdx) const
    {
        if (scvf_.hasHalfParallelNeighbor(localSubFaceIdx))
        {
            return fvGeometry_.gridGeometry().element(scvf_.outsideScvIdx());
        }
        else if (scvf_.hasCornerParallelNeighbor(localSubFaceIdx))
        {
            const SubControlVolumeFace& lateralFace = fvGeometry_.scvf(scvf_.insideScvIdx(), scvf_.pairData(localSubFaceIdx).localLateralFaceIdx);
            const SubControlVolumeFace& parallelFace = fvGeometry_.scvf(lateralFace.outsideScvIdx(), scvf_.localFaceIdx());

            return fvGeometry_.gridGeometry().element(parallelFace.outsideScvIdx());
        }
        else
        {
            DUNE_THROW(Dune::InvalidStateException, "When entering boundaryElement_ scvf_ should have either hasHalfParallelNeighbor or hasCornerParallelNeighbor true. Not the case here.");
        }
    }

   /*!
    * \brief Sets the bools hasDirichletCornerParallelNeighbor and hasDirichletHalfParallelNeighbor.
    *
    * ------------
    * |     xxxx o
    * |     xxxx o
    * |     xxxx o
    * ------------bbbbbbbbbbb
    * |     yyyy o          |
    * |     yyyy o boundary |
    * |     yyyy o  element |
    * -----------------------
    *
    * This function will be entered in such a corner geometry (there is no cell in the upper right, --- and |
    * stand for the grid cells). The scvf_ will be one of the two ones denoted by o (upper one
    * hasCornerParallelNeighbor, lower one hasHalfParallelNeighbor). x and y are the two possible corresponding
    * half-control volumes. In both cases, we check if the face bbb, part of the edge of element boundaryElement,
    * is a Dirichlet boundary.
    */
    bool dirichletParallelNeighbor_(const int localSubFaceIdx) const
    {
        const auto& problem = elemVolVars_.gridVolVars().problem();
        const Element& boundaryElement = boundaryElement_(localSubFaceIdx);
        const SubControlVolumeFace& boundaryScvf = boundaryScvf_(localSubFaceIdx);

        return problem.boundaryTypes(boundaryElement, boundaryScvf).isDirichlet(Indices::velocity(scvf_.directionIndex()));
    }

   /*!
    * \brief Gets the parallel velocity from a corner geometry.
    *
    * ------------
    * |     xxxx o
    * |     xxxx o
    * |     xxxx o
    * -----------*-----------
    * |     yyyy o          |
    * |     yyyy o          |
    * |     yyyy o          |
    * -----------------------
    *
    * This function will be entered in such a corner geometry (there is no cell in the upper right, --- and |
    * stand for the grid cells). The scvf_ will be one of the two ones denoted by o (upper one
    * hasCornerParallelNeighbor, lower one hasHalfParallelNeighbor). x and y are the two possible corresponding
    * half-control volumes. In both cases, the returned velocity is situated in the corner (*).
    */
    Scalar getParallelVelocityFromCorner_(const int localSubFaceIdx) const
    {
        const auto& problem = elemVolVars_.gridVolVars().problem();
        const Element& boundaryElement = boundaryElement_(localSubFaceIdx);
        const SubControlVolumeFace& boundaryScvf = boundaryScvf_(localSubFaceIdx);
        const auto ghostFace = makeStaggeredBoundaryFace(boundaryScvf, scvf_.pairData(localSubFaceIdx).lateralStaggeredFaceCenter);

        return problem.dirichlet(boundaryElement, ghostFace)[Indices::velocity(scvf_.directionIndex())];
    }

    const Element& element_;
    const FVElementGeometry& fvGeometry_;
    const SubControlVolumeFace& scvf_;
    const ElementFaceVariables& elemFaceVars_;
    const FaceVariables& faceVars_;
    const ElementVolumeVariables& elemVolVars_;
    const UpwindScheme& upwindScheme_;
};

} // end namespace Dumux

#endif
