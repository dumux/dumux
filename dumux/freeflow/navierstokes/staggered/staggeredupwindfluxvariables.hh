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
 * \copydoc Dumux::StaggeredUpwindFluxVariables
 */
#ifndef DUMUX_NAVIERSTOKES_STAGGERED_UPWINDVARIABLES_HH
#define DUMUX_NAVIERSTOKES_STAGGERED_UPWINDVARIABLES_HH

#include <array>
#include <optional>

#include <dumux/common/math.hh>
#include <dumux/common/exceptions.hh>
#include <dumux/common/parameters.hh>
#include <dumux/common/properties.hh>

#include <dumux/discretization/method.hh>
#include <dumux/freeflow/staggeredupwindmethods.hh>
#include "velocitygradients.hh"

namespace Dumux {

/*!
 * \ingroup NavierStokesModel
 * \brief The upwinding variables class for the Navier-Stokes model using the staggered grid discretization.
 */
template<class TypeTag, int upwindSchemeOrder>
class StaggeredUpwindFluxVariables
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
    using SubControlVolume = typename FVElementGeometry::SubControlVolume;
    using SubControlVolumeFace = typename FVElementGeometry::SubControlVolumeFace;
    using CellCenterPrimaryVariables = GetPropType<TypeTag, Properties::CellCenterPrimaryVariables>;
    using FacePrimaryVariables = GetPropType<TypeTag, Properties::FacePrimaryVariables>;
    using BoundaryTypes = GetPropType<TypeTag, Properties::BoundaryTypes>;
    using VelocityGradients = StaggeredVelocityGradients<Scalar, GridGeometry, BoundaryTypes, Indices>;

    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;
    static constexpr bool useHigherOrder = upwindSchemeOrder > 1;

public:
    /*!
     * \brief Returns the momentum in the frontal directon.
     *
     *        Checks if the model has higher order methods enabled and if the scvf in
     *        question is far enough from the boundary such that higher order methods can be employed.
     *        Then the corresponding set of momenta are collected and the prescribed
     *        upwinding method is used to calculate the momentum.
     */
    static FacePrimaryVariables computeUpwindedFrontalMomentum(const SubControlVolumeFace& scvf,
                                                               const ElementFaceVariables& elemFaceVars,
                                                               const ElementVolumeVariables& elemVolVars,
                                                               const GridFluxVariablesCache& gridFluxVarsCache,
                                                               const Scalar transportingVelocity)
    {
        const bool canHigherOrder = canFrontalSecondOrder_(scvf, transportingVelocity);
        const auto upwindingMomenta = getFrontalUpwindingMomenta_(scvf, elemFaceVars, elemVolVars[scvf.insideScvIdx()].density(),
                                                                  transportingVelocity, canHigherOrder);
        return doFrontalMomentumUpwinding_(scvf, upwindingMomenta, transportingVelocity, gridFluxVarsCache, canHigherOrder);
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
    static FacePrimaryVariables computeUpwindedLateralMomentum(const Problem& problem,
                                                               const FVElementGeometry& fvGeometry,
                                                               const Element& element,
                                                               const SubControlVolumeFace& scvf,
                                                               const ElementVolumeVariables& elemVolVars,
                                                               const FaceVariables& faceVars,
                                                               const GridFluxVariablesCache& gridFluxVarsCache,
                                                               const int localSubFaceIdx,
                                                               const std::optional<BoundaryTypes>& currentScvfBoundaryTypes,
                                                               const std::optional<BoundaryTypes>& lateralFaceBoundaryTypes)
    {
        // Check whether the own or the neighboring element is upstream.
        const auto eIdx = scvf.insideScvIdx();
        const auto& lateralFace = fvGeometry.scvf(eIdx, scvf.pairData(localSubFaceIdx).localLateralFaceIdx);
        // Get the transporting velocity, located at the scvf perpendicular to the current scvf where the dof
        // of interest is located.
        const Scalar transportingVelocity = faceVars.velocityLateralInside(localSubFaceIdx);

        // Check whether the own or the neighboring element is upstream.
        const bool selfIsUpstream = ( lateralFace.directionSign() == sign(transportingVelocity) );
        const bool canHigherOrder = canLateralSecondOrder_(scvf, selfIsUpstream, localSubFaceIdx);
        const auto parallelUpwindingMomenta = getLateralUpwindingMomenta_(problem, fvGeometry, element, scvf, elemVolVars, faceVars,
                                                                          transportingVelocity, localSubFaceIdx, currentScvfBoundaryTypes,
                                                                          lateralFaceBoundaryTypes, canHigherOrder);
        return doLateralMomentumUpwinding_(fvGeometry, scvf, parallelUpwindingMomenta, transportingVelocity,
                                           localSubFaceIdx, gridFluxVarsCache, canHigherOrder);
    }

private:

    /*!
     * \brief Returns whether or not the face in question is far enough from the wall to handle higher order methods.
     *
     *        Evaluates which face is upstream.
     *        If the face is upstream, and the scvf has a forward neighbor, higher order methods are possible.
     *        If the face is downstream, and the scvf has a backwards neighbor, higher order methods are possible.
     *        Otherwise, higher order methods are not possible.
     *        If higher order methods are not prescribed, this function will return false.
     *
     * \param ownScvf the sub control volume face
     * \param transportingVelocity The average of the self and opposite velocities.
     */
    static bool canFrontalSecondOrder_([[maybe_unused]] const SubControlVolumeFace& ownScvf,
                                       [[maybe_unused]] const Scalar transportingVelocity)
    {
        if constexpr (useHigherOrder)
        {
            // Depending on selfIsUpstream I have to check if I have a forward or a backward neighbor to retrieve
            const bool selfIsUpstream = ownScvf.directionSign() != sign(transportingVelocity);
            return selfIsUpstream ? ownScvf.hasForwardNeighbor(0) : ownScvf.hasBackwardNeighbor(0);
        }
        else
            return false;
    }

    /*!
     * \brief Returns an array of the three momenta needed for higher order upwinding methods.
     *
     * \param scvf The sub control volume face
     * \param elemFaceVars The element face variables
     * \param density The given density \f$\mathrm{[kg/m^3]}\f$
     * \param transportingVelocity The average of the self and opposite velocities.
     */
    static auto getFrontalUpwindingMomenta_(const SubControlVolumeFace& scvf,
                                            const ElementFaceVariables& elemFaceVars,
                                            const Scalar density,
                                            const Scalar transportingVelocity,
                                            [[maybe_unused]] const bool canHigherOrder)
    {
        const bool selfIsUpstream = scvf.directionSign() != sign(transportingVelocity);
        const Scalar momentumSelf = elemFaceVars[scvf].velocitySelf() * density;
        const Scalar momentumOpposite = elemFaceVars[scvf].velocityOpposite() * density;

        if constexpr (useHigherOrder)
        {
            if (canHigherOrder)
                return selfIsUpstream ? std::array<Scalar, 3>{momentumOpposite, momentumSelf,
                                                              elemFaceVars[scvf].velocityForward(0)*density}
                                      : std::array<Scalar, 3>{momentumSelf, momentumOpposite,
                                                              elemFaceVars[scvf].velocityBackward(0)*density};
            else
                return selfIsUpstream ? std::array<Scalar, 3>{momentumOpposite, momentumSelf}
                                      : std::array<Scalar, 3>{momentumSelf, momentumOpposite};
        }
        else
            return selfIsUpstream ? std::array<Scalar, 2>{momentumOpposite, momentumSelf}
                                  : std::array<Scalar, 2>{momentumSelf, momentumOpposite};
    }

    /*!
     * \brief Returns the upwinded momentum for higher order upwind schemes
     *
     * \param scvf The sub control volume face
     * \param momenta The momenta to be upwinded
     * \param transportingVelocity The average of the self and opposite velocities.
     * \param gridFluxVarsCache The grid flux variables cache
     */
    template<class MomentaArray>
    static Scalar doFrontalMomentumUpwinding_([[maybe_unused]] const SubControlVolumeFace& scvf,
                                              const MomentaArray& momenta,
                                              [[maybe_unused]] const Scalar transportingVelocity,
                                              const GridFluxVariablesCache& gridFluxVarsCache,
                                              [[maybe_unused]] const bool canHigherOrder)
    {
        const auto& upwindScheme = gridFluxVarsCache.staggeredUpwindMethods();
        if constexpr (useHigherOrder)
        {
            if (canHigherOrder)
            {
                const bool selfIsUpstream = scvf.directionSign() != sign(transportingVelocity);
                const std::array<Scalar,3> distances = getFrontalDistances_(scvf, selfIsUpstream);
                return upwindScheme.tvd(momenta, distances, selfIsUpstream, upwindScheme.tvdApproach());
            }
            else
                return upwindScheme.upwind(momenta[0], momenta[1]);
        }
        else
            return upwindScheme.upwind(momenta[0], momenta[1]);
    }

    /*!
     * \brief Returns an array of distances needed for non-uniform higher order upwind schemes
     *
     * \param ownScvf The sub control volume face
     * \param selfIsUpstream bool describing upstream face.
     */
    static std::array<Scalar, 3> getFrontalDistances_(const SubControlVolumeFace& ownScvf,
                                                      const bool selfIsUpstream)
    {
        static_assert(useHigherOrder, "Should only be reached if higher order methods are enabled");
        // Depending on selfIsUpstream the downstream and the (up)upstream distances are saved.
        // distances {upstream to downstream distance, up-upstream to upstream distance, downstream staggered cell size}
        std::array<Scalar, 3> distances;

        if (selfIsUpstream)
        {
            distances[0] = ownScvf.selfToOppositeDistance();
            distances[1] = ownScvf.axisData().inAxisForwardDistances[0];
            distances[2] = 0.5 * (ownScvf.axisData().selfToOppositeDistance + ownScvf.axisData().inAxisBackwardDistances[0]);
        }
        else
        {
            distances[0] = ownScvf.selfToOppositeDistance();
            distances[1] = ownScvf.axisData().inAxisBackwardDistances[0];
            distances[2] = 0.5 * (ownScvf.axisData().selfToOppositeDistance + ownScvf.axisData().inAxisForwardDistances[0]);
        }
        return distances;
    }

    /*!
     * \brief Check if a second order approximation for the lateral part of the advective term can be used
     *
     * This helper function checks if the scvf of interest is not too near to the
     * boundary so that a dof upstream with respect to the upstream dof is available.
     *
     * If higher order methods are not prescribed, false is returned.
     *
     * \param ownScvf The SubControlVolumeFace we are considering
     * \param selfIsUpstream @c true if the velocity at ownScvf is upstream wrt the transporting velocity
     * \param localSubFaceIdx The local subface index
     */
    static bool canLateralSecondOrder_(const SubControlVolumeFace& ownScvf,
                                       [[maybe_unused]] const bool selfIsUpstream,
                                       [[maybe_unused]] const int localSubFaceIdx)
    {
        if constexpr (useHigherOrder)
        {
            if (selfIsUpstream)
            {
                // The self velocity is upstream. The downstream velocity can be assigned or retrieved
                // from the boundary, even if there is no parallel neighbor.
                return true;
            }
            else
            {
                // The self velocity is downstream. If there is no parallel neighbor I cannot use a second order approximation.
                return ownScvf.hasParallelNeighbor(localSubFaceIdx, 0);
            }
        }
        else
            return false;
    }

    /*!
     * \brief Returns an array of momenta needed for higher order or calls a function to return an array for basic upwinding methods.
     */
    static auto getLateralUpwindingMomenta_(const Problem& problem,
                                            const FVElementGeometry& fvGeometry,
                                            const Element& element,
                                            const SubControlVolumeFace& ownScvf,
                                            const ElementVolumeVariables& elemVolVars,
                                            const FaceVariables& faceVars,
                                            const Scalar transportingVelocity,
                                            const int localSubFaceIdx,
                                            const std::optional<BoundaryTypes>& currentScvfBoundaryTypes,
                                            const std::optional<BoundaryTypes>& lateralFaceBoundaryTypes,
                                            [[maybe_unused]] const bool canHigherOrder)
    {
        // Check whether the own or the neighboring element is upstream.
        const SubControlVolumeFace& lateralFace = fvGeometry.scvf(ownScvf.insideScvIdx(), ownScvf.pairData(localSubFaceIdx).localLateralFaceIdx);

        // Get the volume variables of the own and the neighboring element
        const auto& insideVolVars = elemVolVars[lateralFace.insideScvIdx()];
        const auto& outsideVolVars = elemVolVars[lateralFace.outsideScvIdx()];

        // Check whether the own or the neighboring element is upstream.
        const bool selfIsUpstream = lateralFace.directionSign() == sign(transportingVelocity);

        if constexpr (useHigherOrder)
        {
            if (canHigherOrder)
            {
                // If the lateral face lies on a boundary, we assume that the parallel velocity on the boundary is actually known,
                // thus we always use this value for the computation of the transported momentum.
                if (!ownScvf.hasParallelNeighbor(localSubFaceIdx, 0))
                {
                    const Scalar boundaryMomentum = getParallelVelocityFromBoundary_(problem, element, fvGeometry, ownScvf,
                                                                                    faceVars, currentScvfBoundaryTypes, lateralFaceBoundaryTypes,
                                                                                    localSubFaceIdx) * insideVolVars.density();

                    return std::array<Scalar, 3>{boundaryMomentum, boundaryMomentum, boundaryMomentum};
                }

                std::array<Scalar, 3> momenta;
                if (selfIsUpstream)
                {
                    momenta[1] = faceVars.velocitySelf() * insideVolVars.density();

                    if (ownScvf.hasParallelNeighbor(localSubFaceIdx, 0))
                        momenta[0] = faceVars.velocityParallel(localSubFaceIdx, 0) * insideVolVars.density();
                    else
                        momenta[0] = getParallelVelocityFromBoundary_(problem, element, fvGeometry, ownScvf,
                                                                    faceVars, currentScvfBoundaryTypes, lateralFaceBoundaryTypes,
                                                                    localSubFaceIdx) * insideVolVars.density();

                    // The local index of the faces that is opposite to localSubFaceIdx
                    const int oppositeSubFaceIdx = localSubFaceIdx % 2 ? localSubFaceIdx - 1 : localSubFaceIdx + 1;

                    // The "upstream-upstream" velocity is retrieved from the other parallel neighbor or from the boundary.
                    if (ownScvf.hasParallelNeighbor(oppositeSubFaceIdx, 0))
                        momenta[2] = faceVars.velocityParallel(oppositeSubFaceIdx, 0) * insideVolVars.density();
                    else
                        momenta[2] = getParallelVelocityFromOppositeBoundary_(problem, element, fvGeometry, ownScvf,
                                                                            faceVars, currentScvfBoundaryTypes,
                                                                            oppositeSubFaceIdx) * insideVolVars.density();
                }
                else
                {
                    momenta[0] = faceVars.velocitySelf() * outsideVolVars.density();
                    momenta[1] = faceVars.velocityParallel(localSubFaceIdx, 0) * outsideVolVars.density();

                    // If there is another parallel neighbor I can assign the "upstream-upstream" velocity, otherwise I retrieve it from the boundary.
                    if (ownScvf.hasParallelNeighbor(localSubFaceIdx, 1))
                        momenta[2] = faceVars.velocityParallel(localSubFaceIdx, 1) * outsideVolVars.density();
                    else
                    {
                        const Element& elementParallel = fvGeometry.gridGeometry().element(lateralFace.outsideScvIdx());
                        const SubControlVolumeFace& firstParallelScvf = fvGeometry.scvf(lateralFace.outsideScvIdx(), ownScvf.localFaceIdx());

                        momenta[2] = getParallelVelocityFromOppositeBoundary_(problem, elementParallel, fvGeometry, firstParallelScvf,
                                                                            faceVars, problem.boundaryTypes(elementParallel, firstParallelScvf),
                                                                            localSubFaceIdx) * outsideVolVars.density();
                    }
                }
                return momenta;
            }
            else
                return getFirstOrderLateralUpwindingMomenta_(problem, element, fvGeometry, ownScvf, faceVars,
                                                             currentScvfBoundaryTypes, lateralFaceBoundaryTypes,
                                                             localSubFaceIdx, selfIsUpstream, insideVolVars, outsideVolVars);
        }
        else
            return getFirstOrderLateralUpwindingMomenta_(problem, element, fvGeometry, ownScvf, faceVars,
                                                         currentScvfBoundaryTypes, lateralFaceBoundaryTypes,
                                                         localSubFaceIdx, selfIsUpstream, insideVolVars, outsideVolVars);
    }

    /*!
     * \brief Returns an array of momenta needed for basic upwinding methods.
     */
    static auto getFirstOrderLateralUpwindingMomenta_(const Problem& problem,
                                                      const Element& element,
                                                      const FVElementGeometry& fvGeometry,
                                                      const SubControlVolumeFace& ownScvf,
                                                      const FaceVariables& faceVars,
                                                      const std::optional<BoundaryTypes>& currentScvfBoundaryTypes,
                                                      const std::optional<BoundaryTypes>& lateralFaceBoundaryTypes,
                                                      const int localSubFaceIdx,
                                                      const bool selfIsUpstream,
                                                      const VolumeVariables& insideVolVars,
                                                      const VolumeVariables& outsideVolVars)
    {
        const Scalar momentumParallel = ownScvf.hasParallelNeighbor(localSubFaceIdx, 0)
                                      ? faceVars.velocityParallel(localSubFaceIdx, 0) * outsideVolVars.density()
                                      : (getParallelVelocityFromBoundary_(problem, element, fvGeometry, ownScvf,
                                                                          faceVars, currentScvfBoundaryTypes, lateralFaceBoundaryTypes,
                                                                          localSubFaceIdx) * insideVolVars.density());
        // If the lateral face lies on a boundary, we assume that the parallel velocity on the boundary is actually known,
        // thus we always use this value for the computation of the transported momentum.
        if (!ownScvf.hasParallelNeighbor(localSubFaceIdx, 0))
        {
            if constexpr (useHigherOrder)
                return std::array<Scalar, 3>{momentumParallel, momentumParallel};
            else
                return std::array<Scalar, 2>{momentumParallel, momentumParallel};
        }

        const Scalar momentumSelf = faceVars.velocitySelf() * insideVolVars.density();
        if constexpr (useHigherOrder)
            return selfIsUpstream ? std::array<Scalar, 3>{momentumParallel, momentumSelf}
                                  : std::array<Scalar, 3>{momentumSelf, momentumParallel};
        else
            return selfIsUpstream ? std::array<Scalar, 2>{momentumParallel, momentumSelf}
                                  : std::array<Scalar, 2>{momentumSelf, momentumParallel};

    }

    /*!
     * \brief Returns the upwinded momentum for higher order or basic upwind schemes
     *
     * If higher order methods are not enabled, this is fowarded to the Frontal Momentum Upwinding method
     *
     * \param scvf The sub control volume face
     * \param momenta The momenta to be upwinded
     * \param transportingVelocity The average of the self and opposite velocities.
     * \param localSubFaceIdx  The local index of the subface
     * \param gridFluxVarsCache The grid flux variables cache
     */
    template<class MomentaArray>
    static Scalar doLateralMomentumUpwinding_([[maybe_unused]] const FVElementGeometry& fvGeometry,
                                              const SubControlVolumeFace& scvf,
                                              const MomentaArray& momenta,
                                              [[maybe_unused]] const Scalar transportingVelocity,
                                              [[maybe_unused]] const int localSubFaceIdx,
                                              const GridFluxVariablesCache& gridFluxVarsCache,
                                              [[maybe_unused]] const bool canHigherOrder)
    {
        const auto& upwindScheme = gridFluxVarsCache.staggeredUpwindMethods();
        if constexpr (useHigherOrder)
        {
            if (canHigherOrder)
            {
                const auto eIdx = scvf.insideScvIdx();
                const auto& lateralFace = fvGeometry.scvf(eIdx, scvf.pairData(localSubFaceIdx).localLateralFaceIdx);
                const bool selfIsUpstream = ( lateralFace.directionSign() == sign(transportingVelocity) );
                const std::array<Scalar,3> distances = getLateralDistances_(scvf, localSubFaceIdx, selfIsUpstream);
                return upwindScheme.tvd(momenta, distances, selfIsUpstream, upwindScheme.tvdApproach());
            }
            else
                return upwindScheme.upwind(momenta[0], momenta[1]);
        }
        else
            return upwindScheme.upwind(momenta[0], momenta[1]);

    }

    /*!
     * \brief Returns an array of distances needed for non-uniform higher order upwind schemes
     *
     * \param ownScvf The sub control volume face
     * \param localSubFaceIdx  The local index of the subface
     * \param selfIsUpstream bool describing upstream face.
     */
    static std::array<Scalar, 3> getLateralDistances_(const SubControlVolumeFace& ownScvf,
                                                      const int localSubFaceIdx,
                                                      const bool selfIsUpstream)
    {
        static_assert(useHigherOrder, "Should only be reached if higher order methods are enabled");
        // distances {upstream to downstream distance, up-upstream to upstream distance, downstream staggered cell size}
        std::array<Scalar, 3> distances;

        if (selfIsUpstream)
        {
            // The local index of the faces that is opposite to localSubFaceIdx
            const int oppositeSubFaceIdx = localSubFaceIdx % 2 ? localSubFaceIdx - 1 : localSubFaceIdx + 1;

            distances[0] = ownScvf.parallelDofsDistance(localSubFaceIdx, 0);
            distances[1] = ownScvf.parallelDofsDistance(oppositeSubFaceIdx, 0);
            if (ownScvf.hasParallelNeighbor(localSubFaceIdx, 0))
                distances[2] = ownScvf.pairData(localSubFaceIdx).parallelCellWidths[0];
            else
                distances[2] = ownScvf.area() / 2.0;
        }
        else
        {
            distances[0] = ownScvf.parallelDofsDistance(localSubFaceIdx, 0);
            distances[1] = ownScvf.parallelDofsDistance(localSubFaceIdx, 1);
            distances[2] = ownScvf.pairData(localSubFaceIdx).parallelCellWidths[0];
        }

        return distances;
    }

    /*!
     * \brief Return the outer parallel velocity for normal faces that are on the boundary and therefore have no neighbor.
     *
     * Calls the problem to retrieve a fixed value set on the boundary.
     *
     * \param problem The problem
     * \param scvf The SubControlVolumeFace that is normal to the boundary
     * \param velocitySelf the velocity at scvf
     * \param localSubFaceIdx The local index of the face that is on the boundary
     * \param element The element that is on the boundary
     * \param lateralFaceBoundaryTypes stores the type of boundary at the lateral face
     */
    static Scalar getParallelVelocityFromBoundary_(const Problem& problem,
                                                   const Element& element,
                                                   const FVElementGeometry& fvGeometry,
                                                   const SubControlVolumeFace& scvf,
                                                   const FaceVariables& faceVars,
                                                   const std::optional<BoundaryTypes>& currentScvfBoundaryTypes,
                                                   const std::optional<BoundaryTypes>& lateralFaceBoundaryTypes,
                                                   const int localSubFaceIdx)
    {
        // Find out what boundary type is set on the lateral face
        const bool lateralFaceHasBJS = lateralFaceBoundaryTypes && lateralFaceBoundaryTypes->isBeaversJoseph(Indices::velocity(scvf.directionIndex()));
        const bool lateralFaceHasDirichletVelocity = lateralFaceBoundaryTypes && lateralFaceBoundaryTypes->isDirichlet(Indices::velocity(scvf.directionIndex()));

        // If there is a Dirichlet condition for the pressure we assume zero gradient for the velocity,
        // so the velocity at the boundary equal to that on the scvf.
        if (!lateralFaceHasBJS && !lateralFaceHasDirichletVelocity)
            return faceVars.velocitySelf();

        if (lateralFaceHasBJS)
            return VelocityGradients::beaversJosephVelocityAtLateralScvf(problem, element, fvGeometry, scvf,  faceVars,
                                                                         currentScvfBoundaryTypes, lateralFaceBoundaryTypes, localSubFaceIdx);

        else if (lateralFaceHasDirichletVelocity)
        {
            //     ________________
            //     ---------------o                 || frontal face of staggered half-control-volume
            //     |      ||      # current scvf    #  ghostFace of interest, lies on boundary
            //     |      ||      #                 x  dof position
            //     |      ||      x~~~~> vel.Self   -- element boundaries
            //     |      ||      #                 __ domain boundary
            //     |      ||      #                 o  position at which the boundary conditions will be evaluated
            //     ---------------#

            const auto eIdx = fvGeometry.gridGeometry().elementMapper().index(element);
            const auto& lateralFace = fvGeometry.scvf(eIdx, scvf.pairData(localSubFaceIdx).localLateralFaceIdx);

            const auto ghostFace = lateralFace.makeBoundaryFace(scvf.pairData(localSubFaceIdx).lateralStaggeredFaceCenter);
            return problem.dirichlet(element, ghostFace)[Indices::velocity(scvf.directionIndex())];
        }
        else
        {
            // Neumann conditions are not well implemented
            DUNE_THROW(Dune::InvalidStateException, "Something went wrong with the boundary conditions for the momentum equations at global position " << fvGeometry.scvf(scvf.insideScvIdx(), scvf.pairData(localSubFaceIdx).localLateralFaceIdx).center());
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
    static Scalar getParallelVelocityFromOppositeBoundary_(const Problem& problem,
                                                           const Element& element,
                                                           const FVElementGeometry& fvGeometry,
                                                           const SubControlVolumeFace& scvf,
                                                           const FaceVariables& faceVars,
                                                           const std::optional<BoundaryTypes>& currentScvfBoundaryTypes,
                                                           const int localOppositeSubFaceIdx)
    {
        // A ghost subface at the boundary is created, featuring the location of the sub face's center
        const SubControlVolumeFace& lateralOppositeScvf = fvGeometry.scvf(scvf.insideScvIdx(), scvf.pairData(localOppositeSubFaceIdx).localLateralFaceIdx);
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
        const auto lateralOppositeFaceBoundaryTypes = problem.boundaryTypes(element, lateralOppositeScvf.makeBoundaryFace(center));

        return getParallelVelocityFromBoundary_(problem, element, fvGeometry, scvf,
                                                faceVars, currentScvfBoundaryTypes, lateralOppositeFaceBoundaryTypes,
                                                localOppositeSubFaceIdx);
    }
};

} // end namespace Dumux

#endif // DUMUX_NAVIERSTOKES_STAGGERED_UPWINDVARIABLES_HH
