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
 * \copydoc Dumux::MomentumUpwindSchemeHelper
 */
#ifndef DUMUX_MOMENTUM_UPWIND_SCHEME_HELPER_HH
#define DUMUX_MOMENTUM_UPWIND_SCHEME_HELPER_HH

#include <array>
#include <optional>

#include <dumux/common/math.hh>
#include <dumux/common/exceptions.hh>
#include <dumux/common/parameters.hh>
#include <dumux/common/properties.hh>
#include <dumux/common/typetraits/problem.hh>

#include <dumux/discretization/method.hh>
#include <dumux/discretization/facecentered/staggered/upwindgeometryhelper.hh>
#include <dumux/freeflow/staggeredupwindmethods.hh>
#include "velocitygradients.hh"

namespace Dumux {

/*!
 * \ingroup NavierStokesModel
 * \brief The upwinding variables class for the Navier-Stokes model using the staggered grid discretization.
 */
template<class FVElementGeometry, class ElementVolumeVariables>
class FaceCenteredStaggeredUpwindHelper
{
    using Problem = std::decay_t<decltype(std::declval<ElementVolumeVariables>().gridVolVars().problem())>;
    using GridGeometry = typename FVElementGeometry::GridGeometry;
    using GridView = typename GridGeometry::GridView;
    using Element = typename GridView::template Codim<0>::Entity;
    using SubControlVolume = typename FVElementGeometry::SubControlVolume;
    using SubControlVolumeFace = typename FVElementGeometry::SubControlVolumeFace;
    using Scalar = double; // TODO
    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;

    static constexpr int upwindSchemeOrder = GridGeometry::upwindSchemeOrder;
    static_assert(upwindSchemeOrder <= 2, "Not implemented: Order higher than 2!");
    static constexpr bool useHigherOrder = upwindSchemeOrder > 1;
    using UpwindScheme = StaggeredUpwindMethods<Scalar, upwindSchemeOrder>;

public:
    FaceCenteredStaggeredUpwindHelper(const FVElementGeometry& fvGeometry,
                                      const SubControlVolumeFace& scvf,
                                      const ElementVolumeVariables& elemVolVars)
    : element_(fvGeometry.element())
    , fvGeometry_(fvGeometry)
    , problem_(elemVolVars.gridVolVars().problem())
    , scvf_(scvf)
    , elemVolVars_(elemVolVars)
    , upwindScheme_(fvGeometry.upwindMethods())
    , upwindgeometryHelper_(fvGeometry)
    {}

    /*!
     * \brief Returns the momentum in the frontal directon.
     *
     *        Checks if the model has higher order methods enabled and if the scvf in
     *        question is far enough from the boundary such that higher order methods can be employed.
     *        Then the corresponding set of momenta are collected and the prescribed
     *        upwinding method is used to calculate the momentum.
     */
    Scalar upwindFrontalMomentum(const bool selfIsUpstream) const
    {
        const auto density = problem_.density(element_, fvGeometry_, scvf_);

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
    Scalar upwindLateralMomentum(const bool selfIsUpstream) const
    {
        getLateralDistances_(selfIsUpstream);
        // Check whether the own or the neighboring element is upstream.
        // Get the transporting velocity, located at the scvf perpendicular to the current scvf where the dof
        // of interest is located.

        // collect the inside and outside densities
        const auto densityPair = problem_.insideAndOutsideDensity(element_, fvGeometry_, scvf_);

        // for higher order schemes do higher order upwind reconstruction
        if constexpr (useHigherOrder)
        {
            // only second order is implemented so far
            if (canDoLateralSecondOrder_(selfIsUpstream))
            {
                const auto distances = getLateralDistances_(selfIsUpstream);
                const auto upwindMomenta = getLateralSecondOrderUpwindMomenta_(densityPair, selfIsUpstream);
                return upwindScheme_.tvd(upwindMomenta, distances, selfIsUpstream, upwindScheme_.tvdApproach());
            }
        }

        // otherwise apply first order upwind scheme
        const auto upwindMomenta = getLateralFirstOrderUpwindMomenta_(densityPair, selfIsUpstream);
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
        static_assert(useHigherOrder, "Should only be reached if higher order methods are enabled");
        // Depending on selfIsUpstream I have to check if I have a forward or a backward neighbor to retrieve
        return selfIsUpstream ? fvGeometry_.hasForwardNeighbor(scvf_) : fvGeometry_.hasBackwardNeighbor(scvf_);
    }

    /*!
     * \brief Returns an array of the three momenta needed for second order upwinding methods.
     */
    std::array<Scalar, 3> getFrontalSecondOrderUpwindMomenta_(const Scalar density, bool selfIsUpstream) const
    {
        static_assert(useHigherOrder, "Should only be reached if higher order methods are enabled");
        const auto momentumSelf = elemVolVars_[scvf_.insideScvIdx()].velocity() * density;
        const auto momentumOpposite = elemVolVars_[scvf_.outsideScvIdx()].velocity() * density;
        if (selfIsUpstream)
        {
            const auto momentumForward =  elemVolVars_[fvGeometry_.forwardScvIdx(scvf_)].velocity() * density;
            return {momentumOpposite, momentumSelf, momentumForward};
        }
        else
        {
            const auto momentumBackward = elemVolVars_[fvGeometry_.backwardScvIdx(scvf_)].velocity() * density;
            return {momentumSelf, momentumOpposite, momentumBackward};
        }
    }

    /*!
     * \brief Returns an array of the two momenta needed for first order upwinding method
     */
    std::array<Scalar, 2> getFrontalFirstOrderUpwindMomenta_(const Scalar density, bool selfIsUpstream) const
    {
        const auto momentumSelf = elemVolVars_[scvf_.insideScvIdx()].velocity() * density;
        const auto momentumOpposite = elemVolVars_[scvf_.outsideScvIdx()].velocity() * density;
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
            distances[0] = upwindgeometryHelper_.selfToOppositeDistance(scvf_);
            distances[1] = upwindgeometryHelper_.selfToForwardDistance(scvf_);
            distances[2] = 0.5 * (upwindgeometryHelper_.selfToOppositeDistance(scvf_) + upwindgeometryHelper_.oppositeToBackwardDistance(scvf_));
            return distances;
        }
        else
        {
            std::array<Scalar, 3> distances;
            distances[0] = upwindgeometryHelper_.selfToOppositeDistance(scvf_);
            distances[1] = upwindgeometryHelper_.oppositeToBackwardDistance(scvf_);
            distances[2] = 0.5 * (upwindgeometryHelper_.selfToOppositeDistance(scvf_) + upwindgeometryHelper_.selfToForwardDistance(scvf_));
            return distances;
        }
    }

    /*!
     * \brief Check if a second order approximation for the lateral part of the advective term can be used
     *
     * This helper function checks if the scvf of interest is not too near to the
     * boundary so that a dof upstream with respect to the upstream dof is available.
     */
    bool canDoLateralSecondOrder_(const bool selfIsUpstream) const
    {
        static_assert(useHigherOrder, "Should only be reached if higher order methods are enabled");
        // If the self velocity is upstream, the downstream velocity can be assigned or retrieved
        // from the boundary, even if there is no parallel neighbor.
        // If the self velocity is downstream and there is no parallel neighbor I cannot use a second order approximation.
        return selfIsUpstream || fvGeometry_.hasParallelNeighbor(scvf_);
    }

    /*!
     * \brief Returns an array of the three momenta needed for second order upwinding methods. (Lateral)
     */
    template<class DensityPair>
    std::array<Scalar, 3> getLateralSecondOrderUpwindMomenta_(const DensityPair densityPair,
                                                              const bool selfIsUpstream) const
    {
        assert(scvf_.isLateral());
        static_assert(useHigherOrder, "Should only be reached if higher order methods are enabled");

        // If the lateral face lies on a boundary, we assume that the parallel velocity on the boundary is actually known,
        // thus we always use this value for the computation of the transported momentum.
        if ( fvGeometry_.lateralFaceContactsBoundary(scvf_) )
        {
            const Scalar boundaryMomentum = getParallelVelocityFromBoundary_(scvf_)  * densityPair.first;
            return {boundaryMomentum, boundaryMomentum, boundaryMomentum};
        }

        std::array<Scalar, 3> momenta;
        if (selfIsUpstream)
        {
            momenta[0] = elemVolVars_[fvGeometry_.parallelScvIdx(scvf_)].velocity() * densityPair.first;
            momenta[1] = elemVolVars_[scvf_.insideScvIdx()].velocity() * densityPair.first;

            // The "upstream-upstream" velocity is either retrieved from the opposite parallel neighbor or from the boundary.
            const auto& oppositeLateralScvf = fvGeometry_.oppositeLateralScvf(scvf_);
            if (fvGeometry_.hasParallelNeighbor(oppositeLateralScvf) && !fvGeometry_.lateralFaceContactsBoundary(oppositeLateralScvf))
                momenta[2] = elemVolVars_[fvGeometry_.parallelScvIdx(oppositeLateralScvf)].velocity() * densityPair.first;
            else
                momenta[2] = getParallelVelocityFromOppositeBoundary_(scvf_) * densityPair.first;

            return momenta;
        }
        else
        {
            momenta[0] = elemVolVars_[scvf_.insideScvIdx()].velocity() * densityPair.second;
            momenta[1] = elemVolVars_[fvGeometry_.parallelScvIdx(scvf_)].velocity() * densityPair.second;

            // If there is another parallel neighbor I can assign the "upstream-upstream" velocity, otherwise I retrieve it from the boundary.
            if (fvGeometry_.hasSecondParallelNeighbor(scvf_) && !fvGeometry_.lateralFaceContactsBoundary(fvGeometry_.outerParallelLateralScvf(scvf_)) )
                momenta[2] = elemVolVars_[fvGeometry_.secondParallelScvIdx(scvf_)].velocity() * densityPair.second;
            else
                momenta[2] = getParallelVelocityFromOutsideParallelBoundary_(scvf_) * densityPair.second;

            return momenta;
        }
    }

   /*!
    * \brief Returns an array of momenta needed for basic upwinding methods.
    */
    template<class DensityPair>
    std::array<Scalar, 2> getLateralFirstOrderUpwindMomenta_(const DensityPair densityPair,
                                                             const bool selfIsUpstream) const
    {
        assert(scvf_.isLateral());

        // use the Dirichlet velocity as for transported momentum if the lateral face is on a Dirichlet boundary
        if (fvGeometry_.lateralFaceContactsBoundary(scvf_))
        {
                const Scalar boundaryMomentum =  getParallelVelocityFromBoundary_(scvf_) * densityPair.first;
                return {boundaryMomentum, boundaryMomentum};
        }

        const Scalar momentumSelf = elemVolVars_[scvf_.insideScvIdx()].velocity() * densityPair.first;
        const Scalar momentumParallel =  elemVolVars_[scvf_.outsideScvIdx()].velocity()* densityPair.second;
        if (selfIsUpstream)
            return {momentumParallel, momentumSelf};
        else
            return {momentumSelf, momentumParallel};
    }

   /*!
    * \brief Returns an array of distances needed for non-uniform higher order upwind schemes
    * computes lateral distances {upstream to downstream distance, up-upstream to upstream distance, downstream staggered cell size}
    */
    std::array<Scalar, 3> getLateralDistances_(const bool selfIsUpstream) const
    {
        // static_assert(useHigherOrder, "Should only be reached if higher order methods are enabled");
        assert(scvf_.isLateral());

        std::array<Scalar, 3> distances;
        if (selfIsUpstream)
        {
            distances[0] = upwindgeometryHelper_.selfToParallelDistance(scvf_);
            distances[1] = upwindgeometryHelper_.selfToParallelDistance(fvGeometry_.oppositeLateralScvf(scvf_));

            if (fvGeometry_.hasParallelNeighbor(scvf_))
                distances[2] = fvGeometry_.outsideScvLateralLength(scvf_);
            else
                distances[2] = fvGeometry_.insideScvLateralLength(scvf_) / 2.0;

            return distances;
        }
        else
        {
            distances[0] = fvGeometry_.selfToParallelDistance(scvf_);
            distances[1] = fvGeometry_.paralleltoSecondParallelDistance(scvf_);
            distances[2] = fvGeometry_.outsideScvLateralLength(scvf_);
            return distances;
        }
    }

   /*!
    * \brief Return the parallel velocity for lateral faces on the boundary
    *
    *   =====================
    *   ----------LLLLLLL (*)      /  (scv) frontal face
    *   -       / ssssssss -       s  current Scv
    *   -       / ssssssss -       L  current lateral face
    *   -       / ssssssss x       -- element boundaries
    *   -       / ssssssss -       x  dof position
    *   -       / ssssssss -       == domain boundary
    *   --------------------       (*) queried velocity's location
    *
    *   Calls the problem to retrieve a fixed value set on the boundary.
    */
    Scalar getParallelVelocityFromBoundary_(const SubControlVolumeFace& lateralScvf) const
    {
        assert(lateralScvf.boundary() || fvGeometry_.hasHalfParallelNeighbor(lateralScvf) || fvGeometry_.hasCornerParallelNeighbor(lateralScvf));

        const auto& scv = fvGeometry_.scv(lateralScvf.insideScvIdx());
        const auto& element = fvGeometry_.gridGeometry().element(scv.elementIndex());
        const auto& lateralFaceBcTypes = problem_.boundaryTypes(element, lateralScvf);

        // In the case there is a neumann condition for momentum at the face in question, we force a zero gradient condition
        // This would include cases such as a "symmetry" condition, or a dirichlet pressure condition
        // In the case of a zero gradient condition, we return the velocity at the scv dof, forcing a zero gradient
        const bool useZeroGradient = lateralFaceBcTypes.hasNeumann();
        if (useZeroGradient)
            return elemVolVars_[lateralScvf.insideScvIdx()].velocity();

// TODO: Adapt to new BJS interface (Pseudo code)
//         const bool hasBJS = lateralFaceBcTypes.hasBeaversJoesephSlipVelocityCondition();
//         if (hasBJS)
//             return object.beaversJoesephSlipVelocity(x,y,z);

        if (lateralFaceBcTypes.isDirichlet(scv.dofAxis()))
        {
            const auto& tempCornerFace = makeCornerFace_(lateralScvf);
            return problem_.dirichlet(element, tempCornerFace)[scv.dofAxis()];
        }

        // Neumann conditions are not implemented well, throw this error in the case this other cases are skipped
        DUNE_THROW(Dune::InvalidStateException,
                   "Something went wrong with the boundary conditions for the momentum equations at global position " << lateralScvf.center());
    }

    SubControlVolumeFace makeCornerFace_(const SubControlVolumeFace& lateralScvf) const
    {
        return lateralScvf; // TODO!!
        // assert(lateralScvf.boundary() || fvGeometry_.hasHalfParallelNeighbor(lateralScvf) || fvGeometry_.hasCornerParallelNeighbor(lateralScvf));
        // const GlobalPosition& cornerPosition = fvGeometry_.cornerBoundaryPosition(lateralScvf);
        //
        //
        // return SubControlVolumeFace(element.geometry(), )
        //
        //
        //
        // return makeFaceCenteredSCVF(lateralScvf, cornerPosition);
    }

   /*!
    * \brief Return a velocity value from the opposite face's boundary a boundary for which the boundary conditions have to be checked.
    */
    Scalar getParallelVelocityFromOppositeBoundary_(const SubControlVolumeFace& lateralScvf) const
    {
        //   ----------LLLLLLLL--       /  (scv) frontal face
        //   -       / ssssssss -       s  current Scv
        //   -       / ssssssss -       L  current lateral face
        //   -       / ssssssss x       -- element boundaries
        //   -       / ssssssss -       x  dof position
        //   -       / ssssssss -       O  opposite lateral face
        //   ----------OOOOOOOO--       == domain boundary
        //    ================ (*)      (*) queried velocity's location

        // A ghost subface at the boundary is created, featuring the location of the sub face's center
        const auto& oppositeLateralScvf = fvGeometry_.oppositeLateralScvf(scvf_);
        return getParallelVelocityFromBoundary_(oppositeLateralScvf);
    }

    /*!
     * \brief Return a velocity value from the outer parallel cell's boundary for which the boundary conditions have to be checked.
     */
    Scalar getParallelVelocityFromOutsideParallelBoundary_(const SubControlVolumeFace& lateralScvf) const
    {
        //                \=
        //   \==========\ (*) \==   /    (scv) frontal face
        //   ---------------        s    current Scv
        //   -       ppppp -        p    parallel Scv
        //   -       ppppp x        L    current lateral face
        //   -       ppppp -        x    dof position
        //   --------LLLLL--        \==\ possible domain boundary
        //   -     / sssss -        --   element boundaries
        //   -     / sssss x        (*)  queried velocity's location
        //   -     / sssss -
        //   ---------------

        // Collect the outer lateral face and check to make sure it is on a boundary
        const auto& outerLateralFace = fvGeometry_.outerParallelLateralScvf(lateralScvf);
        assert(outerLateralFace.boundary() || fvGeometry_.hasHalfParallelNeighbor(outerLateralFace) || fvGeometry_.hasCornerParallelNeighbor(outerLateralFace));

        const auto& parallelScv = fvGeometry_.scv(outerLateralFace.insideScvIdx());
        const auto& parallelElement = fvGeometry_.gridGeometry().element(parallelScv.elementIndex());
        const auto& boundaryTypes = problem_.boundaryTypes(parallelElement, outerLateralFace);

        // if there is a neumann boundary condition, we assume a zero gradient, and take the
        if (boundaryTypes.hasNeumann())
            return elemVolVars_[parallelScv.dofIndex()].velocity();

// TODO: Adapt to new BJS interface (Pseudo code)
//         const bool hasBJS = lateralFaceBcTypes.hasBeaversJoesephSlipVelocityCondition();
//         if (hasBJS)
//             return object.beaversJoesephSlipVelocity(x,y,z);

        if (boundaryTypes.isDirichlet(parallelScv.dofAxis()))
        {
            const auto& tempCornerFace = makeCornerFace_(outerLateralFace);
            return problem_.dirichlet(parallelElement, tempCornerFace)[parallelScv.dofAxis()];
        }

        // Neumann conditions are not implemented well, throw this error in the case this other cases are skipped
        DUNE_THROW(Dune::InvalidStateException,
                   "Something went wrong with the boundary conditions for the momentum equations at global position " << outerLateralFace.center());
    }

    const Element& element_;
    const FVElementGeometry& fvGeometry_;
    const Problem& problem_;
    const SubControlVolumeFace& scvf_;
    const ElementVolumeVariables& elemVolVars_;
    const UpwindScheme& upwindScheme_;
    FaceCenteredStaggeredUpwindGeometryHelper<FVElementGeometry> upwindgeometryHelper_;
};

} // end namespace Dumux

#endif
