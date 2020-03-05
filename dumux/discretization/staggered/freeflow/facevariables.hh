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
 * \ingroup StaggeredDiscretization
 * \copydoc Dumux::StaggeredFaceVariables
 */
#ifndef DUMUX_DISCRETIZATION_STAGGERED_FREEFLOW_FACEVARIABLES_HH
#define DUMUX_DISCRETIZATION_STAGGERED_FREEFLOW_FACEVARIABLES_HH

#include <array>
#include <vector>
#include <type_traits>

namespace Dumux {

namespace Detail {

template<class Scalar, int upwindSchemeOrder>
struct InAxisVelocities
{
    Scalar self = 0.0;
    Scalar opposite = 0.0;
    std::array<Scalar, upwindSchemeOrder-1> forward{};
    std::array<Scalar, upwindSchemeOrder-1> backward{};
};

template<class Scalar>
struct InAxisVelocities<Scalar, 1>
{
    Scalar self = 0.0;
    Scalar opposite = 0.0;
};

} // end namespace Detail

/*!
 * \ingroup StaggeredDiscretization
 * \brief The face variables class for free flow staggered grid models.
 *        Contains all relevant velocities for the assembly of the momentum balance.
 *        When the upwindSchemeOrder is set to 2, additional velocities located at Dofs
 *        further from the central stencil will be added and used when calculating the
 *        advective term. When the order remains at 1, these velocities will not be provided.
 */
template<class FacePrimaryVariables, int dim, int upwindSchemeOrder>
class StaggeredFaceVariables
{
    static constexpr int numPairs = (dim == 2) ? 2 : 4;
    static constexpr bool useHigherOrder = upwindSchemeOrder > 1;

    using Scalar = typename FacePrimaryVariables::block_type;
    using InAxisVelocities = Detail::InAxisVelocities<Scalar, upwindSchemeOrder>;

public:

    /*!
    * \brief Partial update of the face variables. Only the face itself is considered.
    *
    * \param priVars The face-specific primary variales
    */
    void updateOwnFaceOnly(const FacePrimaryVariables& priVars)
    {
        inAxisVelocities_.self = priVars[0];
    }

    /*!
    * \brief Complete update of the face variables (i.e. velocities for free flow)
    *        for a given face
    *
    * \param faceSol The face-specific solution vector
    * \param problem The problem
    * \param element The element
    * \param fvGeometry The finite-volume geometry
    * \param scvf The sub-control volume face of interest
    */
    template<class FaceSolution, class Problem, class Element,
             class FVElementGeometry, class SubControlVolumeFace>
    void update(const FaceSolution& faceSol,
                const Problem& problem,
                const Element& element,
                const FVElementGeometry& fvGeometry,
                const SubControlVolumeFace& scvf)
    {
        static_assert(std::decay_t<decltype(faceSol[0])>::dimension == 1,
                      "\n\n\nVelocity primary variable must be a scalar value. \n\n Make sure to use\n\n ffSol = partial(sol, ffFaceIdx, ffCellCenterIdx);\n\n");

        inAxisVelocities_.self = faceSol[scvf.dofIndex()];
        inAxisVelocities_.opposite = faceSol[scvf.dofIndexOpposingFace()];

        addHigherOrderInAxisVelocities_(faceSol, scvf, std::integral_constant<bool, useHigherOrder>{});

        // handle all sub faces
        for (int i = 0; i < velocityParallel_.size(); ++i)
            velocityParallel_[i].fill(0.0);

        for (int i = 0; i < scvf.pairData().size(); ++i)
        {
            const auto& subFaceData = scvf.pairData(i);

            // treat the velocities normal to the face
            velocityLateralInside_[i] = faceSol[subFaceData.lateralPair.first];

            if (scvf.hasOuterLateral(i))
                velocityLateralOutside_[i] = faceSol[subFaceData.lateralPair.second];

            // treat the velocities parallel to the self face
            for (int j = 0; j < upwindSchemeOrder; j++)
            {
                if (scvf.hasParallelNeighbor(i,j))
                    velocityParallel_[i][j] = faceSol[subFaceData.parallelDofs[j]];
            }
        }
    }

    /*!
    * \brief Returns the velocity at the face itself
    */
    Scalar velocitySelf() const
    {
        return inAxisVelocities_.self;
    }

    /*!
    * \brief Returns the velocity at the opposing face
    */
    Scalar velocityOpposite() const
    {
        return inAxisVelocities_.opposite;
    }

    /*!
    * \brief Returns the velocity at a backward face
    *
    * \param backwardIdx The index describing how many faces backward this dof is from the opposite face
    */
    template<bool enable = useHigherOrder, std::enable_if_t<enable, int> = 0>
    Scalar velocityBackward(const int backwardIdx) const
    {
        return inAxisVelocities_.backward[backwardIdx];
    }

    /*!
    * \brief Returns the velocity at a forward face
    *
    * \param forwardIdx The index describing how many faces forward this dof is of the self face
    */
    template<bool enable = useHigherOrder, std::enable_if_t<enable, int> = 0>
    Scalar velocityForward(const int forwardIdx) const
    {
        return inAxisVelocities_.forward[forwardIdx];
    }

    /*!
    * \brief Returns the velocity at a parallel face
    *
    * \param localSubFaceIdx The local index of the subface
    * \param parallelDegreeIdx The index describing how many faces parallel this dof is of the parallel face
    */
    Scalar velocityParallel(const int localSubFaceIdx, const int parallelDegreeIdx = 0) const
    {
        return velocityParallel_[localSubFaceIdx][parallelDegreeIdx];
    }

    /*!
    * \brief Returns the velocity at the inner normal face
    *
    * \param localSubFaceIdx The local index of the subface
    */
    Scalar velocityLateralInside(const int localSubFaceIdx) const
    {
        return velocityLateralInside_[localSubFaceIdx];
    }

    /*!
    * \brief Returns the velocity at the outer normal face
    *
    * \param localSubFaceIdx The local index of the subface
    */
    Scalar velocityLateralOutside(const int localSubFaceIdx) const
    {
        return velocityLateralOutside_[localSubFaceIdx];
    }

private:

    template<class SolVector, class SubControlVolumeFace>
    void addHigherOrderInAxisVelocities_(const SolVector& faceSol, const SubControlVolumeFace& scvf, std::false_type) {}

    template<class SolVector, class SubControlVolumeFace>
    void addHigherOrderInAxisVelocities_(const SolVector& faceSol, const SubControlVolumeFace& scvf, std::true_type)
    {

        // treat the velocity forward of the self face i.e. the face that is
        // forward wrt the self face by degree i
        for (int i = 0; i < scvf.axisData().inAxisForwardDofs.size(); i++)
        {
             if (scvf.hasForwardNeighbor(i))
                 inAxisVelocities_.forward[i]= faceSol[scvf.axisData().inAxisForwardDofs[i]];
        }

        // treat the velocity at the first backward face i.e. the face that is
        // behind the opposite face by degree i
        for (int i = 0; i < scvf.axisData().inAxisBackwardDofs.size(); i++)
        {
             if (scvf.hasBackwardNeighbor(i))
                 inAxisVelocities_.backward[i] = faceSol[scvf.axisData().inAxisBackwardDofs[i]];
        }
    }

    InAxisVelocities inAxisVelocities_;
    std::array<std::array<Scalar, upwindSchemeOrder>, numPairs> velocityParallel_;
    std::array<Scalar, numPairs>  velocityLateralInside_;
    std::array<Scalar, numPairs>  velocityLateralOutside_;

};

} // end namespace Dumux

#endif
