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
 * \ingroup StaggeredDiscretization
 * \copydoc Dumux::StaggeredFaceVariables
 */
#ifndef DUMUX_DISCRETIZATION_STAGGERED_FREEFLOW_FACEVARIABLES_HH
#define DUMUX_DISCRETIZATION_STAGGERED_FREEFLOW_FACEVARIABLES_HH

#include <dune/common/fvector.hh>
#include <dumux/common/properties.hh>

namespace Dumux
{

/*!
 * \ingroup StaggeredDiscretization
 * \brief The face variables class for free flow staggered grid models.
 *        Contains all relevant velocities for the assembly of the momentum balance.
 */
template<class TypeTag>
class StaggeredFaceVariables
{
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using Problem = typename GET_PROP_TYPE(TypeTag, Problem);
    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
    using FacePrimaryVariables = typename GET_PROP_TYPE(TypeTag, FacePrimaryVariables);
    using FVElementGeometry = typename GET_PROP_TYPE(TypeTag, FVElementGeometry);
    using SubControlVolumeFace = typename FVElementGeometry::SubControlVolumeFace;
    using Indices = typename GET_PROP_TYPE(TypeTag, Indices);

    static constexpr int dimWorld = GridView::dimensionworld;
    static constexpr int numPairs = (dimWorld == 2) ? 2 : 4;

    using GlobalPosition = Dune::FieldVector<Scalar, dimWorld>;
    using Element = typename GridView::template Codim<0>::Entity;

    using DofTypeIndices = typename GET_PROP(TypeTag, DofTypeIndices);
    typename DofTypeIndices::CellCenterIdx cellCenterIdx;
    typename DofTypeIndices::FaceIdx faceIdx;

public:

    /*!
    * \brief Partial update of the face variables. Only the face itself is considered.
    *
    * \param priVars The face-specific primary variales
    */
    void updateOwnFaceOnly(const FacePrimaryVariables& priVars)
    {
        velocitySelf_ = priVars[0];
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
    template<class SolVector>
    void update(const SolVector& faceSol,
                const Problem& problem,
                const Element& element,
                const FVElementGeometry& fvGeometry,
                const SubControlVolumeFace& scvf)
    {
        velocitySelf_ = faceSol[scvf.dofIndex()];
        velocityOpposite_ = faceSol[scvf.dofIndexOpposingFace()];

        // lambda to conveniently create a ghost face which is outside the domain, parallel to the scvf of interest
        auto makeGhostFace = [](const auto& pos)
        {
            return SubControlVolumeFace(pos, std::vector<unsigned int>{0,0});
        };

        // lambda to check whether there is a parallel face neighbor
        auto hasParallelNeighbor = [](const auto& subFaceData)
        {
            return subFaceData.outerParallelFaceDofIdx >= 0;
        };

        // lambda to check whether there is a normal face neighbor
        auto hasNormalNeighbor = [](const auto& subFaceData)
        {
            return subFaceData.normalPair.second >= 0;
        };

        // handle all sub faces
        for(int i = 0; i < scvf.pairData().size(); ++i)
        {
            const auto& subFaceData = scvf.pairData(i);

            // treat the velocities normal to the face
            velocityNormalInside_[i] = faceSol[subFaceData.normalPair.first];

            if(hasNormalNeighbor(subFaceData))
            {
                velocityNormalOutside_[i] = faceSol[subFaceData.normalPair.second];
            }
            else
            {
                const auto& normalFace = fvGeometry.scvf(scvf.insideScvIdx(), subFaceData.localNormalFaceIdx);
                const auto normalDirIdx = normalFace.directionIndex();
                velocityNormalOutside_[i] = problem.dirichlet(element, makeGhostFace(subFaceData.virtualOuterNormalFaceDofPos))[Indices::velocity(normalDirIdx)];
            }

            // treat the velocity parallel to the face
            velocityParallel_[i] = hasParallelNeighbor(subFaceData) ?
                                   velocityParallel_[i] = faceSol[subFaceData.outerParallelFaceDofIdx] :
                                   problem.dirichlet(element, makeGhostFace(subFaceData.virtualOuterParallelFaceDofPos))[Indices::velocity(scvf.directionIndex())];
        }
    }

    /*!
    * \brief Returns the velocity at the face itself
    */
    Scalar velocitySelf() const
    {
        return velocitySelf_;
    }

    /*!
    * \brief Returns the velocity at the opposing face
    */
    Scalar velocityOpposite() const
    {
        return velocityOpposite_;
    }

    /*!
    * \brief Returns the velocity at the parallel face
    *
    * \param localSubFaceIdx The local index of the subface
    */
    Scalar velocityParallel(const int localSubFaceIdx) const
    {
        return velocityParallel_[localSubFaceIdx];
    }

    /*!
    * \brief Returns the velocity at the inner normal face
    *
    * \param localSubFaceIdx The local index of the subface
    */
    Scalar velocityNormalInside(const int localSubFaceIdx) const
    {
        return velocityNormalInside_[localSubFaceIdx];
    }

    /*!
    * \brief Returns the velocity at the outer normal face
    *
    * \param localSubFaceIdx The local index of the subface
    */
    Scalar velocityNormalOutside(const int localSubFaceIdx) const
    {
        return velocityNormalOutside_[localSubFaceIdx];
    }

private:

    Scalar velocitySelf_;
    Scalar velocityOpposite_;
    std::array<Scalar, numPairs> velocityParallel_;
    std::array<Scalar, numPairs> velocityNormalInside_;
    std::array<Scalar, numPairs> velocityNormalOutside_;
};

} // end namespace

#endif
