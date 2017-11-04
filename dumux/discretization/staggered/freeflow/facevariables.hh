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
 * \brief The face variables class for free flow staggered grid models
 */
#ifndef DUMUX_DISCRETIZATION_STAGGERED_FREEFLOW_FACEVARIABLES_HH
#define DUMUX_DISCRETIZATION_STAGGERED_FREEFLOW_FACEVARIABLES_HH

#include <dumux/common/basicproperties.hh>

namespace Dumux
{

namespace Properties
{
    NEW_PROP_TAG(StaggeredFaceSolution);
}

template<class TypeTag>
class StaggeredFaceVariables
{
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
    using FacePrimaryVariables = typename GET_PROP_TYPE(TypeTag, FacePrimaryVariables);
    using SubControlVolumeFace = typename GET_PROP_TYPE(TypeTag, SubControlVolumeFace);

    static constexpr int dimWorld = GridView::dimensionworld;
    static constexpr int numPairs = (dimWorld == 2) ? 2 : 4;

public:

    void update(const FacePrimaryVariables &facePrivars)
    {
        velocitySelf_ = facePrivars;
    }

    template<class SolVector>
    void update(const SubControlVolumeFace& scvf, const SolVector& faceSol)
    {
        velocitySelf_ = faceSol[scvf.dofIndex()];
        velocityOpposite_ = faceSol[scvf.dofIndexOpposingFace()];

        for(int i = 0; i < scvf.pairData().size(); ++i)
        {
            velocityNormalInside_[i] = faceSol[scvf.pairData(i).normalPair.first];

            if(scvf.pairData(i).normalPair.second >= 0)
                velocityNormalOutside_[i] = faceSol[scvf.pairData(i).normalPair.second];

            if(scvf.pairData(i).outerParallelFaceDofIdx >= 0)
                velocityParallel_[i] = faceSol[scvf.pairData(i).outerParallelFaceDofIdx];
        }
    }

    Scalar velocitySelf() const
    {
        return velocitySelf_;
    }

    Scalar velocityOpposite() const
    {
        return velocityOpposite_;
    }

    Scalar velocityParallel(const int localSubFaceIdx) const
    {
        return velocityParallel_[localSubFaceIdx];
    }

    Scalar velocityNormalInside(const int localSubFaceIdx) const
    {
        return velocityNormalInside_[localSubFaceIdx];
    }

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
