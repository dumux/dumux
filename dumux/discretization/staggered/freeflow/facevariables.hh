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
    using FaceSolutionVector = typename GET_PROP_TYPE(TypeTag, FaceSolutionVector);
    using FacePrimaryVariables = typename GET_PROP_TYPE(TypeTag, FacePrimaryVariables);
    using SubControlVolumeFace = typename GET_PROP_TYPE(TypeTag, SubControlVolumeFace);
    using GlobalFaceVars = typename GET_PROP_TYPE(TypeTag, GlobalFaceVars);
    using FaceSolution = typename GET_PROP_TYPE(TypeTag, StaggeredFaceSolution);

    struct SubFaceData
    {
        Scalar velocityNormalInside;
        Scalar velocityNormalOutside;
        Scalar velocityParallelInside;
        Scalar velocityParallelOutside;
    };

public:
    StaggeredFaceVariables() = default;

    StaggeredFaceVariables(const GlobalFaceVars& globalFacesVars) : globalFaceVarsPtr_(&globalFacesVars) {}

    const StaggeredFaceVariables& operator [](const SubControlVolumeFace& scvf) const
    { return globalFaceVars().faceVars(scvf.index()); }

    // // operator for the access with an index
    // // needed for cc methods for the access to the boundary volume variables
    // const StaggeredFaceVariables& operator [](const IndexType scvIdx) const
    // { return globalVolVars().volVars(scvIdx); }


    void update(const FacePrimaryVariables &facePrivars)
    {
        velocitySelf_ = facePrivars;
    }

    void update(const SubControlVolumeFace& scvf, const FaceSolutionVector& sol)
    {
        velocitySelf_ = sol[scvf.dofIndex()];
        velocityOpposite_ = sol[scvf.dofIndexOpposingFace()];

        subFaceVelocities_.resize(scvf.pairData().size());


        for(int i = 0; i < scvf.pairData().size(); ++i)
        {
            subFaceVelocities_[i].velocityNormalInside = sol[scvf.pairData(i).normalPair.first];

            if(scvf.pairData(i).normalPair.second >= 0)
                subFaceVelocities_[i].velocityNormalOutside = sol[scvf.pairData(i).normalPair.second];

            subFaceVelocities_[i].velocityParallelInside = sol[scvf.dofIndex()];

            if(scvf.pairData(i).outerParallelFaceDofIdx >= 0)
                subFaceVelocities_[i].velocityParallelOutside = sol[scvf.pairData(i).outerParallelFaceDofIdx];
        }

    }

    void update(const SubControlVolumeFace& scvf, const FaceSolution& faceSol)
    {
        velocitySelf_ = faceSol[scvf.dofIndex()];
        velocityOpposite_ = faceSol[scvf.dofIndexOpposingFace()];

        subFaceVelocities_.resize(scvf.pairData().size());

        for(int i = 0; i < scvf.pairData().size(); ++i)
        {
            subFaceVelocities_[i].velocityNormalInside = faceSol[scvf.pairData(i).normalPair.first];

            if(scvf.pairData(i).normalPair.second >= 0)
                subFaceVelocities_[i].velocityNormalOutside = faceSol[scvf.pairData(i).normalPair.second];

            subFaceVelocities_[i].velocityParallelInside = faceSol[scvf.dofIndex()];

            if(scvf.pairData(i).outerParallelFaceDofIdx >= 0)
                subFaceVelocities_[i].velocityParallelOutside = faceSol[scvf.pairData(i).outerParallelFaceDofIdx];
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

    auto& subFaceData() const
    {
        return subFaceVelocities_;
    }

    auto& subFaceData(const int localSubFaceIdx) const
    {
        return subFaceVelocities_[localSubFaceIdx];
    }

    //! The global volume variables object we are a restriction of
    const GlobalFaceVars& globalFaceVars() const
    { return *globalFaceVarsPtr_; }


private:
    const GlobalFaceVars* globalFaceVarsPtr_;
    Scalar velocity_;
    Scalar velocitySelf_;
    Scalar velocityOpposite_;
    std::vector<SubFaceData> subFaceVelocities_;
};

} // end namespace

#endif
