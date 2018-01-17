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
 * \copydoc Dumux::StaggeredFaceSolution
 */
#ifndef DUMUX_DISCRETIZATION_STAGGERED_FACE_SOLUTION_HH
#define DUMUX_DISCRETIZATION_STAGGERED_FACE_SOLUTION_HH

#include <vector>
#include <algorithm>
#include <dumux/common/properties.hh>
#include <dumux/discretization/staggered/elementvolumevariables.hh>

namespace Dumux
{

/*!
 * \ingroup StaggeredDiscretization
 * \brief The global face variables class for staggered models
 */
template<class TypeTag>
class StaggeredFaceSolution
{
    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
    using PrimaryVariables = typename GET_PROP_TYPE(TypeTag, PrimaryVariables);
    using Element = typename GridView::template Codim<0>::Entity;
    using FaceSolutionVector = typename GET_PROP_TYPE(TypeTag, FaceSolutionVector);
    using FacePrimaryVariables = typename GET_PROP_TYPE(TypeTag, FacePrimaryVariables);
    using FVGridGeometry = typename GET_PROP_TYPE(TypeTag, FVGridGeometry);
    using FVElementGeometry = typename GET_PROP_TYPE(TypeTag, FVGridGeometry)::LocalView;
    using SubControlVolumeFace = typename FVElementGeometry::SubControlVolumeFace;

    using DofTypeIndices = typename GET_PROP(TypeTag, DofTypeIndices);
    typename DofTypeIndices::CellCenterIdx cellCenterIdx;
    typename DofTypeIndices::FaceIdx faceIdx;

public:

    StaggeredFaceSolution(const SubControlVolumeFace& scvf, const FaceSolutionVector& sol,
                          const FVGridGeometry& fvGridGeometry)
    {
        const auto& connectivityMap = fvGridGeometry.connectivityMap();
        const auto& stencil = connectivityMap(faceIdx, faceIdx, scvf.index());

        facePriVars_.reserve(stencil.size());
        map_.reserve(stencil.size());

        for(const auto dofJ : stencil)
        {
            map_.push_back(dofJ);
            facePriVars_.push_back(sol[dofJ]);
        }
    }

    //! bracket operator const access
    template<typename IndexType>
    const FacePrimaryVariables& operator [](IndexType globalFaceDofIdx) const
    {
        const auto pos = std::find(map_.begin(), map_.end(), globalFaceDofIdx);
        assert (pos != map_.end());
        return facePriVars_[pos - map_.begin()];
    }

    //! bracket operator
    template<typename IndexType>
    FacePrimaryVariables& operator [](IndexType globalFaceDofIdx)
    {
        const auto pos = std::find(map_.begin(), map_.end(), globalFaceDofIdx);
        assert (pos != map_.end());
        return facePriVars_[pos - map_.begin()];
    }


private:

    std::vector<FacePrimaryVariables> facePriVars_;
    std::vector<unsigned int> map_;
};

} // end namespace

#endif
