// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup StaggeredDiscretization
 * \copydoc Dumux::StaggeredFaceSolution
 */
#ifndef DUMUX_DISCRETIZATION_STAGGERED_FACE_SOLUTION_HH
#define DUMUX_DISCRETIZATION_STAGGERED_FACE_SOLUTION_HH

#include <algorithm>
#include <cassert>
#include <type_traits>
#include <vector>

namespace Dumux {

/*!
 * \ingroup StaggeredDiscretization
 * \brief The global face variables class for staggered models
 */
template<class FaceSolutionVector>
class StaggeredFaceSolution
{
    using FacePrimaryVariables = std::decay_t<decltype(std::declval<FaceSolutionVector>()[0])>;

public:

    template<class SubControlVolumeFace, class GridGeometry>
    StaggeredFaceSolution(const SubControlVolumeFace& scvf, const FaceSolutionVector& sol,
                          const GridGeometry& gridGeometry)
    {

        const auto& connectivityMap = gridGeometry.connectivityMap();
        const auto& stencil = connectivityMap(GridGeometry::faceIdx(), GridGeometry::faceIdx(), scvf.index());

        facePriVars_.reserve(stencil.size()+1);
        map_.reserve(stencil.size()+1);

        map_.push_back(scvf.dofIndex());
        facePriVars_.push_back(sol[scvf.dofIndex()]);
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

} // end namespace Dumux

#endif
