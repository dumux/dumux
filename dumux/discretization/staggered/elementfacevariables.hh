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
 * \copydoc Dumux::StaggeredElementFaceVariables
 */
#ifndef DUMUX_DISCRETIZATION_STAGGERED_ELEMENTFACEVARIABLES_HH
#define DUMUX_DISCRETIZATION_STAGGERED_ELEMENTFACEVARIABLES_HH

#include <algorithm>
#include <cassert>
#include <vector>

namespace Dumux {

/*!
 * \ingroup StaggeredDiscretization
 * \brief Base class for the face variables vector
 */
template<class GridFaceVariables, bool cachingEnabled>
class StaggeredElementFaceVariables
{};

/*!
 * \ingroup StaggeredDiscretization
 * \brief Class for the face variables vector. Specialization for the case of storing the face variables globally.
 */
template<class GFV>
class StaggeredElementFaceVariables<GFV, /*cachingEnabled*/true>
{
public:
    //! export type of the grid volume variables
    using GridFaceVariables = GFV;

    //! export type of the volume variables
    using FaceVariables = typename GridFaceVariables::FaceVariables;

    StaggeredElementFaceVariables(const GridFaceVariables& gridFaceVariables) : gridFaceVariablesPtr_(&gridFaceVariables) {}

    //! operator for the access with an scvf
    template<class SubControlVolumeFace, typename std::enable_if_t<!std::is_integral<SubControlVolumeFace>::value, int> = 0>
    const FaceVariables& operator [](const SubControlVolumeFace& scvf) const
    { return gridFaceVariables().faceVars(scvf.index()); }

    //! operator for the access with an index
    //! needed for cc methods for the access to the boundary volume variables
    const FaceVariables& operator [](const std::size_t scvfIdx) const
    { return gridFaceVariables().faceVars(scvfIdx); }


    //! For compatibility reasons with the case of not storing the face vars.
    //! function to be called before assembling an element, preparing the vol vars within the stencil
    template<class FVElementGeometry, class SolutionVector>
    void bind(const typename FVElementGeometry::GridGeometry::GridView::template Codim<0>::Entity& element,
              const FVElementGeometry& fvGeometry,
              const SolutionVector& sol)
    {}

    //! Binding of an element, prepares only the face variables of the element
    //! specialization for Staggered models
    template<class FVElementGeometry, class SolutionVector>
    void bindElement(const typename FVElementGeometry::GridGeometry::GridView::template Codim<0>::Entity& element,
                     const FVElementGeometry& fvGeometry,
                     const SolutionVector& sol)
    {}


    //! The global volume variables object we are a restriction of
    const GridFaceVariables& gridFaceVariables() const
    { return *gridFaceVariablesPtr_; }


private:
    const GridFaceVariables* gridFaceVariablesPtr_;
};

/*!
 * \ingroup StaggeredDiscretization
 * \brief Class for the face variables vector. Specialization for the case of not storing the face variables globally.
 */
template<class GFV>
class StaggeredElementFaceVariables<GFV, /*cachingEnabled*/false>
{
public:
    //! export type of the grid volume variables
    using GridFaceVariables = GFV;

    //! export type of the volume variables
    using FaceVariables = typename GridFaceVariables::FaceVariables;

    StaggeredElementFaceVariables(const GridFaceVariables& globalFacesVars) : gridFaceVariablesPtr_(&globalFacesVars) {}

    //! const operator for the access with an scvf
    template<class SubControlVolumeFace, typename std::enable_if_t<!std::is_integral<SubControlVolumeFace>::value, int> = 0>
    const FaceVariables& operator [](const SubControlVolumeFace& scvf) const
    { return faceVariables_[scvf.localFaceIdx()]; }

    //! const operator for the access with an index
    const FaceVariables& operator [](const std::size_t scvfIdx) const
    { return faceVariables_[getLocalIdx_(scvfIdx)]; }

    //! operator for the access with an scvf
    template<class SubControlVolumeFace, typename std::enable_if_t<!std::is_integral<SubControlVolumeFace>::value, int> = 0>
    FaceVariables& operator [](const SubControlVolumeFace& scvf)
    { return faceVariables_[scvf.localFaceIdx()]; }

    // operator for the access with an index
    FaceVariables& operator [](const std::size_t scvfIdx)
    { return faceVariables_[getLocalIdx_(scvfIdx)]; }

    //! For compatibility reasons with the case of not storing the vol vars.
    //! function to be called before assembling an element, preparing the vol vars within the stencil
    template<class FVElementGeometry, class SolutionVector>
    void bind(const typename FVElementGeometry::GridGeometry::GridView::template Codim<0>::Entity& element,
              const FVElementGeometry& fvGeometry,
              const SolutionVector& sol)
    {
        faceVariables_.resize(fvGeometry.numScvf());
        faceVarIndices_.resize(fvGeometry.numScvf());

        constexpr auto faceIdx = FVElementGeometry::GridGeometry::faceIdx();

        for(auto&& scvf : scvfs(fvGeometry))
        {
            faceVariables_[scvf.localFaceIdx()].update(sol[faceIdx], gridFaceVariables().problem(), element, fvGeometry, scvf);
            faceVarIndices_[scvf.localFaceIdx()] = scvf.index();
        }
    }

    //! Binding of an element, prepares only the face variables of the element
    //! specialization for Staggered models
    template<class FVElementGeometry, class SolutionVector>
    void bindElement(const typename FVElementGeometry::GridGeometry::GridView::template Codim<0>::Entity& element,
                     const FVElementGeometry& fvGeometry,
                     const SolutionVector& sol)
    {
        faceVariables_.resize(fvGeometry.numScvf());
        faceVarIndices_.resize(fvGeometry.numScvf());

        constexpr auto faceIdx = FVElementGeometry::GridGeometry::faceIdx();

        for(auto&& scvf : scvfs(fvGeometry))
        {
            faceVariables_[scvf.localFaceIdx()].updateOwnFaceOnly(sol[faceIdx][scvf.dofIndex()]);
            faceVarIndices_[scvf.localFaceIdx()] = scvf.index();
        }
    }

    //! The global volume variables object we are a restriction of
    const GridFaceVariables& gridFaceVariables() const
    { return *gridFaceVariablesPtr_; }

private:

    int getLocalIdx_(const int scvfIdx) const
    {
        auto it = std::find(faceVarIndices_.begin(), faceVarIndices_.end(), scvfIdx);
        assert(it != faceVarIndices_.end() && "Could not find the current face variables for scvfIdx!");
        return std::distance(faceVarIndices_.begin(), it);
    }

    const GridFaceVariables* gridFaceVariablesPtr_;
    std::vector<std::size_t> faceVarIndices_;
    std::vector<FaceVariables> faceVariables_;
};

} // end namespace

#endif
