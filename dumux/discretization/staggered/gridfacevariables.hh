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
  * \copydoc Dumux::StaggeredGridFaceVariables
  */
#ifndef DUMUX_DISCRETIZATION_STAGGERED_GRID_FACEVARIABLES_HH
#define DUMUX_DISCRETIZATION_STAGGERED_GRID_FACEVARIABLES_HH

//! make the local view function available whenever we use this class
#include <dumux/discretization/localview.hh>

namespace Dumux {

/*!
 * \ingroup StaggeredDiscretization
 * \brief Face variables cache class for staggered models
 */
template<class FVGridGeometry, class Traits, bool enableGlobalFaceVarsCache>
class StaggeredGridFaceVariables;

/*!
 * \ingroup StaggeredDiscretization
 * \brief Face variables cache class for staggered models.
          Specialization in case of storing the face variables.
 */
template<class FVGridGeometry, class Traits>
class StaggeredGridFaceVariables<FVGridGeometry, Traits, /*enableGlobalFaceVarsCache*/true>
{
    using ThisType = StaggeredGridFaceVariables<FVGridGeometry, Traits, true>;
    using IndexType = typename FVGridGeometry::GridView::IndexSet::IndexType;
    using Problem = typename Traits::Problem;

public:
    //! export the type of the local view
    using LocalView = typename Traits::template LocalView<FVGridGeometry, ThisType, true>;
    //! export the type of the face variables
    using FaceVariables = typename Traits::FaceVariables;

    StaggeredGridFaceVariables(const Problem& problem) : problemPtr_(&problem) {}

    //! Update all face variables
    template<class SolutionVector>
    void update(const FVGridGeometry& fvGridGeometry, const SolutionVector& faceSol)
    {

        faceVariables_.resize(fvGridGeometry.numScvf());

        for(auto&& element : elements(fvGridGeometry.gridView()))
        {
            auto fvGeometry = localView(fvGridGeometry);
            fvGeometry.bindElement(element);

            for(auto&& scvf : scvfs(fvGeometry))
            {
                faceVariables_[scvf.index()].update(faceSol, problem(), element, fvGeometry, scvf);
            }
        }
    }

    const FaceVariables& faceVars(const IndexType facetIdx) const
    { return faceVariables_[facetIdx]; }

    FaceVariables& faceVars(const IndexType facetIdx)
    { return faceVariables_[facetIdx]; }

    const Problem& problem() const
    { return *problemPtr_; }

private:

    const Problem* problemPtr_;

    std::vector<FaceVariables> faceVariables_;
};

/*!
 * \ingroup StaggeredDiscretization
 * \brief Face variables cache class for staggered models.
          Specialization in case of not storing the face variables.
 */
template<class FVGridGeometry, class Traits>
class StaggeredGridFaceVariables<FVGridGeometry, Traits, /*enableGlobalFaceVarsCache*/false>
{
    using ThisType = StaggeredGridFaceVariables<FVGridGeometry, Traits, false>;
    using Problem = typename Traits::Problem;

public:
    //! export the type of the local view
    using LocalView = typename Traits::template LocalView<FVGridGeometry, ThisType, false>;
    //! export the type of the face variables
    using FaceVariables = typename Traits::FaceVariables;

    StaggeredGridFaceVariables(const Problem& problem) : problemPtr_(&problem) {}

    //! Do nothing here.
    template<class SolutionVector>
    void update(const FVGridGeometry& fvGridGeometry, const SolutionVector& sol) {}

    const Problem& problem() const
    { return *problemPtr_; }

private:
    const Problem* problemPtr_;
};

} // end namespace Dumux

#endif
