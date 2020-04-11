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
  * \copydoc Dumux::StaggeredGridFaceVariables
  */
#ifndef DUMUX_DISCRETIZATION_STAGGERED_GRID_FACEVARIABLES_HH
#define DUMUX_DISCRETIZATION_STAGGERED_GRID_FACEVARIABLES_HH

// make the local view function available whenever we use this class
#include <dumux/discretization/localview.hh>
#include <dumux/discretization/staggered/elementfacevariables.hh>
#include <dumux/discretization/staggered/facesolution.hh>

namespace Dumux {

/*!
 * \ingroup StaggeredDiscretization
 * \brief Traits class to be used for the StaggeredGridFaceVariables.
 *
 * \tparam P The problem type
 * \tparam FV The face variables type
 */
template<class P, class FV>
struct StaggeredDefaultGridFaceVariablesTraits
{
    template<class GridFaceVariables, bool enableCache>
    using LocalView = StaggeredElementFaceVariables<GridFaceVariables, enableCache>;

    using FaceVariables = FV;
    using Problem = P;
};

/*!
 * \ingroup StaggeredDiscretization
 * \brief Face variables cache class for staggered models
 */
template<class Problem,
         class FaceVariables,
         bool cachingEnabled = false,
         class Traits = StaggeredDefaultGridFaceVariablesTraits<Problem, FaceVariables> >
class StaggeredGridFaceVariables;

/*!
 * \ingroup StaggeredDiscretization
 * \brief Face variables cache class for staggered models.
          Specialization in case of storing the face variables.
 */
template<class P, class FV, class Traits>
class StaggeredGridFaceVariables<P, FV, /*cachingEnabled*/true, Traits>
{
    using ThisType = StaggeredGridFaceVariables<P, FV, /*cachingEnabled*/true, Traits>;
    using Problem = typename Traits::Problem;

public:
    //! make it possible to query if caching is enabled
    static constexpr bool cachingEnabled = true;

    //! export the type of the local view
    using LocalView = typename Traits::template LocalView<ThisType, cachingEnabled>;

    //! export the type of the face variables
    using FaceVariables = typename Traits::FaceVariables;

    StaggeredGridFaceVariables(const Problem& problem) : problemPtr_(&problem) {}

    //! Update all face variables
    template<class GridGeometry, class SolutionVector>
    void update(const GridGeometry& gridGeometry, const SolutionVector& sol)
    {
        faceVariables_.resize(gridGeometry.numScvf());

        for (auto&& element : elements(gridGeometry.gridView()))
        {
            auto fvGeometry = localView(gridGeometry);
            fvGeometry.bindElement(element);

            for (auto&& scvf : scvfs(fvGeometry))
            {
                const auto faceSol = StaggeredFaceSolution(scvf, sol, gridGeometry);
                faceVariables_[scvf.index()].update(faceSol, problem(), element, fvGeometry, scvf);
            }
        }
    }

    const FaceVariables& faceVars(const std::size_t facetIdx) const
    { return faceVariables_[facetIdx]; }

    FaceVariables& faceVars(const std::size_t facetIdx)
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
template<class P, class FV, class Traits>
class StaggeredGridFaceVariables<P, FV, /*cachingEnabled*/false, Traits>
{
    using ThisType = StaggeredGridFaceVariables<P, FV, /*cachingEnabled*/false, Traits>;
    using Problem = typename Traits::Problem;

public:
    //! make it possible to query if caching is enabled
    static constexpr bool cachingEnabled = false;

    //! export the type of the local view
    using LocalView = typename Traits::template LocalView<ThisType, cachingEnabled>;

    //! export the type of the face variables
    using FaceVariables = typename Traits::FaceVariables;

    StaggeredGridFaceVariables(const Problem& problem) : problemPtr_(&problem) {}

    //! Do nothing here.
    template<class GridGeometry, class SolutionVector>
    void update(const GridGeometry& gridGeometry, const SolutionVector& sol) {}

    const Problem& problem() const
    { return *problemPtr_; }

private:
    const Problem* problemPtr_;
};

} // end namespace Dumux

#endif
