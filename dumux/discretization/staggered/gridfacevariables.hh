// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
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
        auto fvGeometry = localView(gridGeometry);
        for (auto&& element : elements(gridGeometry.gridView()))
        {
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
