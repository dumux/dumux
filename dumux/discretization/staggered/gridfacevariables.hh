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

#include <dumux/common/properties.hh>
#include <dumux/discretization/staggered/facesolution.hh>

//! make the local view function available whenever we use this class
#include <dumux/discretization/localview.hh>

namespace Dumux {

/*!
 * \ingroup StaggeredDiscretization
 * \brief Face variables cache class for staggered models
 */
template<class TypeTag, bool enableGlobalFaceVarsCache>
class StaggeredGridFaceVariables;

/*!
 * \ingroup StaggeredDiscretization
 * \brief Face variables cache class for staggered models.
          Specialization in case of storing the face variables.
 */
template<class TypeTag>
class StaggeredGridFaceVariables<TypeTag, /*enableGlobalFaceVarsCache*/true>
{
    using Problem = typename GET_PROP_TYPE(TypeTag, Problem);
    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
    using SolutionVector = typename GET_PROP_TYPE(TypeTag, SolutionVector);
    using FaceVariables = typename GET_PROP_TYPE(TypeTag, FaceVariables);
    using FVGridGeometry = typename GET_PROP_TYPE(TypeTag, FVGridGeometry);
    using IndexType = typename GridView::IndexSet::IndexType;

    using DofTypeIndices = typename GET_PROP(TypeTag, DofTypeIndices);
    typename DofTypeIndices::FaceIdx faceIdx;

public:
    //! export the type of the local view
    using LocalView = typename GET_PROP_TYPE(TypeTag, ElementFaceVariables);

    StaggeredGridFaceVariables(const Problem& problem) : problemPtr_(&problem) {}

    //! Update all face variables
    void update(const FVGridGeometry& fvGridGeometry, const SolutionVector& sol)
    {
        const auto& faceSol = sol[faceIdx];
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
template<class TypeTag>
class StaggeredGridFaceVariables<TypeTag, /*enableGlobalFaceVarsCache*/false>
{
    using Problem = typename GET_PROP_TYPE(TypeTag, Problem);
    using SolutionVector = typename GET_PROP_TYPE(TypeTag, SolutionVector);
    using FVGridGeometry = typename GET_PROP_TYPE(TypeTag, FVGridGeometry);

public:
    //! export the type of the local view
    using LocalView = typename GET_PROP_TYPE(TypeTag, ElementFaceVariables);

    StaggeredGridFaceVariables(const Problem& problem) : problemPtr_(&problem) {}

    //! Do nothing here.
    void update(const FVGridGeometry& fvGridGeometry, const SolutionVector& sol) {}

    const Problem& problem() const
    { return *problemPtr_; }

private:
    const Problem* problemPtr_;
};

} // end namespace Dumux

#endif
