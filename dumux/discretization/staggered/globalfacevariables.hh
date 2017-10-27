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
 * \brief The global face variables class for staggered grid models
 */
#ifndef DUMUX_DISCRETIZATION_STAGGERED_GLOBAL_FACEVARIABLES_HH
#define DUMUX_DISCRETIZATION_STAGGERED_GLOBAL_FACEVARIABLES_HH

#include <dumux/implicit/properties.hh>

namespace Dumux
{

template<class TypeTag>
class StaggeredGlobalFaceVariables
{
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using Problem = typename GET_PROP_TYPE(TypeTag, Problem);
    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
    using SolutionVector = typename GET_PROP_TYPE(TypeTag, SolutionVector);
    using FaceVariables = typename GET_PROP_TYPE(TypeTag, FaceVariables);
    using FVGridGeometry = typename GET_PROP_TYPE(TypeTag, FVGridGeometry);
    using IndexType = typename GridView::IndexSet::IndexType;

    using DofTypeIndices = typename GET_PROP(TypeTag, DofTypeIndices);
    typename DofTypeIndices::CellCenterIdx cellCenterIdx;
    typename DofTypeIndices::FaceIdx faceIdx;

public:
    StaggeredGlobalFaceVariables(const Problem& problem) : problemPtr_(&problem) {}

    void update(const FVGridGeometry& fvGridGeometry, const SolutionVector& sol)
    {

        const auto& faceSol = sol[faceIdx];

        const auto numFaceDofs = fvGridGeometry.gridView().size(1);

        faceVariables_.resize(numFaceDofs);
        assert(faceVariables_.size() == faceSol.size());

        for(int i = 0; i < numFaceDofs; ++i)
        {
            faceVariables_[i].update(faceSol[i]);
        }
    }

    const FaceVariables& faceVars(const IndexType facetIdx) const
    { return faceVariables_[facetIdx]; }

    FaceVariables& faceVars(const IndexType facetIdx)
    { return faceVariables_[facetIdx]; }
private:
    const Problem& problem_() const
    { return *problemPtr_; }

    const Problem* problemPtr_;

    std::vector<FaceVariables> faceVariables_;
};


} // end namespace

#endif
