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
// TODO: put in separate class
template<class TypeTag>
class StaggeredFaceVariables
{
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using FacePrimaryVariables = typename GET_PROP_TYPE(TypeTag, FacePrimaryVariables);
public:
    void update(const FacePrimaryVariables &facePrivars)
    {
        velocity_ = facePrivars[0];
    }

    Scalar velocity() const
    {
        return velocity_;
    }


private:
    Scalar velocity_;
};




template<class TypeTag>
class StaggeredModel;

template<class TypeTag>
class StaggeredGlobalFaceVariables;


namespace Properties
{
SET_TYPE_PROP(StaggeredModel, FaceVars, Dumux::StaggeredFaceVariables<TypeTag>);
SET_TYPE_PROP(StaggeredModel, GlobalFaceVars, Dumux::StaggeredGlobalFaceVariables<TypeTag>);
}

template<class TypeTag>
class StaggeredGlobalFaceVariables
{
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using Problem = typename GET_PROP_TYPE(TypeTag, Problem);
    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
    using FaceSolutionVector = typename GET_PROP_TYPE(TypeTag, FaceSolutionVector);
    using FaceVariables = typename GET_PROP_TYPE(TypeTag, FaceVars);
    using IndexType = typename GridView::IndexSet::IndexType;

public:
    void update(Problem& problem, const FaceSolutionVector& sol)
    {
        problemPtr_ = &problem;

        faceVariables_.resize(problem.model().numFaceDofs());
        assert(faceVariables_.size == sol.size());

        for(int i = 0; i < problem.model().numFaceDofs(); ++i)
        {
            faceVariables_[i].update(sol[i]);
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
