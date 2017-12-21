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
#ifndef DUMUX_DISCRETIZATION_STAGGERED_ELEMENTFACEVARIABLES_HH
#define DUMUX_DISCRETIZATION_STAGGERED_ELEMENTFACEVARIABLES_HH

#include <dumux/common/properties.hh>

namespace Dumux
{

/*!
 * \ingroup ImplicitModel
 * \brief Base class for the face variables vector
 */
template<class TypeTag, bool enableGridFaceVariablesCache>
class StaggeredElementFaceVariables
{};


template<class TypeTag>
class StaggeredElementFaceVariables<TypeTag, /*enableGridFaceVariablesCache*/true>
{
    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
    using Element = typename GridView::template Codim<0>::Entity;
    using SolutionVector = typename GET_PROP_TYPE(TypeTag, SolutionVector);
    using SubControlVolumeFace = typename GET_PROP_TYPE(TypeTag, SubControlVolumeFace);
    using GridFaceVariables = typename GET_PROP_TYPE(TypeTag, GridFaceVariables);
    using FaceVariables = typename GET_PROP_TYPE(TypeTag, FaceVariables);
    using FVElementGeometry = typename GET_PROP_TYPE(TypeTag, FVElementGeometry);
    using IndexType = typename GridView::IndexSet::IndexType;

public:

    StaggeredElementFaceVariables(const GridFaceVariables& gridFaceVariables) : gridFaceVariablesPtr_(&gridFaceVariables) {}

    const FaceVariables& operator [](const SubControlVolumeFace& scvf) const
    { return gridFaceVariables().faceVars(scvf.index()); }

    // operator for the access with an index
    // needed for cc methods for the access to the boundary volume variables
    const FaceVariables& operator [](const IndexType scvfIdx) const
    { return gridFaceVariables().faceVars(scvfIdx); }


    //! For compatibility reasons with the case of not storing the face vars.
    //! function to be called before assembling an element, preparing the vol vars within the stencil
    void bind(const Element& element,
              const FVElementGeometry& fvGeometry,
              const SolutionVector& sol)
    {}

    // Binding of an element, prepares only the face variables of the element
    // specialization for Staggered models
    void bindElement(const Element& element,
                     const FVElementGeometry& fvGeometry,
                     const SolutionVector& sol)
    {}


    //! The global volume variables object we are a restriction of
    const GridFaceVariables& gridFaceVariables() const
    { return *gridFaceVariablesPtr_; }


private:
    const GridFaceVariables* gridFaceVariablesPtr_;
};

template<class TypeTag>
class StaggeredElementFaceVariables<TypeTag, /*enableGridFaceVariablesCache*/false>
{
    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
    using Element = typename GridView::template Codim<0>::Entity;
    using SolutionVector = typename GET_PROP_TYPE(TypeTag, SolutionVector);
    using SubControlVolumeFace = typename GET_PROP_TYPE(TypeTag, SubControlVolumeFace);
    using GridFaceVariables = typename GET_PROP_TYPE(TypeTag, GridFaceVariables);
    using FaceVariables = typename GET_PROP_TYPE(TypeTag, FaceVariables);
    using FVElementGeometry = typename GET_PROP_TYPE(TypeTag, FVElementGeometry);
    using IndexType = typename GridView::IndexSet::IndexType;

    using DofTypeIndices = typename GET_PROP(TypeTag, DofTypeIndices);
    typename DofTypeIndices::CellCenterIdx cellCenterIdx;
    typename DofTypeIndices::FaceIdx faceIdx;

public:

    StaggeredElementFaceVariables(const GridFaceVariables& globalFacesVars) : gridFaceVariablesPtr_(&globalFacesVars) {}

    const FaceVariables& operator [](const SubControlVolumeFace& scvf) const
    { return faceVariables_[scvf.localFaceIdx()]; }

    // operator for the access with an index
    const FaceVariables& operator [](const IndexType scvfIdx) const
    { return faceVariables_[getLocalIdx_(scvfIdx)]; }

    FaceVariables& operator [](const SubControlVolumeFace& scvf)
    { return faceVariables_[scvf.localFaceIdx()]; }

    // operator for the access with an index
    FaceVariables& operator [](const IndexType scvfIdx)
    { return faceVariables_[getLocalIdx_(scvfIdx)]; }

    //! For compatibility reasons with the case of not storing the vol vars.
    //! function to be called before assembling an element, preparing the vol vars within the stencil
    void bind(const Element& element,
              const FVElementGeometry& fvGeometry,
              const SolutionVector& sol)
    {
        faceVariables_.resize(fvGeometry.numScvf());
        faceVarIndices_.resize(fvGeometry.numScvf());

        for(auto&& scvf : scvfs(fvGeometry))
        {
            faceVariables_[scvf.localFaceIdx()].update(sol[faceIdx], GridFaceVariables().problem(), element, fvGeometry, scvf);
            faceVarIndices_[scvf.localFaceIdx()] = scvf.index();
        }
    }

    // Binding of an element, prepares only the face variables of the element
    // specialization for Staggered models
    void bindElement(const Element& element,
                     const FVElementGeometry& fvGeometry,
                     const SolutionVector& sol)
    {
        faceVariables_.resize(fvGeometry.numScvf());
        faceVarIndices_.resize(fvGeometry.numScvf());

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

    const int getLocalIdx_(const int scvfIdx) const
    {
        auto it = std::find(faceVarIndices_.begin(), faceVarIndices_.end(), scvfIdx);
        assert(it != faceVarIndices_.end() && "Could not find the current face variables for scvfIdx!");
        return std::distance(faceVarIndices_.begin(), it);
    }

    const GridFaceVariables* gridFaceVariablesPtr_;
    std::vector<IndexType> faceVarIndices_;
    std::vector<FaceVariables> faceVariables_;
};

} // end namespace

#endif
