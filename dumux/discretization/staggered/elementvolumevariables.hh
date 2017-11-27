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
 * \brief The local (stencil) volume variables class for cell centered models
 */
#ifndef DUMUX_DISCRETIZATION_STAGGERED_ELEMENT_VOLUMEVARIABLES_HH
#define DUMUX_DISCRETIZATION_STAGGERED_ELEMENT_VOLUMEVARIABLES_HH

#include <dumux/common/basicproperties.hh>

namespace Dumux
{

/*!
 * \ingroup ImplicitModel
 * \brief Base class for the volume variables vector
 */
template<class TypeTag, bool enableGlobalVolVarsCache>
class StaggeredElementVolumeVariables
{};

// specialization in case of storing the volume variables globally
template<class TypeTag>
class StaggeredElementVolumeVariables<TypeTag, /*enableGlobalVolVarsCache*/true>
{
    using Problem = typename GET_PROP_TYPE(TypeTag, Problem);
    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
    using SolutionVector = typename GET_PROP_TYPE(TypeTag, SolutionVector);
    using FVElementGeometry = typename GET_PROP_TYPE(TypeTag, FVElementGeometry);
    using VolumeVariables = typename GET_PROP_TYPE(TypeTag, VolumeVariables);
    using GlobalVolumeVariables = typename GET_PROP_TYPE(TypeTag, GlobalVolumeVariables);
    using SubControlVolume = typename GET_PROP_TYPE(TypeTag, SubControlVolume);
    using IndexType = typename GridView::IndexSet::IndexType;

    static const int dim = GridView::dimension;
    using Element = typename GridView::template Codim<0>::Entity;

public:
    //! Constructor
    StaggeredElementVolumeVariables(const GlobalVolumeVariables& globalVolVars)
    : globalVolVarsPtr_(&globalVolVars) {}

    const VolumeVariables& operator [](const SubControlVolume& scv) const
    { return globalVolVars().volVars(scv.dofIndex()); }

    // operator for the access with an index
    // needed for Staggered methods for the access to the boundary volume variables
    const VolumeVariables& operator [](const IndexType scvIdx) const
    { return globalVolVars().volVars(scvIdx); }

    // For compatibility reasons with the case of not storing the vol vars.
    // function to be called before assembling an element, preparing the vol vars within the stencil
    void bind(const Element& element,
              const FVElementGeometry& fvGeometry,
              const SolutionVector& sol)
    {}

    // function to prepare the vol vars within the element
    void bindElement(const Element& element,
                     const FVElementGeometry& fvGeometry,
                     const SolutionVector& sol)
    {}

    //! The global volume variables object we are a restriction of
    const GlobalVolumeVariables& globalVolVars() const
    { return *globalVolVarsPtr_; }

private:
    const GlobalVolumeVariables* globalVolVarsPtr_;
};


// Specialization when the current volume variables are not stored
template<class TypeTag>
class StaggeredElementVolumeVariables<TypeTag, /*enableGlobalVolVarsCache*/false>
{
    using Problem = typename GET_PROP_TYPE(TypeTag, Problem);
    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
    using Indices = typename GET_PROP_TYPE(TypeTag, Indices);
    using SolutionVector = typename GET_PROP_TYPE(TypeTag, SolutionVector);
    using ElementSolution = typename GET_PROP_TYPE(TypeTag, ElementSolutionVector);
    using VolumeVariables = typename GET_PROP_TYPE(TypeTag, VolumeVariables);
    using PrimaryVariables = typename GET_PROP_TYPE(TypeTag, PrimaryVariables);
    using GlobalVolumeVariables = typename GET_PROP_TYPE(TypeTag, GlobalVolumeVariables);
    using SubControlVolume = typename GET_PROP_TYPE(TypeTag, SubControlVolume);
    using FVElementGeometry = typename GET_PROP_TYPE(TypeTag, FVElementGeometry);
    using IndexType = typename GridView::IndexSet::IndexType;

    static const int dim = GridView::dimension;
    using Element = typename GridView::template Codim<0>::Entity;

    using DofTypeIndices = typename GET_PROP(TypeTag, DofTypeIndices);
    typename DofTypeIndices::CellCenterIdx cellCenterIdx;
    typename DofTypeIndices::FaceIdx faceIdx;

    static constexpr auto numEqCellCenter = GET_PROP_VALUE(TypeTag, NumEqCellCenter);

public:

    //! Constructor
    StaggeredElementVolumeVariables(const GlobalVolumeVariables& globalVolVars)
    : globalVolVarsPtr_(&globalVolVars) {}

    // Binding of an element, prepares the volume variables within the element stencil
    // called by the local jacobian to prepare element assembly
    void bind(const Element& element,
              const FVElementGeometry& fvGeometry,
              const SolutionVector& sol)
    {
        clear();

        const auto& problem = globalVolVars().problem();
        const auto& fvGridGeometry = fvGeometry.fvGridGeometry();
        const auto globalI = fvGridGeometry.elementMapper().index(element);
        const auto map = fvGridGeometry.connectivityMap();
        const auto& connectivityMapI = map(cellCenterIdx, cellCenterIdx, globalI);
        const auto numDofs = connectivityMapI.size();

        auto&& scvI = fvGeometry.scv(globalI);

        // resize local containers to the required size (for internal elements)
        volumeVariables_.resize(numDofs);
        volVarIndices_.resize(numDofs);
        int localIdx = 0;

        // Update the volume variables of the element at hand and the neighboring elements
        for (auto globalJ : connectivityMapI)
        {
            const auto& elementJ = fvGridGeometry.element(globalJ);
            auto&& scvJ = fvGeometry.scv(globalJ);
            PrimaryVariables priVars(0.0);
            priVars[cellCenterIdx] = sol[cellCenterIdx][globalJ];
            ElementSolution elemSol{std::move(priVars)};
            volumeVariables_[localIdx].update(elemSol,
                                              problem,
                                              elementJ,
                                              scvJ);
            volVarIndices_[localIdx] = scvJ.dofIndex();
            ++localIdx;
        }

        // Update boundary volume variables
        for (auto&& scvf : scvfs(fvGeometry))
        {
            // if we are not on a boundary, skip to the next scvf
            if (!scvf.boundary())
                continue;

            const auto bcTypes = problem.boundaryTypes(element, scvf);

            PrimaryVariables boundaryPriVars(0.0);

            for(int eqIdx = 0; eqIdx < numEqCellCenter; ++eqIdx)
            {
                if(bcTypes.isDirichlet(eqIdx) || bcTypes.isDirichletCell(eqIdx))
                    boundaryPriVars[cellCenterIdx][eqIdx] = problem.dirichlet(element, scvf)[cellCenterIdx][eqIdx];
                else if(bcTypes.isNeumann(eqIdx) || bcTypes.isOutflow(eqIdx) || bcTypes.isSymmetry())
                    boundaryPriVars[cellCenterIdx][eqIdx] = sol[cellCenterIdx][scvf.insideScvIdx()][eqIdx];
                //TODO: this assumes a zero-gradient for e.g. the pressure on the boundary
                // could be made more general by allowing a non-zero-gradient, provided in problem file
                else
                    if(eqIdx == Indices::pressureIdx)
                        DUNE_THROW(Dune::InvalidStateException, "Face at: " << scvf.center() << " has neither Dirichlet nor Neumann BC.");
            }

            volumeVariables_.resize(localIdx+1);
            volVarIndices_.resize(localIdx+1);

            ElementSolution elemSol{std::move(boundaryPriVars)};
            volumeVariables_[localIdx].update(elemSol,
                                              problem,
                                              element,
                                              scvI);
            volVarIndices_[localIdx] = scvf.outsideScvIdx();
             ++localIdx;
        }
    }

    // Binding of an element, prepares only the volume variables of the element
    // specialization for Staggered models
    void bindElement(const Element& element,
                     const FVElementGeometry& fvGeometry,
                     const SolutionVector& sol)
    {
        clear();

        const auto eIdx = fvGeometry.fvGridGeometry().elementMapper().index(element);
        volumeVariables_.resize(1);
        volVarIndices_.resize(1);

        // update the volume variables of the element
        auto&& scv = fvGeometry.scv(eIdx);
        PrimaryVariables priVars(0.0);
        priVars[cellCenterIdx] = sol[cellCenterIdx][eIdx];
        ElementSolution elemSol{std::move(priVars)};
        volumeVariables_[0].update(elemSol,
                                   globalVolVars().problem(),
                                   element,
                                   scv);
        volVarIndices_[0] = scv.dofIndex();
    }

    const VolumeVariables& operator [](const SubControlVolume& scv) const
    { return volumeVariables_[getLocalIdx_(scv.dofIndex())]; }

    VolumeVariables& operator [](const SubControlVolume& scv)
    { return volumeVariables_[getLocalIdx_(scv.dofIndex())]; }

    const VolumeVariables& operator [](IndexType scvIdx) const
    { return volumeVariables_[getLocalIdx_(scvIdx)]; }

    VolumeVariables& operator [](IndexType scvIdx)
    { return volumeVariables_[getLocalIdx_(scvIdx)]; }

    //! The global volume variables object we are a restriction of
    const GlobalVolumeVariables& globalVolVars() const
    { return *globalVolVarsPtr_; }

    //! Clear all local storage
    void clear()
    {
        volVarIndices_.clear();
        volumeVariables_.clear();
    }

private:
    const GlobalVolumeVariables* globalVolVarsPtr_;

    const int getLocalIdx_(const int volVarIdx) const
    {
        auto it = std::find(volVarIndices_.begin(), volVarIndices_.end(), volVarIdx);
        assert(it != volVarIndices_.end() && "Could not find the current volume variables for volVarIdx!");
        return std::distance(volVarIndices_.begin(), it);
    }

    std::vector<IndexType> volVarIndices_;
    std::vector<VolumeVariables> volumeVariables_;
};

} // end namespace

#endif
