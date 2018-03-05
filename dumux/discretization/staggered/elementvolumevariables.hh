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
 * \copydoc Dumux::StaggeredElementVolumeVariables
 */
#ifndef DUMUX_DISCRETIZATION_STAGGERED_ELEMENT_VOLUMEVARIABLES_HH
#define DUMUX_DISCRETIZATION_STAGGERED_ELEMENT_VOLUMEVARIABLES_HH

#include <dune/common/exceptions.hh>
#include <dumux/common/properties.hh>

namespace Dumux
{

/*!
 * \ingroup StaggeredDiscretization
 * \brief Base class for the element volume variables vector for the staggered model
 */
template<class TypeTag, bool enableGridVolVarsCache>
class StaggeredElementVolumeVariables
{};

/*!
 * \ingroup StaggeredDiscretization
 * \brief Class for the element volume variables vector for the staggered model.
          Specialization in case the volume variables are stored globally.
 */
template<class TypeTag>
class StaggeredElementVolumeVariables<TypeTag, /*enableGridVolVarsCache*/true>
{
    using Problem = typename GET_PROP_TYPE(TypeTag, Problem);
    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
    using FVElementGeometry = typename GET_PROP_TYPE(TypeTag, FVGridGeometry)::LocalView;
    using VolumeVariables = typename GET_PROP_TYPE(TypeTag, VolumeVariables);
    using GridVolumeVariables = typename GET_PROP_TYPE(TypeTag, GridVolumeVariables);
    using SubControlVolume = typename FVElementGeometry::SubControlVolume;
    using IndexType = typename GridView::IndexSet::IndexType;

    static const int dim = GridView::dimension;
    using Element = typename GridView::template Codim<0>::Entity;

public:
    //! Export type of the solution vector
    using SolutionVector = typename GET_PROP_TYPE(TypeTag, SolutionVector);

    //! Constructor
    StaggeredElementVolumeVariables(const GridVolumeVariables& gridVolVars)
    : gridVolVarsPtr_(&gridVolVars) {}

    //! operator for the access with an scv
    const VolumeVariables& operator [](const SubControlVolume& scv) const
    { return gridVolVars().volVars(scv.dofIndex()); }

    //! operator for the access with an index
    //! needed for Staggered methods for the access to the boundary volume variables
    const VolumeVariables& operator [](const IndexType scvIdx) const
    { return gridVolVars().volVars(scvIdx); }

    //! For compatibility reasons with the case of not storing the vol vars.
    //! function to be called before assembling an element, preparing the vol vars within the stencil
    void bind(const Element& element,
              const FVElementGeometry& fvGeometry,
              const SolutionVector& sol)
    {}

    //! function to prepare the vol vars within the element
    void bindElement(const Element& element,
                     const FVElementGeometry& fvGeometry,
                     const SolutionVector& sol)
    {}

    //! The global volume variables object we are a restriction of
    const GridVolumeVariables& gridVolVars() const
    { return *gridVolVarsPtr_; }

private:
    const GridVolumeVariables* gridVolVarsPtr_;
};


/*!
 * \ingroup StaggeredDiscretization
 * \brief Class for the element volume variables vector for the staggered model.
          Specialization in case the volume variables are not stored globally.
 */
template<class TypeTag>
class StaggeredElementVolumeVariables<TypeTag, /*enableGridVolVarsCache*/false>
{
    using Problem = typename GET_PROP_TYPE(TypeTag, Problem);
    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
    using Indices = typename GET_PROP_TYPE(TypeTag, Indices);
    using ElementSolution = typename GET_PROP_TYPE(TypeTag, ElementSolutionVector);
    using VolumeVariables = typename GET_PROP_TYPE(TypeTag, VolumeVariables);
    using CellCenterPrimaryVariables = typename GET_PROP_TYPE(TypeTag, CellCenterPrimaryVariables);
    using GridVolumeVariables = typename GET_PROP_TYPE(TypeTag, GridVolumeVariables);
    using FVElementGeometry = typename GET_PROP_TYPE(TypeTag, FVGridGeometry)::LocalView;
    using SubControlVolume = typename FVElementGeometry::SubControlVolume;
    using IndexType = typename GridView::IndexSet::IndexType;

    static const int dim = GridView::dimension;
    using Element = typename GridView::template Codim<0>::Entity;

    using DofTypeIndices = typename GET_PROP(TypeTag, DofTypeIndices);
    typename DofTypeIndices::CellCenterIdx cellCenterIdx;
    typename DofTypeIndices::FaceIdx faceIdx;

    static constexpr auto numEqCellCenter = GET_PROP_VALUE(TypeTag, NumEqCellCenter);

public:
    //! Export type of the solution vector
    using SolutionVector = typename GET_PROP_TYPE(TypeTag, SolutionVector);

    //! Constructor
    StaggeredElementVolumeVariables(const GridVolumeVariables& gridVolVars)
    : gridVolVarsPtr_(&gridVolVars) {}

    //! Binding of an element, prepares the volume variables within the element stencil
    //! called by the local jacobian to prepare element assembly
    void bind(const Element& element,
              const FVElementGeometry& fvGeometry,
              const SolutionVector& sol)
    {
        clear();

        const auto& problem = gridVolVars().problem();
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
            CellCenterPrimaryVariables priVars(0.0);
            priVars = sol[cellCenterIdx][globalJ];
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

            CellCenterPrimaryVariables boundaryPriVars(0.0);

            for(int eqIdx = 0; eqIdx < numEqCellCenter; ++eqIdx)
            {
                if(bcTypes.isDirichlet(eqIdx) || bcTypes.isDirichletCell(eqIdx))
                    boundaryPriVars[eqIdx] = problem.dirichlet(element, scvf)[cellCenterIdx][eqIdx];
                else if(bcTypes.isNeumann(eqIdx) || bcTypes.isOutflow(eqIdx) || bcTypes.isSymmetry())
                    boundaryPriVars[eqIdx] = sol[cellCenterIdx][scvf.insideScvIdx()][eqIdx];
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

    //! Binding of an element, prepares only the volume variables of the element.
    //! Specialization for Staggered models
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
        CellCenterPrimaryVariables priVars(0.0);
        priVars = sol[cellCenterIdx][eIdx];
        ElementSolution elemSol{std::move(priVars)};
        volumeVariables_[0].update(elemSol,
                                   gridVolVars().problem(),
                                   element,
                                   scv);
        volVarIndices_[0] = scv.dofIndex();
    }

    //! const operator for the access with an scv
    const VolumeVariables& operator [](const SubControlVolume& scv) const
    { return volumeVariables_[getLocalIdx_(scv.dofIndex())]; }

    //! operator for the access with an scv
    VolumeVariables& operator [](const SubControlVolume& scv)
    { return volumeVariables_[getLocalIdx_(scv.dofIndex())]; }

    //! const operator for the access with an index
    const VolumeVariables& operator [](IndexType scvIdx) const
    { return volumeVariables_[getLocalIdx_(scvIdx)]; }

    //! operator for the access with an index
    VolumeVariables& operator [](IndexType scvIdx)
    { return volumeVariables_[getLocalIdx_(scvIdx)]; }

    //! The global volume variables object we are a restriction of
    const GridVolumeVariables& gridVolVars() const
    { return *gridVolVarsPtr_; }

    //! Clear all local storage
    void clear()
    {
        volVarIndices_.clear();
        volumeVariables_.clear();
    }

private:
    const GridVolumeVariables* gridVolVarsPtr_;

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
