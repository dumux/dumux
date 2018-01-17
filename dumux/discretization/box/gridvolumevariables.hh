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
 * \brief The grid volume variables class for box models
 */
#ifndef DUMUX_DISCRETIZATION_BOX_GRID_VOLUMEVARIABLES_HH
#define DUMUX_DISCRETIZATION_BOX_GRID_VOLUMEVARIABLES_HH

#include <dumux/common/properties.hh>
#include <dumux/discretization/box/elementvolumevariables.hh>

//! make the local view function available whenever we use this class
#include <dumux/discretization/localview.hh>

namespace Dumux {

/*!
 * \ingroup BoxModel
 * \brief Base class for the grid volume variables
 */
template<class TypeTag, bool enableGridVolVarsCache>
class BoxGridVolumeVariables;

// specialization in case of storing the volume variables
template<class TypeTag>
class BoxGridVolumeVariables<TypeTag,/*enableGlobalVolVarCache*/true>
{
    using Problem = typename GET_PROP_TYPE(TypeTag, Problem);
    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
    using SolutionVector = typename GET_PROP_TYPE(TypeTag, SolutionVector);
    using VolumeVariables = typename GET_PROP_TYPE(TypeTag, VolumeVariables);
    using FVElementGeometry = typename GET_PROP_TYPE(TypeTag, FVGridGeometry)::LocalView;
    using SubControlVolume = typename FVElementGeometry::SubControlVolume;
    using IndexType = typename GridView::IndexSet::IndexType;
    using FVGridGeometry = typename GET_PROP_TYPE(TypeTag, FVGridGeometry);
    using ElementSolutionVector = typename GET_PROP_TYPE(TypeTag, ElementSolutionVector);

    static const int dim = GridView::dimension;
    using Element = typename GridView::template Codim<0>::Entity;

public:
    //! export the type of the local view
    using LocalView = typename GET_PROP_TYPE(TypeTag, ElementVolumeVariables);

    BoxGridVolumeVariables(const Problem& problem) : problemPtr_(&problem) {}

    void update(const FVGridGeometry& fvGridGeometry, const SolutionVector& sol)
    {
        volumeVariables_.resize(fvGridGeometry.gridView().size(0));
        for (const auto& element : elements(fvGridGeometry.gridView()))
        {
            auto eIdx = fvGridGeometry.elementMapper().index(element);

            auto fvGeometry = localView(fvGridGeometry);
            fvGeometry.bindElement(element);

            // get the element solution
            ElementSolutionVector elemSol(element, sol, fvGeometry);

            // update the volvars of the element
            volumeVariables_[eIdx].resize(fvGeometry.numScv());
            for (auto&& scv : scvs(fvGeometry))
                volumeVariables_[eIdx][scv.indexInElement()].update(elemSol, problem(), element, scv);
        }
    }

    const VolumeVariables& volVars(const SubControlVolume& scv) const
    { return volumeVariables_[scv.elementIndex()][scv.indexInElement()]; }

    VolumeVariables& volVars(const SubControlVolume& scv)
    { return volumeVariables_[scv.elementIndex()][scv.indexInElement()]; }

    const VolumeVariables& volVars(const IndexType eIdx, const IndexType scvIdx) const
    { return volumeVariables_[eIdx][scvIdx]; }

    VolumeVariables& volVars(const IndexType eIdx, const IndexType scvIdx)
    { return volumeVariables_[eIdx][scvIdx]; }

    const Problem& problem() const
    { return *problemPtr_; }

private:
    const Problem* problemPtr_;
    std::vector<std::vector<VolumeVariables>> volumeVariables_;
};


// Specialization when the current volume variables are not stored
template<class TypeTag>
class BoxGridVolumeVariables<TypeTag, /*enableGlobalVolVarCache*/false>
{
    using Problem = typename GET_PROP_TYPE(TypeTag, Problem);
    using SolutionVector = typename GET_PROP_TYPE(TypeTag, SolutionVector);
    using VolumeVariables = typename GET_PROP_TYPE(TypeTag, VolumeVariables);
    using FVGridGeometry = typename GET_PROP_TYPE(TypeTag, FVGridGeometry);

public:
    //! export the type of the local view
    using LocalView = typename GET_PROP_TYPE(TypeTag, ElementVolumeVariables);

    BoxGridVolumeVariables(const Problem& problem) : problemPtr_(&problem) {}

    void update(const FVGridGeometry& fvGridGeometry, const SolutionVector& sol) {}

    const Problem& problem() const
    { return *problemPtr_;}

private:
    const Problem* problemPtr_;
};

} // end namespace Dumux

#endif
