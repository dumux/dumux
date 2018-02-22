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
 * \brief The local volume variables class
 */
#ifndef DUMUX_DISCRETIZATION_BOX_ELEMENT_VOLUMEVARIABLES_HH
#define DUMUX_DISCRETIZATION_BOX_ELEMENT_VOLUMEVARIABLES_HH

#include <dumux/common/properties.hh>
#include <dumux/discretization/methods.hh>

namespace Dumux
{

/*!
 * \ingroup ImplicitModel
 * \brief Base class for the volume variables vector
 */
template<class TypeTag, bool enableGridVolVarsCache>
class BoxElementVolumeVariables
{};

// specialization in case of storing the volume variables
template<class TypeTag>
class BoxElementVolumeVariables<TypeTag,/*enableGlobalVolVarCache*/true>
{
    using Problem = typename GET_PROP_TYPE(TypeTag, Problem);
    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
    using SolutionVector = typename GET_PROP_TYPE(TypeTag, SolutionVector);
    using VolumeVariables = typename GET_PROP_TYPE(TypeTag, VolumeVariables);
    using GridVolumeVariables = typename GET_PROP_TYPE(TypeTag, GridVolumeVariables);
    using FVElementGeometry = typename GET_PROP_TYPE(TypeTag, FVGridGeometry)::LocalView;
    using SubControlVolume = typename FVElementGeometry::SubControlVolume;
    using IndexType = typename GridView::IndexSet::IndexType;

    static const int dim = GridView::dimension;
    using Element = typename GridView::template Codim<0>::Entity;

    static constexpr bool isBox = GET_PROP_VALUE(TypeTag, DiscretizationMethod) == DiscretizationMethod::box;

public:
    //! Constructor
    BoxElementVolumeVariables(const GridVolumeVariables& gridVolVars)
    : gridVolVarsPtr_(&gridVolVars) {}

    const VolumeVariables& operator [](IndexType scvIdx) const
    { return gridVolVars().volVars(eIdx_, scvIdx); }

    const VolumeVariables& operator [](const SubControlVolume& scv) const
    { return gridVolVars().volVars(eIdx_, scv.indexInElement()); }

    // For compatibility reasons with the case of not storing the vol vars.
    // function to be called before assembling an element, preparing the vol vars within the stencil
    void bind(const Element& element,
              const FVElementGeometry& fvGeometry,
              const SolutionVector& sol)
    {
        bindElement(element, fvGeometry, sol);
    }

    // function to prepare the vol vars within the element
    void bindElement(const Element& element,
                     const FVElementGeometry& fvGeometry,
                     const SolutionVector& sol)
    {
        eIdx_ = fvGeometry.fvGridGeometry().elementMapper().index(element);
    }

    //! The global volume variables object we are a restriction of
    const GridVolumeVariables& gridVolVars() const
    { return *gridVolVarsPtr_; }

private:
    const GridVolumeVariables* gridVolVarsPtr_;
    IndexType eIdx_;
};


// Specialization when the current volume variables are not stored
template<class TypeTag>
class BoxElementVolumeVariables<TypeTag, /*enableGlobalVolVarCache*/false>
{
    using Problem = typename GET_PROP_TYPE(TypeTag, Problem);
    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
    using SolutionVector = typename GET_PROP_TYPE(TypeTag, SolutionVector);
    using VolumeVariables = typename GET_PROP_TYPE(TypeTag, VolumeVariables);
    using GridVolumeVariables = typename GET_PROP_TYPE(TypeTag, GridVolumeVariables);
    using FVElementGeometry = typename GET_PROP_TYPE(TypeTag, FVGridGeometry)::LocalView;
    using SubControlVolume = typename FVElementGeometry::SubControlVolume;
    using IndexType = typename GridView::IndexSet::IndexType;
    using ElementSolutionVector = typename GET_PROP_TYPE(TypeTag, ElementSolutionVector);

    static const int dim = GridView::dimension;
    using Element = typename GridView::template Codim<0>::Entity;

public:
    //! Constructor
    BoxElementVolumeVariables(const GridVolumeVariables& gridVolVars)
    : gridVolVarsPtr_(&gridVolVars) {}

    // specialization for box models, simply forwards to the bindElement method
    void bind(const Element& element,
              const FVElementGeometry& fvGeometry,
              const SolutionVector& sol)
    {
        bindElement(element, fvGeometry, sol);
    }

    // specialization for box models
    void bindElement(const Element& element,
                     const FVElementGeometry& fvGeometry,
                     const SolutionVector& sol)
    {
        // get the solution at the dofs of the element
        ElementSolutionVector elemSol(element, sol, fvGeometry);

        // resize volume variables to the required size
        volumeVariables_.resize(fvGeometry.numScv());
        for (auto&& scv : scvs(fvGeometry))
        {
            // TODO: INTERFACE SOLVER
            // gridVolVars().problem().model().boxInterfaceConditionSolver().updateScvVolVars(element, scv, sol);
            volumeVariables_[scv.indexInElement()].update(elemSol, gridVolVars().problem(), element, scv);
        }
    }

    const VolumeVariables& operator [](IndexType scvIdx) const
    { return volumeVariables_[scvIdx]; }

    VolumeVariables& operator [](IndexType scvIdx)
    { return volumeVariables_[scvIdx]; }

    const VolumeVariables& operator [](const SubControlVolume& scv) const
    { return volumeVariables_[scv.indexInElement()]; }

    VolumeVariables& operator [](const SubControlVolume& scv)
    { return volumeVariables_[scv.indexInElement()]; }

    //! The global volume variables object we are a restriction of
    const GridVolumeVariables& gridVolVars() const
    { return *gridVolVarsPtr_; }

private:
    const GridVolumeVariables* gridVolVarsPtr_;
    std::vector<VolumeVariables> volumeVariables_;
};

} // end namespace

#endif
