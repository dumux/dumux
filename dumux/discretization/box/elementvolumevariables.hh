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

#include <dumux/implicit/properties.hh>
#include <dumux/implicit/model.hh>

namespace Dumux
{

/*!
 * \ingroup ImplicitModel
 * \brief Base class for the volume variables vector
 */
template<class TypeTag, bool enableGlobalVolVarsCache>
class BoxElementVolumeVariables
{};

// specialization in case of storing the volume variables
template<class TypeTag>
class BoxElementVolumeVariables<TypeTag,/*enableGlobalVolVarCache*/true>
{
    friend ImplicitModel<TypeTag>;
    using Problem = typename GET_PROP_TYPE(TypeTag, Problem);
    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
    using SolutionVector = typename GET_PROP_TYPE(TypeTag, SolutionVector);
    using VolumeVariables = typename GET_PROP_TYPE(TypeTag, VolumeVariables);
    using GlobalVolumeVariables = typename GET_PROP_TYPE(TypeTag, GlobalVolumeVariables);
    using SubControlVolume = typename GET_PROP_TYPE(TypeTag, SubControlVolume);
    using IndexType = typename GridView::IndexSet::IndexType;
    using FVElementGeometry = typename GET_PROP_TYPE(TypeTag, FVElementGeometry);

    static const int dim = GridView::dimension;
    using Element = typename GridView::template Codim<0>::Entity;

    enum{ isBox = GET_PROP_VALUE(TypeTag, ImplicitIsBox) };

public:
    //! Constructor
    BoxElementVolumeVariables(const GlobalVolumeVariables& globalVolVars)
    : globalVolVarsPtr_(&globalVolVars) {}

    const VolumeVariables& operator [](IndexType scvIdx) const
    { return globalVolVars().volVars(eIdx_, scvIdx); }

    const VolumeVariables& operator [](const SubControlVolume& scv) const
    { return globalVolVars().volVars(eIdx_, scv.indexInElement()); }

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
        eIdx_ = globalVolVars().problem_().elementMapper().index(element);
    }

    //! The global volume variables object we are a restriction of
    const GlobalVolumeVariables& globalVolVars() const
    { return *globalVolVarsPtr_; }

private:
    const GlobalVolumeVariables* globalVolVarsPtr_;
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
    using GlobalVolumeVariables = typename GET_PROP_TYPE(TypeTag, GlobalVolumeVariables);
    using SubControlVolume = typename GET_PROP_TYPE(TypeTag, SubControlVolume);
    using IndexType = typename GridView::IndexSet::IndexType;
    using FVElementGeometry = typename GET_PROP_TYPE(TypeTag, FVElementGeometry);

    static const int dim = GridView::dimension;
    using Element = typename GridView::template Codim<0>::Entity;

public:
    //! Constructor
    BoxElementVolumeVariables(const GlobalVolumeVariables& globalVolVars)
    : globalVolVarsPtr_(&globalVolVars) {}

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
        auto elemSol = globalVolVars().problem_().model().elementSolution(element, sol);

        // resize volume variables to the required size
        volumeVariables_.resize(fvGeometry.numScv());
        for (auto&& scv : scvs(fvGeometry))
        {
            // TODO: INTERFACE SOLVER
            // globalVolVars().problem_().model().boxInterfaceConditionSolver().updateScvVolVars(element, scv, sol);
            volumeVariables_[scv.indexInElement()].update(elemSol, globalVolVars().problem_(), element, scv);
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
    const GlobalVolumeVariables& globalVolVars() const
    { return *globalVolVarsPtr_; }

private:
    const GlobalVolumeVariables* globalVolVarsPtr_;
    std::vector<VolumeVariables> volumeVariables_;
};

} // end namespace

#endif
