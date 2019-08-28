// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   See the file COPYING for full copying permissions.                      *
 *                                                                           *
 *   This program is free software: you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation, either version 3 of the License, or       *
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
 * \ingroup BoxDiscretization
 * \brief The local volume variables class
 */
#ifndef DUMUX_DISCRETIZATION_BOX_ELEMENT_VOLUMEVARIABLES_HH
#define DUMUX_DISCRETIZATION_BOX_ELEMENT_VOLUMEVARIABLES_HH

#include <type_traits>

#include <dumux/discretization/box/elementsolution.hh>

namespace Dumux {

/*!
 * \ingroup BoxDiscretization
 * \brief The local (stencil) volume variables class for box models
 * \note The class is specilized for versions with and without caching
 * \tparam GVV the grid volume variables type
 * \tparam cachingEnabled if the cache is enabled
 */
template<class GVV, bool cachingEnabled>
class BoxElementVolumeVariables;

/*!
 * \ingroup BoxDiscretization
 * \brief The local (stencil) volume variables class for box models with caching
 * \note the volume variables are stored for the whole grid view in the corresponding GridVolumeVariables class
 */
template<class GVV>
class BoxElementVolumeVariables<GVV, /*cachingEnabled*/true>
{
public:
    //! export type of the grid volume variables
    using GridVolumeVariables = GVV;

    //! export type of the volume variables
    using VolumeVariables = typename GridVolumeVariables::VolumeVariables;

    //! Constructor
    BoxElementVolumeVariables(const GridVolumeVariables& gridVolVars)
    : gridVolVarsPtr_(&gridVolVars) {}

    const VolumeVariables& operator [](std::size_t scvIdx) const
    { return gridVolVars().volVars(eIdx_, scvIdx); }

    template<class SubControlVolume, typename std::enable_if_t<!std::is_integral<SubControlVolume>::value, int> = 0>
    const VolumeVariables& operator [](const SubControlVolume& scv) const
    { return gridVolVars().volVars(eIdx_, scv.indexInElement()); }

    // For compatibility reasons with the case of not storing the vol vars.
    // function to be called before assembling an element, preparing the vol vars within the stencil
    template<class FVElementGeometry, class SolutionVector>
    void bind(const typename FVElementGeometry::GridGeometry::GridView::template Codim<0>::Entity& element,
              const FVElementGeometry& fvGeometry,
              const SolutionVector& sol)
    {
        bindElement(element, fvGeometry, sol);
    }

    // function to prepare the vol vars within the element
    template<class FVElementGeometry, class SolutionVector>
    void bindElement(const typename FVElementGeometry::GridGeometry::GridView::template Codim<0>::Entity& element,
                     const FVElementGeometry& fvGeometry,
                     const SolutionVector& sol)
    {
        eIdx_ = fvGeometry.gridGeometry().elementMapper().index(element);
    }

    //! The global volume variables object we are a restriction of
    const GridVolumeVariables& gridVolVars() const
    { return *gridVolVarsPtr_; }

private:
    const GridVolumeVariables* gridVolVarsPtr_;
    std::size_t eIdx_;
};


/*!
 * \ingroup BoxDiscretization
 * \brief The local (stencil) volume variables class for box models without caching
 */
template<class GVV>
class BoxElementVolumeVariables<GVV, /*cachingEnabled*/false>
{
public:
    //! export type of the grid volume variables
    using GridVolumeVariables = GVV;

    //! export type of the volume variables
    using VolumeVariables = typename GridVolumeVariables::VolumeVariables;

    //! Constructor
    BoxElementVolumeVariables(const GridVolumeVariables& gridVolVars)
    : gridVolVarsPtr_(&gridVolVars) {}

    // specialization for box models, simply forwards to the bindElement method
    template<class FVElementGeometry, class SolutionVector>
    void bind(const typename FVElementGeometry::GridGeometry::GridView::template Codim<0>::Entity& element,
              const FVElementGeometry& fvGeometry,
              const SolutionVector& sol)
    {
        bindElement(element, fvGeometry, sol);
    }

    // specialization for box models
    template<class FVElementGeometry, class SolutionVector>
    void bindElement(const typename FVElementGeometry::GridGeometry::GridView::template Codim<0>::Entity& element,
                     const FVElementGeometry& fvGeometry,
                     const SolutionVector& sol)
    {
        // get the solution at the dofs of the element
        auto elemSol = elementSolution(element, sol, fvGeometry.gridGeometry());

        // resize volume variables to the required size
        volumeVariables_.resize(fvGeometry.numScv());
        for (auto&& scv : scvs(fvGeometry))
            volumeVariables_[scv.indexInElement()].update(elemSol, gridVolVars().problem(), element, scv);
    }

    const VolumeVariables& operator [](std::size_t scvIdx) const
    { return volumeVariables_[scvIdx]; }

    VolumeVariables& operator [](std::size_t scvIdx)
    { return volumeVariables_[scvIdx]; }

    template<class SubControlVolume, typename std::enable_if_t<!std::is_integral<SubControlVolume>::value, int> = 0>
    const VolumeVariables& operator [](const SubControlVolume& scv) const
    { return volumeVariables_[scv.indexInElement()]; }

    template<class SubControlVolume, typename std::enable_if_t<!std::is_integral<SubControlVolume>::value, int> = 0>
    VolumeVariables& operator [](const SubControlVolume& scv)
    { return volumeVariables_[scv.indexInElement()]; }

    //! The global volume variables object we are a restriction of
    const GridVolumeVariables& gridVolVars() const
    { return *gridVolVarsPtr_; }

private:
    const GridVolumeVariables* gridVolVarsPtr_;
    std::vector<VolumeVariables> volumeVariables_;
};

} // end namespace Dumux

#endif
