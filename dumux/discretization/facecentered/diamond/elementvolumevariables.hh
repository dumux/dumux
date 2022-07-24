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
 * \ingroup DiamondDiscretization
 * \copydoc Dumux::FaceCenteredDiamondElementVolumeVariables
 */
#ifndef DUMUX_DISCRETIZATION_FACECENTERED_DIAMOND_ELEMENTVOLUMEVARIABLES_HH
#define DUMUX_DISCRETIZATION_FACECENTERED_DIAMOND_ELEMENTVOLUMEVARIABLES_HH

#include <algorithm>
#include <cassert>
#include <vector>

namespace Dumux {

/*!
 * \ingroup DiamondDiscretization
 * \brief Base class for the face variables vector
 */
template<class GridVolumeVariables, bool cachingEnabled>
class FaceCenteredDiamondElementVolumeVariables
{};

/*!
 * \ingroup DiamondDiscretization
 * \brief Class for the face variables vector. Specialization for the case of storing the face variables globally.
 */
template<class GFV>
class FaceCenteredDiamondElementVolumeVariables<GFV, /*cachingEnabled*/true>
{
public:
    //! export type of the grid volume variables
    using GridVolumeVariables = GFV;

    //! export type of the volume variables
    using VolumeVariables = typename GridVolumeVariables::VolumeVariables;

    FaceCenteredDiamondElementVolumeVariables(const GridVolumeVariables& gridVolumeVariables)
    : gridVolumeVariablesPtr_(&gridVolumeVariables)
    {}

    template<class SubControlVolume, typename std::enable_if_t<!std::is_integral<SubControlVolume>::value, int> = 0>
    const VolumeVariables& operator [](const SubControlVolume& scv) const
    { return gridVolVars().volVars(eIdx_, scv.indexInElement()); }

    const VolumeVariables& operator [](const std::size_t scvIdx) const
    { return gridVolVars().volVars(eIdx_, scvIdx); }

    template<class FVElementGeometry, class SolutionVector>
    void bind(const typename FVElementGeometry::GridGeometry::GridView::template Codim<0>::Entity& element,
              const FVElementGeometry& fvGeometry,
              const SolutionVector& sol) &
    {
        bindElement(element, fvGeometry, sol);
    }

    template<class FVElementGeometry, class SolutionVector>
    void bindElement(const typename FVElementGeometry::GridGeometry::GridView::template Codim<0>::Entity& element,
                     const FVElementGeometry& fvGeometry,
                     const SolutionVector& sol) &
    {
        eIdx_ = fvGeometry.elementIndex();
    }

    template<class FVElementGeometry, class SolutionVector>
    FaceCenteredDiamondElementVolumeVariables
    bind(const typename FVElementGeometry::GridGeometry::GridView::template Codim<0>::Entity& element,
              const FVElementGeometry& fvGeometry,
              const SolutionVector& sol) &&
    {
        this->bind(element, fvGeometry, sol);
        return std::move(*this);
    }

    template<class FVElementGeometry, class SolutionVector>
    FaceCenteredDiamondElementVolumeVariables
    bindElement(const typename FVElementGeometry::GridGeometry::GridView::template Codim<0>::Entity& element,
                     const FVElementGeometry& fvGeometry,
                     const SolutionVector& sol) &&
    {
        this->bindElement(element, fvGeometry, sol);
        return std::move(*this);
    }

    //! The global volume variables object we are a restriction of
    const GridVolumeVariables& gridVolVars() const
    { return *gridVolumeVariablesPtr_; }

private:
    const GridVolumeVariables* gridVolumeVariablesPtr_;
    std::size_t eIdx_;
};

} // end namespace Dumux

#endif
