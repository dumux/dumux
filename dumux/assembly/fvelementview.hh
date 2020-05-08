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
 * \ingroup Assembly
 * \brief Element-centered views on the local context needed for assembly routines
 */
#ifndef DUMUX_ASSEMBLY_FV_ELEMENT_VIEW_HH
#define DUMUX_ASSEMBLY_FV_ELEMENT_VIEW_HH

#include <dumux/discretization/localview.hh>

namespace Dumux {

template<class EG, class EVV, class EFC>
class FVElementView
{
public:
    using Element = typename EG::Element;
    using FVElementGeometry = EG;
    using ElementVolumeVariables = EVV;
    using ElementFluxVariablesCache = EFC;
    using Problem = typename ElementVolumeVariables::GridVolumeVariables::Problem;

    FVElementView(const Element& element,
                  FVElementGeometry& fvGeometry,
                  ElementVolumeVariables& elemVolVars,
                  ElementFluxVariablesCache& elemFluxVarCache)
    : element_(element)
    , fvGeometry_(fvGeometry)
    , elemVolVars_(elemVolVars)
    , elemFluxVarCache_(elemFluxVarCache)
    {}

    template<class SolutionVector>
    void bind(const SolutionVector& sol)
    {
        fvGeometry_.bind(element_);
        elemVolVars_.bind(element_, fvGeometry_, sol);
        elemFluxVarCache_.bind(element_, fvGeometry_, elemVolVars_);
    }

    FVElementGeometry& fvGeometry()
    { return fvGeometry_; }

    ElementVolumeVariables& elemVolVars()
    { return elemVolVars_; }

    ElementFluxVariablesCache& elemFluxVarsCache()
    { return elemFluxVarCache_; }

    const Problem& problem() const
    { return elemVolVars_.gridVolVars().problem(); }

    const Element& element() const
    { return element_; }

private:
    const Element& element_;
    FVElementGeometry& fvGeometry_;
    ElementVolumeVariables& elemVolVars_;
    ElementFluxVariablesCache& elemFluxVarCache_;
};

template<class Element, class GridGeometry, class FVGridVariables>
auto makeElementView(const Element& element, GridGeometry& gg, FVGridVariables& gridVars)
{
    return FVElementView{element, localView(gg), localView(gridVars.curGridVolVars()), localView(gridVars.fluxVarsCache())};
}

} // end namespace Dumux

#endif
