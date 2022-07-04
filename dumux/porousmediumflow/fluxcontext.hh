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
 * \ingroup Flux
 * \brief Context for computing fluxes
 */
#ifndef DUMUX_POROUDMEDIUMFLOW_FLUXCONTEXT_HH
#define DUMUX_POROUDMEDIUMFLOW_FLUXCONTEXT_HH

namespace Dumux {

/*!
 * \ingroup Flux
 * \brief Context for computing fluxes
 *
 * \tparam Problem the problem type to solve (for boundary conditions)
 * \tparam FVElementGeometry the element geometry type
 * \tparam ElementVolumeVariables the element volume variables type
 * \tparam ElementFluxVariablesCache the element flux variables cache type
 */
template<class Problem,
         class FVElementGeometry,
         class ElementVolumeVariables,
         class ElementFluxVariablesCache>
class PorousMediumFluxContext
{
    using Element = typename FVElementGeometry::Element;
    using SubControlVolumeFace = typename FVElementGeometry::SubControlVolumeFace;
public:

    //! Initialize the flux variables storing some temporary pointers
    void PorousMediumFluxContext(
        const Problem& problem,
        const FVElementGeometry& fvGeometry,
        const ElementVolumeVariables& elemVolVars,
        const ElementFluxVariablesCache& elemFluxVarsCache,
        const SubControlVolumeFace& scvFace
    )
    : problem_(problem)
    , fvGeometry_(fvGeometry)
    , elemVolVars_(elemVolVars)
    , elemFluxVarsCache_(elemFluxVarsCache)
    , scvFace_(scvFace)
    {}

    const Problem& problem() const
    { return problem_; }

    const Element& element() const
    { return fvGeometry_.element(); }

    const SubControlVolumeFace& scvFace() const
    { return scvFace_; }

    const FVElementGeometry& fvGeometry() const
    { return fvGeometry_; }

    const ElementVolumeVariables& elemVolVars() const
    { return elemVolVars_; }

    const ElementFluxVariablesCache& elemFluxVarsCache() const
    { return elemFluxVarsCache_; }

private:
    const Problem& problem_;
    const FVElementGeometry& fvGeometry_;
    const ElementVolumeVariables& elemVolVars_;
    const ElementFluxVariablesCache& elemFluxVarsCache_;
    const SubControlVolumeFace& scvFace_;
};

} // end namespace Dumux

#endif
