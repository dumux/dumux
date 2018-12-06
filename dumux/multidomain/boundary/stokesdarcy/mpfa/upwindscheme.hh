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
 * \ingroup Discretization
 * \brief TODO
 */
#ifndef DUMUX_STOKESDARCY_CCMPFA_UPWINDSCHEME_HH
#define DUMUX_STOKESDARCY_CCMPFA_UPWINDSCHEME_HH

#include <dumux/common/parameters.hh>

namespace Dumux {

//! Upwind scheme for the cell-centered tpfa scheme
template<class FVGridGeometry>
class CCMpfaStokesDarcyUpwindScheme
{

public:

    // For grids with dim == dimWorld we use a simple upwinding scheme
    template<class FluxVariables, class UpwindTermFunction, class Scalar>
    static Scalar
    apply(const FluxVariables& fluxVars,
          const UpwindTermFunction& upwindTerm,
          Scalar flux, int phaseIdx)
    {
        static const Scalar upwindWeight = getParam<Scalar>("Implicit.UpwindWeight");

        const auto& scvf = fluxVars.scvFace();
        const auto& elemVolVars = fluxVars.elemVolVars();
        const auto& insideVolVars = elemVolVars[scvf.insideScvIdx()];

        // redefine upwind term to work with navier stokes vol vars (have no mobility() function)
        auto usedUpwindTerm = [phaseIdx] (const auto& volVars) { return volVars.density(phaseIdx)/volVars.viscosity(phaseIdx); };

        using CouplingManager = std::decay_t<decltype(fluxVars.problem().couplingManager())>;
        if (fluxVars.problem().couplingManager().isCoupledEntity(CouplingManager::darcyIdx, scvf))
        {
            const auto& outsideVolVars = fluxVars.problem().couplingManager().darcyCouplingContext(scvf).volVars;

            if (std::signbit(flux)) // if sign of flux is negative
                return flux*(upwindWeight*usedUpwindTerm(outsideVolVars)
                             + (1.0 - upwindWeight)*usedUpwindTerm(insideVolVars));
            else
                return flux*(upwindWeight*usedUpwindTerm(insideVolVars)
                             + (1.0 - upwindWeight)*usedUpwindTerm(outsideVolVars));
        }
        else
        {
            const auto& outsideVolVars = elemVolVars[scvf.outsideScvIdx()];

            if (std::signbit(flux)) // if sign of flux is negative
                return flux*(upwindWeight*usedUpwindTerm(outsideVolVars)
                             + (1.0 - upwindWeight)*usedUpwindTerm(insideVolVars));
            else
                return flux*(upwindWeight*usedUpwindTerm(insideVolVars)
                             + (1.0 - upwindWeight)*usedUpwindTerm(outsideVolVars));
        }
    }
};

} // end namespace Dumux

#endif
