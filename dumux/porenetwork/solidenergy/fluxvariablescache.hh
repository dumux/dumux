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
 * \brief Flux variables cache for 1p PNM
 */
#ifndef DUMUX_PNM_SOLID_ENERGY_FLUXVARIABLESCACHE_HH
#define DUMUX_PNM_SOLID_ENERGY_FLUXVARIABLESCACHE_HH

namespace Dumux::PoreNetwork {

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
//! The cache is dependent on the active physical processes (advection, diffusion, heat conduction)
//! For each type of process there is a base cache storing the data required to compute the respective fluxes
//! Specializations of the overall cache are provided for combinations of processes
///////////////////////////////////////////////////////////////////////////////////////////////////////////////


/*!
 * \ingroup ImplicitModel
 * \brief The flux variables cache classes for porous media.
 *        Store data required for flux calculation. For each type of physical process (advection, diffusion, heat conduction)
 *        there is a base cache storing the data required to compute the respective fluxes. Specializations of the overall
 *        cache class are provided for different combinations of processes.
 */

//! We only store discretization-related quantities for the box method.
//! Thus, we need no physics-dependent specialization.
template<class Scalar>
class SolidEnergyFluxVariablesCache
{
public:

    template<class Problem, class Element, class FVElementGeometry,
             class ElementVolumeVariables, class SubControlVolumeFace>
    void update(const Problem& problem,
                const Element& element,
                const FVElementGeometry& fvGeometry,
                const ElementVolumeVariables& elemVolVars,
                const SubControlVolumeFace& scvf)
    {
        grainContactArea_ = problem.spatialParams().throatCrossSectionalArea(element, elemVolVars);
        throatLength_ = problem.spatialParams().throatLength(element, elemVolVars);
        throatInscribedRadius_ = problem.spatialParams().throatInscribedRadius(element, elemVolVars);
    }

    Scalar grainContactArea() const
    { return grainContactArea_; }

    Scalar throatLength() const
    { return throatLength_; }

    Scalar throatInscribedRadius() const
    { return throatInscribedRadius_; }

private:
    Scalar grainContactArea_;
    Scalar throatLength_;
    Scalar throatInscribedRadius_;
};

} // end namespace Dumux::PoreNetwork

#endif
