// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup PNMSolidEnergyModel
 * \brief Flux variables cache for PNM solid-energy model
 */
#ifndef DUMUX_PNM_SOLID_ENERGY_FLUXVARIABLESCACHE_HH
#define DUMUX_PNM_SOLID_ENERGY_FLUXVARIABLESCACHE_HH

namespace Dumux::PoreNetwork {

template<class Scalar>
class SolidEnergyFluxVariablesCache
{
public:
    //! whether the cache needs an update when the solution changes
    static bool constexpr isSolDependent = false;

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
