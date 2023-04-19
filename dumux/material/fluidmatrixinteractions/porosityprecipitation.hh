// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup Fluidmatrixinteractions
 * \brief Class for the evaluation of the porosity subject to precipitation.
 */
#ifndef DUMUX_POROSITY_PRECIPITATION_HH
#define DUMUX_POROSITY_PRECIPITATION_HH

#include <dumux/discretization/evalsolution.hh>

namespace Dumux {

/*!
 * \ingroup Fluidmatrixinteractions
 * \brief Calculates the porosity depending on the volume fractions of precipitated minerals.
 *
 * \tparam Scalar The type used for scalar values
 * \param numComp The number of components in the fluid phases
 * \param numSolidPhases The number of precipitating solid phases
 */
template<class Scalar, int numComp, int numSolidPhases>
class PorosityPrecipitation
{
public:
    /*!
     * \brief Calculates the porosity in a sub-control volume
     * \param element Element
     * \param elemSol The element solution
     * \param scv Sub control volume
     * \param refPoro The solid matrix porosity without precipitates
     * \param minPoro A minimum porosity value
     */
    template<class Element, class SubControlVolume, class ElemSol>
    Scalar evaluatePorosity(const Element& element,
                            const SubControlVolume& scv,
                            const ElemSol& elemSol,
                            Scalar refPoro,
                            Scalar minPoro = 0.0) const
    {
        auto priVars = evalSolution(element, element.geometry(), elemSol, scv.center());

        Scalar sumPrecipitates = 0.0;
        for (unsigned int solidPhaseIdx = 0; solidPhaseIdx < numSolidPhases; ++solidPhaseIdx)
            sumPrecipitates += priVars[numComp + solidPhaseIdx];

        using std::max;
        return max(minPoro, refPoro - sumPrecipitates);
    }
};

} // namespace Dumux

#endif
