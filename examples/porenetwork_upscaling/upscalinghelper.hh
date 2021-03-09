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

#ifndef DUMUX_PNM_ONEP_PERMEABILITY_UPSCALING_HELPER_HH
#define DUMUX_PNM_ONEP_PERMEABILITY_UPSCALING_HELPER_HH

// ## Upscaling helper struct (`upscalinghelper.hh`)
//
// This file contains the __upscaling helper struct__ which considers the volume flux leaving
// the pore network in flow direction in order to find the upscaled Darcy permeability.
// [[content]]
#include <iostream>
#include <ostream>
#include <iomanip>
#include <numeric>
#include <functional>

namespace Dumux {

struct UpscalingHelper
{
    // ### Calculate the intrinsic permeability
    // This function first evaluates the mass flux leaving the network in the direction of the applied pressure gradient.
    // Afterwards, the mass flux is converted into an area specify volume flux from which finally the intrinsic Darcy
    // permeability $`\mathbf{K}`$ [m$`^2`$] can be evaluated.
    // [[codeblock]]
    template<class Problem, class Scalar>
    static Scalar getDarcyPermeability(const Problem& problem, const Scalar totalMassFlux)
    {
        // get the domain side lengths from the problem
        auto sideLengths = problem.sideLengths();

        // create temporary stringstream with fixed scientifc formatting without affecting std::cout
        std::ostream tmp(std::cout.rdbuf());
        tmp << std::fixed << std::scientific;
        static constexpr char dirNames[] = "xyz";

        // convert mass to volume flux
        const auto volumeFlux = totalMassFlux / problem.liquidDensity();;

        sideLengths[problem.direction()] = 1.0;
        const auto outflowArea = std::accumulate(sideLengths.begin(), sideLengths.end(), 1.0, std::multiplies<Scalar>());
        const auto vDarcy = volumeFlux / outflowArea;
        const auto K = vDarcy / problem.pressureGradient() * problem.liquidDynamicViscosity();
        tmp << "\n########################################\n" << std::endl;
        tmp << dirNames[problem.direction()] << "-direction";
        tmp << ": Area = " << outflowArea << " m^2";
        tmp << "; Massflux = " << totalMassFlux << " kg/s";
        tmp << "; v_Darcy = " << vDarcy << " m/s";
        tmp << "; K = " << K << " m^2" << std::endl;
        tmp << "\n########################################\n" << std::endl;

        return K;
    }
    // [[/codeblock]]

    // ### Determine the domain's side lengths automatically based on the bounding box of the network.
    // [[codeblock]]
    template<class GridGeometry>
    static auto getSideLengths(const GridGeometry& gridGeometry)
    {
        using GlobalPosition = typename GridGeometry::GlobalCoordinate;
        GlobalPosition result;

        std::cout << "Automatically determining side lengths of REV based on bounding box of pore network" << std::endl;
        for (int dimIdx = 0; dimIdx < GridGeometry::GridView::dimensionworld; ++dimIdx)
            result[dimIdx] = gridGeometry.bBoxMax()[dimIdx] - gridGeometry.bBoxMin()[dimIdx];

        return result;
    }
    // [[/codeblock]]
};

} // end namespace Dumux
// [[/content]]
#endif
