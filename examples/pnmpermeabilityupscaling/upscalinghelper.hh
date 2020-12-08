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

// ## Upscaling helper class (`upscalinghelper.hh`)
//
#include <vector>
#include <dune/common/exceptions.hh>
// This file contains the __upscaling helper class__ which considers the volume flux leaving
// the pore network in flow direction in order to find the upscaled Darcy permeability-
// [[content]]
// We include the boundary flux helper class in order to conveniently determine the mass and volume flux
// leaving the network.
#include <dumux/porenetworkflow/common/boundaryflux.hh>

// ### A helper class to find the upscaled intrinsic Darcy permeability.
// [[codeblock]]
namespace Dumux {

template<class Assembler>
class UpscalingHelper
{
    using Scalar = typename Assembler::Scalar;
    using GridView = typename Assembler::GridGeometry::GridView;
    static constexpr auto dimWorld = GridView::dimensionworld;
public:

    UpscalingHelper(const Assembler& assembler)
    : assembler_(assembler)
    {
        // The dimensions of the domain must be known in order to calculate the permeability.
        // One can either specify the domain size (e.g., based on the size of the sample used for the CT-scan image) ....
        sideLength_ = getParam<std::vector<Scalar>>("Problem.SideLength", std::vector<Scalar>{});
        if (!sideLength_.empty())
        {
            if (sideLength_.size() != dimWorld)
                DUNE_THROW(Dune::IOError, "Problem.SideLength must have exactly " << dimWorld << " entries");
        }
        // ... or get the size automatically based on the bounding box the pore network.
        else
        {
            std::cout << "Automatically determining side lengths of REV based on bounding box of pore network" << std::endl;
            sideLength_.resize(dimWorld);
            for (int dimIdx = 0; dimIdx < dimWorld; ++dimIdx)
                sideLength_[dimIdx] = assembler.gridGeometry().bBoxMax()[dimIdx] - assembler.gridGeometry().bBoxMin()[dimIdx];
        }
    }
    // [[/codeblock]]

    // #### Calculate the intrinsic permeability
    // This function first evaluates the mass flux leaving the network in the direction of the applied pressure gradient.
    // Afterwards, the mass flux is converted into an area specify volume flux from which finally the intrinsic Darcy
    // permeability K [m^2] can be evaluated.
    // [[codeblock]]
    template<class SolutionVector>
    void doUpscaling(const SolutionVector& x, const int direction) const
    {
        const auto boundaryFlux = PoreNetworkModelBoundaryFlux<Assembler>(assembler_, x);
        const auto outletPoreLabel = 2 + 2*direction;
        const auto massFlux = boundaryFlux.getFlux(std::vector<int>{outletPoreLabel});

        // create temporary stringstream with fixed scientifc formatting without affecting std::cout
        std::ostream tmp(std::cout.rdbuf());
        tmp << std::fixed << std::scientific;
        static constexpr char dirNames[] = "xyz";

        // convert mass to volume flux
        static const Scalar density = getParam<Scalar>("Component.LiquidDensity");
        static const Scalar dynamicViscosity = getParam<Scalar>("Component.LiquidKinematicViscosity") * density;
        static const auto pressureGradient = getParam<Scalar>("Problem.PressureGradient");
        const auto volumeFlux = massFlux / density;

        auto length = sideLength_;
        length[direction] = 1.0;
        const auto outflowArea = std::accumulate(length.begin(), length.end(), 1.0, std::multiplies<Scalar>());
        const auto vDarcy = volumeFlux / outflowArea;
        const auto K = vDarcy / pressureGradient * dynamicViscosity;
        tmp << "\n########################################\n" << std::endl;
        tmp << dirNames[direction] << "-direction";
        tmp << ": Area = " << outflowArea << " m^2";
        tmp << "; Massflux = " << massFlux << " kg/s";
        tmp << "; v_Darcy = " << vDarcy << " m/s";
        tmp << "; K = " << K << " m^2" << std::endl;
        tmp << "\n########################################\n" << std::endl;
    }
    // [[/codeblock]]

private:
    const Assembler& assembler_;
    std::vector<Scalar> sideLength_;
};

} // end namespace Dumux
// [[/content]]
#endif
