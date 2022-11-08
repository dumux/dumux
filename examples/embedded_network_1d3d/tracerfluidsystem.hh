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

#ifndef DUMUX_TRACER_FLUID_SYSTEM_HH
#define DUMUX_TRACER_FLUID_SYSTEM_HH

#include <dumux/material/fluidsystems/base.hh>

namespace Dumux  {

// A simple fluid system with one tracer component
template<class Scalar>
class TracerFluidSystem
: public FluidSystems::Base<Scalar, TracerFluidSystem<Scalar>>
{
public:
    // If the fluid system only contains tracer components
    static constexpr bool isTracerFluidSystem() { return true; }

    // The number of components
    static constexpr int numComponents = 1;

    // Human readable phase name (liquid)
    static std::string phaseName(int phaseIdx = 0) { return "l"; }

    // Human readable component name
    static std::string componentName(int compIdx = 0)
    {
        static const std::string name = getParam<std::string>("Tracer.Name", "tracer");
        return name;
    }

    // Molar mass in kg/mol of the component with index compIdx
    static Scalar molarMass(int compIdx)
    {
        static const Scalar molarMass = getParam<Scalar>("Tracer.MolarMass", 0.018);
        return molarMass;
    }

    // binary diffusion coefficient
    template<class Problem, class Element, class SubControlVolume>
    static Scalar binaryDiffusionCoefficient(int compIdx,
                                             const Problem& problem,
                                             const Element& element,
                                             const SubControlVolume& scv)
    {
        static const Scalar D = getParam<Scalar>("Tracer.DiffusionCoefficient");
        return D;
    }
};

} // end namespace Dumux

#endif
