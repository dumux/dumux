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
/**
 * \file
 * \ingroup TracerTests
 * \brief The test fluidsystem
 */

#ifndef DUMUX_TRACER_MULTIPHASE_TEST_FLUIDSYSTEM_HH
#define DUMUX_TRACER_MULTIPHASE_TEST_FLUIDSYSTEM_HH

#include <dumux/common/parameters.hh>
#include <dumux/material/fluidsystems/base.hh>

namespace Dumux::FluidSystems {

//! A simple fluid system with two independent tracer components
template<class Scalar>
class TracerTest : public Base<Scalar, TracerTest<Scalar>>
{
public:
    static constexpr bool isTracerFluidSystem()
    { return true; }

    //! The number of components
    static constexpr int numComponents = 2;

    //! Human readable component name (index compIdx) (for vtk output)
    static std::string componentName(int compIdx)
    { return "tracer_" + std::to_string(compIdx); }

    //! Human readable phase name (index phaseIdx) (for velocity vtk output)
    static std::string phaseName(int phaseIdx = 0)
    { return "aq"; }

    //! Molar mass in kg/mol of the component with index compIdx
    static Scalar molarMass(unsigned int compIdx)
    { return 0.300; }

    //! Binary diffusion coefficient
    //! (might depend on spatial parameters like pressure / temperature)
    template<class Problem, class Element, class SubControlVolume>
    static Scalar binaryDiffusionCoefficient(unsigned int compIdx,
                                             const Problem& problem,
                                             const Element& element,
                                             const SubControlVolume& scv)
    {
        static const Scalar D = getParam<Scalar>("Problem.D");
        static const Scalar D2 = getParam<Scalar>("Problem.D2");
        if (compIdx == 0)
            return D;
        else
            return D2;
    }

    /*!
     * \copydoc Dumux::FluidSystems::Base::isCompressible
     */
    static constexpr bool isCompressible(int phaseIdx)
    { return false; }

     /*!
     * \copydoc  Dumux::FluidSystems::Base::viscosityIsConstant
     */
    static constexpr bool viscosityIsConstant(int phaseIdx)
    { return true; }
};

} // end namespace Dumux::FluidSystems

#endif
