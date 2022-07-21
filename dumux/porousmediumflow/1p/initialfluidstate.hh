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
 * \ingroup OnePModel
 * \brief Quantities required by the one-phase fully implicit model defined on a vertex.
 */

#ifndef DUMUX_1P_INITIAL_FLUIDSTATE_HH
#define DUMUX_1P_INITIAL_FLUIDSTATE_HH

namespace Dumux{
template<class Scalar, class FluidSystem, class FluidState>
class InitialImmiscibleFluidState
{
    static_assert(std::is_same_v<FluidState, ImmiscibleFluidState<Scalar,FluidSystem>>);
    static constexpr int numPhases = FluidSystem::numPhases;

    using Array = typename std::array<Scalar, numPhases>;
public:
    /*!
     * \brief complete the fluid state from an initial state.
     *
     * \param initFS
     * \return FluidState
     */

    InitialImmiscibleFluidState(const Scalar& p)
    {
        saturation_[0] = 1.0;
        pressure_[0] = p;
    }

    void setTemperature(const std::size_t& phaseIdx, const Scalar& value)
    {
        temperature_[phaseIdx] = value;
    }

    template<class EnergyVolumeVariables>
    void completeFluidState(FluidState& fs)
    {
        for (int phaseIdx = 0; phaseIdx < numPhases; phaseIdx++)
        {
            fs.setPressure(phaseIdx,pressure_[phaseIdx]);
            fs.setSaturation(phaseIdx, saturation_[phaseIdx]);
            fs.setTemperature(phaseIdx, temperature_[phaseIdx]);
        }

        typename FluidSystem::ParameterCache paramCache;
        paramCache.updatePhase(fs, /*phaseIdx=*/0);

        Scalar value = FluidSystem::density(fs, paramCache, /*phaseIdx=*/0);
        fs.setDensity(/*phaseIdx=*/0, value);

        value = FluidSystem::viscosity(fs, paramCache, /*phaseIdx=*/0);
        fs.setViscosity(/*phaseIdx=*/0, value);

        // compute and set the enthalpy
        value = EnergyVolumeVariables::enthalpy(fs, paramCache, /*phaseIdx=*/0);
        fs.setEnthalpy(/*phaseIdx=*/0, value);
    }

private:
    Array pressure_ = {};
    Array saturation_ = {};
    Array temperature_ = {};
};
} //namespace dumux

#endif
