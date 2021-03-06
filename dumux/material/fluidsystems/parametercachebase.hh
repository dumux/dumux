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
 * \ingroup Fluidsystems
 * \brief @copybrief Dumux::ParameterCacheBase
 */
#ifndef DUMUX_PARAMETER_CACHE_BASE_HH
#define DUMUX_PARAMETER_CACHE_BASE_HH

namespace Dumux {

/*!
 * \ingroup Fluidsystems
 * \brief The base class of the parameter cache classes for fluid systems
 */
template <class Implementation>
class ParameterCacheBase
{
public:
    enum ExceptQuantities {
        None = 0,
        Temperature = 1,
        Pressure = 2,
        Composition = 4
    };

    /*!
     * \brief Update all cached quantities for all phases.
     *
     * The <tt>except</tt> argument contains a bit field of the  quantities
     * which have not been modified since the last call to an <tt>update()</tt>
     * method.
     */
    template <class FluidState>
    void updateAll(const FluidState &fs, int exceptQuantities = None)
    {
        for (int phaseIdx = 0; phaseIdx < FluidState::numPhases; ++phaseIdx)
            asImp_().updatePhase(fs, phaseIdx);
    }

    /*!
     * \brief Update all cached quantities which depend on the pressure of any
     * fluid phase.
     */
    template <class FluidState>
    void updateAllPressures(const FluidState &fs)
    {
        for (int phaseIdx = 0; phaseIdx < FluidState::numPhases; ++phaseIdx)
            asImp_().updatePhase(fs, phaseIdx);
    }

    /*!
     * \brief Update all cached quantities which depend on the temperature of any
     * fluid phase.
     */
    template <class FluidState>
    void updateAllTemperatures(const FluidState &fs)
    {
        for (int phaseIdx = 0; phaseIdx < FluidState::numPhases; ++phaseIdx)
            asImp_().updatePhase(fs, phaseIdx);
    }

    /*!
     * \brief Update all cached parameters of a specific fluid phase
     *
     * The quantities specified by the <tt>except</tt> bit field have not been
     * modified since since the last call to an <tt>update()</tt> method.
     */
    template <class FluidState>
    void updatePhase(const FluidState &fs, int phaseIdx, int exceptQuantities = None)
    {}

    /*!
     * \brief Update all cached parameters of a specific fluid phase
     *        which depend on temperature
     *
     * \b Only use this method if only the temperature of a phase
     * changed between two update*() calls. If more changed, call
     * updatePhase()!
     */
    template <class FluidState>
    void updateTemperature(const FluidState &fs, int phaseIdx)
    {
        asImp_().updatePhase(fs, phaseIdx);
    }

    /*!
     * \brief Update all cached parameters of a specific fluid phase
     *        which depend on pressure
     *
     * \b Only use this method if only the pressure of a phase changed
     * between two update*() calls. If more changed, call
     * updatePhase()!
     */
    template <class FluidState>
    void updatePressure(const FluidState &fs, int phaseIdx)
    {
        asImp_().updatePhase(fs, phaseIdx);
    }

    /*!
     * \brief Update all cached parameters of a specific fluid phase
     *        which depend on composition
     *
     * \b Only use this method if neither the pressure nor the
     * temperature of the phase changed between two update*()
     * calls. If more changed, call updatePhase()!
     */
    template <class FluidState>
    void updateComposition(const FluidState &fs, int phaseIdx)
    {
        asImp_().updatePhase(fs, phaseIdx, /*except=*/Temperature | Pressure);
    }

    /*!
     * \brief Update all cached parameters of a specific fluid phase
     *        which depend on the mole fraction of a single component
     *
     * \b Only use this method if just a single component's
     * concentration changed between two update*() calls. If more than
     * one concentration changed, call updatePhaseComposition() of
     * updatePhase()!
     */
    template <class FluidState>
    void updateSingleMoleFraction(const FluidState &fs,
                                  int phaseIdx,
                                  int compIdx)
    {
        asImp_().updateComposition(fs, phaseIdx);
    }

private:
    Implementation &asImp_()
    { return *static_cast<Implementation*>(this); }
};

} // end namespace

#endif
