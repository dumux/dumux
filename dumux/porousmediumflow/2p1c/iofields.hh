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
 * \ingroup TwoPOneCModel
 * \copydoc Dumux::TwoPOneCIOFields
 */

#ifndef DUMUX_TWOP_ONEC_IO_FIELDS_HH
#define DUMUX_TWOP_ONEC_IO_FIELDS_HH

#include <dumux/porousmediumflow/2p/iofields.hh>
#include <dumux/io/name.hh>

namespace Dumux {

/*!
 * \ingroup TwoPOneCModel
 * \brief Adds I/O fields specific to two-phase one-component model.
 */
class TwoPOneCIOFields
{
public:
    template <class OutputModule>
    static void initOutputModule(OutputModule& out)
    {
        // use default fields from the 2p model
        TwoPIOFields::initOutputModule(out);

        // output additional to TwoP output:
        out.addVolumeVariable([](const auto& v){ return v.priVars().state(); },
                              IOName::phasePresence());
    }

    template <class ModelTraits, class FluidSystem, class SolidSystem = void>
    static std::string primaryVariableName(int pvIdx, int state)
    {
        using Indices = typename ModelTraits::Indices;

        if (ModelTraits::priVarFormulation() == TwoPFormulation::p0s1)
            return (pvIdx == 0) ? IOName::pressure<FluidSystem>(FluidSystem::phase0Idx) :
                                  (state == Indices::twoPhases)
                                  ? IOName::saturation<FluidSystem>(FluidSystem::phase1Idx)
                                  : IOName::temperature();
        else
            return (pvIdx == 0) ? IOName::pressure<FluidSystem>(FluidSystem::phase1Idx) :
                                  (state == Indices::twoPhases)
                                  ? IOName::saturation<FluidSystem>(FluidSystem::phase0Idx)
                                  : IOName::temperature();
    }
};

} // end namespace Dumux

#endif
