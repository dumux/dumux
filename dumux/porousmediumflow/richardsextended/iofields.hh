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
 * \ingroup ExtendedRichardsModel
 * \brief Adds I/O fields specific to the extended Richards model.
 */

#ifndef DUMUX_RICHARDSEXTENDED_IO_FIELDS_HH
#define DUMUX_RICHARDSEXTENDED_IO_FIELDS_HH

#include <dumux/io/name.hh>
#include <dumux/porousmediumflow/richards/iofields.hh>

namespace Dumux {

/*!
 * \ingroup ExtendedRichardsModel
 * \brief Adds I/O fields specific to the extended Richards model.
 */
class ExtendedRichardsIOFields : public RichardsIOFields
{
    using ParentType = RichardsIOFields;
public:
    template <class OutputModule>
    static void initOutputModule(OutputModule& out)
    {
        using VV = typename OutputModule::VolumeVariables;
        using FS = typename VV::FluidSystem;

        ParentType::initOutputModule(out);

        out.addVolumeVariable([](const auto& v){ return v.moleFraction(FS::phase1Idx, FS::comp0Idx); },
                                  IOName::moleFraction<FS>(FS::phase1Idx, FS::comp0Idx));
        out.addVolumeVariable([](const auto& v){ return v.priVars().state(); },
                                  IOName::phasePresence());
    }
};

    template<class ModelTraits, class FluidSystem, class SolidSystem = void>
    static std::string primaryVariableName(int pvIdx, int state)
    {
        using Indices = typename ModelTraits::Indices;

        if (state == Indices::gasPhaseOnly)
            return IOName::moleFraction<FluidSystem>(FluidSystem::phase1Idx, FluidSystem::phase0Idx);
        else
            return IOName::pressure<FluidSystem>(FluidSystem::phase0Idx);
    }
} // end namespace Dumux

#endif
