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
 * \ingroup LowReKEpsilonModel
 * \copydoc Dumux::LowReKEpsilonIOFields
 */
#ifndef DUMUX_LOWREKEPSILON_IO_FIELDS_HH
#define DUMUX_LOWREKEPSILON_IO_FIELDS_HH

#include <dumux/freeflow/rans/iofields.hh>

namespace Dumux {

/*!
 * \ingroup LowReKEpsilonModel
 * \brief Adds I/O fields for the low-Re k-epsilon turbulence model
 */
struct LowReKEpsilonIOFields
{
    //! Initialize the LowReKEpsilon specific output fields.
    template <class OutputModule>
    static void initOutputModule(OutputModule& out)
    {
        RANSIOFields::initOutputModule(out);
        out.addVolumeVariable([](const auto& v){ return v.turbulentKineticEnergy(); }, "k");
        out.addVolumeVariable([](const auto& v){ return v.dissipationTilde(); }, "epsilon");
    }

    //! return the names of the primary variables
    template <class ModelTraits, class FluidSystem>
    static std::string primaryVariableName(int pvIdx = 0, int state = 0)
    {
        if (pvIdx < ModelTraits::dim() + ModelTraits::numFluidComponents())
            return RANSIOFields::template primaryVariableName<ModelTraits, FluidSystem>(pvIdx, state);
        else if (pvIdx == ModelTraits::dim() + ModelTraits::numFluidComponents())
            return "k";
        else
            return "epsilon";
    }
};

} // end namespace Dumux

#endif
