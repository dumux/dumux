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
 * \ingroup OneEqModel
 * \copydoc Dumux::OneEqIOFields
 */
#ifndef DUMUX_ONEEQ_IO_FIELDS_HH
#define DUMUX_ONEEQ_IO_FIELDS_HH

#include <dumux/freeflow/rans/iofields.hh>

namespace Dumux {

/*!
 * \ingroup OneEqModel
 * \brief Adds I/O fields for the one-equation turbulence model by Spalart-Allmaras
 */
struct OneEqIOFields
{
    //! Initialize the OneEq specific output fields.
    template <class OutputModule>
    static void initOutputModule(OutputModule& out)
    {
        RANSIOFields::initOutputModule(out);
        out.addVolumeVariable([](const auto& v){ return v.viscosityTilde(); }, "nu_tilde");
    }

    //! return the names of the primary variables
    template <class ModelTraits, class FluidSystem>
    static std::string primaryVariableName(int pvIdx = 0, int state = 0)
    {
        if (pvIdx < ModelTraits::dim() + 1)
            return RANSIOFields::template primaryVariableName<ModelTraits, FluidSystem>(pvIdx, state);
        else
            return "nu_tilde";
    }
};

} // end namespace Dumux

#endif
