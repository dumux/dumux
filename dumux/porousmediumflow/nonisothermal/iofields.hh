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
 * \ingroup NIModel
 * \brief Adds I/O fields specific to non-isothermal models.
 */
#ifndef DUMUX_ENERGY_IO_FIELDS_HH
#define DUMUX_ENERGY_IO_FIELDS_HH

#include <string>
#include <dumux/io/name.hh>

namespace Dumux {

/*!
 * \ingroup NIModel
 * \brief Adds I/O fields specific to non-isothermal models
 * \tparam IsothermalIOFields the isothermal io fields to adapt to non-isothermal io fields
 */
template<class IsothermalIOFields = void>
class EnergyIOFields
{
public:
    template <class OutputModule>
    static void initOutputModule(OutputModule& out)
    {
        IsothermalIOFields::initOutputModule(out);
        out.addVolumeVariable( [](const auto& v){ return v.temperature(); }, IOName::temperature());
    }

    template <class ModelTraits, class FluidSystem = void, class SolidSystem = void>
    static std::string primaryVariableName(int pvIdx, int state = 0)
    {
        using IsothermalTraits = typename ModelTraits::IsothermalTraits;

        if (pvIdx < ModelTraits::numEq() - 1)
            return IsothermalIOFields::template primaryVariableName<IsothermalTraits, FluidSystem, SolidSystem>(pvIdx, state);
        else
            return IOName::temperature();
    }
};

/*!
 * \ingroup NIModel
 * \brief Adds I/O fields specific to non-isothermal models.
 * \note Specialization if this is class used on its own (not as an adapter)
 */
template<>
class EnergyIOFields<void>
{
public:
    template <class OutputModule>
    static void initOutputModule(OutputModule& out)
    {
        out.addVolumeVariable( [](const auto& v){ return v.temperature(); }, IOName::temperature());
    }

    template <class ModelTraits, class FluidSystem = void, class SolidSystem = void>
    static std::string primaryVariableName(int pvIdx, int state = 0)
    {
        return IOName::temperature();
    }
};

} // end namespace Dumux

#endif
