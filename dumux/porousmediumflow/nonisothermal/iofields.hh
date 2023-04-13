// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
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
