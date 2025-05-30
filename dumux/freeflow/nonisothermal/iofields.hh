// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup FreeflowNIModel
 * \copydoc Dumux::FreeflowNonIsothermalIOFields
 */
#ifndef DUMUX_FREEFLOW_NI_IO_FIELDS_HH
#define DUMUX_FREEFLOW_NI_IO_FIELDS_HH


#include <dumux/io/name.hh>

namespace Dumux {

/*!
 * \ingroup FreeflowNIModel
 * \brief Adds I/O fields specific to non-isothermal free-flow models
 */
template<class IsothermalIOFields, bool turbulenceModel = false>
struct FreeflowNonIsothermalIOFields
{

    //! Add the non-isothermal specific output fields.
    template <class OutputModule>
    static void initOutputModule(OutputModule& out)
    {
        IsothermalIOFields::initOutputModule(out);

        out.addVolumeVariable([](const auto& v){ return v.temperature(); }, IOName::temperature());
        out.addVolumeVariable([](const auto& v){ return v.thermalConductivity(); }, "lambda");
        if (turbulenceModel)
            out.addVolumeVariable([](const auto& v){ return v.effectiveThermalConductivity() - v.thermalConductivity(); }, "lambda_t");
    }

    //! return the names of the primary variables
    template<class ModelTraits, class FluidSystem = void>
    static std::string primaryVariableName(int pvIdx, int state = 0)
    {
        if (pvIdx < ModelTraits::numEq() - 1)
            return IsothermalIOFields::template primaryVariableName<ModelTraits, FluidSystem>(pvIdx, state);
        else
            return IOName::temperature();
    }
};

} // end namespace Dumux

#endif
