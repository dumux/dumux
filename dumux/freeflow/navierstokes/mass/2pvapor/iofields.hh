// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
#ifndef DUMUX_NAVIERSTOKES_MASS_2PVAPOR_IO_FIELDS_HH
#define DUMUX_NAVIERSTOKES_MASS_2PVAPOR_IO_FIELDS_HH

#include <dumux/freeflow/navierstokes/mass/2p/iofields.hh>

namespace Dumux {

class NavierStokesMassTwoPVaporIOFields
{
public:
    template <class OutputModule>
    static void initOutputModule(OutputModule& out)
    {
        NavierStokesMassTwoPIOFields::initOutputModule(out);
        out.addVolumeVariable([](const auto& v){ return v.vaporConcentration(); }, "c_v");
    }

    template <class ModelTraits, class FluidSystem = void>
    static std::string primaryVariableName(int pvIdx = 0, int state = 0)
    {
        static constexpr std::array<std::string_view, 4> names = { "p", "phi", "mu", "c_v" };
        return std::string(names[pvIdx]);
    }
};

} // end namespace Dumux

#endif
