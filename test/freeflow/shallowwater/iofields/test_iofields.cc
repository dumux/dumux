// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup ShallowWaterTests
 * \brief Test for the primary variable names used to restart shallow water simulations
 */
#include <config.h>

#include <array>
#include <iostream>
#include <string>

#include <dune/common/exceptions.hh>

#include <dumux/freeflow/shallowwater/iofields.hh>

int main(int argc, char** argv)
{
    const std::array<std::string, 3> expectedNames{"waterDepth", "velocityX", "velocityY"};
    for (int pvIdx = 0; pvIdx < 3; ++pvIdx)
    {
        const auto name = Dumux::ShallowWaterIOFields::primaryVariableName<void>(pvIdx);
        if (name != expectedNames[pvIdx])
            DUNE_THROW(Dune::Exception, "Primary variable " << pvIdx << " is named " << name
                                        << " instead of " << expectedNames[pvIdx]);
    }

    std::cout << "All shallow water primary variable names are correct" << std::endl;
    return 0;
}
