// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup MaterialTests
 * \brief Test the component traits.
 */

#include "config.h"

#include <type_traits>

#include <dumux/material/components/air.hh>
#include <dumux/material/components/componenttraits.hh>

int main(int argc, char *argv[])
{
    using namespace Dumux;

    using Traits = ComponentTraits<Components::Air<double>>;
    static_assert(Traits::hasGasState, "Air component is reported to have no gas state?!");
    static_assert(!Traits::hasSolidState, "Air component is reported to implement a solid state?!");
    static_assert(!Traits::hasLiquidState, "Air component is reported to implement a liquid state?!");
    static_assert(!Traits::isIon, "Air component is reported to be an ion?!");
    static_assert(std::is_same<double, Traits::Scalar>::value, "Scalar type not correctly reported!");
}
