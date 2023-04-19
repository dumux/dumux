// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup MaterialTests
 * \brief This file provides the actual code for the solid systems test.
 *
 * It is not directly in test_solidsystems.cc so that external modules
 * like dumux-devel can use it easily
 */

#ifndef DUMUX_CHECK_SOLIDSYSTEM_HH
#define DUMUX_CHECK_SOLIDSYSTEM_HH

#include <type_traits>
#include <iostream>

#include <dune/common/classname.hh>
#include <dune/common/exceptions.hh>

// include all solid states
#include <dumux/material/solidstates/inertsolidstate.hh>
#include <dumux/material/solidstates/compositionalsolidstate.hh>

namespace Dumux {

template<class Scalar, class SolidSystem>
int checkSolidSystem()
{
    static_assert((SolidSystem::isInert() && SolidSystem::numInertComponents == SolidSystem::numComponents)
                  || !SolidSystem::isInert(), "numInertComponents == numComponents has to be equal to SolidSystem::isInert()");

    // inert solid state
    using Solidstate = std::conditional_t<SolidSystem::isInert(),
                                          InertSolidState<Scalar, SolidSystem>,
                                          CompositionalSolidState<Scalar, SolidSystem>>;
    Solidstate sst;

    int success = 0;
    std::cout << "Testing solid system '" << Dune::className<SolidSystem>() << "'\n";

    // output strings
    std::string collectedErrors;
    std::string collectedWarnings;

    // make sure the solid system provides the number of phases and
    // the number of components
    enum
    {
        numComponents = SolidSystem::numComponents
    };

    // some value to make sure the return values of the solid system
    // are convertible to scalars
    [[maybe_unused]] Scalar val;

    // test for componentName and isCompressible
    for (int phaseIdx = 0; phaseIdx < numComponents; ++phaseIdx)
    {
        [[maybe_unused]] std::string name = SolidSystem::componentName(phaseIdx);
        [[maybe_unused]] bool bVal = SolidSystem::isCompressible(phaseIdx);
        val = SolidSystem::molarMass(phaseIdx);
    }

    // test for componentName
    for (int compIdx = 0; compIdx < numComponents; ++compIdx)
    {
        [[maybe_unused]] std::string name = SolidSystem::componentName(compIdx);
    }

    // test for name
    [[maybe_unused]] std::string name = SolidSystem::name();

    try
    {
        val = SolidSystem::heatCapacity(sst);
    } catch (Dune::NotImplemented&)
    {
        collectedWarnings += "warning: SolidSystem::heatCapacity() is not implemented\n";
    } catch (...)
    {
        collectedErrors += "error: SolidSystem::heatCapacity() throws exception!\n";
    }

    try
    {
        val = SolidSystem::thermalConductivity(sst);
    } catch (Dune::NotImplemented&)
    {
        collectedWarnings += "warning: SolidSystem::thermalConductivity() is not implemented\n";
    } catch (...)
    {
        collectedErrors += "error: SolidSystem::thermalConductivity() throws exception!\n";
    }

    try
    {
        val = SolidSystem::density(sst);
    } catch (Dune::Exception &e)
    {
        collectedErrors += "error: SolidSystem::density() throws exception!\n";
    }

    std::cout << collectedErrors;
    std::cout << collectedWarnings;
    if (collectedErrors.empty()) // success
    {
        std::cout << "... successful" << std::endl;
        std::cout << "----------------------------------" << std::endl;
        return 0;
    }
    else
    {
        std::cout << "... failed" << std::endl;
        std::cout << "----------------------------------" << std::endl;
        return 1;
    }
    std::cout << "----------------------------------\n";
    return success;
}

} // end namespace Dumux

#endif
