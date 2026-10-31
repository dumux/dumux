// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \brief Test that the assembler of the general interface picks the local assembler
 *        belonging to the discretization method and the kind of grid variables
 */
#include <config.h>

#include <iostream>
#include <type_traits>

#include <dumux/assembly/assembler.hh>
#include <dumux/assembly/fvassembler.hh>
#include <dumux/common/concepts/variables_.hh>
#include <dumux/discretization/method.hh>

namespace Dumux {

struct MockFVGridVariables
{
    struct GridVolumeVariables {};
    struct VolumeVariables {};
    struct GridFluxVariablesCache {};
    int curGridVolVars() const { return 0; }
    int prevGridVolVars() const { return 0; }
};

struct MockGridVariables
{
    struct GridVariablesCache {};
    struct Variables {};
    int curGridVars() const { return 0; }
    int prevGridVars() const { return 0; }
};

template<class DM, bool usesVariables>
using Chosen = typename Experimental::Detail::LocalAssemblerChooser<DM, usesVariables>
    ::template type<void, void, DiffMethod::numeric, true>;

} // end namespace Dumux

int main()
{
    using namespace Dumux;
    namespace DMs = DiscretizationMethods;

    static_assert(Concept::FVGridVariables<MockFVGridVariables> && !Concept::GridVariables<MockFVGridVariables>);
    static_assert(Concept::GridVariables<MockGridVariables> && !Concept::FVGridVariables<MockGridVariables>);

    static_assert(std::is_same_v<Chosen<DMs::CCTpfa, false>, CCLocalAssembler<void, void, DiffMethod::numeric, true>>);
    static_assert(std::is_same_v<Chosen<DMs::CCMpfa, false>, CCLocalAssembler<void, void, DiffMethod::numeric, true>>);
    static_assert(std::is_same_v<Chosen<DMs::FCStaggered, false>, FaceCenteredLocalAssembler<void, void, DiffMethod::numeric, true>>);
    static_assert(std::is_same_v<Chosen<DMs::Box, false>, CVFELocalAssembler<void, void, DiffMethod::numeric, true>>);
    static_assert(std::is_same_v<Chosen<DMs::Box, true>, Experimental::CVFELocalAssembler<void, void, DiffMethod::numeric, true>>);

    std::cout << "All local assembler choices are as expected" << std::endl;
    return 0;
}
