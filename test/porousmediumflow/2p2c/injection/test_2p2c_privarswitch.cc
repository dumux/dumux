// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup TwoPTwoCTests
 * \brief Detection of the primary variable switch for the volume variables and for the variables per local dof
 *
 * With `USEVARIABLES`, the box scheme assembles with the variables per local dof, for which the
 * primary variable switch is not implemented. Instantiating the switch adapter of the Newton
 * solver (`INSTANTIATE_ADAPTER`), which evaluates `hasPriVarsSwitch`, then fails to compile.
 */
#include <config.h>

#include <dune/common/std/type_traits.hh>

#include <dumux/common/properties.hh>
#include <dumux/nonlinear/primaryvariableswitchadapter.hh>

#include "properties.hh"

namespace Dumux::Properties {

template<class TypeTag>
struct Transition::UseVariables<TypeTag, TTag::InjectionBox> { static constexpr bool value = USEVARIABLES; };

} // end namespace Dumux::Properties

int main()
{
    using namespace Dumux;
    using GridVariables = GetPropType<Properties::TTag::InjectionBox, Properties::GridVariables>;

#if USEVARIABLES
    static_assert(Dune::Std::is_detected_v<Detail::DetectLocalDofVariablesPVSwitch, GridVariables>);
#else
    static_assert(hasPriVarsSwitch<GridVariables>);
    static_assert(!Dune::Std::is_detected_v<Detail::DetectLocalDofVariablesPVSwitch, GridVariables>);
#endif

#if INSTANTIATE_ADAPTER
    [[maybe_unused]] constexpr auto adapterSize = sizeof(PrimaryVariableSwitchAdapter<GridVariables>);
#endif

    return 0;
}
