// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup Concepts
 * \brief Concepts related to local dofs
 */
#ifndef DUMUX_CONCEPTS_LOCALDOFS__HH
#define DUMUX_CONCEPTS_LOCALDOFS__HH

#include <type_traits>
#include <concepts>

#include <dumux/discretization/cvfe/localdof.hh>

namespace Dumux::Concept {

template<class L>
concept LocalDof = std::is_base_of_v<Dumux::CVFE::LocalDof<typename L::LocalIndexType, typename L::GridIndexType>, L>;

} // end namespace Dumux::Concept

#endif
