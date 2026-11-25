// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \brief The assembler a type tag needs, chosen by the kind of grid variables it declares
 */

#ifndef DUMUX_PYTHON_ASSEMBLY_ASSEMBLER_HH
#define DUMUX_PYTHON_ASSEMBLY_ASSEMBLER_HH

#include <type_traits>

#include <dumux/common/properties.hh>
#include <dumux/common/concepts/variables_.hh>
#include <dumux/assembly/diffmethod.hh>
#include <dumux/assembly/fvassembler.hh>
#include <dumux/assembly/assembler.hh>
#include <dumux/python/assembly/fvassembler.hh>

namespace Dumux::Python {

/*!
 * \brief The finite-volume assembler for a model on finite-volume grid variables, and the general
 *        assembler for a model on the general grid variables, which the former cannot assemble
 */
template<class TypeTag, DiffMethod diffMethod, bool isImplicit = true>
using AssemblerFor = std::conditional_t<
    Dumux::Concept::FVGridVariables<GetPropType<TypeTag, Properties::GridVariables>>,
    Dumux::FVAssembler<TypeTag, diffMethod, isImplicit>,
    Dumux::Experimental::Assembler<TypeTag, diffMethod, isImplicit>
>;

} // namespace Dumux::Python

#endif
