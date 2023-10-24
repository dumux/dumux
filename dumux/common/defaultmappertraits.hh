// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup Core
 * \brief Defines the default element and vertex mapper types
 */
#ifndef DUMUX_DEFAULT_MAPPER_TRAITS_HH
#define DUMUX_DEFAULT_MAPPER_TRAITS_HH

#include <dune/grid/common/mcmgmapper.hh>

namespace Dumux {

template <class GridView,
          class EM = Dune::MultipleCodimMultipleGeomTypeMapper<GridView>,
          class VM = Dune::MultipleCodimMultipleGeomTypeMapper<GridView>>
struct DefaultMapperTraits
{
    using ElementMapper = EM;
    using VertexMapper = VM;
};

} // namespace Dumux

#endif
