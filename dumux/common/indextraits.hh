// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup Core
 * \brief Defines the index types used for grid and local indices.
 */
#ifndef DUMUX_COMMON_INDEX_TRAITS_HH
#define DUMUX_COMMON_INDEX_TRAITS_HH

#include <cstdint>

namespace Dumux {

/*!
 * \ingroup Core
 * \brief Structure to define the index types used for grid and local indices.
 * \tparam GridView The grid view type
 */
template<class GridView>
struct IndexTraits
{
    using GridIndex = typename GridView::IndexSet::IndexType;
    using LocalIndex = unsigned int;
    using SmallLocalIndex = std::uint_least8_t;
};

} // namespace Dumux

#endif
