// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup Grids
 * \brief Grid properties related to periodicity
 */

#ifndef DUMUX_IO_GRID_PERIODIC_GRID_TRAITS_HH
#define DUMUX_IO_GRID_PERIODIC_GRID_TRAITS_HH

#include <type_traits>
#include <dune/common/std/type_traits.hh>

namespace Dumux {

template<typename Grid>
struct PeriodicGridTraits
{
    struct SupportsPeriodicity : public std::false_type {};

    PeriodicGridTraits(const Grid& grid) {};

    bool isPeriodic (const typename Grid::LeafIntersection& intersection) const
    {
        return false;
    }
};

template<class T>
class SupportsPeriodicity
{
    template<class G>
    using SP = typename G::SupportsPeriodicity;
public:
    using type = typename Dune::Std::detected_or<std::false_type, SP, T>::type;
};

template<class T>
static constexpr bool supportsPeriodicity()
{ return typename SupportsPeriodicity<T>::type(); }

} // end namespace Dumux

#endif
