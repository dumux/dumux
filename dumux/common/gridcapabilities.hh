// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup Core
 * \brief dune-grid capabilities compatibility layer
 */
#ifndef DUMUX_COMMON_GRID_CAPABILITIES_HH
#define DUMUX_COMMON_GRID_CAPABILITIES_HH

#include <dune/grid/common/capabilities.hh>

// TODO: The following is a temporary solution to make canCommunicate work.
// Once it is resolved upstream
// (https://gitlab.dune-project.org/core/dune-grid/issues/78),
// it should be guarded by a DUNE_VERSION macro and removed later.

#if HAVE_DUNE_UGGRID
namespace Dune {
template<int dim>
class UGGrid;
} // end namespace Dumux
#endif // HAVE_DUNE_UGGRID

namespace Dumux::Temp::Capabilities {

template<class Grid, int codim>
struct canCommunicate
{
  static const bool v = false;
};

#if HAVE_DUNE_UGGRID
template<int dim, int codim>
struct canCommunicate<Dune::UGGrid<dim>, codim>
{
  static const bool v = true;
};
#endif // HAVE_DUNE_UGGRID

} // namespace Dumux::Temp::Capabilities
// end workaround

namespace Dumux::Detail {

template<class Grid, int dofCodim>
static constexpr bool canCommunicate =
    Dune::Capabilities::canCommunicate<Grid, dofCodim>::v
    || Dumux::Temp::Capabilities::canCommunicate<Grid, dofCodim>::v;

} // namespace Dumux

namespace Dumux::Grid::Capabilities {

// Default implementation
// The grid capability gives an absolute guarantee
// however this does not mean that multithreading is not
// supported at all. It might still work for some special cases.
// This knowledge is encoded in specializations for the different
// grid managers, see dumux/grid/io/gridmanager_*.hh
template<class Grid>
struct MultithreadingSupported
{
    template<class GV>
    static bool eval(const GV&) // default is independent of the grid view
    { return Dune::Capabilities::viewThreadSafe<Grid>::v; }
};

template<class GridView>
inline bool supportsMultithreading(const GridView& gridView)
{ return MultithreadingSupported<typename GridView::Grid>::eval(gridView); }

} // namespace Dumux::Grid::Capabilities

#endif
