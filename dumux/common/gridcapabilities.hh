// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   See the file COPYING for full copying permissions.                      *
 *                                                                           *
 *   This program is free software: you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation, either version 3 of the License, or       *
 *   (at your option) any later version.                                     *
 *                                                                           *
 *   This program is distributed in the hope that it will be useful,         *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of          *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the            *
 *   GNU General Public License for more details.                            *
 *                                                                           *
 *   You should have received a copy of the GNU General Public License       *
 *   along with this program.  If not, see <http://www.gnu.org/licenses/>.   *
 *****************************************************************************/
/*!
 * \file
 * \ingroup Common
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

#endif
