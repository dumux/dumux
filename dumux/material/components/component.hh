// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   See the file COPYING for full copying permissions.                      *
 *                                                                           *
 *   This program is free software: you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation, either version 2 of the License, or       *
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
 * \ingroup Components
 * \brief Abstract base class of a pure chemical species.
 */
#ifndef DUMUX_COMPONENT_HH
#define DUMUX_COMPONENT_HH

#warning "This header is deprecated. Use base.hh/solid.hh/liquid.hh/gas.hh"
#include <dune/common/stdstreams.hh>
#include <dumux/material/components/base.hh>
#include <dumux/material/components/liquid.hh>
#include <dumux/material/components/gas.hh>
#include <dumux/material/components/solid.hh>

namespace Dumux {

/*!
 * \ingroup Components
 * \brief Abstract base class of a pure chemical species.
 *
 * \tparam Scalar The type used for scalar values
 * \tparam Implementation Necessary for static polymorphism
 */
template <class Scalar, class C>
class DUNE_DEPRECATED_MSG("Derive from Base and Liquid and/or Gas and/or Solid directly") Component
: public Components::Base<Scalar, C>
, public Components::Liquid<Scalar, C>
, public Components::Gas<Scalar, C>
, public Components::Solid<Scalar, C>
{ };

} // end namespace Dumux

#endif
