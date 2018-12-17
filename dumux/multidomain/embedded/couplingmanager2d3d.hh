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
 * \ingroup EmbeddedCoupling
 * \brief Coupling manager for embedded fractures
 */

#ifndef DUMUX_MULTIDOMAIN_EMBEDDED_COUPLINGMANAGER_2D3D_HH
#define DUMUX_MULTIDOMAIN_EMBEDDED_COUPLINGMANAGER_2D3D_HH

#include <dumux/multidomain/embedded/couplingmanagerbase.hh>

namespace Dumux {

/*!
 * \ingroup EmbeddedCoupling
 * \brief Coupling manager for embedded fractures
 * \note we just use the default coupling manager
 */
template<class MDTraits>
class EmbeddedCouplingManager2d3d
: public EmbeddedCouplingManagerBase<MDTraits,
                                     EmbeddedCouplingManager2d3d<MDTraits>>
{
    using ThisType = EmbeddedCouplingManager2d3d<MDTraits>;
    using ParentType = EmbeddedCouplingManagerBase<MDTraits, ThisType>;
public:
    using ParentType::ParentType;
};

} // end namespace Dumux

#endif
