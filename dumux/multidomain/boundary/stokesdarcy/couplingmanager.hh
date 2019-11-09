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
 * \ingroup StokesDarcyCoupling
 * \copydoc Dumux::StokesDarcyCouplingManager
 */

#ifndef DUMUX_STOKES_DARCY_COUPLINGMANAGER_HH
#define DUMUX_STOKES_DARCY_COUPLINGMANAGER_HH

#include <dumux/common/properties.hh>

namespace Dumux {

template <class MDTraits, DiscretizationMethod darcyDM>
class StokesDarcyCouplingManagerImplementation;

template <class MDTraits>
class FreeFlowPorousMediumCouplingManagerBase
{
public:
    static constexpr auto freeFlowFaceIdx = typename MDTraits::template SubDomain<0>::Index();
    static constexpr auto freeFlowCellCenterIdx = typename MDTraits::template SubDomain<1>::Index();
    static constexpr auto freeFlowIdx = freeFlowFaceIdx;
    static constexpr auto porousMediumIdx = typename MDTraits::template SubDomain<2>::Index();

    [[deprecated("Use freeFlowFaceIdx instead!")]] static constexpr auto stokesFaceIdx = freeFlowFaceIdx;
    [[deprecated("Use freeFlowCellCenterIdx instead!")]] static constexpr auto stokesCellCenterIdx = freeFlowCellCenterIdx;
    [[deprecated("Use freeFlowIdx instead!")]] static constexpr auto stokesIdx = freeFlowIdx;
    [[deprecated("Use porousMediumIdx instead!")]] static constexpr auto darcyIdx = porousMediumIdx;
};

/*!
 * \ingroup StokesDarcyCoupling
 * \brief Coupling manager for Stokes and Darcy domains with equal dimension.
 */
template <class MDTraits>
using StokesDarcyCouplingManager = StokesDarcyCouplingManagerImplementation<MDTraits, GetPropType<typename MDTraits::template SubDomain<2>::TypeTag, Properties::GridGeometry>::discMethod>;


} // end namespace Dumux

#include  <dumux/multidomain/boundary/stokesdarcy/cellcentered/tpfa/couplingmanager.hh>
#include  <dumux/multidomain/boundary/stokesdarcy/box/couplingmanager.hh>

#endif
