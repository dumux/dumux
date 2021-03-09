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
 * \brief Coupling manager for low-dimensional domains embedded in the bulk domain.
 */

#ifndef DUMUX_MULTIDOMAIN_EMBEDDED_COUPLINGMANAGER_1D3D_HH
#define DUMUX_MULTIDOMAIN_EMBEDDED_COUPLINGMANAGER_1D3D_HH

namespace Dumux {

/*!
 * \ingroup EmbeddedCoupling
 * \brief Manages the coupling between bulk elements and lower dimensional elements
 *        Point sources on each integration point are computed by an AABB tree.
 */
template<class MDTraits, class CouplingMode>
class Embedded1d3dCouplingManager;

/*!
 * \ingroup EmbeddedCoupling
 * \brief The coupling mode
 * \deprecated Will be removed after 3.3
 */
enum class EmbeddedCouplingMode
{ line, average, cylindersources, kernel };

// deprecated and will be removed after 3.3
template<class MDTraits, EmbeddedCouplingMode mode>
struct EmbeddedCouplingManager1d3d;

} // end namespace Dumux

#include <dumux/multidomain/embedded/couplingmanager1d3d_line.hh>
#include <dumux/multidomain/embedded/couplingmanager1d3d_average.hh>
#include <dumux/multidomain/embedded/couplingmanager1d3d_surface.hh>
#include <dumux/multidomain/embedded/couplingmanager1d3d_kernel.hh>

// TODO: remove this part after the release 3.3
namespace Dumux {

template<class MDTraits>
struct [[deprecated("Removed after 3.3. Use Embedded1d3dCouplingManager.")]] EmbeddedCouplingManager1d3d<MDTraits, EmbeddedCouplingMode::line>
: public Embedded1d3dCouplingManager<MDTraits, Embedded1d3dCouplingMode::Line> {
    using ParentType = Embedded1d3dCouplingManager<MDTraits, Embedded1d3dCouplingMode::Line>;
    using ParentType::ParentType;
};

template<class MDTraits>
struct [[deprecated("Removed after 3.3. Use Embedded1d3dCouplingManager.")]] EmbeddedCouplingManager1d3d<MDTraits, EmbeddedCouplingMode::average>
: public Embedded1d3dCouplingManager<MDTraits, Embedded1d3dCouplingMode::Average> {
    using ParentType = Embedded1d3dCouplingManager<MDTraits, Embedded1d3dCouplingMode::Average>;
    using ParentType::ParentType;
};

template<class MDTraits>
struct [[deprecated("Removed after 3.3. Use Embedded1d3dCouplingManager.")]] EmbeddedCouplingManager1d3d<MDTraits, EmbeddedCouplingMode::cylindersources>
: public Embedded1d3dCouplingManager<MDTraits, Embedded1d3dCouplingMode::Surface> {
    using ParentType = Embedded1d3dCouplingManager<MDTraits, Embedded1d3dCouplingMode::Surface>;
    using ParentType::ParentType;
};

template<class MDTraits>
struct [[deprecated("Removed after 3.3. Use Embedded1d3dCouplingManager.")]] EmbeddedCouplingManager1d3d<MDTraits, EmbeddedCouplingMode::kernel>
: public Embedded1d3dCouplingManager<MDTraits, Embedded1d3dCouplingMode::Kernel> {
    using ParentType = Embedded1d3dCouplingManager<MDTraits, Embedded1d3dCouplingMode::Kernel>;
    using ParentType::ParentType;
};

} // end namespace Dumux

#endif
