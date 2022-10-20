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
 * \ingroup CCTpfaDiscretization
 * \copydoc Dumux::CCTpfaGridGeometry
 */
#ifndef DUMUX_DISCRETIZATION_CCTPFA_GRID_GEOMETRY_HH
#define DUMUX_DISCRETIZATION_CCTPFA_GRID_GEOMETRY_HH

#include <dune/common/exceptions.hh>
#include <dumux/common/defaultmappertraits.hh>

#include <dumux/discretization/checkoverlapsize.hh>
#include <dumux/discretization/method.hh>

namespace Dumux {

#ifndef DOXYGEN
namespace CCTpfa::Detail {

template<typename GridView>
void checkOverlapSize(const GridView& gridView)
{
    if (!CheckOverlapSize<DiscretizationMethods::CCTpfa>::isValid(gridView))
        DUNE_THROW(Dune::InvalidStateException,
                    "The cctpfa discretization method needs at least an overlap of 1 for parallel computations. " <<
                    "Set the parameter \"Grid.Overlap\" in the input file.");
}

} // CCTpfa::Detail
#endif // DOXYGEN

/*!
 * \ingroup CCTpfaDiscretization
 * \brief Default traits for `CCTpfaGridGeometry`
 */
template<typename GV>
struct DefaultCCTpfaGridGeometryTraits
: public DefaultMapperTraits<GV>
{
    static constexpr int numFacetsCube = 2*int(GV::dimension);
public:
    //! Per default, we use cubes with one hanging node per facet as maximum
    static constexpr int maxNumScvfsPerElement = 2*numFacetsCube;
    //! Per default, we set an upper limit of 7 branches (8 neighbors) on surface/network grids
    static constexpr int maxNumBranchesPerScvf = int(GV::dimension) < int(GV::dimensionworld) ? 7 : 1;
};

/*!
 * \ingroup CCTpfaDiscretization
 * \brief The finite volume grid geometry for the cell-centered TPFA scheme on a grid view.
 *        Provides an element-local view on the scvs and scvfs of the discretization.
 */
template<typename GV,
         bool cacheGridGeometries = false,
         typename Traits = DefaultCCTpfaGridGeometryTraits<GV>>
class CCTpfaGridGeometry;

// Deduction guide to allow for template argument deduction
template<typename GV>
CCTpfaGridGeometry(const GV& gv) -> CCTpfaGridGeometry<GV>;

} // end namespace Dumux

#include "gridgeometry_caching.hh"
#include "gridgeometry_nocaching.hh"

#endif
