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

#include <numeric>

#include <dumux/common/defaultmappertraits.hh>

#include <dumux/discretization/method.hh>
#include <dumux/discretization/extrusion.hh>
#include <dumux/discretization/checkoverlapsize.hh>
#include <dumux/discretization/basegridgeometry.hh>

#include <dumux/experimental/new_assembly/dumux/discretization/ccgridconnectivity.hh>

#include "detail.hh"

namespace Dumux {

/*!
 * \ingroup CCTpfaDiscretization
 * \brief Default traits for the `CCTpfaGridGeometry`
 */
template<typename GV>
struct DefaultCCTpfaGridGeometryTraits
: public DefaultMapperTraits<GV>
{
    //! Per default, we allow for maximally one hanging node per facet
    static constexpr int maxAdjacentElementLevelDifference = 1;
    //! Per default, we set an upper limit of 8 neighbors on surface/network grids
    static constexpr int maxNumBranchesPerScvf = int(GV::dimension) < int(GV::dimensionworld) ? 8 : 1;
};

/*!
 * \ingroup CCTpfaDiscretization
 * \brief Base class for finite volume grid geometries of the cell-centered TPFA scheme.
 */
template<typename GV, typename Traits = DefaultCCTpfaGridGeometryTraits<GV>>
class CCTpfaGridGeometryBase
: public BaseGridGeometry<GV, Traits>
{
    using ThisType = CCTpfaGridGeometryBase<GV, Traits>;
    using ParentType = BaseGridGeometry<GV, Traits>;
    using Element = typename GV::template Codim<0>::Entity;

    static constexpr bool isNetworkGrid = int(GV::dimension) < int(GV::dimensionworld);
    static constexpr auto maxNumElementFacets = CCTpfa::Detail::maxNumElementFacets<GV::dimension>;
    static constexpr auto maxNumFaceNeighbors = CCTpfa::Detail::maxNumFaceNeighbors<Traits>;
    static constexpr auto maxNumFacesPerFacet = CCTpfa::Detail::maxNumFacesPerFacet<Traits>;

    using Connectivity = CCGridConnectivity<GV, maxNumFaceNeighbors, maxNumFacesPerFacet, maxNumElementFacets>;

public:
    using GridView = GV;
    using Extrusion = Extrusion_t<Traits>;

    using DiscretizationMethod = DiscretizationMethods::CCTpfa;
    static constexpr DiscretizationMethod discMethod{};
    static constexpr auto maxElementStencilSize = CCTpfa::Detail::maxElementStencilSize<GV::dimension, Traits>;

    explicit CCTpfaGridGeometryBase(const GridView& gridView)
    : ParentType(gridView)
    , connectivity_(gridView, this->elementMapper())
    , numScvf_(
        std::accumulate(
            connectivity_.faceSeeds().begin(),
            connectivity_.faceSeeds().end(),
            std::size_t{0},
            [] (std::size_t current, const auto& face) {
                return current + face.numNeighbors();
        })
    )
    {
        if (!CheckOverlapSize<DiscretizationMethod>::isValid(gridView))
            DUNE_THROW(Dune::InvalidStateException,
                       "The cctpfa discretization method needs at least an overlap of 1 for parallel computations. " <<
                       "Set the parameter \"Grid.Overlap\" in the input file.");

    }

    //! Return the total number of degrees of freedom
    std::size_t numDofs() const
    { return this->gridView().size(0); }

    //! Return the total number of faces of the discretization
    std::size_t numFaces() const
    { return connectivity_.numFaces(); }

    //! Return the total number of sub control volumes
    std::size_t numScv() const
    { return this->gridView().size(0); }

    //! Return the total number of sub control volume faces
    std::size_t numScvf() const
    { return numScvf_; }

protected:
    using FaceSeed = typename Connectivity::FaceSeed;

    //! Return the face seeds associated with the given facet
    decltype(auto) faceSeeds_(std::size_t eIdx, unsigned int facetIndex) const
    {
        using GridIndex = typename Connectivity::Facet::GridIndex;
        using LocalIndex = typename Connectivity::Facet::LocalIndex;
        return connectivity_.faceSeeds({
            static_cast<GridIndex>(eIdx),
            static_cast<LocalIndex>(facetIndex)
        });
    }

private:
    Connectivity connectivity_;
    std::size_t numScvf_;
};

/*!
 * \ingroup CCTpfaDiscretization
 * \brief The finite volume grid geometry for the cell-centered TPFA scheme on a grid view.
 *        Provides an element-local view on the scvs and scvfs of the discretization.
 */
template<typename GV,
         bool cacheGridGeometries = true,
         typename Traits = DefaultCCTpfaGridGeometryTraits<GV>>
class CCTpfaGridGeometry;

// Deduction guide to allow for template argument deduction
template<typename GV>
CCTpfaGridGeometry(const GV& gv) -> CCTpfaGridGeometry<GV>;

} // end namespace Dumux

#include "gridgeometry_caching.hh"
#include "gridgeometry_nocaching.hh"

#endif
