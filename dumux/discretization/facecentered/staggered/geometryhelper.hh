// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup FaceCenteredStaggeredDiscretization
 * \brief Geometry helper for face-centered staggered scheme.
 */
#ifndef DUMUX_DISCRETIZATION_FACECENTERED_STAGGERED_GEOMETRY_HELPER_HH
#define DUMUX_DISCRETIZATION_FACECENTERED_STAGGERED_GEOMETRY_HELPER_HH

#include <array>
#include <cassert>
#include <dune/common/exceptions.hh>
#include <dumux/common/indextraits.hh>
#include <dumux/discretization/facecentered/staggered/gridsupportsconcavecorners.hh>

namespace Dumux {
/*!
 * \ingroup FaceCenteredStaggeredDiscretization
 * \brief Face centered staggered geometry helper
 */
template<class GridView>
class FaceCenteredStaggeredGeometryHelper
{
    using GridIndexType = typename IndexTraits<GridView>::GridIndex;
    using SmallLocalIndexType = typename IndexTraits<GridView>::SmallLocalIndex;
    using Element = typename GridView::template Codim<0>::Entity;
    using Facet = typename GridView::template Codim<1>::Entity;

public:
    static constexpr auto dim = GridView::Grid::dimension;
    static constexpr auto numElementFaces = dim * 2;
    static constexpr auto numLateralFacesPerScv = 2 * (dim - 1);

    FaceCenteredStaggeredGeometryHelper(const GridView& gridView) : gridView_(gridView) {}

    //! Returns the local index of the opposing face.
    static constexpr SmallLocalIndexType localOppositeIdx(const SmallLocalIndexType ownLocalFaceIndex)
    {
        return isOdd_(ownLocalFaceIndex) ? (ownLocalFaceIndex - 1) : (ownLocalFaceIndex + 1);
    }

    //! Return the local index of a lateral orthogonal scvf
    static constexpr int lateralOrthogonalScvfLocalIndex(const SmallLocalIndexType ownLocalScvfIndex)
    {
        if constexpr(GridView::Grid::dimension == 1)
        {
            assert(false && "There are no lateral scvfs in 1D");
            return -1;
        }

        if constexpr (GridView::Grid::dimension == 2)
        {
            switch (ownLocalScvfIndex)
            {
                case 1: return 7;
                case 7: return 1;
                case 2: return 10;
                case 10: return 2;
                case 4: return 8;
                case 8: return 4;
                case 5: return 11;
                case 11: return 5;
                default:
                {
                    assert(false && "No lateral orthogonal scvf found");
                    return -1;
                }
            }
        }
        else
        {
            switch (ownLocalScvfIndex)
            {
                case 1: return 11;
                case 11: return 1;
                case 2: return 16;
                case 16: return 2;
                case 3: return 21;
                case 21: return 3;
                case 4: return 26;
                case 26: return 4;
                case 6: return 12;
                case 12: return 6;
                case 7: return 17;
                case 17: return 7;
                case 8: return 22;
                case 22: return 8;
                case 9: return 27;
                case 27: return 9;
                case 13: return 23;
                case 23: return 13;
                case 14: return 28;
                case 28: return 14;
                case 18: return 24;
                case 24: return 18;
                case 19: return 29;
                case 29: return 19;

                default:
                {
                    assert(false && "No lateral orthogonal scvf found");
                    return -1;
                }
            }
        }
    }

     //! Returns the local indices of the faces lateral to the own one.
    static constexpr auto localLaterFaceIndices(const SmallLocalIndexType ownLocalFaceIndex)
    {
        constexpr auto table = []
        {
            using Table = std::array<std::array<SmallLocalIndexType, numLateralFacesPerScv>, numElementFaces>;
            if constexpr (dim == 1)
                return Table{};
            else if constexpr (dim == 2)
                return Table {{ {2,3}, {2,3}, {0,1}, {0,1} }};
            else
                return Table {{ {2,3,4,5}, {2,3,4,5}, {0,1,4,5}, {0,1,4,5}, {0,1,2,3}, {0,1,2,3} }};
        }();

        return table[ownLocalFaceIndex];
    }

    //! Returns an element's facet based on the local facet index.
    static Facet facet(const SmallLocalIndexType localFacetIdx, const Element& element)
    {
        return element.template subEntity <1> (localFacetIdx);
    }

    //! Returns an element's intersection based on the local facet index.
    auto intersection(const SmallLocalIndexType localFacetIdx, const Element& element) const
    {
        for (const auto& intersection : intersections(gridView(), element))
        {
            if (intersection.indexInInside() == localFacetIdx)
                return intersection;
        }
        DUNE_THROW(Dune::InvalidStateException, "localFacetIdx " << localFacetIdx << " out of range");
    }

     //! Returns true if the IP of an scvf lies on a concave corner
    template<class FVElementGeometry, class SubControlVolumeFace>
    static bool scvfIntegrationPointInConcaveCorner(const FVElementGeometry& fvGeometry, const SubControlVolumeFace& scvf)
    {
        assert (scvf.isLateral());

        using Grid = typename FVElementGeometry::GridGeometry::Grid;
        if constexpr (!GridSupportsConcaveCorners<Grid>::value)
            return false;
        else
        {
            if (scvf.boundary())
                return false;

            const auto& insideScv = fvGeometry.scv(scvf.insideScvIdx());
            const auto& outsideScv = fvGeometry.scv(scvf.outsideScvIdx());
            int onBoundaryCounter = 0;
            onBoundaryCounter += static_cast<int>(insideScv.boundary());
            onBoundaryCounter += static_cast<int>(outsideScv.boundary());
            return onBoundaryCounter == 1;
        }
    }

    template<class SubControlVolumeFace, class SubControlVolume>
    static SmallLocalIndexType localIndexOutsideScvfWithSameIntegrationPoint(const SubControlVolumeFace& scvf,
                                                                             const SubControlVolume& scv)
    {
        // In 2D there are 3 non-boundary faces per scv. In 3D, there are 5.
        // This number of scvfs per scv is used as an offset to find the indexes in the outside half-scv.
        const SmallLocalIndexType offset = (dim == 2) ? 3 : 5;
        // For half-scvs with odd indexes, the outside half-scv has scvf local indexes with + offset.
        // For half-scvs with even indexes, the outside half-scv has scvf local indexes have a - offset.
        return isOdd_(scv.indexInElement()) ? scvf.localIndex() - offset : scvf.localIndex() + offset;
    }

    const GridView& gridView() const
    { return gridView_; }

private:

    static constexpr bool isOdd_(int number)
    { return number % 2; }

    GridView gridView_;
};

} // end namespace Dumux

#endif
