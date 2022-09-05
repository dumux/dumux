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
 * \brief TODO: Doc me
 */
#ifndef DUMUX_DISCRETIZATION_CCTPFA_GRID_GEOMETRY_LOCAL_VIEW_HH
#define DUMUX_DISCRETIZATION_CCTPFA_GRID_GEOMETRY_LOCAL_VIEW_HH

#include <cstdint>
#include <ranges>

#include <dumux/experimental/new_assembly/dumux/discretization/cctpfa/concepts.hh>
#include <dumux/experimental/new_assembly/dumux/discretization/cctpfa/gridgeometrycache.hh>

namespace Dumux::CCTpfa {

template<int codim, typename Id>
struct SubControlEntity
{ Id id; };

template<typename Id>
using SubControlVolume = SubControlEntity<0, Id>;

template<typename Id>
using SubControlVolumeFace = SubControlEntity<1, Id>;

template<typename GG, typename Traits, bool cacheGlobally>
class GridGeometryLocalView
{
    using GridGeometryCache = Detail::GridGeometryCache<GG, Traits, cacheGlobally>;
    using LocalCache = typename GridGeometryCache::LocalCache;

    using GridView = typename GG::GridView;
    using Element = typename GridView::template Codim<0>::Entity;

public:
    using GridGeometry = GG;
    using SubControlVolume = CCTpfa::SubControlVolume<std::size_t>;
    using SubControlVolumeFace = CCTpfa::SubControlVolumeFace<std::size_t>;

    using Index = typename GridView::IndexSet::IndexType;
    using Coordinate = typename Element::Geometry::GlobalCoordinate;
    using ctype = typename GridView::ctype;

    GridGeometryLocalView(const GridGeometry& gg,
                          const GridGeometryCache& cache)
    : gridGeometry_(gg)
    , geometryCache_(cache)
    {}

    std::size_t numScv() const { return std::ranges::size(localCache_.insideScvIndices()); }
    std::size_t numScvf() const { return std::ranges::size(localCache_.insideScvfIndices()); }

    //! prepare the geometries inside the given element
    void bindElement(const Element& element) &
    { geometryCache_.bindElement(localCache_, gridGeometry_, element); }

    //! prepare all geometries in the stencil of the given element
    void bind(const Element& element) &
    { geometryCache_.bindStencil(localCache_, gridGeometry_, element); }

    // overloads for rvalue references
    GridGeometryLocalView bind(const Element& e) && { bind(e); return *this; }
    GridGeometryLocalView bindElement(const Element& e) && { bindElement(e); return *this; }

    //! return the dof index associated with a sub-control volume
    Index dofIndex(const SubControlVolume& scv) const
    { return cache_(scv).dofIndex; }

    //! return the volume of a sub-control volume
    ctype volume(const SubControlVolume& scv) const
    { return cache_(scv).volume; }

    //! return the area of a sub-control volume face
    ctype area(const SubControlVolumeFace& scvf) const
    { return face_(scvf).area; }

    //! return the center point of a sub-control volume
    const Coordinate& center(const SubControlVolume& scv) const
    { return cache_(scv).center; }

    //! return the center point of a sub-control volume face
    const Coordinate& center(const SubControlVolumeFace& scvf) const
    { return face_(scvf).center; }

    //! return the outer normal vector of a sub-control volume face
    const Coordinate& unitOuterNormal(const SubControlVolumeFace& scvf) const
    { return cache_(scvf).normal; }

    //! return the geometry of a sub-control volume
    auto geometry(const SubControlVolume& scv) const
    { return gridGeometry_.element(dofIndex(scv)).geometry(); }

    //! return the geometry of a sub-control volume face
    auto geometry(const SubControlVolumeFace& scvf) const
    {
        const auto& facet = *face_(scvf).facet;
        const auto& element = gridGeometry_.element(facet.elementIndex);
        return element.template subEntity<1>(facet.facetIndex).geometry();
    }

    //! return the scv on the "inside" of the given sub-control volume face
    SubControlVolume insideScv(const SubControlVolumeFace& scvf) const
    {
        const auto nIdx = cache_(scvf).indexInFaceNeighbors;
        return SubControlVolume{face_(scvf).neighborScvCacheIndices[nIdx]};
    }

    //! return the i-th scv on the "outside" of the given sub-control volume face
    SubControlVolume outsideScv(const SubControlVolumeFace& scvf, unsigned int i = 0) const
    {
        assert(i < numOutsideScvs(scvf));
        if (i >= cache_(scvf).indexInFaceNeighbors)
            i++;
        return SubControlVolume{face_(scvf).neighborScvCacheIndices[i]};
    }

    //! return the number of "outside" scvs adjacent to the given sub-control volume face
    std::size_t numOutsideScvs(const SubControlVolumeFace& scvf) const
    { return face_(scvf).neighborScvCacheIndices.size() - 1; }

    // //! return the i-th "outside" scvf that this sub-control volume face coincides with
    // SubControlVolumeFace flipScvf(const SubControlVolumeFace& scvf, unsigned int i = 0) const
    // {
    //     assert(i < numFlipScvf(scvf));
    //     if (i >= cache_(scvf).indexInFaceNeighbors)
    //         i++;
    //     return SubControlVolumeFace{face_(scvf).neighborScvfCacheIndices[i]};
    // }

    // //! return the number of "outside" scvfs that this scvf coincides with
    // std::size_t numFlipScvf(const SubControlVolumeFace& scvf) const
    // { return numOutsideScvs(scvf); }

    //! return true if the given sub-control volume face is on the boundary
    bool onBoundary(const SubControlVolumeFace& scvf) const
    { return numOutsideScvs(scvf) == 0; }

    //! iterator range over all scvs inside the bound element
    friend auto scvs(const GridGeometryLocalView& lv)
    {
        return std::views::transform(
            lv.localCache_.insideScvIndices(),
            [&] (std::size_t idx) { return SubControlVolume{idx}; }
        );
    }

    //! iterator range over all scvfs inside the bound element
    friend auto scvfs(const GridGeometryLocalView& lv)
    {
        return std::views::transform(
            lv.localCache_.insideScvfIndices(),
            [&] (std::size_t idx) { return SubControlVolumeFace{idx}; }
        );
    }

private:
    const auto& cache_(const SubControlVolume& scv) const
    { return localCache_.scvCache(scv.id); }

    const auto& cache_(const SubControlVolumeFace& scvf) const
    { return localCache_.scvfCache(scvf.id); }

    const auto& face_(const SubControlVolumeFace& scvf) const
    { return localCache_.faceCache(localCache_.scvfCache(scvf.id)); }

    const GridGeometry& gridGeometry_;
    const GridGeometryCache& geometryCache_;
    LocalCache localCache_;
};

} // end namespace Dumux::CCTpfa

#endif
