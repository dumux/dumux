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
#ifndef DUMUX_DISCRETIZATION_CCTPFA_GRID_GEOMETRY_CACHE_HH
#define DUMUX_DISCRETIZATION_CCTPFA_GRID_GEOMETRY_CACHE_HH
#ifndef DOXYGEN

#include <vector>
#include <ranges>
#include <cstdint>
#include <type_traits>
#include <concepts>
#include <limits>

#include <dune/common/reservedvector.hh>

#include <dumux/common/enumerate.hh>
#include <dumux/geometry/volume.hh>

#include <dumux/experimental/new_assembly/dumux/common/typetraits.hh>
#include <dumux/experimental/new_assembly/dumux/common/storage.hh>
#include <dumux/experimental/new_assembly/dumux/common/concepts.hh>
#include <dumux/experimental/new_assembly/dumux/geometry/normal.hh>

#include <dumux/experimental/new_assembly/dumux/discretization/cctpfa/concepts.hh>

namespace Dumux::CCTpfa::Detail {

template<std::size_t dim>
inline constexpr std::size_t maxNumElementFacets = 2*dim; // cubes

template<Concepts::CCTpfaGridGeometryTraits Traits>
inline constexpr auto maxNumScvfsPerFacet = Size::pow<2, Traits::maxAdjacentElementLevelDifference>();

template<std::size_t dim, Concepts::CCTpfaGridGeometryTraits Traits>
inline constexpr auto maxNumElementScvfs = maxNumElementFacets<dim>*maxNumScvfsPerFacet<Traits>;

template<std::size_t dim, Concepts::CCTpfaGridGeometryTraits Traits>
inline constexpr auto maxNumLocalScvfs = maxNumElementScvfs<dim, Traits>*Traits::maxNumBranchesPerScvf;

template<std::size_t dim, Concepts::CCTpfaGridGeometryTraits Traits>
inline constexpr auto maxElementStencilSize = maxNumLocalScvfs<dim, Traits> + 1;

template<Concepts::CCTpfaGridGeometryTraits Traits>
inline constexpr auto maxNumScvfNeighbors = Traits::maxNumBranchesPerScvf + 1;

template<typename GG>
using GGElement = typename GG::GridView::template Codim<0>::Entity;

template<typename GV>
using GridViewCType = typename GV::ctype;

template<typename GV>
using GridViewIndex = typename GV::IndexSet::IndexType;

template<typename GV>
using GridViewCoordinate = typename GV::template Codim<0>::Geometry::GlobalCoordinate;

template<typename GV>
inline constexpr bool gridViewIsNetwork = int(GV::dimension) < int(GV::dimensionworld);

template<typename GV,
         typename Facet,
         std::integral LocalIndex,
         auto size = Size::dynamic>
struct FaceCache
{
    const Facet* facet{nullptr};
    DefaultStorage<LocalIndex, size> neighborScvCacheIndices;
    // ScvfNeighborDataStorage<gridViewIsNetwork<GV>, LocalIndex> neighborScvfCacheIndices;
    GridViewCType<GV> area;
    GridViewCoordinate<GV> center;
};

template<typename GV>
struct ScvCache
{
    GridViewCType<GV> volume;
    GridViewCoordinate<GV> center;
    GridViewIndex<GV> dofIndex;
};

template<typename GV, std::integral LocalIndex>
struct ScvfCache
{
    LocalIndex faceCacheIndex;
    std::uint_least8_t indexInFaceNeighbors;
    GridViewCoordinate<GV> normal;
};

template<Concepts::Indexable ScvCacheStorage,
         Concepts::Indexable ScvfCacheStorage,
         Concepts::Indexable FaceCacheStorage>
struct GridGeometryCacheStorage
{
private:
    template<typename Storage>
    using StoredType = std::decay_t<decltype(std::declval<Storage>()[0])>;

public:
    using ScvCache = StoredType<ScvCacheStorage>;
    using ScvfCache = StoredType<ScvfCacheStorage>;
    using FaceCache = StoredType<FaceCacheStorage>;

    ScvCacheStorage scvCaches;
    ScvfCacheStorage scvfCaches;
    FaceCacheStorage faceCaches;
};


template<typename GridGeometry, typename Traits, bool cacheGlobally>
class GridGeometryCache;

template<typename GridGeometry, typename Traits>
class GridGeometryCache<GridGeometry, Traits, false>
{
    using GridView = typename GridGeometry::GridView;
    using Element = typename GridView::template Codim<0>::Entity;
    using ElementGeometry = typename Element::Geometry;
    using Coordinate = typename ElementGeometry::GlobalCoordinate;

    using GridIndex = typename GridView::IndexSet::IndexType;
    using LocalIndex = std::uint_least8_t;

    using ScvfSeed = typename GridGeometry::ScvfSeed;
    using Facet = typename ScvfSeed::Facet;

    using FaceCache = Detail::FaceCache<GridView, Facet, LocalIndex, maxNumScvfNeighbors<Traits>>;
    using ScvCache = Detail::ScvCache<GridView>;
    using ScvfCache = Detail::ScvfCache<GridView, LocalIndex>;

    using CacheStorage = GridGeometryCacheStorage<
        DefaultStorage<ScvCache, maxElementStencilSize<GridView::dimension, Traits>>,
        DefaultStorage<ScvfCache, maxNumLocalScvfs<GridView::dimension, Traits>>,
        DefaultStorage<FaceCache, maxNumElementScvfs<GridView::dimension, Traits>>
    >;

    struct BindContext
    {
        const GridGeometry& gridGeometry;
        const Element& element;
        const ElementGeometry elementGeometry;
        const GridIndex elementIndex;

        BindContext(const GridGeometry& gg, const Element& e)
        : gridGeometry(gg)
        , element(e)
        , elementGeometry(e.geometry())
        , elementIndex(gg.elementMapper().index(e))
        {}
    };

public:
    class LocalCache
    {
    public:
        const ScvCache& scvCache(std::size_t i) const { return cacheStorage_.scvCaches[i]; }
        const ScvfCache& scvfCache(std::size_t i) const { return cacheStorage_.scvfCaches[i]; }

        constexpr auto insideScvIndices() const { return std::views::single(std::size_t{0}); }
        auto insideScvfIndices() const { return std::views::iota(LocalIndex{0}, numInsideScvfs_); }

        const FaceCache& faceCache(const ScvfCache& cache) const
        { return cacheStorage_.faceCaches[static_cast<std::size_t>(cache.faceCacheIndex)]; }

    private:
        friend class GridGeometryCache;
        CacheStorage cacheStorage_;
        LocalIndex numInsideScvfs_;
    };

    void update(const GridGeometry&)
    {}

    void bindElement(LocalCache& localCache,
                     const GridGeometry& gg,
                     const Element& element) const
    { updateElementCaches_(localCache, BindContext{gg, element}); }

    void bindStencil(LocalCache& localCache,
                     const GridGeometry& gg,
                     const Element& element) const
    { updateStencilCaches_(localCache, BindContext{gg, element}); }

private:
    template<std::integral I>
    LocalIndex makeLocalIndex_(const I& i) const
    { return static_cast<LocalIndex>(i); }

    std::ranges::range auto scvfSeeds_(const BindContext& context,
                                       unsigned int facetIdx) const
    { return context.gridGeometry.scvfSeeds(context.element, facetIdx); }

    void clearCaches_(LocalCache& localCache) const
    {
        localCache.cacheStorage_.scvCaches.clear();
        localCache.cacheStorage_.scvfCaches.clear();
        localCache.cacheStorage_.faceCaches.clear();
    }

    void updateStencilCaches_(LocalCache& localCache, const BindContext& context) const
    {
        updateElementCaches_(localCache, context);
        makeOutsideEntities_(localCache, context);
    }

    void updateElementCaches_(LocalCache& localCache, const BindContext& context) const
    {
        clearCaches_(localCache);
        addInsideScv_(localCache, context);
        makeInsideScvfs_(localCache, context);
    }

    void addInsideScv_(LocalCache& localCache, const BindContext& context) const
    {
        localCache.cacheStorage_.scvCaches.push_back(ScvCache{
            convexPolytopeVolume(context.elementGeometry),
            context.elementGeometry.center(),
            context.elementIndex
        });
    }

    void makeInsideScvfs_(LocalCache& localCache, const BindContext& context) const
    {
        const auto facetIndices = std::views::iota(unsigned{0}, context.element.subEntities(1));
        std::ranges::for_each(facetIndices, [&] (unsigned int facetIdx) {
            std::ranges::for_each(scvfSeeds_(context, facetIdx), [&] (const ScvfSeed& seed) {
                addInsideScvf_(localCache, context, facetIdx, seed);
            });
        });
        localCache.numInsideScvfs_ = localCache.cacheStorage_.scvfCaches.size();
    }

    void addInsideScvf_(LocalCache& localCache,
                        const BindContext& context,
                        const unsigned int facetIndex,
                        const ScvfSeed& seed) const
    {
        const auto& seedFacet = seed.facet();
        const Element& element = context.gridGeometry.element(seedFacet.elementIndex);
        const auto& faceGeo = element.template subEntity<1>(seedFacet.facetIndex).geometry();
        const auto localFaceIdx = addFace_(localCache, faceGeo, seedFacet, LocalIndex{0});
        auto outerNormal = getNormal(context.elementGeometry, facetIndex);
        addScvf_(localCache, std::move(outerNormal), localFaceIdx, LocalIndex{0});
    }

    template<typename FaceGeometry>
    LocalIndex addFace_(LocalCache& localCache,
                        const FaceGeometry& faceGeo,
                        const Facet& facet,
                        const LocalIndex insideScvCacheIndex) const
    {
        const auto nextFaceIdx = localCache.cacheStorage_.faceCaches.size();
        localCache.cacheStorage_.faceCaches.push_back(FaceCache{
            &facet,
            {insideScvCacheIndex},
            convexPolytopeVolume(faceGeo),
            faceGeo.center()
        });
        return nextFaceIdx;
    }

    LocalIndex addScvf_(LocalCache& localCache,
                        Coordinate&& normal,
                        const LocalIndex faceIndex,
                        const LocalIndex idxInFaceNeighbors) const
    {
        const auto nextScvfIdx = localCache.cacheStorage_.scvfCaches.size();
        localCache.cacheStorage_.scvfCaches.push_back(ScvfCache{
            faceIndex,
            idxInFaceNeighbors,
            std::move(normal)
        });
        return makeLocalIndex_(nextScvfIdx);
    }

    void makeOutsideEntities_(LocalCache& localCache, const BindContext& context) const
    {
        LocalIndex curFaceIndex = 0;
        const auto facetIndices = std::views::iota(unsigned{0}, context.element.subEntities(1));
        std::ranges::for_each(facetIndices, [&] (unsigned int facetIdx) {
            std::ranges::for_each(scvfSeeds_(context, facetIdx), [&] (const ScvfSeed& seed) {
                makeOutsideScvfsAndScvs_(localCache, context, curFaceIndex, seed);
                curFaceIndex++;
            });
        });
        assert(curFaceIndex == localCache.numInsideScvfs_);
    }

    void makeOutsideScvfsAndScvs_(LocalCache& localCache,
                                  const BindContext& context,
                                  const LocalIndex faceIndex,
                                  const ScvfSeed& seed) const
    {
        auto& faceCacheStorage = localCache.cacheStorage_.faceCaches[faceIndex];
        auto facetsToVisit = std::views::filter(seed.outsideFacets(), [&] (const Facet& facet) {
            return facet.elementIndex != context.elementIndex;
        });
        std::ranges::for_each(facetsToVisit, [&] (const Facet& facet) {
            const auto& element = context.gridGeometry.element(facet.elementIndex);
            const auto& elementGeometry = element.geometry();
            const auto elementIndex = context.gridGeometry.elementMapper().index(element);

            const auto localScvIdx = addScv_(localCache, elementGeometry, elementIndex);
            const auto idxInNeighbors = faceCacheStorage.neighborScvCacheIndices.size();
            faceCacheStorage.neighborScvCacheIndices.push_back(localScvIdx);
            addScvf_(
                localCache,
                getNormal(elementGeometry, facet.facetIndex),
                faceIndex,
                makeLocalIndex_(idxInNeighbors));
        });
    }

    LocalIndex addScv_(LocalCache& localCache,
                       const ElementGeometry& elemGeo,
                       const GridIndex elemIdx) const
    {
        const auto nextScvIdx = localCache.cacheStorage_.scvCaches.size();
        localCache.cacheStorage_.scvCaches.push_back(ScvCache{
            convexPolytopeVolume(elemGeo),
            elemGeo.center(),
            elemIdx
        });
        return makeLocalIndex_(nextScvIdx);
    }
};

template<typename GridGeometry, typename Traits>
class GridGeometryCache<GridGeometry, Traits, true>
{
    using GridView = typename GridGeometry::GridView;
    using Element = typename GridView::template Codim<0>::Entity;

    using GridIndex = typename GridView::IndexSet::IndexType;
    using LocalIndex = std::size_t;
    using ScvfSeed = typename GridGeometry::ScvfSeed;
    using Facet = typename ScvfSeed::Facet;

    using FaceCache = Detail::FaceCache<GridView, Facet, LocalIndex, maxNumScvfNeighbors<Traits>>;
    using ScvCache = Detail::ScvCache<GridView>;
    using ScvfCache = Detail::ScvfCache<GridView, LocalIndex>;

    using CacheStorage = GridGeometryCacheStorage<
        DefaultStorage<ScvCache, Size::dynamic>,
        DefaultStorage<ScvfCache, Size::dynamic>,
        DefaultStorage<FaceCache, Size::dynamic>
    >;

    static constexpr int dim = GridView::dimension;
    using FacetScvfIndices = DefaultStorage<std::size_t, maxNumScvfsPerFacet<Traits>>;
    using FacetToScvfIndices = DefaultStorage<FacetScvfIndices, maxNumElementFacets<dim>>;

public:
    class LocalCache
    {
    public:
        const ScvCache& scvCache(std::size_t i) const
        { return gridCachePtr_->caches_.scvCaches[i]; }

        const ScvfCache& scvfCache(std::size_t i) const
        { return gridCachePtr_->caches_.scvfCaches[i]; }

        constexpr auto insideScvIndices() const
        { return std::views::single(*elementIndex_); }

        auto insideScvfIndices() const
        {
            return std::views::join(std::views::transform(
                std::views::iota(unsigned{0}, element_->subEntities(1)), [&] (auto facetIdx) {
                    return std::views::all(
                        gridCachePtr_->elementFacetToScvfIndices_[*elementIndex_][facetIdx]
                    );
            }));
        }

        const FaceCache& faceCache(const ScvfCache& cache) const
        { return gridCachePtr_->caches_.faceCaches[cache.faceCacheIndex]; }

    private:
        friend class GridGeometryCache;
        const GridGeometryCache* gridCachePtr_;
        std::optional<Element> element_;
        std::optional<GridIndex> elementIndex_;
    };

    void update(const GridGeometry& gridGeometry)
    {
        resizeCaches_(gridGeometry);
        for (const auto& element : elements(gridGeometry.gridView()))
        {
            makeScv_(gridGeometry, element);
            makeFaces_(gridGeometry, element);
        }

        assert(std::ranges::all_of(elementFacetToScvfIndices_, [] (const auto& elemFacetScvfs) {
            return std::ranges::none_of(elemFacetScvfs, [] (const auto& facetScvfs) {
                return facetScvfs.empty();
            });
        }));
    }

    void bindElement(LocalCache& localCache,
                     const GridGeometry& gg,
                     const Element& element) const
    { bind_(localCache, gg, element); }

    void bindStencil(LocalCache& localCache,
                     const GridGeometry& gg,
                     const Element& element) const
    { bind_(localCache, gg, element); }

private:
    friend class LocalCache;

    void resizeCaches_(const GridGeometry& gridGeometry)
    {
        caches_.scvfCaches.clear();
        caches_.scvfCaches.reserve(gridGeometry.numScvf());

        caches_.faceCaches.clear();
        caches_.faceCaches.reserve(gridGeometry.gridView().size(1));

        caches_.scvCaches.assign(gridGeometry.numScv(), {});

        elementFacetToScvfIndices_.resize(gridGeometry.gridView().size(0));
        for (const auto& e : elements(gridGeometry.gridView()))
            elementFacetToScvfIndices_[gridGeometry.elementMapper().index(e)].resize(
                e.subEntities(1)
            );
    }

    void makeScv_(const GridGeometry& gridGeometry, const Element& element)
    {
        const auto eIdx = gridGeometry.elementMapper().index(element);
        const auto& geo = element.geometry();
        caches_.scvCaches[eIdx] = ScvCache{
            convexPolytopeVolume(geo),
            geo.center(),
            eIdx
        };
    }

    void makeFaces_(const GridGeometry& gridGeometry, const Element& element)
    {
        const auto eIdx = gridGeometry.elementMapper().index(element);
        const auto facetIndices = std::views::iota(unsigned{0}, element.subEntities(1));
        std::ranges::for_each(facetIndices, [&] (unsigned int facetIdx) {
            const auto& facetScvfs = gridGeometry.scvfSeeds(element, facetIdx);
            std::ranges::for_each(facetScvfs, [&] (const ScvfSeed& seed) {
                if (seed.facet().elementIndex == eIdx && seed.facet().facetIndex == facetIdx)
                {
                    const auto curFaceIdx = addFace_(gridGeometry, element, seed.facet());
                    addFaceScvfs_(curFaceIdx, seed, gridGeometry);
                }
            });
        });
    }

    std::size_t addFace_(const GridGeometry& gridGeometry,
                         const Element& element,
                         const Facet& facet)
    {
        const auto nextFaceIndex = caches_.faceCaches.size();
        const auto& facetGeom = element.template subEntity<1>(facet.facetIndex).geometry();
        caches_.faceCaches.push_back(FaceCache{
            &facet,
            {},
            // {},
            convexPolytopeVolume(facetGeom),
            facetGeom.center()
        });
        return nextFaceIndex;
    }

    // TODO: stream line argument order!?
    void addFaceScvfs_(std::size_t faceIndex,
                       const ScvfSeed& seed,
                       const GridGeometry& gridGeometry)
    {
        std::uint_least8_t idxInNeighbors = 0;
        for (const auto& facet : seed.neighboringFacets())
        {
            const auto element = gridGeometry.element(facet.elementIndex);
            const auto nextScvfIndex = caches_.scvfCaches.size();
            caches_.scvfCaches.push_back(ScvfCache{
                faceIndex,
                idxInNeighbors++,
                getNormal(element.geometry(), facet.facetIndex)
            });
            caches_.faceCaches[faceIndex].neighborScvCacheIndices.push_back(facet.elementIndex);
            elementFacetToScvfIndices_[facet.elementIndex][facet.facetIndex].push_back(nextScvfIndex);
        }
    }

    void bind_(LocalCache& localCache,
               const GridGeometry& gridGeom,
               const Element& element) const
    {
        localCache.gridCachePtr_ = this;
        localCache.element_ = element;
        localCache.elementIndex_ = gridGeom.elementMapper().index(element);
    }

    CacheStorage caches_;
    std::vector<FacetToScvfIndices> elementFacetToScvfIndices_;
};

} // namespace Dumux::CCTpfa::Detail

#endif // DOXYGEN
#endif
