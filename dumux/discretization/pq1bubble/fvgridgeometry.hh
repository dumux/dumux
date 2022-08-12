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
 * \ingroup PQ1BubbleDiscretization
 * \brief Base class for the finite volume geometry vector for the pq1bubble method
 *        This builds up the sub control volumes and sub control volume faces
 *        for each element of the grid partition.
 */
#ifndef DUMUX_DISCRETIZATION_PQ1BUBBLE_GRID_GEOMETRY_HH
#define DUMUX_DISCRETIZATION_PQ1BUBBLE_GRID_GEOMETRY_HH

#include <utility>
#include <unordered_map>

#include <dune/grid/common/mcmgmapper.hh>

#include <dumux/discretization/method.hh>
#include <dumux/common/indextraits.hh>
#include <dumux/common/defaultmappertraits.hh>
#include <dumux/discretization/basegridgeometry.hh>
#include <dumux/discretization/pq1bubble/geometryhelper.hh>
#include <dumux/discretization/pq1bubble/fvelementgeometry.hh>
#include <dumux/discretization/pq1bubble/subcontrolvolume.hh>
#include <dumux/discretization/pq1bubble/subcontrolvolumeface.hh>
#include <dumux/discretization/pq1bubble/pq1bubblefecache.hh>
#include <dumux/discretization/extrusion.hh>

namespace Dumux {

template <class GridView>
struct PQ1BubbleMapperTraits :public DefaultMapperTraits<GridView>
{
    using DofMapper = Dune::MultipleCodimMultipleGeomTypeMapper<GridView>;

    /**
     * \brief layout for elements and vertices
     *
     */
    static Dune::MCMGLayout layout()
    {
        return [](Dune::GeometryType gt, int dimgrid) {
            return (gt.dim() == dimgrid) || (gt.dim() == 0);
        };
    }
};

/*!
 * \ingroup PQ1BubbleDiscretization
 * \brief The default traits for the pq1bubble finite volume grid geometry
 *        Defines the scv and scvf types and the mapper types
 * \tparam the grid view type
 */
template<class GridView, class MapperTraits = PQ1BubbleMapperTraits<GridView>>
struct PQ1BubbleDefaultGridGeometryTraits
: public MapperTraits
{
    using SubControlVolume = PQ1BubbleSubControlVolume<GridView>;
    using SubControlVolumeFace = PQ1BubbleSubControlVolumeFace<GridView>;

    template<class GridGeometry, bool enableCache>
    using LocalView = PQ1BubbleFVElementGeometry<GridGeometry, enableCache>;
};

/*!
 * \ingroup PQ1BubbleDiscretization
 * \brief Base class for the finite volume geometry vector for pq1bubble schemes
 *        This builds up the sub control volumes and sub control volume faces
 * \note For caching enabled we store the fv geometries for the whole grid view which is memory intensive but faster
 */
template<class Scalar,
         class GV,
         bool enableCaching = true,
         class Traits = PQ1BubbleDefaultGridGeometryTraits<GV>>
class PQ1BubbleFVGridGeometry
: public BaseGridGeometry<GV, Traits>
{
    using ThisType = PQ1BubbleFVGridGeometry<Scalar, GV, enableCaching, Traits>;
    using ParentType = BaseGridGeometry<GV, Traits>;
    using GridIndexType = typename IndexTraits<GV>::GridIndex;
    using LocalIndexType = typename IndexTraits<GV>::LocalIndex;

    using Element = typename GV::template Codim<0>::Entity;
    using CoordScalar = typename GV::ctype;
    static const int dim = GV::dimension;
    static const int dimWorld = GV::dimensionworld;

    using GeometryHelper = PQ1BubbleGeometryHelper<
        GV, typename Traits::SubControlVolume, typename Traits::SubControlVolumeFace
    >;

    static_assert(dim > 1, "Only implemented for dim > 1");

public:
    //! export the discretization method this geometry belongs to
    using DiscretizationMethod = DiscretizationMethods::PQ1Bubble;
    static constexpr DiscretizationMethod discMethod{};

    //! export the type of the fv element geometry (the local view type)
    using LocalView = typename Traits::template LocalView<ThisType, true>;
    //! export the type of sub control volume
    using SubControlVolume = typename Traits::SubControlVolume;
    //! export the type of sub control volume
    using SubControlVolumeFace = typename Traits::SubControlVolumeFace;
    //! export the type of extrusion
    using Extrusion = Extrusion_t<Traits>;
    //! export dof mapper type
    using DofMapper = typename Traits::DofMapper;
    //! export the finite element cache type
    using FeCache = Dumux::PQ1BubbleFECache<CoordScalar, Scalar, dim>;
    //! export the grid view type
    using GridView = GV;

    //! Constructor
    PQ1BubbleFVGridGeometry(const GridView gridView)
    : ParentType(gridView), dofMapper_(gridView, Traits::layout())
    {
        update_();
    }

    //! The dofMapper
    const DofMapper& dofMapper() const
    { return dofMapper_; }

    //! The total number of sub control volumes
    std::size_t numScv() const
    { return numScv_; }

    //! The total number of sub control volume faces
    std::size_t numScvf() const
    { return numScvf_; }

    //! The total number of boundary sub control volume faces
    std::size_t numBoundaryScvf() const
    { return numBoundaryScvf_; }

    //! The total number of degrees of freedom
    std::size_t numDofs() const
    { return this->dofMapper().size(); }

    //! update all geometries (call this after grid adaption)
    void update(const GridView& gridView)
    {
        ParentType::update(gridView);
        update_();
    }

    //! update all geometries (call this after grid adaption)
    void update(GridView&& gridView)
    {
        ParentType::update(std::move(gridView));
        update_();
    }

    //! The finite element cache for creating local FE bases
    const FeCache& feCache() const
    { return feCache_; }

    //! Get the local scvs for an element
    const std::vector<SubControlVolume>& scvs(GridIndexType eIdx) const
    { return scvs_[eIdx]; }

    //! Get the local scvfs for an element
    const std::vector<SubControlVolumeFace>& scvfs(GridIndexType eIdx) const
    { return scvfs_[eIdx]; }

    //! If a vertex / d.o.f. is on the boundary
    bool dofOnBoundary(GridIndexType dofIdx) const
    { return boundaryDofIndices_[dofIdx]; }

    //! If a vertex / d.o.f. is on a periodic boundary
    bool dofOnPeriodicBoundary(GridIndexType dofIdx) const
    { return periodicVertexMap_.count(dofIdx); }

    //! The index of the vertex / d.o.f. on the other side of the periodic boundary
    GridIndexType periodicallyMappedDof(GridIndexType dofIdx) const
    { return periodicVertexMap_.at(dofIdx); }

    //! Returns the map between dofs across periodic boundaries
    const std::unordered_map<GridIndexType, GridIndexType>& periodicVertexMap() const
    { return periodicVertexMap_; }

    //! Returns whether one of the geometry's scvfs lies on a boundary
    bool hasBoundaryScvf(GridIndexType eIdx) const
    { return hasBoundaryScvf_[eIdx]; }

private:
    void update_()
    {
        scvs_.clear();
        scvfs_.clear();
        dofMapper_.update(this->gridView());

        auto numElements = this->gridView().size(0);
        scvs_.resize(numElements);
        scvfs_.resize(numElements);
        hasBoundaryScvf_.resize(numElements, false);

        boundaryDofIndices_.assign(numDofs(), false);

        numScv_ = 0;
        numScvf_ = 0;
        numBoundaryScvf_ = 0;

        // Build the scvs and scv faces
        for (const auto& element : elements(this->gridView()))
        {
            auto eIdx = this->elementMapper().index(element);

            // get the element geometry
            auto elementGeometry = element.geometry();
            const auto refElement = referenceElement(elementGeometry);

            // instantiate the geometry helper
            GeometryHelper geometryHelper(elementGeometry);

            numScv_ += geometryHelper.numScv();
            // construct the sub control volumes
            scvs_[eIdx].resize(geometryHelper.numScv());
            for (LocalIndexType scvLocalIdx = 0; scvLocalIdx < geometryHelper.numScv(); ++scvLocalIdx)
            {
                auto corners = geometryHelper.getScvCorners(scvLocalIdx);
                scvs_[eIdx][scvLocalIdx] = SubControlVolume(geometryHelper.scvVolume(scvLocalIdx, corners),
                                                            geometryHelper.dofPosition(scvLocalIdx),
                                                            geometryHelper.scvCenter(scvLocalIdx, corners),
                                                            scvLocalIdx,
                                                            eIdx,
                                                            geometryHelper.dofIndex(this->dofMapper(), element, scvLocalIdx)
                                                            );
            }

            // construct the sub control volume faces
            numScvf_ += geometryHelper.numInteriorScvf();
            scvfs_[eIdx].resize(geometryHelper.numInteriorScvf());
            LocalIndexType scvfLocalIdx = 0;
            for (; scvfLocalIdx < geometryHelper.numInteriorScvf(); ++scvfLocalIdx)
            {
                auto scvPair = geometryHelper.getScvPairForScvf(scvfLocalIdx);
                auto corners = geometryHelper.getScvfCorners(scvfLocalIdx);

                // the sub control volume face
                typename SubControlVolumeFace::Traits::Geometry scvfGeo{
                    geometryHelper.getScvfGeometryType(scvfLocalIdx), corners
                };

                scvfs_[eIdx][scvfLocalIdx] = SubControlVolumeFace(scvfGeo.center(),
                                                                  scvfGeo.volume(),
                                                                  geometryHelper.normal(corners, scvPair),
                                                                  std::move(scvPair),
                                                                  scvfLocalIdx,
                                                                  geometryHelper.isOverlapping(scvfLocalIdx));
            }

            // construct the sub control volume faces on the domain boundary
            for (const auto& intersection : intersections(this->gridView(), element))
            {
                if (intersection.boundary() && !intersection.neighbor())
                {
                    hasBoundaryScvf_[eIdx] = true;

                    numScvf_ += geometryHelper.numBoundaryIsScvf(intersection);
                    numBoundaryScvf_ += geometryHelper.numBoundaryIsScvf(intersection);

                    for (unsigned int isScvfLocalIdx = 0; isScvfLocalIdx < geometryHelper.numBoundaryIsScvf(intersection); ++isScvfLocalIdx)
                    {
                        // find the scvs this scvf is belonging to
                        auto scvPair = geometryHelper.getScvPairForBoundaryScvf(intersection, isScvfLocalIdx);
                        auto corners = geometryHelper.getBoundaryScvfCorners(intersection.indexInInside(), isScvfLocalIdx);

                        // the sub control volume face
                        typename SubControlVolumeFace::Traits::Geometry scvfGeo{
                            geometryHelper.getBoundaryScvfGeometryType(scvfLocalIdx), corners
                        };

                        scvfs_[eIdx].emplace_back(scvfGeo.center(),
                                                  scvfGeo.volume(),
                                                  intersection.centerUnitOuterNormal(),
                                                  std::move(scvPair),
                                                  scvfLocalIdx,
                                                  typename SubControlVolumeFace::Traits::BoundaryFlag{ intersection });

                        // increment local counter
                        scvfLocalIdx++;
                    }

                    // TODO also move this to helper class

                    // add all vertices on the intersection to the set of
                    // boundary vertices
                    const auto fIdx = intersection.indexInInside();
                    const auto numFaceVerts = refElement.size(fIdx, 1, dim);
                    for (int localVIdx = 0; localVIdx < numFaceVerts; ++localVIdx)
                    {
                        const auto vIdx = refElement.subEntity(fIdx, 1, localVIdx, dim);
                        const auto vIdxGlobal = this->dofMapper().subIndex(element, vIdx, dim);
                        boundaryDofIndices_[vIdxGlobal] = true;
                    }
                }

                // inform the grid geometry if we have periodic boundaries
                else if (intersection.boundary() && intersection.neighbor())
                {
                    this->setPeriodic();

                    // find the mapped periodic vertex of all vertices on periodic boundaries
                    const auto fIdx = intersection.indexInInside();
                    const auto numFaceVerts = refElement.size(fIdx, 1, dim);
                    const auto eps = 1e-7*(elementGeometry.corner(1) - elementGeometry.corner(0)).two_norm();
                    for (int localVIdx = 0; localVIdx < numFaceVerts; ++localVIdx)
                    {
                        const auto vIdx = refElement.subEntity(fIdx, 1, localVIdx, dim);
                        const auto vIdxGlobal = this->dofMapper().subIndex(element, vIdx, dim);
                        const auto vPos = elementGeometry.corner(vIdx);

                        const auto& outside = intersection.outside();
                        const auto outsideGeometry = outside.geometry();
                        for (const auto& isOutside : intersections(this->gridView(), outside))
                        {
                            // only check periodic vertices of the periodic neighbor
                            if (isOutside.boundary() && isOutside.neighbor())
                            {
                                const auto fIdxOutside = isOutside.indexInInside();
                                const auto numFaceVertsOutside = refElement.size(fIdxOutside, 1, dim);
                                for (int localVIdxOutside = 0; localVIdxOutside < numFaceVertsOutside; ++localVIdxOutside)
                                {
                                    const auto vIdxOutside = refElement.subEntity(fIdxOutside, 1, localVIdxOutside, dim);
                                    const auto vPosOutside = outsideGeometry.corner(vIdxOutside);
                                    const auto shift = std::abs((this->bBoxMax()-this->bBoxMin())*intersection.centerUnitOuterNormal());
                                    if (std::abs((vPosOutside-vPos).two_norm() - shift) < eps)
                                        periodicVertexMap_[vIdxGlobal] = this->dofMapper().subIndex(outside, vIdxOutside, dim);
                                }
                            }
                        }
                    }
                }
            }
        }

        // error check: periodic boundaries currently don't work for cvfe in parallel
        if (this->isPeriodic() && this->gridView().comm().size() > 1)
            DUNE_THROW(Dune::NotImplemented, "Periodic boundaries for cvfe method for parallel simulations!");
    }

    DofMapper dofMapper_;

    const FeCache feCache_;

    std::vector<std::vector<SubControlVolume>> scvs_;
    std::vector<std::vector<SubControlVolumeFace>> scvfs_;

    std::size_t numScv_;
    std::size_t numScvf_;
    std::size_t numBoundaryScvf_;

    // vertices on the boundary
    std::vector<bool> boundaryDofIndices_;
    std::vector<bool> hasBoundaryScvf_;

    // a map for periodic boundary vertices
    std::unordered_map<GridIndexType, GridIndexType> periodicVertexMap_;
};

} // end namespace Dumux

#endif
