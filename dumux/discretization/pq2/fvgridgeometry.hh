// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup PQ2Discretization
 * \brief Base class for the finite volume geometry vector for the pq2 method
 *        This builds up the sub control volumes and sub control volume faces
 *        for each element of the grid partition.
 */
#ifndef DUMUX_DISCRETIZATION_PQ2_GRID_GEOMETRY_HH
#define DUMUX_DISCRETIZATION_PQ2_GRID_GEOMETRY_HH

#include <array>
#include <cstddef>
#include <memory>
#include <vector>
#include <utility>
#include <unordered_map>
#include <type_traits>

#include <dune/grid/common/mcmgmapper.hh>
#include <dune/geometry/type.hh>

#include <dumux/discretization/method.hh>
#include <dumux/discretization/pq2/pq2hierarchicalfecache.hh>

#include <dumux/common/indextraits.hh>
#include <dumux/common/defaultmappertraits.hh>

#include <dumux/geometry/center.hh>
#include <dumux/geometry/volume.hh>

#include <dumux/discretization/basegridgeometry.hh>
#include <dumux/discretization/pq1bubble/fvgridgeometry.hh>
#include <dumux/discretization/pq2/geometryhelper.hh>
#include <dumux/discretization/pq2/fvelementgeometry.hh>
#include <dumux/discretization/pq2/subcontrolvolume.hh>
#include <dumux/discretization/pq2/subcontrolvolumeface.hh>
#include <dumux/discretization/extrusion.hh>

#include <dumux/io/grid/periodicgridtraits.hh>

namespace Dumux {

namespace Detail {

template<class GV, class T>
using PQ2GeometryHelper_t = Dune::Std::detected_or_t<
    std::conditional_t<enablesHybridCVFE<T>,
        Dumux::HybridPQ2GeometryHelper<GV, typename T::SubControlVolume, typename T::SubControlVolumeFace>,
        void // we currently only support hybrid schemes
    >,
    SpecifiesGeometryHelper,
    T
>;
} // end namespace Detail

template <class GridView>
struct PQ2MapperTraits :public DefaultMapperTraits<GridView>
{
    using DofMapper = Dune::MultipleCodimMultipleGeomTypeMapper<GridView>;

    /**
     * \brief layout for vertices and edges
     *
     */
    static Dune::MCMGLayout layout()
    {
        return [](Dune::GeometryType gt, int dimgrid) {
            return (gt.dim() == 0) || (gt.dim() == 1);
        };
    }
};

/*!
 * \ingroup PQ2Discretization
 * \brief Quadrature rule traits for PQ2 discretization
 */
template<class GridView,
         class ScvRule = Dumux::QuadratureRules::MidpointQuadrature,
         class ScvfRule = Dumux::QuadratureRules::DuneQuadrature<2>,
         class ElementRule = Dumux::QuadratureRules::DuneQuadrature<4>,
         class IntersectionRule = Dumux::QuadratureRules::DuneQuadrature<4>>
using PQ2QuadratureTraits = CVFE::DefaultQuadratureTraits<GridView, ScvRule, ScvfRule, ElementRule, IntersectionRule>;

/*!
 * \ingroup PQ2Discretization
 * \brief The default traits for the pq2 finite volume grid geometry
 *        Defines the scv and scvf types and the mapper types
 * \tparam the grid view type
 */
template<class GridView, class MapperTraits = PQ2MapperTraits<GridView>, class QuadratureTraits = PQ2QuadratureTraits<GridView>>
struct PQ2DefaultGridGeometryTraits
: public MapperTraits, public QuadratureTraits
{
    using SubControlVolume = PQ2SubControlVolume<GridView>;
    using SubControlVolumeFace = PQ2SubControlVolumeFace<GridView>;

    template<class GridGeometry, bool enableCache>
    using LocalView = PQ2FVElementGeometry<GridGeometry, enableCache>;

    using EnableHybridCVFE = std::true_type;

    // The maximum values correspond to cubes
    // This can be overwritten by the user when knowing the grid entity types
    static constexpr std::size_t maxNumElementDofs = []()
        {
            if constexpr (GridView::dimension == 1)
                return 3;
            else if constexpr (GridView::dimension == 2)
                return 9;
            else if constexpr (GridView::dimension == 3)
                return 27;
        }();
};

/*!
 * \ingroup PQ2Discretization
 * \brief Base class for the finite volume geometry vector for pq2 schemes
 *        This builds up the sub control volumes and sub control volume faces
 * \note For caching enabled we store the fv geometries for the whole grid view which is memory intensive but faster
 */
template<class Scalar,
         class GV,
         bool enableCaching = true,
         class Traits = PQ2DefaultGridGeometryTraits<GV>>
class PQ2FVGridGeometry
: public BaseGridGeometry<GV, Traits>
{
    using ThisType = PQ2FVGridGeometry<Scalar, GV, enableCaching, Traits>;
    using ParentType = BaseGridGeometry<GV, Traits>;
    using GridIndexType = typename IndexTraits<GV>::GridIndex;
    using LocalIndexType = typename IndexTraits<GV>::LocalIndex;

    using Element = typename GV::template Codim<0>::Entity;
    using CoordScalar = typename GV::ctype;
    static const int dim = GV::dimension;
    static const int dimWorld = GV::dimensionworld;

    static_assert(dim > 1, "Only implemented for dim > 1");

public:
    //! export the discretization method this geometry belongs to
    using DiscretizationMethod = DiscretizationMethods::PQ2;
    static constexpr DiscretizationMethod discMethod{};

    static constexpr bool enableHybridCVFE = Detail::enablesHybridCVFE<Traits>;
    static_assert(enableHybridCVFE, "Only hybrid scheme implemented for pq2");

    static constexpr std::size_t maxNumElementDofs = Traits::maxNumElementDofs;

    //! export basic grid geometry type for the alternative constructor
    using BasicGridGeometry = BasicGridGeometry_t<GV, Traits>;
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
    using FeCache = PQ2HierarchicalFECache<CoordScalar, Scalar, dim>;
    //! export the grid view type
    using GridView = GV;
    //! export whether the grid(geometry) supports periodicity
    using SupportsPeriodicity = typename PeriodicGridTraits<typename GV::Grid>::SupportsPeriodicity;
    //! the quadrature rule type for scvs
    using ScvQuadratureRule = typename Traits::ScvQuadratureRule;
    //! the quadrature rule type for scvfs
    using ScvfQuadratureRule = typename Traits::ScvfQuadratureRule;
    //! the quadrature rule type for elements
    using ElementQuadratureRule = typename Traits::ElementQuadratureRule;
    //! the quadrature rule type for intersections
    using IntersectionQuadratureRule = typename Traits::IntersectionQuadratureRule;

    //! Constructor with basic grid geometry used to share state with another grid geometry on the same grid view
    PQ2FVGridGeometry(std::shared_ptr<BasicGridGeometry> gg)
    : ParentType(std::move(gg))
    , dofMapper_(this->gridView(), Traits::layout())
    , cache_(*this)
    , periodicGridTraits_(this->gridView().grid())
    {
        update_();
    }

    //! Constructor
    PQ2FVGridGeometry(const GridView& gridView)
    : PQ2FVGridGeometry(std::make_shared<BasicGridGeometry>(gridView))
    {}

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

    //! If a vertex / d.o.f. is on the boundary
    bool dofOnBoundary(GridIndexType dofIdx) const
    { return boundaryDofIndices_[dofIdx]; }

    //! If a vertex / d.o.f. is on a periodic boundary
    bool dofOnPeriodicBoundary(GridIndexType dofIdx) const
    { return periodicDofMap_.count(dofIdx); }

    //! The index of the vertex / d.o.f. on the other side of the periodic boundary
    GridIndexType periodicallyMappedDof(GridIndexType dofIdx) const
    { return periodicDofMap_.at(dofIdx); }

    //! Returns the map between dofs across periodic boundaries
    const std::unordered_map<GridIndexType, GridIndexType>& periodicDofMap() const
    { return periodicDofMap_; }

    //! local view of this object (constructed with the internal cache)
    friend inline LocalView localView(const PQ2FVGridGeometry& gg)
    { return { gg.cache_ }; }

private:

    class PQ2GridGeometryCache
    {
        friend class PQ2FVGridGeometry;
    public:
        //! export the geometry helper type
        using GeometryHelper = Detail::PQ2GeometryHelper_t<GV, Traits>;

        explicit PQ2GridGeometryCache(const PQ2FVGridGeometry& gg)
        : gridGeometry_(&gg)
        {}

        const PQ2FVGridGeometry& gridGeometry() const
        { return *gridGeometry_; }

        //! Get the global sub control volume indices of an element
        const std::vector<SubControlVolume>& scvs(GridIndexType eIdx) const
        { return scvs_[eIdx]; }

        //! Get the global sub control volume face indices of an element
        const std::vector<SubControlVolumeFace>& scvfs(GridIndexType eIdx) const
        { return scvfs_[eIdx]; }

        //! Returns whether one of the geometry's scvfs lies on a boundary
        bool hasBoundaryScvf(GridIndexType eIdx) const
        { return hasBoundaryScvf_[eIdx]; }

        //! Local mappings necessary to construct geometries of scvfs
        const std::vector<std::array<LocalIndexType, 2>>& scvfBoundaryGeometryKeys(GridIndexType eIdx) const
        { return scvfBoundaryGeometryKeys_.at(eIdx); }

    private:
        void clear_()
        {
            scvs_.clear();
            scvfs_.clear();
            hasBoundaryScvf_.clear();
            scvfBoundaryGeometryKeys_.clear();
        }

        std::vector<std::vector<SubControlVolume>> scvs_;
        std::vector<std::vector<SubControlVolumeFace>> scvfs_;
        std::vector<bool> hasBoundaryScvf_;
        std::unordered_map<GridIndexType, std::vector<std::array<LocalIndexType, 2>>> scvfBoundaryGeometryKeys_;

        const PQ2FVGridGeometry* gridGeometry_;
    };

public:
    //! the cache type (only the caching implementation has this)
    //! this alias should only be used by the local view implementation
    using Cache = PQ2GridGeometryCache;

private:
    using GeometryHelper = typename Cache::GeometryHelper;

    void update_()
    {
        cache_.clear_();
        dofMapper_.update(this->gridView());

        auto numElements = this->gridView().size(0);
        cache_.scvs_.resize(numElements);
        cache_.scvfs_.resize(numElements);
        cache_.hasBoundaryScvf_.resize(numElements, false);

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

            const auto& localCoefficients = this->feCache().get(element.type()).localCoefficients();

            // instantiate the geometry helper
            GeometryHelper geometryHelper(elementGeometry);

            numScv_ += geometryHelper.numScv();
            // construct the sub control volumes
            cache_.scvs_[eIdx].resize(geometryHelper.numScv());

            for (LocalIndexType keyIdx = 0; keyIdx < localCoefficients.size(); ++keyIdx)
            {
                const auto& localKey = localCoefficients.localKey(keyIdx);
                // If the dof is a vertex, we construct scvs
                if(localKey.codim() == dim)
                {
                    const auto localIdx = localKey.subEntity();
                    // With the new localIdx, scvs can be constructed as for the Box method
                    auto corners = geometryHelper.getScvCorners(localIdx);
                    cache_.scvs_[eIdx][localIdx] = SubControlVolume(
                        geometryHelper.scvVolume(localIdx, corners),
                        geometryHelper.dofPosition(localKey),
                        Dumux::center(corners),
                        localIdx,
                        keyIdx,
                        eIdx,
                        geometryHelper.dofIndex(this->dofMapper(), element, localKey),
                        false
                    );
                }
            }

            // construct the sub control volume faces, this is the same as for the Box method
            const auto numInteriorScvfs = GeometryHelper::numInteriorScvf(elementGeometry.type());
            numScvf_ += numInteriorScvfs;
            cache_.scvfs_[eIdx].resize(numInteriorScvfs);
            LocalIndexType scvfLocalIdx = 0;
            for (; scvfLocalIdx < numInteriorScvfs; ++scvfLocalIdx)
            {
                const auto scvPair = geometryHelper.getScvPairForScvf(scvfLocalIdx);
                const auto corners = geometryHelper.getScvfCorners(scvfLocalIdx);
                const auto area = Dumux::convexPolytopeVolume<dim-1>(
                    geometryHelper.getInteriorScvfGeometryType(scvfLocalIdx),
                    [&](unsigned int i){ return corners[i]; }
                );

                cache_.scvfs_[eIdx][scvfLocalIdx] = SubControlVolumeFace(
                    Dumux::center(corners),
                    area,
                    geometryHelper.normal(corners, scvPair),
                    std::move(scvPair),
                    scvfLocalIdx,
                    geometryHelper.isOverlappingScvf(scvfLocalIdx)
                );
            }

            // construct the sub control volume faces on the domain boundary
            for (const auto& intersection : intersections(this->gridView(), element))
            {
                if (intersection.boundary() && !intersection.neighbor())
                {
                    cache_.hasBoundaryScvf_[eIdx] = true;

                    const auto localFacetIndex = intersection.indexInInside();
                    const auto numBoundaryScvf = GeometryHelper::numBoundaryScvf(elementGeometry.type(), localFacetIndex);
                    numScvf_ += numBoundaryScvf;
                    numBoundaryScvf_ += numBoundaryScvf;

                    for (unsigned int isScvfLocalIdx = 0; isScvfLocalIdx < numBoundaryScvf; ++isScvfLocalIdx)
                    {
                        // find the scvs this scvf is belonging to
                        const auto scvPair = geometryHelper.getScvPairForBoundaryScvf(localFacetIndex, isScvfLocalIdx);
                        const auto corners = geometryHelper.getBoundaryScvfCorners(localFacetIndex, isScvfLocalIdx);
                        const auto area = Dumux::convexPolytopeVolume<dim-1>(
                            geometryHelper.getBoundaryScvfGeometryType(isScvfLocalIdx),
                            [&](unsigned int i){ return corners[i]; }
                        );
                        cache_.scvfs_[eIdx].emplace_back(
                            Dumux::center(corners),
                            area,
                            intersection.centerUnitOuterNormal(),
                            std::move(scvPair),
                            scvfLocalIdx,
                            typename SubControlVolumeFace::Traits::BoundaryFlag{ intersection },
                            geometryHelper.isOverlappingBoundaryScvf(localFacetIndex)
                        );

                        // store look-up map to construct boundary scvf geometries
                        cache_.scvfBoundaryGeometryKeys_[eIdx].emplace_back(std::array<LocalIndexType, 2>{{
                            static_cast<LocalIndexType>(localFacetIndex),
                            static_cast<LocalIndexType>(isScvfLocalIdx)
                        }});

                        // increment local counter
                        scvfLocalIdx++;
                    }

                    for (LocalIndexType keyIdx = 0; keyIdx < localCoefficients.size(); ++keyIdx)
                    {
                        if(GeometryHelper::localDofOnIntersection(elementGeometry.type(), intersection.indexInInside(), localCoefficients.localKey(keyIdx)))
                        {
                            const auto dofIdxGlobal = GeometryHelper::dofIndex(this->dofMapper(), element, localCoefficients.localKey(keyIdx));
                            boundaryDofIndices_[dofIdxGlobal] = true;
                        }
                    }
                }

                // inform the grid geometry if we have periodic boundaries
                else if (periodicGridTraits_.isPeriodic(intersection))
                {
                    this->setPeriodic();

                    // find the mapped periodic vertex of all vertices on periodic boundaries
                    const auto eps = 1e-7*(elementGeometry.corner(1) - elementGeometry.corner(0)).two_norm();
                    for (int localDofIdx = 0; localDofIdx < localCoefficients.size(); ++localDofIdx)
                    {
                        if(!GeometryHelper::localDofOnIntersection(elementGeometry.type(), intersection.indexInInside(), localCoefficients.localKey(localDofIdx)))
                            continue;

                        const auto dofIdxGlobal = GeometryHelper::dofIndex(this->dofMapper(), element, localCoefficients.localKey(localDofIdx));
                        const auto dofPos = geometryHelper.dofPosition(localCoefficients.localKey(localDofIdx));

                        const auto& outside = intersection.outside();
                        const auto outsideGeometry = outside.geometry();
                        const auto& localCoefficientsOut = this->feCache().get(outsideGeometry.type()).localCoefficients();
                        for (const auto& isOutside : intersections(this->gridView(), outside))
                        {
                            // only check periodic vertices of the periodic neighbor
                            if (isOutside.boundary() && isOutside.neighbor())
                            {
                                for (int localDofIdxOut = 0; localDofIdxOut < localCoefficientsOut.size(); ++localDofIdxOut)
                                {
                                    const auto& localKeyOut = localCoefficientsOut.localKey(localDofIdxOut);
                                    if(!GeometryHelper::localDofOnIntersection(outsideGeometry.type(), isOutside.indexInInside(), localKeyOut))
                                        continue;

                                    const auto dofIdxGlobalOut = GeometryHelper::dofIndex(this->dofMapper(), outside, localKeyOut);
                                    const auto dofPosOutside =  GeometryHelper::dofPosition(outsideGeometry, localKeyOut);
                                    const auto shift = std::abs((this->bBoxMax()-this->bBoxMin())*intersection.centerUnitOuterNormal());
                                    if (std::abs((dofPosOutside-dofPos).two_norm() - shift) < eps)
                                        periodicDofMap_[dofIdxGlobal] = dofIdxGlobalOut;
                                }
                            }
                        }
                    }
                }
            }
        }

        // error check: periodic boundaries currently don't work for pq2 in parallel
        if (this->isPeriodic() && this->gridView().comm().size() > 1)
            DUNE_THROW(Dune::NotImplemented, "Periodic boundaries for pq2 method for parallel simulations!");
    }

    DofMapper dofMapper_;

    const FeCache feCache_;

    std::size_t numScv_;
    std::size_t numScvf_;
    std::size_t numBoundaryScvf_;

    // dofs on the boundary
    std::vector<bool> boundaryDofIndices_;

    // a map for periodic boundary dofs
    std::unordered_map<GridIndexType, GridIndexType> periodicDofMap_;

    Cache cache_;

    PeriodicGridTraits<typename GridView::Grid> periodicGridTraits_;
};

} // end namespace Dumux

#endif
