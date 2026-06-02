// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup PQ3Discretization
 * \brief Base class for the finite volume geometry vector for the pq3 method
 *        This builds up the sub control volumes and sub control volume faces
 *        for each element of the grid partition.
 */
#ifndef DUMUX_DISCRETIZATION_PQ3_GRID_GEOMETRY_HH
#define DUMUX_DISCRETIZATION_PQ3_GRID_GEOMETRY_HH

#include <array>
#include <cstddef>
#include <memory>
#include <vector>
#include <utility>
#include <unordered_map>
#include <type_traits>
#include <span>

#include <dune/grid/common/mcmgmapper.hh>
#include <dune/geometry/type.hh>
#include <dune/localfunctions/lagrange/lagrangelfecache.hh>

#include <dumux/discretization/method.hh>
#include <dumux/common/indextraits.hh>
#include <dumux/common/defaultmappertraits.hh>
#include <dumux/geometry/center.hh>
#include <dumux/geometry/volume.hh>
#include <dumux/discretization/basegridgeometry.hh>
#include <dumux/discretization/pq1bubble/fvgridgeometry.hh>
#include <dumux/discretization/pq3/geometryhelper.hh>
// Reuse PQ2 SCV/SCVF/FVElementGeometry (all templated, work for any GG)
#include <dumux/discretization/pq2/subcontrolvolume.hh>
#include <dumux/discretization/pq2/subcontrolvolumeface.hh>
#include <dumux/discretization/pq2/fvelementgeometry.hh>
#include <dumux/discretization/boundaryface.hh>
#include <dumux/discretization/extrusion.hh>
#include <dumux/io/grid/periodicgridtraits.hh>

namespace Dumux {

// Reuse PQ2 SCV/SCVF types under PQ3 aliases
template<class GV>
using PQ3SubControlVolume = PQ2SubControlVolume<GV>;

template<class GV>
using PQ3SubControlVolumeFace = PQ2SubControlVolumeFace<GV>;

// FVElementGeometry reused from PQ2 (it's fully templated on the grid geometry)
template<class GG, bool enableCache>
using PQ3FVElementGeometry = PQ2FVElementGeometry<GG, enableCache>;

namespace Detail {

template<class GV, class T>
using PQ3GeometryHelper_t = Dune::Std::detected_or_t<
    std::conditional_t<enablesHybridCVFE<T>,
        Dumux::HybridPQ3GeometryHelper<GV, typename T::SubControlVolume, typename T::SubControlVolumeFace>,
        void
    >,
    SpecifiesGeometryHelper,
    T
>;

} // end namespace Detail

/*!
 * \ingroup PQ3Discretization
 * \brief Mapper traits for PQ3: vertices get 1 DOF, edges get 2 DOFs,
 *        quad faces/elements get 4 DOFs, simplex faces/elements get 1 DOF,
 *        hex elements get 8 DOFs.
 */
template<class GV>
struct PQ3MapperTraits : public DefaultMapperTraits<GV>
{
    using DofMapper = Dune::MultipleCodimMultipleGeomTypeMapper<GV>;

    static Dune::MCMGLayout layout()
    {
        return [](Dune::GeometryType gt, int /*dimgrid*/) -> std::size_t {
            if (gt.dim() == 0) return 1;  // vertices: 1 DOF each
            if (gt.dim() == 1) return 2;  // edges: 2 interior DOFs (order 3)
            if (gt.dim() == 2) {
                if (gt == Dune::GeometryTypes::cube(2))    return 4;  // quad: 4 interior DOFs
                if (gt == Dune::GeometryTypes::simplex(2)) return 1;  // triangle: 1 interior DOF
            }
            if (gt.dim() == 3) {
                if (gt == Dune::GeometryTypes::cube(3))    return 8;  // hex: 8 interior DOFs
                // simplex (tet): 0 interior DOFs for P3
            }
            return 0;
        };
    }
};

/*!
 * \ingroup PQ3Discretization
 * \brief Quadrature rule traits for PQ3 discretization
 */
template<class GridView,
         class ScvRule = Dumux::QuadratureRules::DuneQuadrature<3>,
         class ScvfRule = Dumux::QuadratureRules::DuneQuadrature<3>,
         class ElementRule = Dumux::QuadratureRules::DuneQuadrature<6>,
         class IntersectionRule = Dumux::QuadratureRules::DuneQuadrature<6>,
         class BoundaryFaceRule = Dumux::QuadratureRules::DuneQuadrature<6>>
using PQ3QuadratureTraits = CVFE::DefaultQuadratureTraits<GridView, ScvRule, ScvfRule, ElementRule, IntersectionRule, BoundaryFaceRule>;

/*!
 * \ingroup PQ3Discretization
 * \brief Default traits for the pq3 finite volume grid geometry
 */
template<class GridView,
         class MapperTraits = PQ3MapperTraits<GridView>,
         class QuadratureTraits = PQ3QuadratureTraits<GridView>>
struct PQ3DefaultGridGeometryTraits
: public MapperTraits, public QuadratureTraits
{
    using SubControlVolume = PQ3SubControlVolume<GridView>;
    using SubControlVolumeFace = PQ3SubControlVolumeFace<GridView>;

    template<class GridGeometry, bool enableCache>
    using LocalView = PQ3FVElementGeometry<GridGeometry, enableCache>;

    using EnableHybridCVFE = std::true_type;

    // Maximum number of element-local DOFs (for cubes, order 3)
    static constexpr std::size_t maxNumElementDofs = []()
        {
            if constexpr (GridView::dimension == 1)
                return 4;   // Q3 interval: 2 vertices + 2 edge interior
            else if constexpr (GridView::dimension == 2)
                return 16;  // Q3 quad: 4 vertices + 8 edge + 4 interior
            else if constexpr (GridView::dimension == 3)
                return 64;  // Q3 hex: 8 vertices + 24 edge + 24 face + 8 interior
        }();
};

/*!
 * \ingroup PQ3Discretization
 * \brief Finite volume geometry for the pq3 hybrid CVFE scheme (order-3 Lagrange elements).
 *
 * Control volumes are defined only for vertex DOFs (box dual mesh).
 * Edge, face, and element interior DOFs are non-CV ("hybrid") DOFs.
 */
template<class Scalar,
         class GV,
         bool enableCaching = true,
         class Traits = PQ3DefaultGridGeometryTraits<GV>>
class PQ3FVGridGeometry
: public BaseGridGeometry<GV, Traits>
{
    using ThisType = PQ3FVGridGeometry<Scalar, GV, enableCaching, Traits>;
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
    using DiscretizationMethod = DiscretizationMethods::PQ3;
    static constexpr DiscretizationMethod discMethod{};

    static constexpr bool enableHybridCVFE = Detail::enablesHybridCVFE<Traits>;
    static_assert(enableHybridCVFE, "Only hybrid scheme implemented for pq3");

    static constexpr std::size_t maxNumElementDofs = Traits::maxNumElementDofs;

    //! export basic grid geometry type
    using BasicGridGeometry = BasicGridGeometry_t<GV, Traits>;
    //! export the type of the fv element geometry (the local view type)
    using LocalView = typename Traits::template LocalView<ThisType, true>;
    //! export the type of sub control volume
    using SubControlVolume = typename Traits::SubControlVolume;
    //! export the type of sub control volume face
    using SubControlVolumeFace = typename Traits::SubControlVolumeFace;
    //! export the boundary face type
    using BoundaryFace = Experimental::BoundaryFace<GV>;
    //! export the type of extrusion
    using Extrusion = Extrusion_t<Traits>;
    //! export dof mapper type
    using DofMapper = typename Traits::DofMapper;
    //! export the finite element cache type (order 3)
    using FeCache = Dune::LagrangeLocalFiniteElementCache<CoordScalar, Scalar, dim, 3>;
    //! export the grid view type
    using GridView = GV;
    //! export whether the grid supports periodicity
    using SupportsPeriodicity = typename PeriodicGridTraits<typename GV::Grid>::SupportsPeriodicity;
    //! quadrature rule types
    using ScvQuadratureRule = typename Traits::ScvQuadratureRule;
    using ScvfQuadratureRule = typename Traits::ScvfQuadratureRule;
    using ElementQuadratureRule = typename Traits::ElementQuadratureRule;
    using IntersectionQuadratureRule = typename Traits::IntersectionQuadratureRule;
    using BoundaryFaceQuadratureRule = typename Traits::BoundaryFaceQuadratureRule;

    //! Constructor with shared basic grid geometry
    PQ3FVGridGeometry(std::shared_ptr<BasicGridGeometry> gg)
    : ParentType(std::move(gg))
    , dofMapper_(this->gridView(), Traits::layout())
    , cache_(*this)
    , periodicGridTraits_(this->gridView().grid())
    {
        update_();
    }

    //! Constructor
    PQ3FVGridGeometry(const GridView& gridView)
    : PQ3FVGridGeometry(std::make_shared<BasicGridGeometry>(gridView))
    {}

    const DofMapper& dofMapper() const
    { return dofMapper_; }

    std::size_t numScv() const
    { return numScv_; }

    std::size_t numScvf() const
    { return numScvf_; }

    std::size_t numBoundaryScvf() const
    { return numBoundaryScvf_; }

    std::size_t numDofs() const
    { return this->dofMapper().size(); }

    void update(const GridView& gridView)
    {
        ParentType::update(gridView);
        update_();
    }

    void update(GridView&& gridView)
    {
        ParentType::update(std::move(gridView));
        update_();
    }

    const FeCache& feCache() const
    { return feCache_; }

    //! Flip-aware global DOF index for a local key.
    //! Used by CVFEElementSolution::update to correctly extract the solution
    //! for simplex edge DOFs where Dune::LagrangeLocalFiniteElementCache returns
    //! unflipped keys but the globally-consistent ordering requires a flip.
    template<class LocalKey>
    auto dofIndex(const Element& element, const LocalKey& localKey) const
    { return Detail::PQ3GeometryHelper_t<GV, Traits>::DofHelper::dofIndex(dofMapper_, element, localKey, this->gridView().grid().globalIdSet()); }

    bool dofOnBoundary(GridIndexType dofIdx) const
    { return boundaryDofIndices_[dofIdx]; }

    bool dofOnPeriodicBoundary(GridIndexType dofIdx) const
    { return periodicDofMap_.count(dofIdx); }

    GridIndexType periodicallyMappedDof(GridIndexType dofIdx) const
    { return periodicDofMap_.at(dofIdx); }

    const std::unordered_map<GridIndexType, GridIndexType>& periodicDofMap() const
    { return periodicDofMap_; }

    friend inline LocalView localView(const PQ3FVGridGeometry& gg)
    { return { gg.cache_ }; }

private:

    class PQ3GridGeometryCache
    {
        friend class PQ3FVGridGeometry;
    public:
        using GeometryHelper = Detail::PQ3GeometryHelper_t<GV, Traits>;
        using DofHelper = GeometryHelper::DofHelper;

        explicit PQ3GridGeometryCache(const PQ3FVGridGeometry& gg)
        : gridGeometry_(&gg)
        {}

        const PQ3FVGridGeometry& gridGeometry() const
        { return *gridGeometry_; }

        const std::vector<SubControlVolume>& scvs(GridIndexType eIdx) const
        { return scvs_[eIdx]; }

        const std::vector<SubControlVolumeFace>& scvfs(GridIndexType eIdx) const
        { return scvfs_[eIdx]; }

        bool hasBoundaryScvf(GridIndexType eIdx) const
        { return hasBoundaryScvf_[eIdx]; }

        const std::vector<std::array<LocalIndexType, 2>>& scvfBoundaryGeometryKeys(GridIndexType eIdx) const
        { return scvfBoundaryGeometryKeys_.at(eIdx); }

        auto boundaryFaces(GridIndexType eIdx) const -> std::span<const BoundaryFace>
        {
            if (auto it = boundaryFaces_.find(eIdx); it != boundaryFaces_.end())
                return {it->second};
            return {};
        }

        const auto& boundaryFaceScvfRanges(GridIndexType eIdx) const
        { return boundaryFaceScvfRanges_.at(eIdx); }

    private:
        void clear_()
        {
            scvs_.clear();
            scvfs_.clear();
            hasBoundaryScvf_.clear();
            scvfBoundaryGeometryKeys_.clear();
            boundaryFaces_.clear();
            boundaryFaceScvfRanges_.clear();
        }

        std::vector<std::vector<SubControlVolume>> scvs_;
        std::vector<std::vector<SubControlVolumeFace>> scvfs_;
        std::vector<bool> hasBoundaryScvf_;
        std::unordered_map<GridIndexType, std::vector<std::array<LocalIndexType, 2>>> scvfBoundaryGeometryKeys_;
        std::unordered_map<GridIndexType, Dune::ReservedVector<typename PQ3FVGridGeometry::BoundaryFace, 2*dim>> boundaryFaces_;
        std::unordered_map<GridIndexType, Dune::ReservedVector<std::array<LocalIndexType, 2>, 2*dim>> boundaryFaceScvfRanges_;

        const PQ3FVGridGeometry* gridGeometry_;
    };

public:
    using Cache = PQ3GridGeometryCache;

private:
    using GeometryHelper = typename Cache::GeometryHelper;
    using DofHelper = typename GeometryHelper::DofHelper;

    void update_()
    {
        cache_.clear_();
        dofMapper_.update(this->gridView());

        const auto& gidSet = this->gridView().grid().globalIdSet();

        auto numElements = this->gridView().size(0);
        cache_.scvs_.resize(numElements);
        cache_.scvfs_.resize(numElements);
        cache_.hasBoundaryScvf_.resize(numElements, false);

        boundaryDofIndices_.assign(numDofs(), false);

        numScv_ = 0;
        numScvf_ = 0;
        numBoundaryScvf_ = 0;

        for (const auto& element : elements(this->gridView()))
        {
            auto eIdx = this->elementMapper().index(element);
            auto elementGeometry = element.geometry();
            const auto& localCoefficients = this->feCache().get(element.type()).localCoefficients();

            GeometryHelper geometryHelper(elementGeometry);

            numScv_ += geometryHelper.numScv();
            cache_.scvs_[eIdx].resize(geometryHelper.numScv());

            for (LocalIndexType keyIdx = 0; keyIdx < localCoefficients.size(); ++keyIdx)
            {
                const auto& localKey = localCoefficients.localKey(keyIdx);
                if (localKey.codim() == dim)
                {
                    const auto localIdx = localKey.subEntity();
                    auto corners = geometryHelper.getScvCorners(localIdx);
                    cache_.scvs_[eIdx][localIdx] = SubControlVolume(
                        geometryHelper.scvVolume(localIdx, corners),
                        DofHelper::dofPosition(elementGeometry, localKey),
                        Dumux::center(corners),
                        localIdx,
                        keyIdx,
                        eIdx,
                        DofHelper::dofIndex(this->dofMapper(), element, localKey, gidSet),
                        false
                    );
                }
            }

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

            LocalIndexType numBoundaryFaces = 0;
            for (const auto& intersection : intersections(this->gridView(), element))
            {
                if (intersection.boundary() && !intersection.neighbor())
                {
                    cache_.hasBoundaryScvf_[eIdx] = true;

                    const auto isGeometry = intersection.geometry();
                    cache_.boundaryFaces_[eIdx].push_back(BoundaryFace{
                        isGeometry.center(),
                        isGeometry.volume(),
                        intersection.centerUnitOuterNormal(),
                        numBoundaryFaces++,
                        static_cast<LocalIndexType>(intersection.indexInInside()),
                        typename BoundaryFace::Traits::BoundaryFlag{intersection}
                    });

                    const auto localFacetIndex = intersection.indexInInside();
                    const auto numBoundaryScvf = GeometryHelper::numBoundaryScvf(elementGeometry.type(), localFacetIndex);
                    numScvf_ += numBoundaryScvf;
                    numBoundaryScvf_ += numBoundaryScvf;

                    cache_.boundaryFaceScvfRanges_[eIdx].push_back(std::array<LocalIndexType, 2>{{
                        scvfLocalIdx,
                        static_cast<LocalIndexType>(numBoundaryScvf)
                    }});

                    for (unsigned int isScvfLocalIdx = 0; isScvfLocalIdx < numBoundaryScvf; ++isScvfLocalIdx)
                    {
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

                        cache_.scvfBoundaryGeometryKeys_[eIdx].emplace_back(std::array<LocalIndexType, 2>{{
                            static_cast<LocalIndexType>(localFacetIndex),
                            static_cast<LocalIndexType>(isScvfLocalIdx)
                        }});

                        scvfLocalIdx++;
                    }

                    for (LocalIndexType keyIdx = 0; keyIdx < localCoefficients.size(); ++keyIdx)
                    {
                        if (DofHelper::localDofOnIntersection(elementGeometry.type(), intersection.indexInInside(), localCoefficients.localKey(keyIdx)))
                        {
                            const auto dofIdxGlobal = DofHelper::dofIndex(this->dofMapper(), element, localCoefficients.localKey(keyIdx), gidSet);
                            boundaryDofIndices_[dofIdxGlobal] = true;
                        }
                    }
                }

                else if (periodicGridTraits_.isPeriodic(intersection))
                {
                    this->setPeriodic();

                    const auto eps = 1e-7*(elementGeometry.corner(1) - elementGeometry.corner(0)).two_norm();
                    for (int localDofIdx = 0; localDofIdx < localCoefficients.size(); ++localDofIdx)
                    {
                        if (!DofHelper::localDofOnIntersection(elementGeometry.type(), intersection.indexInInside(), localCoefficients.localKey(localDofIdx)))
                            continue;

                        const auto dofIdxGlobal = DofHelper::dofIndex(this->dofMapper(), element, localCoefficients.localKey(localDofIdx), gidSet);
                        const auto dofPos = DofHelper::dofPosition(elementGeometry, localCoefficients.localKey(localDofIdx));

                        const auto& outside = intersection.outside();
                        const auto outsideGeometry = outside.geometry();
                        const auto& localCoefficientsOut = this->feCache().get(outsideGeometry.type()).localCoefficients();
                        for (const auto& isOutside : intersections(this->gridView(), outside))
                        {
                            if (isOutside.boundary() && isOutside.neighbor())
                            {
                                for (int localDofIdxOut = 0; localDofIdxOut < localCoefficientsOut.size(); ++localDofIdxOut)
                                {
                                    const auto& localKeyOut = localCoefficientsOut.localKey(localDofIdxOut);
                                    if (!DofHelper::localDofOnIntersection(outsideGeometry.type(), isOutside.indexInInside(), localKeyOut))
                                        continue;

                                    const auto dofIdxGlobalOut = DofHelper::dofIndex(this->dofMapper(), outside, localKeyOut, gidSet);
                                    const auto dofPosOutside = DofHelper::dofPosition(outsideGeometry, localKeyOut);
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

        if (this->isPeriodic() && this->gridView().comm().size() > 1)
            DUNE_THROW(Dune::NotImplemented, "Periodic boundaries for pq3 method for parallel simulations!");
    }

    DofMapper dofMapper_;
    const FeCache feCache_;
    std::size_t numScv_;
    std::size_t numScvf_;
    std::size_t numBoundaryScvf_;
    std::vector<bool> boundaryDofIndices_;
    std::unordered_map<GridIndexType, GridIndexType> periodicDofMap_;
    Cache cache_;
    PeriodicGridTraits<typename GridView::Grid> periodicGridTraits_;
};

} // end namespace Dumux

#endif
