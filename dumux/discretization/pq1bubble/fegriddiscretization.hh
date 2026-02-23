// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup PQ1BubbleDiscretization
 * \brief Base class for the finite element grid discretization for the pq1bubble discretization
 */
#ifndef DUMUX_DISCRETIZATION_PQ1BUBBLE_FE_GRID_DISCRETIZATION_HH
#define DUMUX_DISCRETIZATION_PQ1BUBBLE_FE_GRID_DISCRETIZATION_HH

#include <span>
#include <utility>
#include <unordered_map>

#include <dumux/discretization/method.hh>

#include <dumux/common/indextraits.hh>
#include <dumux/common/defaultmappertraits.hh>

#include <dumux/discretization/basegridgeometry.hh>
#include <dumux/discretization/boundaryface.hh>
#include <dumux/discretization/pq1bubble/fvgridgeometry.hh>
#include <dumux/discretization/pq1bubble/geometryhelper.hh>
#include <dumux/discretization/pq1bubble/feelementdiscretization.hh>
#include <dumux/discretization/pq1bubble/pq1bubblefecache.hh>
#include <dumux/discretization/extrusion.hh>
#include <dumux/discretization/cvfe/quadraturerules.hh>

#include <dumux/io/grid/periodicgridtraits.hh>

namespace Dumux {

/*!
 * \ingroup PQ1BubbleDiscretization
 * \brief Quadrature rule traits for PQ1Bubble fe discretization
 */
template<class GridView,
         class ElementRule = Dumux::QuadratureRules::DuneQuadrature<4>,
         class IntersectionRule = Dumux::QuadratureRules::MidpointQuadrature,
         class BoundaryFaceRule = Dumux::QuadratureRules::MidpointQuadrature>
struct PQ1BubbleFEQuadratureTraits
{
    using ElementQuadratureRule = ElementRule;
    using IntersectionQuadratureRule = IntersectionRule;
    using BoundaryFaceQuadratureRule = BoundaryFaceRule;
};

/*!
 * \ingroup PQ1BubbleDiscretization
 * \brief The default traits for the pq1bubble finite element grid discretization
 * \tparam the grid view type
 */
template<class GridView, class MapperTraits = PQ1BubbleMapperTraits<GridView, 2>, class QuadratureTraits = PQ1BubbleFEQuadratureTraits<GridView>>
struct PQ1BubbleFEDefaultGridDiscretizationTraits
: public MapperTraits, public QuadratureTraits
{
    template<class GridDiscretization, bool enableCache>
    using LocalView = PQ1BubbleFEElementDiscretization<GridDiscretization, enableCache>;

    static constexpr std::size_t maxNumElementDofs = (1<<GridView::dimension)
                                                    + MapperTraits::numCubeBubbleDofs;
};

/*!
 * \ingroup PQ1BubbleDiscretization
 * \brief Base class for the finite element grid discretization for the pq1bubble discretization
 */
template<class Scalar,
         class GV,
         bool enableCaching = true,
         class Traits = PQ1BubbleFEDefaultGridDiscretizationTraits<GV>>
class PQ1BubbleFEGridDiscretization
: public BaseGridGeometry<GV, Traits>
{
    using ThisType = PQ1BubbleFEGridDiscretization<Scalar, GV, enableCaching, Traits>;
    using ParentType = BaseGridGeometry<GV, Traits>;
    using GridIndexType = typename IndexTraits<GV>::GridIndex;
    using LocalIndexType = typename IndexTraits<GV>::LocalIndex;
    using CoordScalar = typename GV::ctype;
    static const int dim = GV::dimension;

    static_assert(dim > 1, "Only implemented for dim > 1");

public:
    //! export the discretization method this geometry belongs to
    using DiscretizationMethod = DiscretizationMethods::PQ1Bubble;
    static constexpr DiscretizationMethod discMethod{};

    // ToDo: Rename
    static constexpr bool enableHybridCVFE = true;

    static constexpr std::size_t maxNumElementDofs = Traits::maxNumElementDofs;

    //! export the boundary face type
    using BoundaryFace = Experimental::BoundaryFace<GV>;

    //! export basic grid geometry type for the alternative constructor
    using BasicGridGeometry = BasicGridGeometry_t<GV, Traits>;
    //! export the type of the element discretization (the local view type)
    using LocalView = typename Traits::template LocalView<ThisType, true>;
    //! export the type of extrusion
    using Extrusion = Extrusion_t<Traits>;
    //! export dof mapper type
    using DofMapper = typename Traits::DofMapper;
    //! export the finite element cache type
    using FeCache = Dumux::PQ1BubbleFECache<CoordScalar, Scalar, dim, Traits::numCubeBubbleDofs>;
    //! export the grid view type
    using GridView = GV;
    //! export whether the grid supports periodicity
    using SupportsPeriodicity = typename PeriodicGridTraits<typename GV::Grid>::SupportsPeriodicity;
    //! the quadrature rule type for elements
    using ElementQuadratureRule = typename Traits::ElementQuadratureRule;
    //! the quadrature rule type for intersections
    using IntersectionQuadratureRule = typename Traits::IntersectionQuadratureRule;
    //! the quadrature rule type for boundary faces
    using BoundaryFaceQuadratureRule = typename Traits::BoundaryFaceQuadratureRule;

    //! Constructor with basic grid geometry used to share state with another grid discretization on the same grid view
    PQ1BubbleFEGridDiscretization(std::shared_ptr<BasicGridGeometry> gg)
    : ParentType(std::move(gg))
    , dofMapper_(this->gridView(), Traits::layout())
    , cache_(*this)
    , periodicGridTraits_(this->gridView().grid())
    {
        update_();
    }

    //! Constructor
    PQ1BubbleFEGridDiscretization(const GridView& gridView)
    : PQ1BubbleFEGridDiscretization(std::make_shared<BasicGridGeometry>(gridView))
    {}

    //! The dofMapper
    const DofMapper& dofMapper() const
    { return dofMapper_; }

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
    friend inline LocalView localView(const PQ1BubbleFEGridDiscretization& gg)
    { return { gg.cache_ }; }

private:
    class PQ1BubbleFEGridDiscretizationCache
    {
        friend class PQ1BubbleFEGridDiscretization;
    public:
        //! export the geometry helper type
        using GeometryHelper = FEPQ1BubbleGeometryHelper<GridView, Traits::numCubeBubbleDofs>;

        explicit PQ1BubbleFEGridDiscretizationCache(const PQ1BubbleFEGridDiscretization& gg)
        : gridDiscretization_(&gg)
        {}

        const PQ1BubbleFEGridDiscretization& gridDiscretization() const
        { return *gridDiscretization_; }

        //! Returns the boundary faces of an element
        auto boundaryFaces(GridIndexType eIdx) const -> std::span<const BoundaryFace>
        {
            if (auto it = boundaryFaces_.find(eIdx); it != boundaryFaces_.end())
                return {it->second};
            return {};
        }

    private:
        void clear_()
        { boundaryFaces_.clear(); }

        std::unordered_map<GridIndexType, Dune::ReservedVector<BoundaryFace, 2*dim>> boundaryFaces_;

        const PQ1BubbleFEGridDiscretization* gridDiscretization_;
    };

public:
    //! the cache type (only the caching implementation has this)
    //! this alias should only be used by the local view implementation
    using Cache = PQ1BubbleFEGridDiscretizationCache;
private:
    using GeometryHelper = typename Cache::GeometryHelper;

    void update_()
    {
        cache_.clear_();
        dofMapper_.update(this->gridView());

        boundaryDofIndices_.assign(numDofs(), false);

        for (const auto& element : elements(this->gridView()))
        {
            // get the element geometry
            auto elementGeometry = element.geometry();
            const auto refElement = referenceElement(elementGeometry);
            const auto& localCoefficients = this->feCache().get(element.type()).localCoefficients();
            const auto eIdx = this->elementMapper().index(element);

            // set boundary dofs and periodic dofs
            LocalIndexType numBoundaryFaces = 0;
            for (const auto& intersection : intersections(this->gridView(), element))
            {
                if (intersection.boundary() && !intersection.neighbor())
                {
                    const auto localFacetIndex = intersection.indexInInside();

                    // add one boundary face per boundary intersection
                    const auto isGeometry = intersection.geometry();
                    cache_.boundaryFaces_[eIdx].push_back(BoundaryFace{
                        isGeometry.center(),
                        isGeometry.volume(),
                        intersection.centerUnitOuterNormal(),
                        numBoundaryFaces++,
                        static_cast<LocalIndexType>(intersection.indexInInside()),
                        typename BoundaryFace::Traits::BoundaryFlag{intersection}
                    });

                    const auto numDofsIntersection = GeometryHelper::numLocalDofsIntersection(elementGeometry.type(), localFacetIndex);
                    // TODO also move this to helper class
                    // add all dofs on the intersection to the set of boundary dofs
                    for (int ilocalDofIdx = 0; ilocalDofIdx < numDofsIntersection; ++ilocalDofIdx)
                    {
                        auto localDofIdx = GeometryHelper::localDofIndexIntersection(elementGeometry.type(), localFacetIndex, ilocalDofIdx);
                        const auto vIdxGlobal = GeometryHelper::dofIndex(this->dofMapper(), element, localCoefficients.localKey(localDofIdx));
                        boundaryDofIndices_[vIdxGlobal] = true;
                    }
                }

                // inform the grid discretization if we have periodic boundaries
                else if (periodicGridTraits_.isPeriodic(intersection))
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
                            if (periodicGridTraits_.isPeriodic(isOutside))
                            {
                                const auto fIdxOutside = isOutside.indexInInside();
                                const auto numFaceVertsOutside = refElement.size(fIdxOutside, 1, dim);
                                for (int localVIdxOutside = 0; localVIdxOutside < numFaceVertsOutside; ++localVIdxOutside)
                                {
                                    const auto vIdxOutside = refElement.subEntity(fIdxOutside, 1, localVIdxOutside, dim);
                                    const auto vPosOutside = outsideGeometry.corner(vIdxOutside);
                                    const auto shift = std::abs((this->bBoxMax()-this->bBoxMin())*intersection.centerUnitOuterNormal());
                                    if (std::abs((vPosOutside-vPos).two_norm() - shift) < eps)
                                        periodicDofMap_[vIdxGlobal] = this->dofMapper().subIndex(outside, vIdxOutside, dim);
                                }
                            }
                        }
                    }
                }
            }
        }

        // error check: periodic boundaries currently don't work for pq1bubble in parallel
        if (this->isPeriodic() && this->gridView().comm().size() > 1)
            DUNE_THROW(Dune::NotImplemented, "Periodic boundaries for pq1bubble method for parallel simulations!");
    }

    DofMapper dofMapper_;

    const FeCache feCache_;

    // vertices on the boundary
    std::vector<bool> boundaryDofIndices_;

    // a map for periodic boundary vertices
    std::unordered_map<GridIndexType, GridIndexType> periodicDofMap_;

    Cache cache_;

    PeriodicGridTraits<typename GridView::Grid> periodicGridTraits_;
};

} // end namespace Dumux

#endif
