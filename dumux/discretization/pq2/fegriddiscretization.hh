// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup PQ2Discretization
 * \brief Base class for the finite element grid discretization for the pq2 discretization
 */
#ifndef DUMUX_DISCRETIZATION_PQ2_FE_GRID_DISCRETIZATION_HH
#define DUMUX_DISCRETIZATION_PQ2_FE_GRID_DISCRETIZATION_HH

#include <span>
#include <utility>
#include <unordered_map>

#include <dune/common/reservedvector.hh>
#include <dune/localfunctions/lagrange/lagrangelfecache.hh>

#include <dumux/discretization/method.hh>

#include <dumux/common/indextraits.hh>
#include <dumux/common/defaultmappertraits.hh>

#include <dumux/discretization/basegridgeometry.hh>
#include <dumux/discretization/boundaryface.hh>
#include <dumux/discretization/pq2/fvgridgeometry.hh>
#include <dumux/discretization/pq2/geometryhelper.hh>
#include <dumux/discretization/pq2/feelementdiscretization.hh>
#include <dumux/discretization/extrusion.hh>
#include <dumux/discretization/cvfe/quadraturerules.hh>

#include <dumux/io/grid/periodicgridtraits.hh>

namespace Dumux {

/*!
 * \ingroup PQ2Discretization
 * \brief Quadrature rule traits for PQ2 fe discretization
 */
template<class GridView,
         class ElementRule = Dumux::QuadratureRules::DuneQuadrature<4>,
         class IntersectionRule = Dumux::QuadratureRules::DuneQuadrature<4>,
         class BoundaryFaceRule = Dumux::QuadratureRules::DuneQuadrature<4>>
struct PQ2FEQuadratureTraits
{
    using ElementQuadratureRule = ElementRule;
    using IntersectionQuadratureRule = IntersectionRule;
    using BoundaryFaceQuadratureRule = BoundaryFaceRule;
};

/*!
 * \ingroup PQ2Discretization
 * \brief The default traits for the pq2 finite element grid discretization
 * \tparam the grid view type
 */
template<class GridView,
         class MapperTraits = PQ2MapperTraits<GridView>,
         class QuadratureTraits = PQ2FEQuadratureTraits<GridView>>
struct PQ2FEDefaultGridDiscretizationTraits
: public MapperTraits, public QuadratureTraits
{
    template<class GridDiscretization, bool enableCache>
    using LocalView = PQ2FEElementDiscretization<GridDiscretization, enableCache>;

    // The maximum values correspond to cubes
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
 * \brief Base class for the finite element grid discretization for the pq2 discretization
 */
template<class Scalar,
         class GV,
         bool enableCaching = true,
         class Traits = PQ2FEDefaultGridDiscretizationTraits<GV>>
class PQ2FEGridDiscretization
: public BaseGridGeometry<GV, Traits>
{
    using ThisType = PQ2FEGridDiscretization<Scalar, GV, enableCaching, Traits>;
    using ParentType = BaseGridGeometry<GV, Traits>;
    using GridIndexType = typename IndexTraits<GV>::GridIndex;
    using LocalIndexType = typename IndexTraits<GV>::LocalIndex;
    using CoordScalar = typename GV::ctype;
    static const int dim = GV::dimension;

    static_assert(dim > 1, "Only implemented for dim > 1");

public:
    //! export the discretization method this geometry belongs to
    using DiscretizationMethod = DiscretizationMethods::PQ2;
    static constexpr DiscretizationMethod discMethod{};

    // FE discretization: treat all dofs as FE dofs (no CV split)
    static constexpr bool enableHybridCVFE = false;

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
    using FeCache = Dune::LagrangeLocalFiniteElementCache<CoordScalar, Scalar, dim, 2>;
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
    PQ2FEGridDiscretization(std::shared_ptr<BasicGridGeometry> gg)
    : ParentType(std::move(gg))
    , dofMapper_(this->gridView(), Traits::layout())
    , cache_(*this)
    , periodicGridTraits_(this->gridView().grid())
    {
        update_();
    }

    //! Constructor
    PQ2FEGridDiscretization(const GridView& gridView)
    : PQ2FEGridDiscretization(std::make_shared<BasicGridGeometry>(gridView))
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
    friend inline LocalView localView(const PQ2FEGridDiscretization& gg)
    { return { gg.cache_ }; }

private:
    class PQ2FEGridDiscretizationCache
    {
        friend class PQ2FEGridDiscretization;
    public:
        //! export the geometry helper type
        using GeometryHelper = FEPQ2GeometryHelper<GV>;

        explicit PQ2FEGridDiscretizationCache(const PQ2FEGridDiscretization& gg)
        : gridDiscretization_(&gg)
        {}

        const PQ2FEGridDiscretization& gridDiscretization() const
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

        const PQ2FEGridDiscretization* gridDiscretization_;
    };

public:
    //! the cache type (only the caching implementation has this)
    //! this alias should only be used by the local view implementation
    using Cache = PQ2FEGridDiscretizationCache;
private:
    using GeometryHelper = typename Cache::GeometryHelper;

    void update_()
    {
        cache_.clear_();
        dofMapper_.update(this->gridView());

        boundaryDofIndices_.assign(numDofs(), false);

        for (const auto& element : elements(this->gridView()))
        {
            auto elementGeometry = element.geometry();
            const auto& localCoefficients = this->feCache().get(element.type()).localCoefficients();
            const auto eIdx = this->elementMapper().index(element);

            LocalIndexType numBoundaryFaces = 0;
            for (const auto& intersection : intersections(this->gridView(), element))
            {
                if (intersection.boundary() && !intersection.neighbor())
                {
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

                    // mark all dofs on this intersection as boundary dofs
                    for (LocalIndexType keyIdx = 0; keyIdx < localCoefficients.size(); ++keyIdx)
                    {
                        if (GeometryHelper::localDofOnIntersection(elementGeometry.type(),
                                                                    intersection.indexInInside(),
                                                                    localCoefficients.localKey(keyIdx)))
                        {
                            const auto dofIdxGlobal = GeometryHelper::dofIndex(
                                this->dofMapper(), element, localCoefficients.localKey(keyIdx));
                            boundaryDofIndices_[dofIdxGlobal] = true;
                        }
                    }
                }

                // inform the grid discretization if we have periodic boundaries
                else if (periodicGridTraits_.isPeriodic(intersection))
                {
                    this->setPeriodic();

                    const auto eps = 1e-7*(elementGeometry.corner(1) - elementGeometry.corner(0)).two_norm();
                    for (int localDofIdx = 0; localDofIdx < localCoefficients.size(); ++localDofIdx)
                    {
                        if (!GeometryHelper::localDofOnIntersection(elementGeometry.type(),
                                                                     intersection.indexInInside(),
                                                                     localCoefficients.localKey(localDofIdx)))
                            continue;

                        const auto dofIdxGlobal = GeometryHelper::dofIndex(
                            this->dofMapper(), element, localCoefficients.localKey(localDofIdx));
                        const auto dofPos = GeometryHelper::dofPosition(elementGeometry, localCoefficients.localKey(localDofIdx));

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
                                    if (!GeometryHelper::localDofOnIntersection(outsideGeometry.type(),
                                                                                 isOutside.indexInInside(),
                                                                                 localKeyOut))
                                        continue;

                                    const auto dofIdxGlobalOut = GeometryHelper::dofIndex(this->dofMapper(), outside, localKeyOut);
                                    const auto dofPosOutside = GeometryHelper::dofPosition(outsideGeometry, localKeyOut);
                                    const auto shift = std::abs((this->bBoxMax()-this->bBoxMin())*intersection.centerUnitOuterNormal());
                                    if (std::abs((dofPosOutside - dofPos).two_norm() - shift) < eps)
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
            DUNE_THROW(Dune::NotImplemented, "Periodic boundaries for pq2 FE method for parallel simulations!");
    }

    DofMapper dofMapper_;

    const FeCache feCache_;

    // dofs on the boundary
    std::vector<bool> boundaryDofIndices_;

    // a map for periodic boundary dofs
    std::unordered_map<GridIndexType, GridIndexType> periodicDofMap_;

    Cache cache_;

    PeriodicGridTraits<typename GridView::Grid> periodicGridTraits_;
};

} // end namespace Dumux

#endif
