// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup PQ3Discretization
 * \brief Finite element grid discretization for the PQ3 method.
 */
#ifndef DUMUX_DISCRETIZATION_PQ3_FE_GRID_DISCRETIZATION_HH
#define DUMUX_DISCRETIZATION_PQ3_FE_GRID_DISCRETIZATION_HH

#include <dune/localfunctions/lagrange/lagrangelfecache.hh>

#include <dumux/discretization/method.hh>
#include <dumux/common/defaultmappertraits.hh>
#include <dumux/discretization/cvfe/quadraturerules.hh>
#include <dumux/discretization/fem/feelementdiscretization.hh>
#include <dumux/discretization/fem/fegriddiscretization.hh>
#include <dumux/discretization/pq3/fvgridgeometry.hh>

#include "dofhelper.hh"

namespace Dumux::Experimental {

/*!
 * \ingroup PQ3Discretization
 * \brief Quadrature rule traits for PQ3 FE discretization
 */
template<class GridView,
         class ElementRule = Dumux::QuadratureRules::DuneQuadrature<6>,
         class IntersectionRule = Dumux::QuadratureRules::DuneQuadrature<6>,
         class BoundaryFaceRule = Dumux::QuadratureRules::DuneQuadrature<6>>
struct PQ3FEQuadratureTraits
{
    using ElementQuadratureRule = ElementRule;
    using IntersectionQuadratureRule = IntersectionRule;
    using BoundaryFaceQuadratureRule = BoundaryFaceRule;
};

/*!
 * \ingroup PQ3Discretization
 * \brief The default traits for the PQ3 finite element grid discretization
 * \tparam GridView the grid view type
 * \tparam Scalar the scalar type (used for the FE cache)
 */
template<class GridView,
         class Scalar,
         class MapperTraits = Dumux::PQ3MapperTraits<GridView>,
         class QuadratureTraits = PQ3FEQuadratureTraits<GridView>>
struct PQ3FEDefaultGridDiscretizationTraits
: public MapperTraits, public QuadratureTraits
{
    using FeCache = Dune::LagrangeLocalFiniteElementCache<typename GridView::ctype, Scalar, GridView::dimension, 3>;
    using DofHelper = PQ3LagrangeDofHelper<GridView>;

    template<class GridDiscretization, bool enableCache>
    using LocalView = FEElementDiscretization<GridDiscretization, enableCache>;

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
 * \brief Class for pq3 finite element grid discretizations.
 */
template<class Scalar,
         class GV,
         bool enableCaching = true,
         class Traits = PQ3FEDefaultGridDiscretizationTraits<GV, Scalar>>
class PQ3FEGridDiscretization
: public BaseGridGeometry<GV, Traits>
{
    using ThisType = PQ3FEGridDiscretization<Scalar, GV, enableCaching, Traits>;
    using ParentType = BaseGridGeometry<GV, Traits>;
    using GridIndexType = typename IndexTraits<GV>::GridIndex;
    using LocalIndexType = typename IndexTraits<GV>::LocalIndex;
    static const int dim = GV::dimension;

public:
    //! export the discretization method this geometry belongs to
    using DiscretizationMethod = typename DiscretizationMethods::PQ3;
    static constexpr DiscretizationMethod discMethod{};

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
    using FeCache = typename Traits::FeCache;
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

    //! Constructor
    PQ3FEGridDiscretization(std::shared_ptr<BasicGridGeometry> gg)
    : ParentType(std::move(gg))
    , dofMapper_(this->gridView(), Traits::layout())
    , cache_(*this)
    , periodicGridTraits_(this->gridView().grid())
    {
        update_();
    }

    //! Constructor
    PQ3FEGridDiscretization(const GridView& gridView)
    : PQ3FEGridDiscretization(std::make_shared<BasicGridGeometry>(gridView))
    {}

    //! The dof mapper
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

    //! If a d.o.f. is on the boundary
    bool dofOnBoundary(GridIndexType dofIdx) const
    { return boundaryDofIndices_[dofIdx]; }

    //! If a d.o.f. is on a periodic boundary
    bool dofOnPeriodicBoundary(GridIndexType dofIdx) const
    { return periodicDofMap_.count(dofIdx); }

    //! The index of the d.o.f. on the other side of the periodic boundary
    GridIndexType periodicallyMappedDof(GridIndexType dofIdx) const
    { return periodicDofMap_.at(dofIdx); }

    //! Returns the map between dofs across periodic boundaries
    const std::unordered_map<GridIndexType, GridIndexType>& periodicDofMap() const
    { return periodicDofMap_; }

    //! local view of this object (constructed with the internal cache)
    friend inline LocalView localView(const PQ3FEGridDiscretization& gg)
    { return { gg.cache_ }; }

private:
    class PQ3FEGridDiscretizationCache
    {
        friend class PQ3FEGridDiscretization;
    public:
        //! export the dof helper type
        using DofHelper = typename Traits::DofHelper;

        explicit PQ3FEGridDiscretizationCache(const PQ3FEGridDiscretization& gg)
        : gridDiscretization_(&gg)
        {}

        const PQ3FEGridDiscretization& gridDiscretization() const
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

        const PQ3FEGridDiscretization* gridDiscretization_;
    };

public:
    //! the cache type (only the caching implementation has this)
    //! this alias should only be used by the local view implementation
    using Cache = PQ3FEGridDiscretizationCache;

private:
    //! Default update: fills cache_ using DofHelper::localDofOnIntersection.
    void update_()
    {
        using DofHelper = typename Cache::DofHelper;

        dofMapper_.update(this->gridView());
        cache_.clear_();
        boundaryDofIndices_.assign(numDofs(), false);
        const auto& idSet = this->gridView().grid().globalIdSet();

        for (const auto& element : elements(this->gridView()))
        {
            auto elementGeometry = element.geometry();
            const auto& localCoefficients = this->feCache().get(element.type()).localCoefficients();
            const auto eIdx = this->elementMapper().index(element);

            LocalIndexType numBoundaryFaces = 0;
            for (const auto& intersection : intersections(this->gridView(), element))
            {
                if (intersection.boundary() && !periodicGridTraits_.isPeriodic(intersection))
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
                        if (DofHelper::localDofOnIntersection(elementGeometry.type(),
                                                              intersection.indexInInside(),
                                                              localCoefficients.localKey(keyIdx)))
                        {
                            const auto dofIdxGlobal = DofHelper::dofIndex(
                                this->dofMapper(), element, localCoefficients.localKey(keyIdx), idSet);
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
                        if (!DofHelper::localDofOnIntersection(elementGeometry.type(),
                                                               intersection.indexInInside(),
                                                               localCoefficients.localKey(localDofIdx)))
                            continue;

                        const auto dofIdxGlobal = DofHelper::dofIndex(
                            this->dofMapper(), element, localCoefficients.localKey(localDofIdx), idSet);
                        const auto dofPos = DofHelper::dofPosition(elementGeometry, localCoefficients.localKey(localDofIdx));

                        const auto& outside = intersection.outside();
                        const auto outsideGeometry = outside.geometry();
                        const auto& localCoefficientsOut = this->feCache().get(outsideGeometry.type()).localCoefficients();
                        for (const auto& isOutside : intersections(this->gridView(), outside))
                        {
                            if (periodicGridTraits_.isPeriodic(isOutside))
                            {
                                for (int localDofIdxOut = 0; localDofIdxOut < localCoefficientsOut.size(); ++localDofIdxOut)
                                {
                                    const auto& localKeyOut = localCoefficientsOut.localKey(localDofIdxOut);
                                    if (!DofHelper::localDofOnIntersection(outsideGeometry.type(),
                                                                           isOutside.indexInInside(),
                                                                           localKeyOut))
                                        continue;

                                    const auto dofIdxGlobalOut = DofHelper::dofIndex(this->dofMapper(), outside, localKeyOut, idSet);
                                    const auto dofPosOutside = DofHelper::dofPosition(outsideGeometry, localKeyOut);
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

        // error check: periodic boundaries currently don't work in parallel
        if (this->isPeriodic() && this->gridView().comm().size() > 1)
            DUNE_THROW(Dune::NotImplemented, "Periodic boundaries for FE PQ3 method for parallel simulations!");
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

} // end namespace Dumux::Experimental

#endif
