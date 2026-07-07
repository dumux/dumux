// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup FEDiscretization
 * \brief Default class for finite element grid discretizations.
 */
#ifndef DUMUX_DISCRETIZATION_FE_GRID_DISCRETIZATION_HH
#define DUMUX_DISCRETIZATION_FE_GRID_DISCRETIZATION_HH

#include <span>
#include <utility>
#include <vector>
#include <unordered_map>

#include <dune/common/reservedvector.hh>

#include <dumux/common/indextraits.hh>

#include <dumux/discretization/basegridgeometry.hh>
#include <dumux/discretization/boundaryface.hh>
#include <dumux/discretization/extrusion.hh>
#include <dumux/discretization/cvfe/quadraturerules.hh>

#include <dumux/io/grid/periodicgridtraits.hh>

namespace Dumux::Experimental {

/*!
 * \ingroup FEDiscretization
 * \brief Default class for finite element grid discretizations.
 * \tparam GV the grid view type
 * \tparam Traits the grid geometry traits; must provide FeCache and DofHelper
 */
template<class GV, class Traits>
class FEGridDiscretization
: public BaseGridGeometry<GV, Traits>
{
    using ThisType = FEGridDiscretization<GV, Traits>;
    using ParentType = BaseGridGeometry<GV, Traits>;
    using GridIndexType = typename IndexTraits<GV>::GridIndex;
    using LocalIndexType = typename IndexTraits<GV>::LocalIndex;
    static const int dim = GV::dimension;

public:
    //! export the discretization method this geometry belongs to
    using DiscretizationMethod = typename Traits::DiscretizationMethod;
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
    FEGridDiscretization(std::shared_ptr<BasicGridGeometry> gg)
    : ParentType(std::move(gg))
    , dofMapper_(this->gridView(), Traits::layout())
    , cache_(*this)
    , periodicGridTraits_(this->gridView().grid())
    {
        update_();
    }

    //! Constructor
    FEGridDiscretization(const GridView& gridView)
    : FEGridDiscretization(std::make_shared<BasicGridGeometry>(gridView))
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
    friend inline LocalView localView(const FEGridDiscretization& gg)
    { return { gg.cache_ }; }

private:
    class FEGridDiscretizationCache
    {
        friend class FEGridDiscretization;
    public:
        //! export the dof helper type
        using DofHelper = typename Traits::DofHelper;

        explicit FEGridDiscretizationCache(const FEGridDiscretization& gg)
        : gridDiscretization_(&gg)
        {}

        const FEGridDiscretization& gridDiscretization() const
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

        const FEGridDiscretization* gridDiscretization_;
    };

public:
    //! the cache type (only the caching implementation has this)
    //! this alias should only be used by the local view implementation
    using Cache = FEGridDiscretizationCache;

private:
    //! Default update: fills cache_ using DofHelper::localDofOnIntersection.
    void update_()
    {
        using DofHelper = typename Cache::DofHelper;

        dofMapper_.update(this->gridView());
        cache_.clear_();
        boundaryDofIndices_.assign(numDofs(), false);

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
                        if (!DofHelper::localDofOnIntersection(elementGeometry.type(),
                                                               intersection.indexInInside(),
                                                               localCoefficients.localKey(localDofIdx)))
                            continue;

                        const auto dofIdxGlobal = DofHelper::dofIndex(
                            this->dofMapper(), element, localCoefficients.localKey(localDofIdx));
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

                                    const auto dofIdxGlobalOut = DofHelper::dofIndex(this->dofMapper(), outside, localKeyOut);
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
            DUNE_THROW(Dune::NotImplemented, "Periodic boundaries for FE method for parallel simulations!");
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
