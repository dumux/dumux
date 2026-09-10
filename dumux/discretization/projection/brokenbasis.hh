// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup Discretization
 * \brief A broken (discontinuous) Lagrange basis over an arbitrary entity set
 * \note Exposes the function space basis interface the projections expect, which follows the
 *       global basis interface of dune-functions: size and localView, with the local view
 *       providing bind, tree().finiteElement() and index.
 */
#ifndef DUMUX_DISCRETIZATION_PROJECTION_BROKEN_BASIS_HH
#define DUMUX_DISCRETIZATION_PROJECTION_BROKEN_BASIS_HH

#include <cstddef>
#include <memory>
#include <utility>
#include <vector>

#include <dune/localfunctions/lagrange/lagrangelfecache.hh>

namespace Dumux {

/*!
 * \ingroup Discretization
 * \brief A discontinuous (broken) Lagrange basis of the given order over an entity set.
 *
 * Each entity owns a disjoint block of degrees of freedom, sized by the local finite
 * element for its own geometry type. Because no degrees of freedom are shared, the
 * basis needs no adjacency information and can therefore be defined over a bare list
 * of geometries as readily as over a grid.
 *
 * \tparam ES The entity set type
 * \tparam order The polynomial order
 * \tparam RangeField The type used for shape function values
 */
template<class ES, int order, class RangeField = double>
class BrokenLagrangeBasis
{
    using Entity = typename ES::Entity;
    using Geometry = typename Entity::Geometry;
    using ctype = typename ES::ctype;

    using FECache = Dune::LagrangeLocalFiniteElementCache<ctype, RangeField, Geometry::mydimension, order>;

public:
    using EntitySet = ES;
    using FiniteElement = typename FECache::FiniteElementType;

    static constexpr int dimension = Geometry::mydimension;

    struct LocalTree
    {
        using FiniteElement = typename BrokenLagrangeBasis::FiniteElement;

        const FiniteElement& finiteElement() const
        { return *fe_; }

        const FiniteElement* fe_ = nullptr;
    };

    class LocalView
    {
    public:
        using Tree = LocalTree;

        explicit LocalView(const BrokenLagrangeBasis& basis) : basis_(basis) {}

        void bind(const Entity& entity)
        {
            entityIndex_ = basis_.entitySet_->index(entity);
            tree_.fe_ = &basis_.feCache_.get(entity.geometry().type());
        }

        const Tree& tree() const
        { return tree_; }

        std::size_t index(std::size_t i) const
        { return basis_.offset_[entityIndex_] + i; }

    private:
        const BrokenLagrangeBasis& basis_;
        std::size_t entityIndex_{};
        Tree tree_{};
    };

    explicit BrokenLagrangeBasis(std::shared_ptr<const EntitySet> entitySet)
    : entitySet_(std::move(entitySet))
    {
        offset_.assign(entitySet_->size() + 1, 0);
        for (const auto& entity : *entitySet_)
            offset_[entitySet_->index(entity) + 1] = feCache_.get(entity.geometry().type()).size();
        for (std::size_t i = 1; i < offset_.size(); ++i)
            offset_[i] += offset_[i-1];
    }

    std::size_t size() const
    { return offset_.back(); }

    LocalView localView() const
    { return LocalView(*this); }

    const EntitySet& entitySet() const
    { return *entitySet_; }

private:
    std::shared_ptr<const EntitySet> entitySet_;
    FECache feCache_;
    //! First global index of each entity, with the total count in the last slot
    std::vector<std::size_t> offset_;
};

} // end namespace Dumux

#endif
