// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \brief Dynamic (post-refinement) load-balancing/rebalancing for the box/CVFE Boussinesq
 *        model: migrates (current, old) solution data as elements move between ranks.
 *
 * This is a separate, self-contained path from griddatatransfer.hh's h-adaptive
 * store()/reconstruct() (which uses the model's real, variable-size, discretization-level
 * ElementSolution type -- fine for h-adaptation, where entities never leave their owning
 * rank). Migrating data across ranks during loadBalance() instead needs a payload the
 * Dune::CommDataHandleIF machinery can move as raw bytes between processes; ElementSolution
 * is an opaque type from dumux/discretization/cvfe with no guaranteed "build one from n raw
 * values on an arbitrary rank" constructor, so this file defines its own small, fixed-size,
 * trivially-copyable payload (a padded array of PrimaryVariables, one slot per element
 * vertex, capped at 8 -- enough for any 2D/3D box element in these tests) purely for the
 * migration step, and converts to/from it directly against the solution vectors.
 */
#ifndef DUMUX_BOUSSINESQ_ADAPTIVE_BOXREBALANCE_HH
#define DUMUX_BOUSSINESQ_ADAPTIVE_BOXREBALANCE_HH

#include <array>
#include <cstddef>

#include <dune/grid/common/datahandleif.hh>
#include <dune/grid/common/rangegenerators.hh>
#include <dune/grid/utility/persistentcontainer.hh>

#include <dumux/common/properties.hh>
#include <dumux/discretization/elementsolution.hh>
#include <dumux/discretization/extrusion.hh>

#include "griddatacommunication.hh"
#include "loadbalancer.hh"

namespace Dumux::BoussinesqAdaptive {

//! Maximum number of vertices of any element type used in these 2D/3D box tests
//! (quad = 4, hex = 8; triangle/tet use fewer slots, left unused/zero).
inline constexpr int boxMaxVerticesPerElement = 8;

/*!
 * \brief Fixed-size, trivially-copyable per-element payload used only for load-balance
 *        migration (see file comment for why this isn't the same type griddatatransfer.hh
 *        uses for h-adaptation).
 */
template<class PrimaryVariables>
struct BoxMigrationValue
{
    std::array<PrimaryVariables, boxMaxVerticesPerElement> uCur{};
    std::array<PrimaryVariables, boxMaxVerticesPerElement> uOld{};
    unsigned int numVertices = 0;
};

/*!
 * \brief Dune::CommDataHandleIF that migrates a Dune::PersistentContainer<Grid,
 *        BoxMigrationValue> as elements move rank during loadBalance().
 */
template<class Grid, class PersistentContainer>
class BoxMigrationDataHandle
: public Dune::CommDataHandleIF<BoxMigrationDataHandle<Grid, PersistentContainer>,
                                 typename PersistentContainer::Value>
{
public:
    using DataType = typename PersistentContainer::Value;

    explicit BoxMigrationDataHandle(PersistentContainer& container) : container_(container) {}

    bool contains(int /*dim*/, int codim) const { return codim == 0; }
    bool fixedSize(int /*dim*/, int /*codim*/) const { return true; }

    template<class Entity>
    std::size_t size(const Entity&) const { return 1; }

    template<class MessageBuffer, class Entity>
    void gather(MessageBuffer& buff, const Entity& entity) const
    { buff.write(container_[entity]); }

    template<class MessageBuffer, class Entity>
    void scatter(MessageBuffer& buff, const Entity& entity, std::size_t /*n*/)
    {
        DataType value;
        buff.read(value);
        // A PersistentContainer sized before loadBalance() has no slot yet for an entity
        // that migrates in *during* the call -- resize before every scatter so it grows to
        // fit, exactly like dune-alugrid's own load-balance example does (see
        // examples/loadbalance/datamap.hh's DataHandle::scatter()). Skipping this corrupts
        // the heap (confirmed: reproduced as "malloc(): invalid next size" without it).
        container_.resize();
        container_[entity] = value;
    }

private:
    PersistentContainer& container_;
};

/*!
 * \brief Rebalance the grid, migrating (x, xOld) along with it.
 *
 * \param gridGeometry the box grid geometry (updated in place to match the new partition)
 * \param grid the grid to rebalance
 * \param x current solution, resized/repopulated in place
 * \param xOld previous-timestep solution, resized/repopulated in place
 * \param verticalStrips partition into weighted vertical strips instead of the grid's
 *        default method (ALUGrid only, see loadbalancer.hh's VerticalStripDestinations)
 * \return true if the partitioning actually changed
 */
template<class TypeTag, class Grid, class GridGeometry, class SolutionVector>
bool rebalanceBox(std::shared_ptr<GridGeometry> gridGeometry, Grid& grid,
                  SolutionVector& x, SolutionVector& xOld, bool verticalStrips = false)
{
    using PrimaryVariables = GetPropType<TypeTag, Properties::PrimaryVariables>;
    using MigrationValue = BoxMigrationValue<PrimaryVariables>;
    using PersistentContainer = Dune::PersistentContainer<Grid, MigrationValue>;

    PersistentContainer container(grid, 0);
    container.resize();

    // 1. snapshot the current per-element, per-vertex (cur, old) values into the container
    {
        auto fvGeometry = localView(*gridGeometry);
        for (const auto& element : elements(gridGeometry->gridView()))
        {
            fvGeometry.bindElement(element);
            auto& value = container[element];
            value.numVertices = static_cast<unsigned int>(fvGeometry.numScv());
            for (const auto& scv : scvs(fvGeometry))
            {
                value.uCur[scv.localDofIndex()] = x[scv.dofIndex()];
                value.uOld[scv.localDofIndex()] = xOld[scv.dofIndex()];
            }
        }
    }

    // 2. migrate: elements (and their container entries) move rank as needed
    BoxMigrationDataHandle<Grid, PersistentContainer> dataHandle(container);
    const bool changed = GridLoadBalancer<Grid>::apply(grid, dataHandle, verticalStrips);
    if (!changed)
        return false;

    // 3. indices changed (elements arrived/left), rebuild the grid geometry
    gridGeometry->update(grid.leafGridView());
    container.resize();

    // 4. rebuild (x, xOld) from the (now-migrated) container: same volume-weighted-average
    // reconstruction as griddatatransfer.hh's reconstruct_(), just sourced from the fixed-
    // size migration payload instead of ElementSolution.
    const auto numDofs = gridGeometry->numDofs();
    x.resize(numDofs);
    xOld.resize(numDofs);

    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    std::vector<Scalar> weight(numDofs, 0.0);
    SolutionVector accCur(numDofs);
    SolutionVector accOld(numDofs);
    accCur = 0.0;
    accOld = 0.0;

    auto fvGeometry = localView(*gridGeometry);
    for (const auto& element : elements(gridGeometry->gridView()))
    {
        fvGeometry.bindElement(element);
        const auto& value = container[element];
        for (const auto& scv : scvs(fvGeometry))
        {
            const auto dofIdx = scv.dofIndex();
            const auto w = Dumux::Extrusion_t<GridGeometry>::volume(fvGeometry, scv);
            weight[dofIdx] += w;
            accCur[dofIdx] += value.uCur[scv.localDofIndex()] * w;
            accOld[dofIdx] += value.uOld[scv.localDofIndex()] * w;
        }
    }

    for (std::size_t dofIdx = 0; dofIdx < numDofs; ++dofIdx)
    {
        if (weight[dofIdx] > 0.0)
        {
            x[dofIdx] = accCur[dofIdx] / weight[dofIdx];
            xOld[dofIdx] = accOld[dofIdx] / weight[dofIdx];
        }
    }

    const auto& gridView = gridGeometry->gridView();
    static constexpr int dim = GridGeometry::GridView::dimension;
    syncOwnerToGhost<dim>(gridView, gridGeometry->dofMapper(), x);
    syncOwnerToGhost<dim>(gridView, gridGeometry->dofMapper(), xOld);

    return true;
}

} // end namespace Dumux::BoussinesqAdaptive

#endif
