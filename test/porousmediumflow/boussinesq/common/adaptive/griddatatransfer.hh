// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \brief Transfers the (current, old) solution across h-adaptation for the box/CVFE
 *        Boussinesq models.
 *
 * Follows the same store()/reconstruct() pattern as
 * dumux/porousmediumflow/2p/griddatatransfer.hh (Dune::PersistentContainer indexed by
 * codim-0 elements, each entry holding the element-local, per-vertex solution blob), but is
 * much simpler because there is no phase mass to conserve here -- box's vertex DOFs are
 * plain nodal values of (vector potential, concentration), so reconstruction is a volume-
 * weighted average/interpolation, not a mass-conservative redistribution.
 *
 * Two things this local copy fixes relative to the 2p version:
 *  1. It also carries the *previous* timestep's solution (xOld) through adaptation, not just
 *     the current iterate. 2p's version only has one SolutionVector; a transient model like
 *     Boussinesq needs xOld for the next timestep's time-derivative term, and silently
 *     dropping it after an adaptation step would corrupt the following residual evaluation
 *     without any visible error.
 *  2. reconstruct() finishes with an explicit owner->ghost sync of the finished solution
 *     (see griddatacommunication.hh), so a shared border vertex ends up with the exact same
 *     value on every rank that has a copy of it. 2p's dead code tried to do something in
 *     this spirit for cell data via a `VectorExchange` class that no longer exists in the
 *     tree and was never replaced; the sync used here is the same
 *     dumux/parallel/vectorcommdatahandle.hh primitive that already backs the rest of
 *     DuMux's parallel infrastructure.
 */
#ifndef DUMUX_BOUSSINESQ_ADAPTIVE_GRIDDATATRANSFER_HH
#define DUMUX_BOUSSINESQ_ADAPTIVE_GRIDDATATRANSFER_HH

#include <memory>
#include <vector>

#include <dune/grid/common/rangegenerators.hh>
#include <dune/grid/utility/persistentcontainer.hh>

#include <dumux/common/properties.hh>
#include <dumux/discretization/elementsolution.hh>
#include <dumux/discretization/evalsolution.hh>
#include <dumux/discretization/extrusion.hh>
#include <dumux/adaptive/griddatatransfer.hh>

#include "griddatacommunication.hh"

namespace Dumux {

/*!
 * \brief Store/reconstruct (current, old) solution across h-adaptation for the box/CVFE
 *        Boussinesq model.
 */
template<class TypeTag>
class BoussinesqBoxGridDataTransfer : public GridDataTransfer<GetPropType<TypeTag, Properties::Grid>>
{
    using Grid = GetPropType<TypeTag, Properties::Grid>;
    using ParentType = GridDataTransfer<Grid>;
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using SolutionVector = GetPropType<TypeTag, Properties::SolutionVector>;
    using Element = typename Grid::template Codim<0>::Entity;
    using ElementSolution = std::decay_t<decltype(elementSolution(std::declval<Element>(),
                                                                   std::declval<SolutionVector>(),
                                                                   std::declval<GridGeometry>()))>;

    struct AdaptedValues
    {
        ElementSolution uCur;
        ElementSolution uOld;
        bool wasLeaf = false;
    };

    using PersistentContainer = Dune::PersistentContainer<Grid, AdaptedValues>;

public:
    BoussinesqBoxGridDataTransfer(std::shared_ptr<GridGeometry> gridGeometry,
                                  SolutionVector& x,
                                  SolutionVector& xOld)
    : ParentType()
    , gridGeometry_(gridGeometry)
    , x_(x)
    , xOld_(xOld)
    , adaptionMap_(gridGeometry->gridView().grid(), 0)
    {}

    void store(const Grid& grid) override
    {
        adaptionMap_.resize();

        for (auto level = grid.maxLevel(); level >= 0; level--)
        {
            for (const auto& element : elements(grid.levelGridView(level)))
            {
                auto& adaptedValues = adaptionMap_[element];

                if (element.isLeaf())
                {
                    adaptedValues.uCur = elementSolution(element, x_, *gridGeometry_);
                    adaptedValues.uOld = elementSolution(element, xOld_, *gridGeometry_);
                    adaptedValues.wasLeaf = true;
                }
                else
                {
                    // Non-leaf elements: box DOFs live on vertices that already exist on the
                    // leaf grid, so the coarse element's own nodal values can be read
                    // directly from the current/old solution -- no averaging needed.
                    adaptedValues.uCur = elementSolution(element, x_, *gridGeometry_);
                    adaptedValues.uOld = elementSolution(element, xOld_, *gridGeometry_);
                    adaptedValues.wasLeaf = false;
                }
            }
        }
    }

    void reconstruct(const Grid& grid) override
    {
        gridGeometry_->update(grid.leafGridView());
        reconstruct_();
    }

private:
    void reconstruct_()
    {
        adaptionMap_.resize();

        const auto numDofs = gridGeometry_->numDofs();
        x_.resize(numDofs);
        xOld_.resize(numDofs);

        std::vector<Scalar> weight(numDofs, 0.0);
        SolutionVector accCur(numDofs);
        SolutionVector accOld(numDofs);
        accCur = 0.0;
        accOld = 0.0;

        auto fvGeometry = localView(*gridGeometry_);
        for (const auto& element : elements(gridGeometry_->gridView()))
        {
            fvGeometry.bindElement(element);

            if (!element.isNew())
            {
                const auto& adaptedValues = adaptionMap_[element];
                for (const auto& scv : scvs(fvGeometry))
                {
                    const auto dofIdx = scv.dofIndex();
                    const auto w = Dumux::Extrusion_t<GridGeometry>::volume(fvGeometry, scv);
                    weight[dofIdx] += w;
                    accCur[dofIdx] += adaptedValues.uCur[scv.localDofIndex()] * w;
                    accOld[dofIdx] += adaptedValues.uOld[scv.localDofIndex()] * w;
                }
            }
            else
            {
                // find the closest ancestor that existed on the old grid
                assert(element.hasFather() && "new element does not have a father element!");
                auto fatherElement = element.father();
                while (fatherElement.isNew() && fatherElement.level() > 0)
                    fatherElement = fatherElement.father();

                const auto& adaptedValuesFather = adaptionMap_[fatherElement];
                const auto fatherGeometry = fatherElement.geometry();

                for (const auto& scv : scvs(fvGeometry))
                {
                    const auto dofIdx = scv.dofIndex();
                    const auto w = Dumux::Extrusion_t<GridGeometry>::volume(fvGeometry, scv);
                    weight[dofIdx] += w;
                    accCur[dofIdx] += evalSolution(fatherElement, fatherGeometry, *gridGeometry_,
                                                    adaptedValuesFather.uCur, scv.dofPosition()) * w;
                    accOld[dofIdx] += evalSolution(fatherElement, fatherGeometry, *gridGeometry_,
                                                    adaptedValuesFather.uOld, scv.dofPosition()) * w;
                }
            }
        }

        for (std::size_t dofIdx = 0; dofIdx < numDofs; ++dofIdx)
        {
            // weight[dofIdx] can be zero for a vertex that only exists on a locally-owned
            // ghost element with no corresponding entry (shouldn't happen for dofs relevant
            // to this rank's interior+border set, but guard against div-by-zero regardless)
            if (weight[dofIdx] > 0.0)
            {
                x_[dofIdx] = accCur[dofIdx] / weight[dofIdx];
                xOld_[dofIdx] = accOld[dofIdx] / weight[dofIdx];
            }
        }

        adaptionMap_.resize(typename PersistentContainer::Value());
        adaptionMap_.shrinkToFit();
        adaptionMap_.fill(typename PersistentContainer::Value());

        // Make sure every rank's copy of a shared (border/overlap/ghost) vertex holds
        // exactly the owning rank's reconstructed value -- see griddatacommunication.hh.
        const auto& gridView = gridGeometry_->gridView();
        static constexpr int dim = GridGeometry::GridView::dimension;
        BoussinesqAdaptive::syncOwnerToGhost<dim>(gridView, gridGeometry_->dofMapper(), x_);
        BoussinesqAdaptive::syncOwnerToGhost<dim>(gridView, gridGeometry_->dofMapper(), xOld_);
    }

    std::shared_ptr<GridGeometry> gridGeometry_;
    SolutionVector& x_;
    SolutionVector& xOld_;
    PersistentContainer adaptionMap_;
};

} // end namespace Dumux

#endif
