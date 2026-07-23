// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \brief Parallel communication helpers for the Boussinesq adaptive-grid machinery.
 *
 * This is the local, corrected replacement for the dead `VectorExchange`-based
 * communication that was never finished in dumux/porousmediumflow/2p/{gridadaptindicator,
 * griddatatransfer}.hh (see the "// TODO: fix adaptive simulations in parallel" blocks
 * there, which reference a class that no longer exists). It is built on top of
 * dumux/parallel/vectorcommdatahandle.hh, which already provides a correct, modern
 * Dune::CommDataHandleIF implementation for exchanging border/overlap/ghost entries of a
 * mapper-indexed vector -- nothing about that primitive itself needed fixing, it was just
 * never wired up here.
 *
 * Two things are needed that were missing before:
 *  1. Border/overlap-consistent *indicator* values (per-element scalar deltas), so every
 *     rank computes the same refine/coarsen decision near a partition boundary.
 *  2. Border/overlap-consistent *solution* values before an adaptation step, so the
 *     per-element data stored for h-adaptive reconstruction (see griddatatransfer.hh) is
 *     valid on ghost/overlap copies, not just on interior elements.
 */
#ifndef DUMUX_BOUSSINESQ_ADAPTIVE_GRIDDATACOMMUNICATION_HH
#define DUMUX_BOUSSINESQ_ADAPTIVE_GRIDDATACOMMUNICATION_HH

#include <vector>

#include <dune/grid/common/gridenums.hh>
#include <dumux/parallel/vectorcommdatahandle.hh>

namespace Dumux::BoussinesqAdaptive {

/*!
 * \brief Synchronize a mapper-indexed vector across interior/border/overlap/ghost so that
 *        every rank's local copy of a shared entity holds the owning rank's value.
 *
 * \tparam entityCodim codim of the mapped entities (0 for cctpfa/elements, dim for box/vertices)
 * \param gridView the (parallel) grid view the mapper is defined on
 * \param mapper   entity mapper matching the vector's indexing
 * \param vector   the vector to update in place; entries on interior/owned entities are
 *                 unchanged, entries on border/overlap/ghost entities are overwritten with
 *                 the value gathered from the owning rank
 */
template<int entityCodim, class GridView, class Mapper, class Vector>
void syncOwnerToGhost(const GridView& gridView, const Mapper& mapper, Vector& vector)
{
    if (gridView.comm().size() <= 1)
        return;

    VectorCommDataHandleEqual<Mapper, Vector, entityCodim> handle(mapper, vector);
    gridView.communicate(handle, Dune::InteriorBorder_All_Interface, Dune::ForwardCommunication);
}

/*!
 * \brief Element-wise max-reduce a scalar vector across interior/border so that overlap and
 *        ghost copies of a shared element see the same (maximum) indicator value as its owner.
 *
 * Used by the concentration-gradient indicator to keep refine/coarsen decisions consistent
 * across a partition boundary -- this, together with the missing gridView.comm().max/min
 * reduction of the global refine/coarsen thresholds (done directly in the indicator, no
 * data handle needed for that part), is exactly what the dead 2p code never finished.
 */
template<class GridView, class ElementMapper, class Vector>
void syncElementDataMax(const GridView& gridView, const ElementMapper& elementMapper, Vector& vector)
{
    if (gridView.comm().size() <= 1)
        return;

    VectorCommDataHandleMax<ElementMapper, Vector, 0> handle(elementMapper, vector);
    gridView.communicate(handle, Dune::InteriorBorder_All_Interface, Dune::ForwardCommunication);
}

} // end namespace Dumux::BoussinesqAdaptive

#endif
