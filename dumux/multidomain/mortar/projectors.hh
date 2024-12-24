// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup MultiDomain
 * \ingroup MortarCoupling
 * \brief Predefined projector implementations for mortar-coupling models.
 */
#ifndef DUMUX_MULTIDOMAIN_MORTAR_PROJECTORS_HH
#define DUMUX_MULTIDOMAIN_MORTAR_PROJECTORS_HH

#include <config.h>

#include <memory>
#include <utility>
#include <type_traits>

#include <dune/functions/functionspacebases/lagrangebasis.hh>

#include <dumux/geometry/boundingboxtree.hh>
#include <dumux/geometry/geometricentityset.hh>
#include <dumux/geometry/intersectionentityset.hh>

#include <dumux/io/grid/facetgridmanager.hh>
#include <dumux/discretization/projection/projector.hh>
#include <dumux/discretization/functionspacebasis.hh>

#include "projectorinterface.hh"

namespace Dumux::Mortar {

template<std::size_t o> struct MortarTraceOrder {};
template<std::size_t o> struct ResidualTraceOrder {};

/*!
 * \ingroup MultiDomain
 * \ingroup MortarCoupling
 * \brief Default projector for finite-volume schemes and mortars representing primary variables (e.g. pressure mortars).
 *        Projects between a sub-domain trace and the mortar domain.
 */
template<typename SolutionVector>
class FVDefaultProjector : public Projector<SolutionVector>
{
    // TODO: generalize this? (assumes Dune::BlockVector).. But Dumux::Projector hardcodes BlockVector, anyway...
    using Scalar = typename SolutionVector::field_type;
    using L2Projector = Dumux::Projector<Scalar>;

 public:
    // TODO: can ordering of function space basis be different and break the mapping to boundaries later?
    template<typename MortarGridGeometry, typename TraceGridManager, std::size_t M, std::size_t R>
    FVDefaultProjector(const MortarGridGeometry& mortarGridGeometry,
                       const TraceGridManager& traceGridManager,
                       const MortarTraceOrder<M>&,
                       const ResidualTraceOrder<R>&)
    : FVDefaultProjector(
        mortarGridGeometry,
        traceGridManager,
        Dune::Functions::LagrangeBasis<typename TraceGridManager::Grid::LeafGridView, M>{traceGridManager.grid().leafGridView()},
        Dune::Functions::LagrangeBasis<typename TraceGridManager::Grid::LeafGridView, R>{traceGridManager.grid().leafGridView()}
    )
    {}

 private:
    template<typename MortarGridGeometry,
             typename TraceGridManager,
             typename MortarTraceBasis,
             typename ResidualTraceBasis>
    FVDefaultProjector(const MortarGridGeometry& mortarGridGeometry,
                       const TraceGridManager& traceGridManager,
                       const MortarTraceBasis& mortarTraceBasis,
                       const ResidualTraceBasis& residualTraceBasis)
    {
        using MortarEntitySet = typename std::remove_cvref_t<decltype(mortarGridGeometry.boundingBoxTree())>::EntitySet;
        using TraceEntitySet = GridViewGeometricEntitySet<typename TraceGridManager::Grid::LeafGridView>;
        BoundingBoxTree<TraceEntitySet> traceTree{std::make_shared<TraceEntitySet>(traceGridManager.grid().leafGridView())};
        IntersectionEntitySet<MortarEntitySet, TraceEntitySet> glue;
        glue.build(mortarGridGeometry.boundingBoxTree(), traceTree);

        const auto& mortarBasis = getFunctionSpaceBasis(mortarGridGeometry);
        to_ = std::make_unique<L2Projector>(makeProjector(mortarBasis, mortarTraceBasis, glue));
        from_ = std::make_unique<L2Projector>(makeProjectorPair(mortarBasis, residualTraceBasis, glue).second);
    }

    SolutionVector toTrace_(const SolutionVector& x) const override { return to_->project(x); }
    SolutionVector fromTrace_(const SolutionVector& x) const override { return from_->project(x); }

    std::unique_ptr<L2Projector> to_;
    std::unique_ptr<L2Projector> from_;
};

} // end namespace Dumux::Mortar

#endif
