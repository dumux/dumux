// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup PQ2Discretization
 * \brief Lightweight DOF helper for order-2 Lagrange elements.
 *
 * Provides `dofIndex` and `dofPosition` as pure static methods without
 * the scv/scvf machinery of HybridPQ2GeometryHelper.  The uniform
 * 4-argument `dofIndex` signature (with an ignored id-set argument) lets
 * callers use PQ2 and PQ3 helpers interchangeably.
 */
#ifndef DUMUX_DISCRETIZATION_PQ2_DOF_HELPER_HH
#define DUMUX_DISCRETIZATION_PQ2_DOF_HELPER_HH

#include <dune/common/fvector.hh>
#include <dune/geometry/referenceelements.hh>

namespace Dumux {

/*!
 * \ingroup PQ2Discretization
 * \brief DOF index and position helper for order-2 Lagrange discretizations.
 *
 * For order 2, every sub-entity carries exactly one DOF so no orientation
 * permutation is needed: `dofIndex` simply returns `mapper.subIndex(...)`.
 *
 * \tparam GridView  The Dune grid view type.
 */
template<class GridView>
struct PQ2LagrangeDofHelper
{
    using Scalar = typename GridView::ctype;
    static constexpr int dim = GridView::dimension;
    using GlobalPosition = Dune::FieldVector<Scalar, GridView::dimensionworld>;

    //! Global DOF index for the given local key.
    //! The id-set argument is accepted but ignored (no orientation needed for k=2).
    template<class DofMapper, class Element, class LocalKey, class IdSet>
    static std::size_t dofIndex(const DofMapper& m, const Element& e,
                                const LocalKey& lk, const IdSet& /*idSet*/)
    { return m.subIndex(e, lk.subEntity(), lk.codim()); }

    //! Global DOF index for the given local key (3-argument overload, no id-set).
    template<class DofMapper, class Element, class LocalKey>
    static std::size_t dofIndex(const DofMapper& m, const Element& e, const LocalKey& lk)
    { return m.subIndex(e, lk.subEntity(), lk.codim()); }

    //! Physical position of the DOF in global coordinates.
    template<class Geometry, class LocalKey>
    static GlobalPosition dofPosition(const Geometry& geo, const LocalKey& lk)
    {
        if (lk.codim() == dim) return geo.corner(lk.subEntity());
        if (lk.codim() == 0)  return geo.center();
        const auto& ref = Dune::referenceElement<Scalar, dim>(geo.type());
        return geo.global(ref.position(lk.subEntity(), lk.codim()));
    }
};

} // namespace Dumux

#endif // DUMUX_DISCRETIZATION_PQ2_DOF_HELPER_HH
