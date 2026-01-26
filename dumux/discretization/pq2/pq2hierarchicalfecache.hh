// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup PQ2Discretization
 * \brief A finite element cache for P2/Q2 hierarchical function spaces
 */
#ifndef DUMUX_DISCRETIZATION_PQ2_HIERARCHICAL_FECACHE_HH
#define DUMUX_DISCRETIZATION_PQ2_HIERARCHICAL_FECACHE_HH

#include <memory>

#include <dune/common/exceptions.hh>
#include <dune/geometry/type.hh>

#include <dune/localfunctions/common/virtualinterface.hh>
#include <dune/localfunctions/common/virtualwrappers.hh>

#include <dumux/discretization/pq2/pq2hierarchicallocalfiniteelement.hh>

namespace Dumux {

template< class CoordScalar, class Scalar, unsigned int dim>
class PQ2HierarchicalFECache
{
    static_assert(dim == 2 || dim == 3, "P2/Q2 hierarchical FE spaces only implemented for 2D and 3D grids");

    using P2Hierarchical = Dumux::PQ2HierarchicalLocalFiniteElement<CoordScalar, Scalar, dim, Dune::GeometryTypes::simplex(dim).toId()>;
    using Q2Hierarchical = Dumux::PQ2HierarchicalLocalFiniteElement<CoordScalar, Scalar, dim, Dune::GeometryTypes::cube(dim).toId()>;

public:
    using FiniteElementType = Dune::LocalFiniteElementVirtualInterface<typename P2Hierarchical::Traits::LocalBasisType::Traits>;

    PQ2HierarchicalFECache()
    : p2HierarchicalBasis_(std::make_unique<Dune::LocalFiniteElementVirtualImp<P2Hierarchical>>(P2Hierarchical{}))
    , q2HierarchicalBasis_(std::make_unique<Dune::LocalFiniteElementVirtualImp<Q2Hierarchical>>(Q2Hierarchical{}))
    {}

    //! Get local finite element for given GeometryType
    const FiniteElementType& get(const Dune::GeometryType& gt) const
    {
        if (gt.isSimplex())
            return *p2HierarchicalBasis_;
        if (gt.isCube())
            return *q2HierarchicalBasis_;
        else
            DUNE_THROW(Dune::NotImplemented, "P2/Q2 hierarchical local finite element for geometry type " << gt);
    }

private:
    std::unique_ptr<FiniteElementType> p2HierarchicalBasis_;
    std::unique_ptr<FiniteElementType> q2HierarchicalBasis_;
};

} // end namespace Dumux

#endif
