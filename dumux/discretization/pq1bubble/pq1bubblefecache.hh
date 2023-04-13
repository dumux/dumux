// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup PQ1BubbleDiscretization
 * \brief A finite element cache for P1/Q1 function with bubble
 */
#ifndef DUMUX_DISCRETIZATION_PQ1BUBBLE_FECACHE_HH
#define DUMUX_DISCRETIZATION_PQ1BUBBLE_FECACHE_HH

#include <memory>

#include <dune/common/exceptions.hh>
#include <dune/geometry/type.hh>

#include <dune/localfunctions/common/virtualinterface.hh>
#include <dune/localfunctions/common/virtualwrappers.hh>

#include <dumux/discretization/pq1bubble/pq1bubblelocalfiniteelement.hh>

namespace Dumux {

template< class CoordScalar, class Scalar, unsigned int dim>
class PQ1BubbleFECache
{
    static_assert(dim == 2 || dim == 3, "P1/Q1 bubble FE spaces only implemented for 2D and 3D grids");

    // These are so-called non-conforming finite element spaces
    // the local basis is only continuous at given points on the faces
    using P1Bubble = Dumux::PQ1BubbleLocalFiniteElement<CoordScalar, Scalar, dim, Dune::GeometryTypes::simplex(dim).toId()>;
    using Q1Bubble = Dumux::PQ1BubbleLocalFiniteElement<CoordScalar, Scalar, dim, Dune::GeometryTypes::cube(dim).toId()>;

public:
    using FiniteElementType = Dune::LocalFiniteElementVirtualInterface<typename P1Bubble::Traits::LocalBasisType::Traits>;

    PQ1BubbleFECache()
    : p1BubbleBasis_(std::make_unique<Dune::LocalFiniteElementVirtualImp<P1Bubble>>(P1Bubble{}))
    , q1BubbleBasis_(std::make_unique<Dune::LocalFiniteElementVirtualImp<Q1Bubble>>(Q1Bubble{}))
    {}

    //! Get local finite element for given GeometryType
    const FiniteElementType& get(const Dune::GeometryType& gt) const
    {
        if (gt.isSimplex())
            return *p1BubbleBasis_;
        if (gt.isCube())
            return *q1BubbleBasis_;
        else
            DUNE_THROW(Dune::NotImplemented,
                "Lagrange bubble local finite element for geometry type " << gt
            );
    }

private:
    std::unique_ptr<FiniteElementType> p1BubbleBasis_;
    std::unique_ptr<FiniteElementType> q1BubbleBasis_;
};

} // end namespace Dumux

#endif
