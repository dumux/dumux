// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup PQ1BubbleDiscretization
 * \brief Finite element grid discretization for the PQ1Bubble method.
 */
#ifndef DUMUX_DISCRETIZATION_PQ1BUBBLE_FE_GRID_DISCRETIZATION_HH
#define DUMUX_DISCRETIZATION_PQ1BUBBLE_FE_GRID_DISCRETIZATION_HH

#include <dumux/discretization/method.hh>
#include <dumux/common/defaultmappertraits.hh>
#include <dumux/discretization/cvfe/quadraturerules.hh>
#include <dumux/discretization/fem/feelementdiscretization.hh>
#include <dumux/discretization/fem/fegriddiscretization.hh>
#include <dumux/discretization/pq1bubble/geometryhelper.hh>
#include <dumux/discretization/pq1bubble/fvgridgeometry.hh>
#include <dumux/discretization/pq1bubble/pq1bubblefecache.hh>

namespace Dumux::Experimental {

/*!
 * \ingroup PQ1BubbleDiscretization
 * \brief Quadrature rule traits for PQ1Bubble FE discretization
 */
template<class GridView,
         class ElementRule = Dumux::QuadratureRules::DuneQuadrature<4>,
         class IntersectionRule = Dumux::QuadratureRules::MidpointQuadrature,
         class BoundaryFaceRule = Dumux::QuadratureRules::MidpointQuadrature>
struct PQ1BubbleFEQuadratureTraits
{
    using ElementQuadratureRule = ElementRule;
    using IntersectionQuadratureRule = IntersectionRule;
    using BoundaryFaceQuadratureRule = BoundaryFaceRule;
};

/*!
 * \ingroup PQ1BubbleDiscretization
 * \brief The default traits for the PQ1Bubble finite element grid discretization
 * \tparam GridView the grid view type
 * \tparam Scalar the scalar type (used for the FE cache)
 */
template<class GridView,
         class Scalar,
         class MapperTraits = Dumux::PQ1BubbleMapperTraits<GridView, 2>,
         class QuadratureTraits = PQ1BubbleFEQuadratureTraits<GridView>>
struct PQ1BubbleFEDefaultGridDiscretizationTraits
: public MapperTraits, public QuadratureTraits
{
    using DiscretizationMethod = DiscretizationMethods::PQ1Bubble;
    using FeCache = Dumux::PQ1BubbleFECache<typename GridView::ctype, Scalar, GridView::dimension, MapperTraits::numCubeBubbleDofs>;
    using DofHelper = Dumux::PQ1BubbleDofHelper<GridView, MapperTraits::numCubeBubbleDofs>;

    template<class GridDiscretization, bool enableCache>
    using LocalView = FEElementDiscretization<GridDiscretization, enableCache>;

    static constexpr std::size_t maxNumElementDofs = (1 << GridView::dimension)
                                                    + MapperTraits::numCubeBubbleDofs;
};

/*!
 * \ingroup PQ1BubbleDiscretization
 * \brief Finite element grid discretization for the PQ1Bubble method.
 */
template<class Scalar,
         class GV,
         bool enableCaching = true,
         class Traits = PQ1BubbleFEDefaultGridDiscretizationTraits<GV, Scalar>>
using PQ1BubbleFEGridDiscretization = FEGridDiscretization<GV, Traits>;

} // end namespace Dumux::Experimental

#endif
