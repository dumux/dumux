// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup PQ2Discretization
 * \brief Finite element grid discretization for the PQ2 method.
 */
#ifndef DUMUX_DISCRETIZATION_PQ2_FE_GRID_DISCRETIZATION_HH
#define DUMUX_DISCRETIZATION_PQ2_FE_GRID_DISCRETIZATION_HH

#include <dune/localfunctions/lagrange/lagrangelfecache.hh>

#include <dumux/discretization/method.hh>
#include <dumux/common/defaultmappertraits.hh>
#include <dumux/discretization/cvfe/quadraturerules.hh>
#include <dumux/discretization/fem/fedofhelper.hh>
#include <dumux/discretization/fem/feelementdiscretization.hh>
#include <dumux/discretization/fem/fegriddiscretization.hh>
#include <dumux/discretization/pq2/fvgridgeometry.hh>

#include "dofhelper.hh"

namespace Dumux::Experimental {

/*!
 * \ingroup PQ2Discretization
 * \brief Quadrature rule traits for PQ2 FE discretization
 */
template<class GridView,
         class ElementRule = Dumux::QuadratureRules::DuneQuadrature<4>,
         class IntersectionRule = Dumux::QuadratureRules::DuneQuadrature<4>,
         class BoundaryFaceRule = Dumux::QuadratureRules::DuneQuadrature<4>>
struct PQ2FEQuadratureTraits
{
    using ElementQuadratureRule = ElementRule;
    using IntersectionQuadratureRule = IntersectionRule;
    using BoundaryFaceQuadratureRule = BoundaryFaceRule;
};

/*!
 * \ingroup PQ2Discretization
 * \brief The default traits for the PQ2 finite element grid discretization
 * \tparam GridView the grid view type
 * \tparam Scalar the scalar type (used for the FE cache)
 */
template<class GridView,
         class Scalar,
         class MapperTraits = Dumux::PQ2MapperTraits<GridView>,
         class QuadratureTraits = PQ2FEQuadratureTraits<GridView>>
struct PQ2FEDefaultGridDiscretizationTraits
: public MapperTraits, public QuadratureTraits
{
    using DiscretizationMethod = DiscretizationMethods::PQ2;
    using FeCache = Dune::LagrangeLocalFiniteElementCache<typename GridView::ctype, Scalar, GridView::dimension, 2>;
    using DofHelper = PQ2LagrangeDofHelper<GridView>;

    template<class GridDiscretization, bool enableCache>
    using LocalView = FEElementDiscretization<GridDiscretization, enableCache>;

    // The maximum values correspond to cubes
    static constexpr std::size_t maxNumElementDofs = []()
        {
            if constexpr (GridView::dimension == 1)
                return 3;
            else if constexpr (GridView::dimension == 2)
                return 9;
            else if constexpr (GridView::dimension == 3)
                return 27;
        }();
};

/*!
 * \ingroup PQ2Discretization
 * \brief Finite element grid discretization for the PQ2 method.
 */
template<class Scalar,
         class GV,
         bool enableCaching = true,
         class Traits = PQ2FEDefaultGridDiscretizationTraits<GV, Scalar>>
using PQ2FEGridDiscretization = FEGridDiscretization<GV, Traits>;

} // end namespace Dumux::Experimental

#endif
