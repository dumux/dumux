// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup BoxDiscretization
 * \brief Finite element grid discretization for the PQ1 (box) FE method.
 */
#ifndef DUMUX_DISCRETIZATION_PQ1_FE_GRID_DISCRETIZATION_HH
#define DUMUX_DISCRETIZATION_PQ1_FE_GRID_DISCRETIZATION_HH

#include <dune/localfunctions/lagrange/lagrangelfecache.hh>

#include <dumux/discretization/method.hh>
#include <dumux/common/defaultmappertraits.hh>
#include <dumux/discretization/cvfe/quadraturerules.hh>
#include <dumux/discretization/fem/feelementdiscretization.hh>
#include <dumux/discretization/fem/fegriddiscretization.hh>
#include "boxgeometryhelper.hh"

namespace Dumux::Experimental {

/*!
 * \ingroup BoxDiscretization
 * \brief Quadrature rule traits for PQ1 FE discretization
 */
template<class GridView,
         class ElementRule = Dumux::QuadratureRules::DuneQuadrature<2>,
         class IntersectionRule = Dumux::QuadratureRules::DuneQuadrature<2>,
         class BoundaryFaceRule = Dumux::QuadratureRules::DuneQuadrature<2>>
struct PQ1FEQuadratureTraits
{
    using ElementQuadratureRule = ElementRule;
    using IntersectionQuadratureRule = IntersectionRule;
    using BoundaryFaceQuadratureRule = BoundaryFaceRule;
};

/*!
 * \ingroup BoxDiscretization
 * \brief Mapper traits for PQ1 FE discretization (vertex dofs only)
 */
template<class GridView>
struct PQ1FEMapperTraits : public Dumux::DefaultMapperTraits<GridView>
{
    using DofMapper = Dune::MultipleCodimMultipleGeomTypeMapper<GridView>;

    //! layout: only vertices (codim == dim)
    static Dune::MCMGLayout layout()
    {
        return [](Dune::GeometryType gt, int /*dimgrid*/) {
            return gt.dim() == 0;
        };
    }
};

/*!
 * \ingroup BoxDiscretization
 * \brief The default traits for the PQ1 finite element grid discretization
 * \tparam GridView the grid view type
 * \tparam Scalar the scalar type (used for the FE cache)
 */
template<class GridView,
         class Scalar,
         class MapperTraits = PQ1FEMapperTraits<GridView>,
         class QuadratureTraits = PQ1FEQuadratureTraits<GridView>>
struct PQ1FEDefaultGridDiscretizationTraits
: public MapperTraits, public QuadratureTraits
{
    using DiscretizationMethod = DiscretizationMethods::Box;
    using FeCache = Dune::LagrangeLocalFiniteElementCache<typename GridView::ctype, Scalar, GridView::dimension, 1>;
    using DofHelper = Dumux::BoxDofHelper<GridView>;

    template<class GridDiscretization, bool enableCache>
    using LocalView = FEElementDiscretization<GridDiscretization, enableCache>;

    // maximum number of dofs per element: 2^dim for cubes
    static constexpr std::size_t maxNumElementDofs = (1 << GridView::dimension);
};

/*!
 * \ingroup BoxDiscretization
 * \brief Finite element grid discretization for the PQ1 (box) FE method.
 */
template<class Scalar,
         class GV,
         bool enableCaching = true,
         class Traits = PQ1FEDefaultGridDiscretizationTraits<GV, Scalar>>
using PQ1FEGridDiscretization = FEGridDiscretization<GV, Traits>;

} // end namespace Dumux::Experimental

#endif
