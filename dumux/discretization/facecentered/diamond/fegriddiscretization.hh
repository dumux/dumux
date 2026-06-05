// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup PQ1NonconformingDiscretization
 * \brief Finite element grid discretization for the PQ1 non-conforming (face-centered) method.
 */
#ifndef DUMUX_DISCRETIZATION_FACECENTERED_PQ1NONCONFORMING_FE_GRID_DISCRETIZATION_HH
#define DUMUX_DISCRETIZATION_FACECENTERED_PQ1NONCONFORMING_FE_GRID_DISCRETIZATION_HH

#include <dumux/discretization/method.hh>
#include <dumux/common/defaultmappertraits.hh>
#include <dumux/discretization/cvfe/quadraturerules.hh>
#include <dumux/discretization/fem/feelementdiscretization.hh>
#include <dumux/discretization/fem/fegriddiscretization.hh>
#include <dumux/discretization/facecentered/diamond/geometryhelper.hh>
#include <dumux/discretization/nonconformingfecache.hh>

namespace Dumux::Experimental {

/*!
 * \ingroup PQ1NonconformingDiscretization
 * \brief Quadrature rule traits for the PQ1 non-conforming (face-centered) FE discretization
 */
template<class GridView,
         class ElementRule = Dumux::QuadratureRules::DuneQuadrature<2>,
         class IntersectionRule = Dumux::QuadratureRules::DuneQuadrature<2>,
         class BoundaryFaceRule = Dumux::QuadratureRules::DuneQuadrature<2>>
struct PQ1NonconformingFEQuadratureTraits
{
    using ElementQuadratureRule = ElementRule;
    using IntersectionQuadratureRule = IntersectionRule;
    using BoundaryFaceQuadratureRule = BoundaryFaceRule;
};

/*!
 * \ingroup PQ1NonconformingDiscretization
 * \brief Mapper traits for the PQ1 non-conforming (face-centered) FE discretization (face/codim-1 dofs only)
 */
template<class GridView>
struct PQ1NonconformingFEMapperTraits : public Dumux::DefaultMapperTraits<GridView>
{
    using DofMapper = Dune::MultipleCodimMultipleGeomTypeMapper<GridView>;

    //! layout: only faces (codim == 1), matching the PQ1 non-conforming FE spaces
    static Dune::MCMGLayout layout()
    {
        return Dune::mcmgLayout(Dune::Codim<1>{});
    }
};

/*!
 * \ingroup PQ1NonconformingDiscretization
 * \brief The default traits for the PQ1 non-conforming (face-centered) finite element grid discretization
 * \tparam GridView the grid view type
 * \tparam Scalar the scalar type (used for the FE cache)
 */
template<class GridView,
         class Scalar,
         class MapperTraits = PQ1NonconformingFEMapperTraits<GridView>,
         class QuadratureTraits = PQ1NonconformingFEQuadratureTraits<GridView>>
struct PQ1NonconformingFEDefaultGridDiscretizationTraits
: public MapperTraits, public QuadratureTraits
{
    using DiscretizationMethod = DiscretizationMethods::PQ1Nonconforming;
    using FeCache = Dumux::NonconformingFECache<typename GridView::ctype, Scalar, GridView::dimension>;
    using DofHelper = Dumux::DiamondDofHelper<GridView>;

    template<class GridDiscretization, bool enableCache>
    using LocalView = FEElementDiscretization<GridDiscretization, enableCache>;

    // maximum number of dofs per element: 2*dim faces for hypercubes
    static constexpr std::size_t maxNumElementDofs = 2 * GridView::dimension;
};

/*!
 * \ingroup PQ1NonconformingDiscretization
 * \brief Finite element grid discretization for the PQ1 non-conforming (face-centered) method.
 *        Uses Crouzeix-Raviart elements on simplices and Rannacher-Turek elements on cubes.
 */
template<class Scalar,
         class GV,
         bool enableCaching = true,
         class Traits = PQ1NonconformingFEDefaultGridDiscretizationTraits<GV, Scalar>>
using PQ1NonconformingFEGridDiscretization = FEGridDiscretization<GV, Traits>;

} // end namespace Dumux::Experimental

#endif
