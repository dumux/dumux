// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup GeomechanicsTests
 * \brief Properties of Cook's membrane, assembled with local degrees of freedom
 */
#ifndef DUMUX_TEST_ELASTIC_COOKS_MEMBRANE_PROPERTIES_HH
#define DUMUX_TEST_ELASTIC_COOKS_MEMBRANE_PROPERTIES_HH

#include <dune/alugrid/grid.hh>

#include <dumux/discretization/box.hh>
#include <dumux/discretization/pq1.hh>
#include <dumux/discretization/pq1bubble.hh>
#include <dumux/discretization/pq2.hh>
#include <dumux/solidmechanics/elastic/model.hh>

#include "spatialparams.hh"
#include "problem.hh"

namespace Dumux::Properties {

namespace TTag {
struct CooksMembrane
{
    using InheritsFrom = std::tuple<Elastic, DISCRETIZATION_MODEL>;
    using Grid = Dune::ALUGrid<2, 2, Dune::simplex, Dune::conforming>;
    using EnableGridGeometryCache = std::true_type;
};
} // end namespace TTag

template<class TypeTag>
struct Problem<TypeTag, TTag::CooksMembrane>
{ using type = CooksMembraneProblem<TypeTag>; };

template<class TypeTag>
struct SpatialParams<TypeTag, TTag::CooksMembrane>
{
    using type = CooksMembraneSpatialParams<GetPropType<TypeTag, Properties::GridGeometry>,
                                            GetPropType<TypeTag, Properties::Scalar>>;
};

//! Assemble with the variables per local dof
template<class TypeTag>
struct Transition::UseVariables<TypeTag, TTag::CooksMembrane> { static constexpr bool value = true; };

} // end namespace Dumux::Properties

#endif
