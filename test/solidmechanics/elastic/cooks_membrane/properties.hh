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

#include <type_traits>

#include <dune/alugrid/grid.hh>
#include <dune/common/fvector.hh>

#include <dumux/discretization/box.hh>
#include <dumux/discretization/pq1.hh>
#include <dumux/discretization/pq1bubble.hh>
#include <dumux/discretization/pq2.hh>
#include <dumux/discretization/defaultgridvariables.hh>
#include <dumux/discretization/cvfe/variablesadapter.hh>
#include <dumux/solidmechanics/elastic/model.hh>
#include <dumux/solidmechanics/elastic/variables.hh>

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

//! Grid variables for the assembly with local degrees of freedom, with the cache selected depending on the scheme
template<class TypeTag>
struct GridVariables<TypeTag, TTag::CooksMembrane>
{
private:
    static constexpr int dim = GetPropType<TypeTag, Properties::GridGeometry>::GridView::dimension;
    using PV = GetPropType<TypeTag, Properties::PrimaryVariables>;
    using Traits = ElasticVolumeVariablesTraits<PV, Dune::FieldVector<typename PV::value_type, dim>,
                                                GetPropType<TypeTag, Properties::ModelTraits>,
                                                GetPropType<TypeTag, Properties::SolidState>,
                                                GetPropType<TypeTag, Properties::SolidSystem>>;
public:
    using type = Dumux::Experimental::DefaultGridVariables<TypeTag, Dumux::Detail::CVFE::VariablesAdapter<ElasticVariables<Traits>>>;
};

} // end namespace Dumux::Properties

#endif
