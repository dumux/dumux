// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \brief Type tags of the complex-valued facet-coupled Helmholtz test
 */
#ifndef DUMUX_TEST_COMPLEX_FACET_HELMHOLTZ_PROPERTIES_HH
#define DUMUX_TEST_COMPLEX_FACET_HELMHOLTZ_PROPERTIES_HH

#include <type_traits>

#include <dune/alugrid/grid.hh>
#include <dune/foamgrid/foamgrid.hh>

#include <dumux/common/properties.hh>
#include <dumux/discretization/box.hh>
#include <dumux/multidomain/traits.hh>
#include <dumux/multidomain/facet/box/properties.hh>
#include <dumux/multidomain/facet/couplingmapper.hh>
#include <dumux/multidomain/facet/couplingmanager.hh>

#include "model.hh"
#include "problems.hh"

namespace Dumux::Properties {

namespace TTag {
struct FacetHelmholtzBulk { using InheritsFrom = std::tuple<FacetHelmholtzModel, BoxFacetCouplingModel>; };
struct FacetHelmholtzLowDim { using InheritsFrom = std::tuple<FacetHelmholtzModel, BoxModel>; };
} // end namespace TTag

template<class TypeTag>
struct Scalar<TypeTag, TTag::FacetHelmholtzBulk> { using type = double; };
template<class TypeTag>
struct Scalar<TypeTag, TTag::FacetHelmholtzLowDim> { using type = double; };

template<class TypeTag>
struct Grid<TypeTag, TTag::FacetHelmholtzBulk> { using type = Dune::ALUGrid<2, 2, Dune::cube, Dune::nonconforming>; };
template<class TypeTag>
struct Grid<TypeTag, TTag::FacetHelmholtzLowDim> { using type = Dune::FoamGrid<1, 2>; };

template<class TypeTag>
struct Problem<TypeTag, TTag::FacetHelmholtzBulk> { using type = FacetHelmholtzBulkProblem<TypeTag>; };
template<class TypeTag>
struct Problem<TypeTag, TTag::FacetHelmholtzLowDim> { using type = FacetHelmholtzLowDimProblem<TypeTag>; };

template<class BulkTypeTag, class LowDimTypeTag>
struct FacetHelmholtzTestTraits
{
    using MDTraits = Dumux::MultiDomainTraits<BulkTypeTag, LowDimTypeTag>;
    using CouplingMapper = Dumux::FacetCouplingMapper<GetPropType<BulkTypeTag, Properties::GridGeometry>,
                                                      GetPropType<LowDimTypeTag, Properties::GridGeometry>>;
    using CouplingManager = Dumux::FacetCouplingManager<MDTraits, CouplingMapper>;
};

using FacetHelmholtzTraits = FacetHelmholtzTestTraits<TTag::FacetHelmholtzBulk, TTag::FacetHelmholtzLowDim>;

template<class TypeTag>
struct CouplingManager<TypeTag, TTag::FacetHelmholtzBulk> { using type = typename FacetHelmholtzTraits::CouplingManager; };
template<class TypeTag>
struct CouplingManager<TypeTag, TTag::FacetHelmholtzLowDim> { using type = typename FacetHelmholtzTraits::CouplingManager; };

} // end namespace Dumux::Properties

#endif
