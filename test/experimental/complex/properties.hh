// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \brief Type tags for the complex-valued Helmholtz tests
 */
#ifndef DUMUX_TEST_COMPLEX_HELMHOLTZ_PROPERTIES_HH
#define DUMUX_TEST_COMPLEX_HELMHOLTZ_PROPERTIES_HH

#include <type_traits>

#include <dune/grid/yaspgrid.hh>

#include <dumux/common/properties.hh>
#include <dumux/discretization/box.hh>
#include <dumux/discretization/cctpfa.hh>

#include "model.hh"
#include "problem.hh"

namespace Dumux::Properties::TTag {

struct ComplexHelmholtzTest
{
    using InheritsFrom = std::tuple<ComplexHelmholtzModel>;

    using Scalar = double;
    using Grid = Dune::YaspGrid<2>;

    template<class TypeTag>
    using Problem = ComplexHelmholtzManufacturedProblem<TypeTag>;

    using EnableGridVolumeVariablesCache = std::true_type;
    using EnableGridFluxVariablesCache = std::true_type;
    using EnableGridGeometryCache = std::true_type;
};

struct ComplexHelmholtzBox { using InheritsFrom = std::tuple<ComplexHelmholtzTest, BoxModel>; };
struct ComplexHelmholtzTpfa { using InheritsFrom = std::tuple<ComplexHelmholtzTest, CCTpfaModel>; };

struct ComplexHelmholtzSweep
{
    using InheritsFrom = std::tuple<ComplexHelmholtzTest, BoxModel>;

    template<class TypeTag>
    using Problem = ComplexHelmholtzSweepProblem<TypeTag>;
};

} // end namespace Dumux::Properties::TTag

#endif
