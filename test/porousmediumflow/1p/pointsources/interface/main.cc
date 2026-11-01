// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup OnePTests
 * \brief Point sources of problems derived from Experimental::Problem
 *
 * Such problems return their point sources from `pointSources()`. A problem that adds point sources
 * with `addPointSources`, as problems derived from FVProblem do (`ADDPOINTSOURCES` 1 for a function,
 * 2 for a function template), fails to compile.
 */
#include <config.h>

#include <vector>

#include <dune/grid/yaspgrid.hh>

#include <dumux/common/properties.hh>
#include <dumux/common/problem.hh>
#include <dumux/discretization/box.hh>
#include <dumux/porousmediumflow/1p/model.hh>

namespace Dumux {

template<class TypeTag>
class OnePPointSourcesInterfaceProblem : public Experimental::Problem<TypeTag>
{
    using ParentType = Experimental::Problem<TypeTag>;

public:
    using ParentType::ParentType;

#if ADDPOINTSOURCES == 1
    using PointSource = GetPropType<TypeTag, Properties::PointSource>;
    void addPointSources(std::vector<PointSource>& pointSources) const {}
#elif ADDPOINTSOURCES == 2
    template<class PointSource>
    void addPointSources(std::vector<PointSource>& pointSources) const {}
#endif
};

} // end namespace Dumux

namespace Dumux::Properties {

namespace TTag {
struct OnePPointSourcesInterface { using InheritsFrom = std::tuple<OneP, BoxModel>; };
} // end namespace TTag

template<class TypeTag>
struct Grid<TypeTag, TTag::OnePPointSourcesInterface> { using type = Dune::YaspGrid<2>; };

template<class TypeTag>
struct Problem<TypeTag, TTag::OnePPointSourcesInterface> { using type = OnePPointSourcesInterfaceProblem<TypeTag>; };

template<class TypeTag>
struct Transition::UseVariables<TypeTag, TTag::OnePPointSourcesInterface> { static constexpr bool value = true; };

} // end namespace Dumux::Properties

int main()
{
    using namespace Dumux;
    using Problem = GetPropType<Properties::TTag::OnePPointSourcesInterface, Properties::Problem>;

    // the residual evaluates the point sources through this function
    [[maybe_unused]] const auto pointSources = &Problem::pointSources;

    return 0;
}
