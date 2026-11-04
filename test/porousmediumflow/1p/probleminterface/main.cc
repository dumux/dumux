// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup OnePTests
 * \brief Boundary fluxes of problems derived from Experimental::Problem
 *
 * Boundary fluxes are defined with `boundaryFlux` or `boundaryFluxAtPos`. A problem that defines them only
 * with `neumannAtPos` (`INTERFACE` 1) or `neumann` (2), as problems derived from FVProblem do, fails to compile,
 * unless it also defines `boundaryFluxAtPos` (3).
 */
#include <config.h>

#include <memory>

#include <dune/common/exceptions.hh>
#include <dune/grid/yaspgrid.hh>

#include <dumux/common/initialize.hh>
#include <dumux/common/properties.hh>
#include <dumux/common/problem.hh>
#include <dumux/discretization/box.hh>
#include <dumux/porousmediumflow/1p/model.hh>

namespace Dumux {

template<class TypeTag>
class OnePProblemInterfaceProblem : public Experimental::Problem<TypeTag>
{
    using ParentType = Experimental::Problem<TypeTag>;
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using GlobalPosition = typename GridGeometry::GridView::template Codim<0>::Geometry::GlobalCoordinate;
    using PrimaryVariables = GetPropType<TypeTag, Properties::PrimaryVariables>;
    using NumEqVector = Dumux::NumEqVector<PrimaryVariables>;

public:
    using ParentType::ParentType;

#if INTERFACE == 1 || INTERFACE == 3
    NumEqVector neumannAtPos(const GlobalPosition& globalPos) const
    { return NumEqVector(1.0); }
#elif INTERFACE == 2
    using Element = typename GridGeometry::GridView::template Codim<0>::Entity;
    using FVElementGeometry = typename GridGeometry::LocalView;
    using SubControlVolumeFace = typename GridGeometry::SubControlVolumeFace;

    template<class ElementVolumeVariables, class ElementFluxVariablesCache>
    NumEqVector neumann(const Element& element,
                        const FVElementGeometry& fvGeometry,
                        const ElementVolumeVariables& elemVolVars,
                        const ElementFluxVariablesCache& elemFluxVarsCache,
                        const SubControlVolumeFace& scvf) const
    { return NumEqVector(1.0); }
#endif

#if INTERFACE == 3
    NumEqVector boundaryFluxAtPos(const GlobalPosition& globalPos) const
    { return NumEqVector(2.0); }
#endif
};

} // end namespace Dumux

namespace Dumux::Properties {

namespace TTag {
struct OnePProblemInterface { using InheritsFrom = std::tuple<OneP, BoxModel>; };
} // end namespace TTag

template<class TypeTag>
struct Grid<TypeTag, TTag::OnePProblemInterface> { using type = Dune::YaspGrid<2>; };

template<class TypeTag>
struct Problem<TypeTag, TTag::OnePProblemInterface> { using type = OnePProblemInterfaceProblem<TypeTag>; };

} // end namespace Dumux::Properties

int main(int argc, char* argv[])
{
    using namespace Dumux;
    initialize(argc, argv);

    using TypeTag = Properties::TTag::OnePProblemInterface;
    using Grid = GetPropType<TypeTag, Properties::Grid>;
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using Problem = GetPropType<TypeTag, Properties::Problem>;

    Grid grid({1.0, 1.0}, {3, 2});
    const auto problem = std::make_shared<Problem>(std::make_shared<GridGeometry>(grid.leafGridView()));

    // the residual evaluates the boundary fluxes through boundaryFlux, which defaults to this function
    const auto boundaryFlux = problem->boundaryFluxAtPos({0.0, 0.0});
#if INTERFACE == 3
    const double expectedBoundaryFlux = 2.0;
#else
    const double expectedBoundaryFlux = 0.0;
#endif
    if (boundaryFlux[0] != expectedBoundaryFlux)
        DUNE_THROW(Dune::Exception, "boundary flux " << boundaryFlux[0] << ", expected " << expectedBoundaryFlux);

    return 0;
}
