// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup NavierStokesTests
 * \brief Test that the boundary fluxes of the problem do not enter the finite element
 *        continuity residual.
 *
 * A uniform velocity field and a constant density satisfy the continuity equation exactly,
 * so the residual has to vanish although the problem returns the mass flux over all boundaries.
 */
#include <config.h>

#include <algorithm>
#include <cmath>
#include <iostream>
#include <memory>
#include <type_traits>

#include <dune/common/exceptions.hh>
#include <dune/common/fvector.hh>
#include <dune/grid/yaspgrid.hh>

#include <dumux/common/initialize.hh>
#include <dumux/common/parameters.hh>
#include <dumux/common/properties.hh>
#include <dumux/assembly/assembler.hh>
#include <dumux/assembly/boundaryfluxintegral.hh>
#include <dumux/discretization/pq1.hh>
#include <dumux/material/components/constant.hh>
#include <dumux/material/fluidsystems/1pliquid.hh>
#include <dumux/freeflow/navierstokes/mass/1p/model.hh>
#include <dumux/freeflow/navierstokes/mass/problem.hh>

namespace Dumux {

template<class TypeTag>
class UniformFlowProblem : public CVFENavierStokesMassProblem<TypeTag>
{
    using ParentType = CVFENavierStokesMassProblem<TypeTag>;
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using ElementDiscretization = typename GridGeometry::LocalView;
    using GlobalPosition = typename GridGeometry::GridView::template Codim<0>::Entity::Geometry::GlobalCoordinate;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;

public:
    using typename ParentType::BoundaryTypes;
    using typename ParentType::BoundaryFluxes;
    using ParentType::ParentType;

    BoundaryTypes boundaryTypesAtPos(const GlobalPosition& globalPos) const
    {
        BoundaryTypes values;
        values.setAllFluxBoundary();
        return values;
    }

    template<class ElementVariables, class FaceIpData>
    BoundaryFluxes boundaryFlux(const ElementDiscretization& elemDisc,
                                const ElementVariables& elemVars,
                                const FaceIpData& faceIpData) const
    {
        const auto& shapeValues = cache(elemVars, faceIpData).shapeValues();
        Scalar density = 0.0;
        for (const auto& localDof : localDofs(elemDisc))
            density += elemVars[localDof.index()].density()*shapeValues[localDof.index()][0];

        BoundaryFluxes values(0.0);
        values[0] = density*(velocityAtPos(faceIpData.global())*faceIpData.unitOuterNormal());
        return values;
    }

    GlobalPosition velocityAtPos(const GlobalPosition& globalPos) const
    { return {1.0, 0.5}; }
};

namespace Properties {

namespace TTag {
struct UniformFlowFE
{
    using InheritsFrom = std::tuple<NavierStokesMassOnePFE, PQ1FEModel>;
    using EnableGridGeometryCache = std::true_type;
};
} // end namespace TTag

template<class TypeTag>
struct Grid<TypeTag, TTag::UniformFlowFE>
{ using type = Dune::YaspGrid<2>; };

template<class TypeTag>
struct Problem<TypeTag, TTag::UniformFlowFE>
{ using type = UniformFlowProblem<TypeTag>; };

template<class TypeTag>
struct FluidSystem<TypeTag, TTag::UniformFlowFE>
{
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using type = FluidSystems::OnePLiquid<Scalar, Components::Constant<1, Scalar>>;
};

} // end namespace Properties
} // end namespace Dumux

int main(int argc, char** argv)
{
    using namespace Dumux;
    using TypeTag = Properties::TTag::UniformFlowFE;

    initialize(argc, argv);
    Parameters::init([](Dune::ParameterTree& params)
    {
        params["Problem.EnableGravity"] = "false";
        params["Component.LiquidDensity"] = "1000";
        params["Component.LiquidKinematicViscosity"] = "1e-6";
        params["SpatialParams.Temperature"] = "293.15";
    });

    using Grid = GetPropType<TypeTag, Properties::Grid>;
    const Grid grid({1.0, 1.0}, {4, 4});

    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    auto gridGeometry = std::make_shared<GridGeometry>(grid.leafGridView());

    using Problem = GetPropType<TypeTag, Properties::Problem>;
    auto problem = std::make_shared<Problem>(gridGeometry);

    using SolutionVector = GetPropType<TypeTag, Properties::SolutionVector>;
    SolutionVector x(gridGeometry->numDofs());
    x = 1e5;

    using GridVariables = GetPropType<TypeTag, Properties::GridVariables>;
    auto gridVariables = std::make_shared<GridVariables>(problem, gridGeometry);
    gridVariables->init(x);

    using Assembler = Experimental::Assembler<TypeTag, DiffMethod::numeric>;
    Assembler assembler(problem, gridGeometry, gridVariables);
    assembler.assembleResidual(x);

    // what the residual would contain if the boundary fluxes of the problem entered it
    using ElementResidualVector = typename GetPropType<TypeTag, Properties::LocalResidual>::ElementResidualVector;
    SolutionVector boundaryFluxIntegral(gridGeometry->numDofs());
    boundaryFluxIntegral = 0.0;
    auto elemDisc = localView(*gridGeometry);
    for (const auto& element : elements(gridGeometry->gridView()))
    {
        elemDisc.bind(element);
        if (!elemDisc.hasBoundaryFaces())
            continue;

        const auto elemVars = localView(gridVariables->curGridVars()).bind(element, elemDisc, x);
        ElementResidualVector elemResidual(Detail::LocalDofs::numLocalDofs(elemDisc));
        elemResidual = 0.0;
        for (const auto& boundaryFace : boundaryFaces(elemDisc))
            Experimental::addFEBoundaryFluxIntegral(elemResidual, *problem, elemDisc, elemVars,
                                                    boundaryFace, problem->boundaryTypes(elemDisc, boundaryFace));

        for (const auto& localDof : localDofs(elemDisc))
            boundaryFluxIntegral[localDof.dofIndex()] += elemResidual[localDof.index()];
    }

    const auto maxAbs = [](const auto& vector)
    {
        double result = 0.0;
        for (const auto& block : vector)
            result = std::max(result, std::abs(block[0]));
        return result;
    };

    const auto maxBoundaryFluxIntegral = maxAbs(boundaryFluxIntegral);
    const auto maxResidual = maxAbs(assembler.residual());
    std::cout << "max. boundary flux integral: " << maxBoundaryFluxIntegral << "\n"
              << "max. residual: " << maxResidual << std::endl;

    if (maxBoundaryFluxIntegral < 1.0)
        DUNE_THROW(Dune::Exception, "The boundary flux integral of the problem vanishes");

    if (maxResidual > 1e-12*maxBoundaryFluxIntegral)
        DUNE_THROW(Dune::Exception, "The finite element continuity residual contains boundary fluxes");

    return 0;
}
