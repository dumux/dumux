// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \brief Test lazily evaluated flux fields over sub-control volume faces.
 */
#include <config.h>

#include <cmath>
#include <vector>
#include <memory>
#include <iostream>

#include <dune/common/exceptions.hh>
#include <dune/common/fvector.hh>
#include <dune/grid/yaspgrid.hh>

#include <dumux/common/initialize.hh>
#include <dumux/common/parameters.hh>
#include <dumux/common/properties.hh>
#include <dumux/common/boundarytypes.hh>
#include <dumux/common/numeqvector.hh>
#include <dumux/common/concepts/field_.hh>

#include <dumux/discretization/cctpfa.hh>
#include <dumux/porousmediumflow/1p/model.hh>
#include <dumux/porousmediumflow/problem.hh>
#include <dumux/porousmediumflow/fvspatialparams1p.hh>
#include <dumux/material/components/constant.hh>
#include <dumux/material/fluidsystems/1pliquid.hh>
#include <dumux/flux/cctpfa/darcyslaw.hh>
#include <dumux/flux/tracefields.hh>

namespace {

constexpr double permeability = 1e-11;
constexpr double gradientX = 3.0;
constexpr double gradientY = -1.5;
constexpr double offset = 1e5;

double analyticPressure(const auto& globalPos)
{ return offset + gradientX*globalPos[0] + gradientY*globalPos[1]; }

} // end anonymous namespace

namespace Dumux {

template<class GridGeometry, class Scalar>
class TraceFieldSpatialParams
: public FVPorousMediumFlowSpatialParamsOneP<GridGeometry, Scalar, TraceFieldSpatialParams<GridGeometry, Scalar>>
{
    using ParentType = FVPorousMediumFlowSpatialParamsOneP<GridGeometry, Scalar, TraceFieldSpatialParams<GridGeometry, Scalar>>;
    using GlobalPosition = typename GridGeometry::GridView::template Codim<0>::Entity::Geometry::GlobalCoordinate;

public:
    using PermeabilityType = Scalar;

    explicit TraceFieldSpatialParams(std::shared_ptr<const GridGeometry> gridGeometry)
    : ParentType(gridGeometry) {}

    PermeabilityType permeabilityAtPos(const GlobalPosition&) const
    { return permeability; }

    Scalar porosityAtPos(const GlobalPosition&) const
    { return 0.4; }
};

template<class TypeTag>
class TraceFieldProblem : public PorousMediumFlowProblem<TypeTag>
{
    using ParentType = PorousMediumFlowProblem<TypeTag>;
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using PrimaryVariables = GetPropType<TypeTag, Properties::PrimaryVariables>;
    using BoundaryTypes = Dumux::BoundaryTypes<GetPropType<TypeTag, Properties::ModelTraits>::numEq()>;
    using GlobalPosition = typename GridGeometry::GridView::template Codim<0>::Entity::Geometry::GlobalCoordinate;

public:
    using ParentType::ParentType;

    BoundaryTypes boundaryTypesAtPos(const GlobalPosition&) const
    {
        BoundaryTypes values;
        values.setAllDirichlet();
        return values;
    }

    PrimaryVariables dirichletAtPos(const GlobalPosition& globalPos) const
    { return PrimaryVariables(analyticPressure(globalPos)); }

    Scalar temperatureAtPos(const GlobalPosition&) const
    { return 293.15; }
};

namespace Properties {
namespace TTag {
struct TraceFieldTest { using InheritsFrom = std::tuple<OneP, CCTpfaModel>; };
} // end namespace TTag

template<class TypeTag>
struct Grid<TypeTag, TTag::TraceFieldTest> { using type = Dune::YaspGrid<2>; };

template<class TypeTag>
struct Problem<TypeTag, TTag::TraceFieldTest> { using type = TraceFieldProblem<TypeTag>; };

template<class TypeTag>
struct SpatialParams<TypeTag, TTag::TraceFieldTest>
{
    using type = TraceFieldSpatialParams<GetPropType<TypeTag, Properties::GridGeometry>,
                                         GetPropType<TypeTag, Properties::Scalar>>;
};

template<class TypeTag>
struct FluidSystem<TypeTag, TTag::TraceFieldTest>
{
    using type = FluidSystems::OnePLiquid<GetPropType<TypeTag, Properties::Scalar>,
                                          Components::Constant<1, GetPropType<TypeTag, Properties::Scalar>>>;
};
} // end namespace Properties
} // end namespace Dumux

int main(int argc, char** argv)
{
    using namespace Dumux;
    initialize(argc, argv);

    Parameters::init([] (auto& params) {
        params["Problem.EnableGravity"] = "false";
        params["Component.LiquidDensity"] = "1000";
        params["Component.LiquidKinematicViscosity"] = "1e-6";
    });

    using TypeTag = Properties::TTag::TraceFieldTest;
    using Grid = GetPropType<TypeTag, Properties::Grid>;
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using Problem = GetPropType<TypeTag, Properties::Problem>;
    using GridVariables = GetPropType<TypeTag, Properties::GridVariables>;
    using SolutionVector = GetPropType<TypeTag, Properties::SolutionVector>;
    using DarcysLaw = GetPropType<TypeTag, Properties::AdvectionType>;

    Grid grid{{1.0, 1.0}, {4, 4}};
    auto gridGeometry = std::make_shared<GridGeometry>(grid.leafGridView());
    auto problem = std::make_shared<Problem>(gridGeometry);

    auto x = std::make_shared<SolutionVector>(gridGeometry->numDofs());
    for (const auto& element : elements(gridGeometry->gridView()))
    {
        const auto eIdx = gridGeometry->elementMapper().index(element);
        (*x)[eIdx] = analyticPressure(element.geometry().center());
    }

    auto gridVariables = std::make_shared<GridVariables>(problem, gridGeometry);
    gridVariables->init(*x);

    const auto darcyInvoker = [] (const auto& context, const auto& scvf)
    {
        return DarcysLaw::flux(
            context.problem(), context.element(), context.fvGeometry(),
            context.elemVolVars(), scvf, 0, context.elemFluxVarsCache()
        );
    };

    const auto darcy = makeLazyFluxField(
        std::shared_ptr<const Problem>(problem),
        std::shared_ptr<const GridVariables>(gridVariables),
        std::shared_ptr<const SolutionVector>(x),
        darcyInvoker
    );

    using Element = typename GridGeometry::GridView::template Codim<0>::Entity;
    using Scvf = typename GridGeometry::SubControlVolumeFace;
    static_assert(Concept::BindableField<decltype(darcy), Element>);
    static_assert(Concept::FaceField<decltype(darcy), Element, Scvf>);
    static_assert(!Concept::BindableField<Scvf, Element>);
    static_assert(!Concept::FaceField<Scvf, Element, Scvf>);

    // -----------------------------------------------------------------
    // the flux of a linear pressure field is exact: -A n.K grad(p)
    // -----------------------------------------------------------------
    const Dune::FieldVector<double, 2> gradient{gradientX, gradientY};
    double maxError = 0.0;
    std::size_t checked = 0;

    for (const auto& element : elements(gridGeometry->gridView()))
    {
        const auto bound = darcy.bind(element);
        const auto fvGeometry = localView(*gridGeometry).bind(element);
        for (const auto& scvf : scvfs(fvGeometry))
        {
            const auto value = bound(scvf);

            const auto expected = -permeability*(gradient*scvf.unitOuterNormal())*scvf.area();
            using std::abs;
            using std::max;
            maxError = max(maxError, abs(value - expected));
            ++checked;
        }
    }

    if (checked == 0)
        DUNE_THROW(Dune::InvalidStateException, "No faces visited");
    if (maxError > 1e-18)
        DUNE_THROW(Dune::MathError, "Darcy flux deviates from the analytic value by " << maxError);
    std::cout << "Lazy Darcy flux against the analytic value on " << checked
              << " faces: max error = " << maxError << std::endl;

    // -----------------------------------------------------------------
    // a field that does not own what it reads gives the same values
    // -----------------------------------------------------------------
    {
        const auto nonOwning = makeNonOwningFluxField(*problem, *gridVariables, *x, darcyInvoker);
        static_assert(Concept::FaceField<decltype(nonOwning), Element, Scvf>);

        double diff = 0.0;
        for (const auto& element : elements(gridGeometry->gridView()))
        {
            const auto lazyBound = darcy.bind(element);
            const auto nonOwningBound = nonOwning.bind(element);
            const auto fvGeometry = localView(*gridGeometry).bind(element);
            for (const auto& scvf : scvfs(fvGeometry))
            {
                using std::abs;
                using std::max;
                diff = max(diff, abs(lazyBound(scvf) - nonOwningBound(scvf)));
            }
        }
        if (diff != 0.0)
            DUNE_THROW(Dune::MathError, "Non-owning and owning fields differ by " << diff);
        std::cout << "Non-owning field agrees with the owning field exactly" << std::endl;
    }

    // -----------------------------------------------------------------
    // a bound field outlives the field it came from
    // -----------------------------------------------------------------
    {
        const auto element = *elements(gridGeometry->gridView()).begin();
        const auto bound = makeLazyFluxField(
            std::shared_ptr<const Problem>(problem),
            std::shared_ptr<const GridVariables>(gridVariables),
            std::shared_ptr<const SolutionVector>(x),
            darcyInvoker
        ).bind(element);
        // the field was a temporary and is gone; the bound object must still be valid

        const auto fvGeometry = localView(*gridGeometry).bind(element);
        for (const auto& scvf : scvfs(fvGeometry))
        {
            const auto expected = -permeability*(gradient*scvf.unitOuterNormal())*scvf.area();
            using std::abs;
            if (abs(bound(scvf) - expected) > 1e-18)
                DUNE_THROW(Dune::MathError, "Field bound from a temporary gives a wrong value");
        }
        std::cout << "A field bound from temporaries stays valid after the full expression ends" << std::endl;
    }

    std::cout << "All trace field checks passed" << std::endl;
    return 0;
}
