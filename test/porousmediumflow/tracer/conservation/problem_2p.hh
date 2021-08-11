// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   See the file COPYING for full copying permissions.                      *
 *                                                                           *
 *   This program is free software: you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation, either version 3 of the License, or       *
 *   (at your option) any later version.                                     *
 *                                                                           *
 *   This program is distributed in the hope that it will be useful,         *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of          *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the            *
 *   GNU General Public License for more details.                            *
 *                                                                           *
 *   You should have received a copy of the GNU General Public License       *
 *   along with this program.  If not, see <http://www.gnu.org/licenses/>.   *
 *****************************************************************************/
/*!
 * \ingroup TwoPTests
 * \brief The properties for the flow problem of a radioactive 2p test.
 */
#ifndef DUMUX_TWOP_FLOW_TEST_PROBLEM_HH
#define DUMUX_TWOP_FLOW_TEST_PROBLEM_HH

#include <iostream>

#include <dune/grid/yaspgrid.hh>

#include <dumux/common/boundarytypes.hh>
#include <dumux/common/numeqvector.hh>

#include <dumux/discretization/box.hh>
#include <dumux/discretization/cctpfa.hh>

#include <dumux/material/components/simpleh2o.hh>
#include <dumux/material/components/constant.hh>
#include <dumux/material/fluidsystems/h2oair.hh>
#include <dumux/material/fluidsystems/1pliquid.hh>
#include <dumux/material/fluidsystems/1pgas.hh>
#include <dumux/material/fluidsystems/2pimmiscible.hh>

#include <dumux/porousmediumflow/problem.hh>
#include <dumux/porousmediumflow/2p/model.hh>
#include <dumux/porousmediumflow/2p/incompressiblelocalresidual.hh>

#include "spatialparams_2p.hh"

namespace Dumux {
// forward declarations
template<class TypeTag> class TwoPFlowTestProblem;

namespace Properties {
// Create new type tags
namespace TTag {
struct TwoPFlow { using InheritsFrom = std::tuple<TwoP>; };
struct TwoPFlowTpfa { using InheritsFrom = std::tuple<TwoPFlow, CCTpfaModel>; };
} // end namespace TTag

// Set the grid type
template<class TypeTag>
struct Grid<TypeTag, TTag::TwoPFlow> { using type = Dune::YaspGrid<2>; };

// Set the problem type
template<class TypeTag>
struct Problem<TypeTag, TTag::TwoPFlow> { using type = TwoPFlowTestProblem<TypeTag>; };

// Set the fluid system
template<class TypeTag>
struct FluidSystem<TypeTag, TTag::TwoPFlow>
{
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using WettingPhase = FluidSystems::OnePLiquid<Scalar, Components::Constant<0, Scalar> >;
    using NonwettingPhase = FluidSystems::OnePGas<Scalar, Components::Constant<1, Scalar> >;
    using type = FluidSystems::TwoPImmiscible<Scalar, WettingPhase, NonwettingPhase>;
};

// Set the spatial parameters
template<class TypeTag>
struct SpatialParams<TypeTag, TTag::TwoPFlow>
{
private:
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
public:
    using type = TwoPFlowTestSpatialParams<GridGeometry, Scalar>;
};

} // end namespace Properties

/*!
 * \ingroup TwoPTests
 * \brief The properties for the flow problem of a radioactive 2p test.
 */
template<class TypeTag>
class TwoPFlowTestProblem : public PorousMediumFlowProblem<TypeTag>
{
    using ParentType = PorousMediumFlowProblem<TypeTag>;
    using GridView = typename GetPropType<TypeTag, Properties::GridGeometry>::GridView;
    using Element = typename GridView::template Codim<0>::Entity;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using PrimaryVariables = GetPropType<TypeTag, Properties::PrimaryVariables>;
    using NumEqVector = Dumux::NumEqVector<PrimaryVariables>;
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using BoundaryTypes = Dumux::BoundaryTypes<GetPropType<TypeTag, Properties::ModelTraits>::numEq()>;
    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;
    using Indices = typename GetPropType<TypeTag, Properties::ModelTraits>::Indices;

    static constexpr int pressureIdx = Indices::pressureIdx;
    static constexpr int saturationIdx = Indices::saturationIdx;
    static constexpr int conti0EqIdx = Indices::conti0EqIdx;

public:
    TwoPFlowTestProblem(std::shared_ptr<const GridGeometry> gridGeometry)
    : ParentType(gridGeometry) {
         injectionMass_ = getParam<Scalar>("Problem.InjectionMass", -0.0);
    }

    BoundaryTypes boundaryTypesAtPos(const GlobalPosition &globalPos) const
    {
        BoundaryTypes values;

        if (onLowerBoundary_(globalPos))
            values.setAllDirichlet();
        else
            values.setAllNeumann();

        return values;
    }

    NumEqVector neumannAtPos(const GlobalPosition& globalPos) const
    {
        NumEqVector values(0.0);

        if (onUpperBoundary_(globalPos))
        {
            values[conti0EqIdx] = -injectionMass_;
        }

        return values;
    }

    PrimaryVariables dirichletAtPos(const GlobalPosition &globalPos) const
    {
        return initialAtPos(globalPos);
    }

    PrimaryVariables initialAtPos(const GlobalPosition &globalPos) const
    {
        PrimaryVariables values(0.0);

        const Scalar sw = 0.1;

        values[pressureIdx] = 1e5;
        values[saturationIdx] = 1.0 - sw;

        return values;
    }

    Scalar temperature() const
    { return 293.15; }

private:

    bool onLowerBoundary_(const GlobalPosition &globalPos) const
    {
        return globalPos[1] < this->gridGeometry().bBoxMin()[1] + eps_;
    }

    bool onUpperBoundary_(const GlobalPosition &globalPos) const
    {
        return globalPos[1] > this->gridGeometry().bBoxMax()[1] - eps_;
    }

    static constexpr Scalar eps_ = 1e-6;

    Scalar injectionMass_;
};

} // end namespace Dumux

#endif
