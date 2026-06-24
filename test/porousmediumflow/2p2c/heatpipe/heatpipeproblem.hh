// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux-Lecture contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
#ifndef DUMUX_HEATPIPE_PROBLEM_HH
#define DUMUX_HEATPIPE_PROBLEM_HH

#include <dumux/common/properties.hh>
#include <dumux/common/parameters.hh>
#include <dumux/common/boundarytypes.hh>
#include <dumux/common/numeqvector.hh>
#include <dumux/porousmediumflow/problem.hh>

namespace Dumux
{
template <class TypeTag >
class HeatPipeProblem : public PorousMediumFlowProblem<TypeTag>
{
    using ParentType = PorousMediumFlowProblem<TypeTag>;

    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using FluidSystem =GetPropType<TypeTag, Properties::FluidSystem>;
    using Indices = typename GetPropType<TypeTag, Properties::ModelTraits>::Indices;

    using FVGridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using GridView = typename GetPropType<TypeTag, Properties::GridGeometry>::GridView;
    using Element = typename GridView::template Codim<0>::Entity;
    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;

    using PrimaryVariables = GetPropType<TypeTag, Properties::PrimaryVariables>;
    using NumEqVector = Dumux::NumEqVector<PrimaryVariables>;
    using BoundaryTypes = Dumux::BoundaryTypes<GetPropType<TypeTag, Properties::ModelTraits>::numEq()>;

public:
    HeatPipeProblem(std::shared_ptr<const FVGridGeometry> fVGridGeometry)
    : ParentType(fVGridGeometry)
    {
        FluidSystem::init();
    }

    NumEqVector sourceAtPos(const GlobalPosition &globalPos) const
    {
        return NumEqVector(0.0);
    }

    BoundaryTypes boundaryTypesAtPos(const GlobalPosition &globalPos) const
    {
        BoundaryTypes values;

        if(globalPos[0] < eps_)
            values.setAllDirichlet();
        else
            values.setAllNeumann();

        return values;
    }

    PrimaryVariables dirichletAtPos(const GlobalPosition &globalPos) const
    {
        PrimaryVariables values(0.0);

        // left boundary: two-phase conditions, atmospheric pressure, almost full water saturation
        // 68.6 degree C;
        // we could use another phase presence (onlyWater) on the left, here simplified for better
        // convergence behavior
        values[Indices::pressureIdx] = 1.013e5;
        values[Indices::switchIdx] = 0.99;
        values[Indices::temperatureIdx] = 341.75;
        values.setState(Indices::bothPhases);

        return values;
    }

    NumEqVector neumannAtPos( const GlobalPosition &globalPos) const
    {
        NumEqVector values(0.0);

        // negative values for injection
        // right boundary: constant heat-only flux, like from a cooktop
        if (globalPos[0] > (this->gridGeometry().bBoxMax()[0] - eps_))
        {
            values[Indices::energyEqIdx] = getParam<Scalar>("Problem.HeatFlux");
        }

        return values;
    }

    PrimaryVariables initialAtPos(const GlobalPosition &globalPos) const
    {
      return initial_(globalPos);
    }

private:
    PrimaryVariables initial_(const GlobalPosition &globalPos) const
    {
        PrimaryVariables values;

        // since we are interested in a steady-state,
        // the initial conditions serve mainly to illustrate the way towards the
        // steady state, but they may be changed as desired
        values[Indices::pressureIdx] = 1.013e5;
        values[Indices::switchIdx] = 0.5;
        values[Indices::temperatureIdx] = 343.15;
        values.setState(Indices::bothPhases);

        return values;
    }

    static constexpr Scalar eps_ = 1e-6;
};

} // namespace Dumux

#endif
