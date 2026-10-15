// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup TwoPOneCTests
 * \brief Horizontal heat pipe inside the two-phase zone for the 2p1cni model.
 */
#ifndef DUMUX_TEST_2P1C_HEATPIPE_PROBLEM_HH
#define DUMUX_TEST_2P1C_HEATPIPE_PROBLEM_HH

#include <dumux/common/properties.hh>
#include <dumux/common/parameters.hh>
#include <dumux/common/boundarytypes.hh>
#include <dumux/common/numeqvector.hh>
#include <dumux/porousmediumflow/problem.hh>

namespace Dumux {

/*!
 * \ingroup TwoPOneCTests
 * \brief Horizontal heat pipe inside the two-phase zone for the 2p1cni model.
 *
 * The cold end at the left boundary has a fixed gas pressure and liquid saturation. The right boundary is
 * closed for mass and receives the heat flux. In the steady state, vapor flows from the heated end to the
 * cold end and the liquid flows back driven by the capillary pressure gradient.
 */
template <class TypeTag>
class TwoPOneCHeatPipeProblem : public PorousMediumFlowProblem<TypeTag>
{
    using ParentType = PorousMediumFlowProblem<TypeTag>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using FluidSystem = GetPropType<TypeTag, Properties::FluidSystem>;
    using ModelTraits = GetPropType<TypeTag, Properties::ModelTraits>;
    using Indices = typename ModelTraits::Indices;
    using PrimaryVariables = GetPropType<TypeTag, Properties::PrimaryVariables>;
    using NumEqVector = Dumux::NumEqVector<PrimaryVariables>;
    using BoundaryTypes = Dumux::BoundaryTypes<ModelTraits::numEq()>;
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using GlobalPosition = typename GridGeometry::GridView::template Codim<0>::Geometry::GlobalCoordinate;

public:
    TwoPOneCHeatPipeProblem(std::shared_ptr<const GridGeometry> gridGeometry)
    : ParentType(gridGeometry)
    {
        FluidSystem::init();
        heatFlux_ = getParam<Scalar>("Problem.HeatFlux");
        coldEndGasPressure_ = getParam<Scalar>("Problem.ColdEndGasPressure");
        coldEndSaturation_ = getParam<Scalar>("Problem.ColdEndSaturation");
    }

    BoundaryTypes boundaryTypesAtPos(const GlobalPosition& globalPos) const
    {
        BoundaryTypes values;
        if (globalPos[0] < this->gridGeometry().bBoxMin()[0] + eps_)
            values.setAllDirichlet();
        else
            values.setAllNeumann();
        return values;
    }

    PrimaryVariables dirichletAtPos(const GlobalPosition& globalPos) const
    {
        PrimaryVariables values(0.0);
        values.setState(Indices::twoPhases);
        values[Indices::pressureIdx] = coldEndGasPressure_;
        values[Indices::switchIdx] = coldEndSaturation_;
        return values;
    }

    NumEqVector neumannAtPos(const GlobalPosition& globalPos) const
    {
        NumEqVector values(0.0);
        if (globalPos[0] > this->gridGeometry().bBoxMax()[0] - eps_)
            values[Indices::energyEqIdx] = -heatFlux_;
        return values;
    }

    PrimaryVariables initialAtPos(const GlobalPosition& globalPos) const
    { return dirichletAtPos(globalPos); }

private:
    static constexpr Scalar eps_ = 1e-6;
    Scalar heatFlux_;
    Scalar coldEndGasPressure_;
    Scalar coldEndSaturation_;
};

} // end namespace Dumux

#endif
