// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup TwoPTests
 * \brief McWhorter-Sunada test problem for incompressible immiscible two-phase flow.
 */
#ifndef DUMUX_TEST_TWOP_MCWHORTERSUNADA_PROBLEM_HH
#define DUMUX_TEST_TWOP_MCWHORTERSUNADA_PROBLEM_HH

#include <dumux/common/boundarytypes.hh>
#include <dumux/common/numeqvector.hh>
#include <dumux/common/parameters.hh>
#include <dumux/common/properties.hh>

#include <dumux/porousmediumflow/problem.hh>

namespace Dumux {

template<class TypeTag>
class McWhorterSunadaProblem : public PorousMediumFlowProblem<TypeTag>
{
    using ParentType = PorousMediumFlowProblem<TypeTag>;
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using GridView = typename GridGeometry::GridView;
    using Element = typename GridView::template Codim<0>::Entity;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using PrimaryVariables = GetPropType<TypeTag, Properties::PrimaryVariables>;
    using NumEqVector = Dumux::NumEqVector<PrimaryVariables>;
    using BoundaryTypes = Dumux::BoundaryTypes<GetPropType<TypeTag, Properties::ModelTraits>::numEq()>;
    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;
    using Indices = typename GetPropType<TypeTag, Properties::ModelTraits>::Indices;

public:
    McWhorterSunadaProblem(std::shared_ptr<const GridGeometry> gridGeometry)
    : ParentType(gridGeometry)
    , referencePressure_(getParam<Scalar>("Problem.ReferencePressure"))
    , injectionPressureNw_(getParam<Scalar>("Problem.InjectionPressure"))
    {}

    BoundaryTypes boundaryTypesAtPos(const GlobalPosition& globalPos) const
    {
        BoundaryTypes values;

        if (onLeftBoundary_(globalPos))
            values.setAllDirichlet();
        else
            values.setAllNeumann();

        return values;
    }

    PrimaryVariables dirichletAtPos(const GlobalPosition& globalPos) const
    {
        PrimaryVariables values;
        const auto snr = residualNonwettingSaturation_(globalPos);
        const auto swInj = 1.0 - snr;
        const auto pc = this->spatialParams().fluidMatrixInteractionAtPos(globalPos).pc(swInj);
        values[Indices::saturationIdx] = snr;
        values[Indices::pressureIdx] = injectionPressureNw_ - pc;
        return values;
    }

    NumEqVector neumannAtPos(const GlobalPosition& globalPos) const
    {
        NumEqVector values(0.0);
        return values;
    }

    PrimaryVariables initialAtPos(const GlobalPosition& globalPos) const
    {
        PrimaryVariables values;
        const auto swr = residualWettingSaturation_(globalPos);
        const auto snInit = 1.0 - swr;
        const auto pc = this->spatialParams().fluidMatrixInteractionAtPos(globalPos).pc(swr);
        values[Indices::saturationIdx] = snInit;
        values[Indices::pressureIdx] = injectionPressureNw_ - pc;
        return values;
    }

    Scalar referencePressure() const
    { return referencePressure_; }

private:
    Scalar residualWettingSaturation_(const GlobalPosition& globalPos) const
    { return this->spatialParams().fluidMatrixInteractionAtPos(globalPos).pcSwCurve().effToAbsParams().swr(); }

    Scalar residualNonwettingSaturation_(const GlobalPosition& globalPos) const
    { return this->spatialParams().fluidMatrixInteractionAtPos(globalPos).pcSwCurve().effToAbsParams().snr(); }

    bool onLeftBoundary_(const GlobalPosition& globalPos) const
    { return globalPos[0] < this->gridGeometry().bBoxMin()[0] + eps_; }

    Scalar referencePressure_, injectionPressureNw_;

    static constexpr Scalar eps_ = 1e-6;
};

} // end namespace Dumux

#endif
