// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup PoromechanicsTests
 * \brief The flow sub-problem of Mandel's problem
 */
#ifndef DUMUX_TEST_MANDEL_ONEP_SUB_PROBLEM_HH
#define DUMUX_TEST_MANDEL_ONEP_SUB_PROBLEM_HH

#include <memory>
#include <string>

#include <dumux/common/properties.hh>
#include <dumux/common/parameters.hh>
#include <dumux/common/boundarytypes.hh>
#include <dumux/common/numeqvector.hh>
#include <dumux/porousmediumflow/problem.hh>

namespace Dumux {

/*!
 * \ingroup PoromechanicsTests
 * \brief The flow sub-problem of Mandel's problem
 *
 * The right boundary is drained, all other boundaries are impermeable.
 */
template <class TypeTag>
class OnePSubProblem : public PorousMediumFlowProblem<TypeTag>
{
    using ParentType = PorousMediumFlowProblem<TypeTag>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using GlobalPosition = typename GridGeometry::GridView::template Codim<0>::Entity::Geometry::GlobalCoordinate;
    using PrimaryVariables = GetPropType<TypeTag, Properties::PrimaryVariables>;
    using NumEqVector = Dumux::NumEqVector<PrimaryVariables>;
    using BoundaryTypes = Dumux::BoundaryTypes<GetPropType<TypeTag, Properties::ModelTraits>::numEq()>;

public:
    OnePSubProblem(std::shared_ptr<const GridGeometry> gridGeometry,
                   std::shared_ptr<GetPropType<TypeTag, Properties::SpatialParams>> spatialParams,
                   const std::string& paramGroup = "OneP")
    : ParentType(gridGeometry, spatialParams, paramGroup)
    {
        problemName_ = getParam<std::string>("Vtk.OutputName") + "_"
            + getParamFromGroup<std::string>(this->paramGroup(), "Problem.Name");
    }

    const std::string& name() const
    { return problemName_; }

    BoundaryTypes boundaryTypesAtPos(const GlobalPosition& globalPos) const
    {
        BoundaryTypes values;
        if (globalPos[0] > this->gridGeometry().bBoxMax()[0] - eps_)
            values.setAllDirichlet();
        else
            values.setAllNeumann();
        return values;
    }

    PrimaryVariables dirichletAtPos(const GlobalPosition& globalPos) const
    { return PrimaryVariables(0.0); }

    PrimaryVariables initialAtPos(const GlobalPosition& globalPos) const
    { return PrimaryVariables(this->spatialParams().analyticalSolution().initialPressure(globalPos)); }

    NumEqVector sourceAtPos(const GlobalPosition& globalPos) const
    { return NumEqVector(0.0); }

private:
    static constexpr Scalar eps_ = 1e-6;
    std::string problemName_;
};

} // end namespace Dumux

#endif
