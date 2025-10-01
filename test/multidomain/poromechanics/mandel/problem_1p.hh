// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
#ifndef DUMUX_TEST_MD_MANDEL_ONEP_SUB_PROBLEM_HH
#define DUMUX_TEST_MD_MANDEL_ONEP_SUB_PROBLEM_HH

#include <dumux/common/properties.hh>
#include <dumux/common/parameters.hh>
#include <dumux/common/boundarytypes.hh>
#include <dumux/common/numeqvector.hh>

#include <dumux/porousmediumflow/problem.hh>

namespace Dumux {

template <class TypeTag>
class OnePSubProblem : public PorousMediumFlowProblem<TypeTag>
{
    using ParentType = PorousMediumFlowProblem<TypeTag>;

    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using GridView = typename GridGeometry::GridView;
    using Element = typename GridView::template Codim<0>::Entity;
    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;
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

    PrimaryVariables dirichletAtPos(const GlobalPosition& globalPos) const
    { return PrimaryVariables(0.0); }

    PrimaryVariables initialAtPos(const GlobalPosition& globalPos) const
    { return this->spatialParams().analyticalSolution().initialPressure(globalPos); }

    NumEqVector sourceAtPos(const GlobalPosition& globalPos) const
    { return NumEqVector(0.0);}

    BoundaryTypes boundaryTypesAtPos(const GlobalPosition &globalPos) const
    {
        BoundaryTypes values;
        values.setAllNeumann();

        // right boundary
        if (globalPos[0] > this->gridGeometry().bBoxMax()[0]- eps_)
            values.setAllDirichlet();

        return values;
    }

private:
    static constexpr Scalar eps_ = 1.0e-6;
    std::string problemName_;
};

} // end namespace Dumux

#endif
