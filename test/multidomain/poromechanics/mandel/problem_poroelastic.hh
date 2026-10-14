// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup PoromechanicsTests
 * \brief The poroelastic sub-problem of Mandel's problem
 */
#ifndef DUMUX_TEST_MANDEL_POROELASTIC_SUB_PROBLEM_HH
#define DUMUX_TEST_MANDEL_POROELASTIC_SUB_PROBLEM_HH

#include <memory>
#include <string>

#include <dumux/common/properties.hh>
#include <dumux/common/parameters.hh>
#include <dumux/common/boundarytypes.hh>
#include <dumux/common/fvproblemwithspatialparams.hh>

namespace Dumux {

/*!
 * \ingroup PoromechanicsTests
 * \brief The poroelastic sub-problem of Mandel's problem
 *
 * The rigid plate on top prescribes the vertical displacement of the analytical solution.
 * The left and bottom boundaries are symmetry planes and the right boundary is traction-free.
 */
template<class TypeTag>
class PoroElasticSubProblem : public FVProblemWithSpatialParams<TypeTag>
{
    using ParentType = FVProblemWithSpatialParams<TypeTag>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using BoundaryTypes = Dumux::BoundaryTypes<GetPropType<TypeTag, Properties::ModelTraits>::numEq()>;
    using PrimaryVariables = GetPropType<TypeTag, Properties::PrimaryVariables>;
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using GlobalPosition = typename GridGeometry::GridView::template Codim<0>::Entity::Geometry::GlobalCoordinate;

public:
    PoroElasticSubProblem(std::shared_ptr<const GridGeometry> gridGeometry,
                          std::shared_ptr<GetPropType<TypeTag, Properties::SpatialParams>> spatialParams,
                          const std::string& paramGroup = "PoroElastic")
    : ParentType(gridGeometry, spatialParams, paramGroup)
    {
        problemName_ = getParam<std::string>("Vtk.OutputName") + "_"
            + getParamFromGroup<std::string>(this->paramGroup(), "Problem.Name");
    }

    const std::string& name() const
    { return problemName_; }

    //! Set the time at which the boundary conditions are evaluated
    void setTime(Scalar time) const
    { time_ = time; }

    PrimaryVariables initialAtPos(const GlobalPosition& globalPos) const
    { return this->spatialParams().analyticalSolution().initialDisplacement(globalPos); }

    BoundaryTypes boundaryTypesAtPos(const GlobalPosition& globalPos) const
    {
        BoundaryTypes values;
        values.setAllNeumann();

        if (globalPos[0] < this->gridGeometry().bBoxMin()[0] + eps_)
            values.setDirichlet(0);

        if (globalPos[1] < this->gridGeometry().bBoxMin()[1] + eps_
            || globalPos[1] > this->gridGeometry().bBoxMax()[1] - eps_)
            values.setDirichlet(1);

        return values;
    }

    PrimaryVariables dirichletAtPos(const GlobalPosition& globalPos) const
    {
        PrimaryVariables values(0.0);
        if (globalPos[1] > this->gridGeometry().bBoxMax()[1] - eps_)
            values[1] = this->spatialParams().analyticalSolution().displacement(globalPos, time_)[1];
        return values;
    }

    PrimaryVariables neumannAtPos(const GlobalPosition& globalPos) const
    { return PrimaryVariables(0.0); }

private:
    static constexpr Scalar eps_ = 1e-6;
    std::string problemName_;
    mutable Scalar time_ = 0.0;
};

} // end namespace Dumux

#endif
