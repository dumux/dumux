// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup PoromechanicsTests
 * \brief The poro-elastic sub-problem in the el1p coupled problem.
 */

#ifndef DUMUX_TEST_MD_MANDEL_POROELASTIC_SUB_PROBLEM_HH
#define DUMUX_TEST_MD_MANDEL_POROELASTIC_SUB_PROBLEM_HH

#include <dune/common/fmatrix.hh>

#include <dumux/common/properties.hh>
#include <dumux/common/parameters.hh>
#include <dumux/common/boundarytypes.hh>
#include <dumux/common/fvproblemwithspatialparams.hh>

namespace Dumux {

template<class TypeTag>
class PoroElasticSubProblem : public FVProblemWithSpatialParams<TypeTag>
{
    using ParentType = FVProblemWithSpatialParams<TypeTag>;

    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using BoundaryTypes = Dumux::BoundaryTypes<GetPropType<TypeTag, Properties::ModelTraits>::numEq()>;
    using PrimaryVariables = GetPropType<TypeTag, Properties::PrimaryVariables>;

    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using GridView = typename GridGeometry::GridView;
    using Element = typename GridView::template Codim<0>::Entity;
    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;
    using TimeLoopPtr = std::shared_ptr<TimeLoop<Scalar>>;
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

    PrimaryVariables initialAtPos(const GlobalPosition& globalPos) const
    { return this->spatialParams().analyticalSolution().initialDisplacement(globalPos); }

    //! Evaluates the boundary conditions for a Dirichlet boundary segment.
    PrimaryVariables dirichletAtPos(const GlobalPosition& globalPos) const
    {
        PrimaryVariables values(0.0);

        // left boundary
        if (globalPos[0] < this->gridGeometry().bBoxMin()[0] + eps_)
            values[0] = 0.0;

        // bottom boundary
        if (globalPos[1] < this->gridGeometry().bBoxMin()[1] + eps_)
            values[1] = 0.0;

        // top boundary
        if (globalPos[1] > this->gridGeometry().bBoxMax()[1] - eps_)
            values = this->spatialParams().analyticalSolution().displacement(globalPos, timeLoopPtr_->time()  +
            timeLoopPtr_->timeStepSize());

        // right boundary
        if (globalPos[0] > this->gridGeometry().bBoxMax()[0] - eps_)
        values = this->spatialParams().analyticalSolution().displacement(globalPos, timeLoopPtr_->time()  +
        timeLoopPtr_->timeStepSize());
        return values;
    }

    PrimaryVariables neumannAtPos(const GlobalPosition& globalPos) const
    { return PrimaryVariables(0.0); }

    BoundaryTypes boundaryTypesAtPos(const GlobalPosition& globalPos) const
    {
        BoundaryTypes values;
        values.setAllNeumann();

        // left boundary
        if (globalPos[0] < this->gridGeometry().bBoxMin()[0] + eps_)
            values.setDirichlet(0);

        // bottom boundary
        if (globalPos[1] < this->gridGeometry().bBoxMin()[1] + eps_)
            values.setDirichlet(1);

        // top boundary
        if (globalPos[1] > this->gridGeometry().bBoxMax()[1] - eps_)
            values.setDirichlet(1);


        // right boundary
        //if (globalPos[0] > this->gridGeometry().bBoxMax()[0] - eps_)
        //    values.setDirichlet(0);
        return values;
    }

    //! set the time
    void setTimeLoop(TimeLoopPtr timeLoop)
    {
        timeLoopPtr_ = timeLoop;
    }

private:
    TimeLoopPtr timeLoopPtr_;
    static constexpr Scalar eps_ = 3e-6;
    std::string problemName_;
};

} // end namespace Dumux

#endif
