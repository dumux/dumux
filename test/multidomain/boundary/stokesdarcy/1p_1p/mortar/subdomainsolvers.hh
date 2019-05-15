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
 * \file
 * \ingroup BoundaryTests
 * \brief Stokes & Darcy solvers for the mortar-coupling stokes darcy problem.
 */
#include <dumux/assembly/fvassembler.hh>
#include <dumux/assembly/staggeredfvassembler.hh>
#include <dumux/assembly/diffmethod.hh>

#include <dumux/common/parameters.hh>
#include <dumux/common/properties.hh>
#include <dumux/common/timeloop.hh>

#include <dumux/io/grid/gridmanager.hh>
#include <dumux/io/vtkoutputmodule.hh>
#include <dumux/io/staggeredvtkoutputmodule.hh>

#include <dumux/linear/seqsolverbackend.hh>
#include <dumux/nonlinear/newtonsolver.hh>

namespace Dumux {

template<class TypeTag, class LinearSolver = Dumux::UMFPackBackend>
class StokesSolver
{
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using FVGridGeometry = GetPropType<TypeTag, Properties::FVGridGeometry>;
    using Problem = GetPropType<TypeTag, Properties::Problem>;
    using SolutionVector = GetPropType<TypeTag, Properties::SolutionVector>;
    using GridVariables = GetPropType<TypeTag, Properties::GridVariables>;
    using OutputModule = StaggeredVtkOutputModule<GridVariables, SolutionVector>;
    using Assembler = StaggeredFVAssembler<TypeTag, DiffMethod::numeric>;
    using NewtonSolver = Dumux::NewtonSolver<Assembler, LinearSolver>;

public:

    void init(const std::string& paramGroup)
    {
        // try to create a grid (from the given grid file or the input file)
        gridManager_.init(paramGroup);

        // we compute on the leaf grid view
        const auto& leafGridView = gridManager_.grid().leafGridView();

        // create the finite volume grid geometry
        fvGridGeometry_ = std::make_shared<FVGridGeometry>(leafGridView);
        fvGridGeometry_->update();

        // the problem (initial and boundary conditions)
        problem_ = std::make_shared<Problem>(fvGridGeometry_, paramGroup);

        // resize and initialize the given solution vector
        x_ = std::make_shared<SolutionVector>();
        (*x_)[FVGridGeometry::cellCenterIdx()].resize(fvGridGeometry_->numCellCenterDofs());
        (*x_)[FVGridGeometry::faceIdx()].resize(fvGridGeometry_->numFaceDofs());
        problem_->applyInitialSolution(*x_);

        // the grid variables
        gridVariables_ = std::make_shared<GridVariables>(problem_, fvGridGeometry_);
        gridVariables_->init(*x_);

        // initialize the vtk output module
        using IOFields = GetPropType<TypeTag, Properties::IOFields>;
        vtkWriter_ = std::make_unique<OutputModule>(*gridVariables_, *x_, problem_->name());
        IOFields::initOutputModule(*vtkWriter_);

        // the assembler without time loop for stationary problem
        assembler_ = std::make_shared<Assembler>(problem_, fvGridGeometry_, gridVariables_);

        // the linear/non-linear solvers
        auto linearSolver = std::make_shared<LinearSolver>();
        newtonSolver_ = std::make_unique<NewtonSolver>(assembler_, linearSolver);
    }

    //! Solve the system
    void solve()
    {
        newtonSolver_->solve(*x_);
    }

    //! Write current state to disk
    void write(Scalar t)
    {
        vtkWriter_->write(t);
    }

    //! Return a pointer to the grid geometry
    std::shared_ptr<FVGridGeometry> gridGeometryPointer()
    { return fvGridGeometry_; }

    //! Return a pointer to the grid variables
    std::shared_ptr<GridVariables> gridVariablesPointer()
    { return gridVariables_; }

    //! Return a pointer to the problem
    std::shared_ptr<Problem> problemPointer()
    { return problem_; }

    //! Return a pointer to the solution
    std::shared_ptr<SolutionVector> solutionPointer()
    { return x_; }

private:
    GridManager<GetPropType<TypeTag, Properties::Grid>> gridManager_;
    std::shared_ptr<FVGridGeometry> fvGridGeometry_;
    std::shared_ptr<Problem> problem_;
    std::shared_ptr<GridVariables> gridVariables_;
    std::shared_ptr<Assembler> assembler_;

    std::unique_ptr<NewtonSolver> newtonSolver_;
    std::unique_ptr<OutputModule> vtkWriter_;

    std::shared_ptr<SolutionVector> x_;
};



template<class TypeTag, class LinearSolver = Dumux::UMFPackBackend>
class DarcySolver
{
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using FVGridGeometry = GetPropType<TypeTag, Properties::FVGridGeometry>;
    using Problem = GetPropType<TypeTag, Properties::Problem>;
    using SpatialParams = GetPropType<TypeTag, Properties::SpatialParams>;
    using SolutionVector = GetPropType<TypeTag, Properties::SolutionVector>;
    using GridVariables = GetPropType<TypeTag, Properties::GridVariables>;
    using OutputModule = VtkOutputModule<GridVariables, SolutionVector>;
    using Assembler = FVAssembler<TypeTag, DiffMethod::numeric>;
    using NewtonSolver = Dumux::NewtonSolver<Assembler, LinearSolver>;

public:

    void init(const std::string& paramGroup)
    {
        // try to create a grid (from the given grid file or the input file)
        gridManager_.init(paramGroup);

        // we compute on the leaf grid view
        const auto& leafGridView = gridManager_.grid().leafGridView();

        // create the finite volume grid geometry
        fvGridGeometry_ = std::make_shared<FVGridGeometry>(leafGridView);
        fvGridGeometry_->update();

        // the problem (initial and boundary conditions)
        auto params = std::make_shared<SpatialParams>(fvGridGeometry_, paramGroup);
        problem_ = std::make_shared<Problem>(fvGridGeometry_, params, paramGroup);

        // resize and initialize the given solution vector
        x_ = std::make_shared<SolutionVector>();
        x_->resize(fvGridGeometry_->numDofs());
        problem_->applyInitialSolution(*x_);

        // the grid variables
        gridVariables_ = std::make_shared<GridVariables>(problem_, fvGridGeometry_);
        gridVariables_->init(*x_);

        // initialize the vtk output module
        using IOFields = GetPropType<TypeTag, Properties::IOFields>;
        using VelocityOutput = GetPropType<TypeTag, Properties::VelocityOutput>;
        vtkWriter_ = std::make_unique<OutputModule>(*gridVariables_, *x_, problem_->name());
        IOFields::initOutputModule(*vtkWriter_);
        vtkWriter_->addVelocityOutput(std::make_shared<VelocityOutput>(*gridVariables_));

        // the assembler without time loop for stationary problem
        assembler_ = std::make_shared<Assembler>(problem_, fvGridGeometry_, gridVariables_);

        // the linear/non-linear solvers
        auto linearSolver = std::make_shared<LinearSolver>();
        newtonSolver_ = std::make_unique<NewtonSolver>(assembler_, linearSolver);
    }

    //! Solve the system
    void solve()
    {
        newtonSolver_->solve(*x_);
    }

    //! Write current state to disk
    void write(Scalar t)
    {
        vtkWriter_->write(t);
    }

    //! Return a pointer to the grid geometry
    std::shared_ptr<FVGridGeometry> gridGeometryPointer()
    { return fvGridGeometry_; }

    //! Return a pointer to the grid variables
    std::shared_ptr<GridVariables> gridVariablesPointer()
    { return gridVariables_; }

    //! Return a pointer to the problem
    std::shared_ptr<Problem> problemPointer()
    { return problem_; }

    //! Return a pointer to the solution
    std::shared_ptr<SolutionVector> solutionPointer()
    { return x_; }

private:
    GridManager<GetPropType<TypeTag, Properties::Grid>> gridManager_;
    std::shared_ptr<FVGridGeometry> fvGridGeometry_;
    std::shared_ptr<Problem> problem_;
    std::shared_ptr<GridVariables> gridVariables_;
    std::shared_ptr<Assembler> assembler_;

    std::unique_ptr<NewtonSolver> newtonSolver_;
    std::unique_ptr<OutputModule> vtkWriter_;

    std::shared_ptr<SolutionVector> x_;
};

} // end namespace Dumux
