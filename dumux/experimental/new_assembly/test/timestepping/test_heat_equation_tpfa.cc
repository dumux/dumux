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
 * \brief Solves a poisson problem discretized using TPFA.
 *        This exposes the minimal requirements for adding a new model.
 */
#include <config.h>

#include <iostream>
#include <vector>
#include <string>
#include <memory>
#include <cmath>

#include <dune/grid/yaspgrid.hh>
#include <dune/grid/utility/structuredgridfactory.hh>
#include <dune/grid/io/file/vtk/vtkwriter.hh>

#include <dumux/common/initialize.hh>
#include <dumux/discretization/cellcentered/tpfa/fvgridgeometry.hh>

#include <dumux/experimental/new_assembly/dumux/common/multiindex.hh>

#include <dumux/experimental/new_assembly/dumux/assembly/assembler.hh>
#include <dumux/experimental/new_assembly/dumux/assembly/jacobianpattern.hh>

#include <dumux/experimental/new_assembly/dumux/linear/system.hh>
#include <dumux/experimental/new_assembly/dumux/linear/operator.hh>
#include <dumux/experimental/new_assembly/dumux/linear/dune/solvers.hh>
#include <dumux/experimental/new_assembly/dumux/linear/dune/labackend.hh>
#include <dumux/experimental/new_assembly/dumux/nonlinear/newton.hh>

#include <dumux/experimental/new_assembly/dumux/timestepping/multistageresidualfunction.hh>
#include <dumux/experimental/new_assembly/dumux/timestepping/multistagetimestepper.hh>
#include <dumux/experimental/new_assembly/dumux/timestepping/multistagemethods.hh>

#if HAVE_EIGEN3
#include <dumux/experimental/new_assembly/dumux/linear/eigen/labackend.hh>
#include <dumux/experimental/new_assembly/dumux/linear/eigen/solvers.hh>
#endif

#include "heat_equation.hh"

template<typename LABackend>
struct DuneSolverTrait
{
    using type = Dumux::DuneCGSolver<typename LABackend::Vector>;
};

#if HAVE_EIGEN3
template<typename LABackend>
struct EigenSolverTrait
{
    using type = Dumux::EigenCGSolver<typename LABackend::Matrix,
                                      typename LABackend::Vector>;
};
#endif

template<typename Scalar,
         typename LABackend,
         template<typename> typename SolverTrait,
         typename TimeIntegrationMethod>
void runSimulation(unsigned int numCellsPerSide,
                   std::shared_ptr<TimeIntegrationMethod> method,
                   const std::string& filenameSuffix)
{
    using Grid = Dune::YaspGrid<2>;
    using GridFactory = Dune::StructuredGridFactory<Grid>;
    auto grid = GridFactory::createCubeGrid(
            {0.0, 0.0},
            {1.0, 1.0},
            {numCellsPerSide, numCellsPerSide}
    );

    using GridGeometry = Dumux::CCTpfaFVGridGeometry<typename Grid::LeafGridView>;
    auto gridGeometry = std::make_shared<GridGeometry>(grid->leafGridView());

    // select matrix/vector types
    LABackend backend;
    using Matrix = typename LABackend::Matrix;
    using Vector = typename LABackend::Vector;
    using IndexStrategy = typename LABackend::IndexStrategy;
    using GridVariables = Dumux::HeatEquationGridVariables<GridGeometry, Vector, IndexStrategy>;

    auto gridVariables = std::make_shared<GridVariables>(
        gridGeometry,
        std::make_shared<IndexStrategy>(
            backend.makeIndexStrategy(gridGeometry->numDofs())
        ),
        Dumux::makeVector(backend, *gridGeometry)
    );

    using Assembler = Dumux::Assembler<Dumux::HeatEquationLocalAssembler<GridGeometry, GridVariables>>;
    using LinearOperator = Dumux::Linear::OperatorType<LABackend>;

    // For multi-staging, this PDE type must be selected
    using PDE = Dumux::MultiStageResidualFunction<Vector, Matrix, Assembler, LinearOperator>;
    auto pde = std::make_shared<PDE>(std::make_shared<Assembler>(gridGeometry),
                                     Dumux::makeMatrix(backend, *gridGeometry),
                                     Dumux::makeVector(backend, *gridGeometry));

    using LinearSolver = typename SolverTrait<LABackend>::type;
    auto linearSolver = std::make_shared<LinearSolver>();

    using NonLinearSolver = Dumux::NewtonSolver<PDE, LinearSolver>;
    auto newtonSolver = std::make_shared<NonLinearSolver>(pde, linearSolver);

    const Scalar dt = 0.0001;
    const std::size_t numTimeSteps = 100;

    Dumux::MultiStageTimeStepper timeStepper{newtonSolver, method};
    for (std::size_t i = 0; i < numTimeSteps; ++i)
    {
        std::cout << "Performing time integration at t = " << dt*i << "..." << std::endl;
        timeStepper.step(*gridVariables, dt*i, dt);
    }

    Dune::VTKWriter<typename GridGeometry::GridView> vtkWriter{gridGeometry->gridView()};
    vtkWriter.addCellData(gridVariables->dofs(), "u");
    vtkWriter.write("test_heat_equation_" + filenameSuffix);
}

int main (int argc, char *argv[])
{
    Dumux::initialize(argc, argv);

    using Scalar = double;
    static constexpr int numEq = 1;

    std::cout << "Running simulation (implicit Euler) with default blocked linear algebra backend" << std::endl;
    using BlockedBackend = Dumux::DefaultBlockedLinearAlgebraBackend<Scalar, numEq>;
    runSimulation<Scalar, BlockedBackend, DuneSolverTrait>(
        25,
        std::make_shared<Dumux::MultiStage::ImplicitEuler<Scalar>>(),
        "default_backend_imp_euler"
    );

    runSimulation<Scalar, BlockedBackend, DuneSolverTrait>(
        25,
        std::make_shared<Dumux::MultiStage::RungeKuttaExplicitFourthOrder<Scalar>>(),
        "default_backend_rk4"
    );

    std::cout << "\n";
    std::cout << "Testing default flat linear algebra backend" << std::endl;
    using FlatBackend = Dumux::DefaultFlatLinearAlgebraBackend<Scalar, numEq>;
    runSimulation<Scalar, FlatBackend, DuneSolverTrait>(
        25,
        std::make_shared<Dumux::MultiStage::Theta<Scalar>>(0.5),
        "flat_backend_crank_nicholson"
    );

#if HAVE_EIGEN3
    std::cout << "\n";
    std::cout << "Testing eigen linear algebra backend" << std::endl;
    using EigenBackend = Dumux::EigenLinearAlgebraBackend<Scalar, numEq>;
    runSimulation<Scalar, EigenBackend, EigenSolverTrait>(
        25,
        std::make_shared<Dumux::MultiStage::ExplicitEuler<Scalar>>(),
        "eigen_backend_exp_euler"
    );
#endif

    std::cout << "\nAll tests passed" << std::endl;
    return 0;
}
