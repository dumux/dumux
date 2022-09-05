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
 */
#include <config.h>

#include <iostream>
#include <vector>
#include <memory>
#include <cmath>

#include <dune/grid/yaspgrid.hh>
#include <dune/grid/utility/structuredgridfactory.hh>

#include <dumux/common/initialize.hh>
#include <dumux/discretization/cellcentered/tpfa/fvgridgeometry.hh>

#include <dumux/experimental/new_assembly/dumux/common/multiindex.hh>

#include <dumux/experimental/new_assembly/dumux/assembly/assembler.hh>
#include <dumux/experimental/new_assembly/dumux/assembly/residualfunction.hh>
#include <dumux/experimental/new_assembly/dumux/assembly/jacobianpattern.hh>

#include <dumux/experimental/new_assembly/dumux/nonlinear/newton.hh>
#include <dumux/experimental/new_assembly/dumux/linear/system.hh>
#include <dumux/experimental/new_assembly/dumux/linear/operator.hh>

// solvers & backend implementations
#include <dumux/experimental/new_assembly/dumux/linear/dune/solvers.hh>
#include <dumux/experimental/new_assembly/dumux/linear/dune/labackend.hh>

#if HAVE_EIGEN3
#include <dumux/experimental/new_assembly/dumux/linear/eigen/labackend.hh>
#include <dumux/experimental/new_assembly/dumux/linear/eigen/solvers.hh>
#endif

#include "poisson.hh"

template<typename Scalar, typename GridVariables>
Scalar computeError(const GridVariables& gridVariables)
{
    const auto& gridGeometry = gridVariables.gridGeometry();
    const auto& problem = gridVariables.problem();

    Scalar error = 0.0;
    std::vector<Scalar> uExact(gridGeometry.numDofs());
    for (const auto& element : elements(gridGeometry.gridView()))
    {
        auto fvGeometry = localView(gridGeometry).bindElement(element);
        for (const auto& scv : scvs(fvGeometry))
        {
            const auto exact = problem.exactAtPos(scv.center());
            const auto localError = exact - Dumux::LinearSystem::get(
                gridVariables.dofs(),
                gridVariables.getDofIndex(Dumux::MultiIndex{scv.dofIndex(), 0})
            );
            uExact[scv.dofIndex()] = exact;
            error += scv.volume()*localError*localError;
        }
    }

    using std::sqrt;
    return sqrt(error);
}

template<typename LABackend>
struct DefaultSolverTrait
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
         template<typename> typename SolverTrait = DefaultSolverTrait>
Scalar runSimulation(unsigned int numCellsPerSide)
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
    using GridVariables = Dumux::PoissonGridVariables<GridGeometry, Vector, IndexStrategy>;

    auto indexStrategy = std::make_shared<IndexStrategy>(
        backend.makeIndexStrategy(gridGeometry->numDofs())
    );
    auto gridVariables = std::make_shared<GridVariables>(
        gridGeometry, indexStrategy, Dumux::makeVector(backend, *gridGeometry)
    );

    using Assembler = Dumux::Assembler<Dumux::PoissonLocalAssembler<GridGeometry, GridVariables>>;
    using LinearOperator = Dumux::Linear::OperatorType<LABackend>;
    using Function = Dumux::ResidualFunction<Vector, Matrix, Assembler, LinearOperator>;
    auto function = std::make_shared<Function>(
        std::make_shared<const Assembler>(gridGeometry),
        Dumux::makeMatrix(backend, *gridGeometry),
        Dumux::makeVector(backend, *gridGeometry)
    );

    Dumux::NewtonSolver newtonSolver{function, std::make_shared<typename SolverTrait<LABackend>::type>()};
    newtonSolver.solve(*gridVariables);
    return computeError<Scalar>(*gridVariables);
}

template<typename Scalar,
         typename LABackend,
         template<typename> typename SolverTrait = DefaultSolverTrait>
void testConvergenceRates(unsigned int numRefinements = 4)
{
    using std::pow;
    std::vector<double> errors;
    for (int refinement = 0; refinement < numRefinements + 1; ++refinement)
        errors.push_back(
            runSimulation<Scalar, LABackend, SolverTrait>(pow(2, refinement))
        );

    std::cout << "Errors/rates\n";
    std::vector<double> rates;
    for (int i = 0; i < errors.size(); ++i)
    {
        if (i > 0)
            rates.push_back(
                (std::log(errors[i]) - std::log(errors[i-1]))/std::log(2)
            );
        std::cout << errors[i]
                  << ", "
                  << (i == 0 ? "" : std::to_string(rates.back()))
                  << std::endl;
    }

    if (rates.back() > -1.9)
        DUNE_THROW(Dune::InvalidStateException, "Final rate below 1.9");
}

int main (int argc, char *argv[])
{
    Dumux::initialize(argc, argv);

    using Scalar = double;
    static constexpr int numEq = 1;

    std::cout << "Testing default blocked linear algebra backend" << std::endl;
    using BlockedBackend = Dumux::DefaultBlockedLinearAlgebraBackend<Scalar, numEq>;
    testConvergenceRates<Scalar, BlockedBackend>();

    std::cout << "\n";
    std::cout << "Testing default flat linear algebra backend" << std::endl;
    using FlatBackend = Dumux::DefaultFlatLinearAlgebraBackend<Scalar, numEq>;
    testConvergenceRates<Scalar, FlatBackend>();

#if HAVE_EIGEN3
    std::cout << "\n";
    std::cout << "Testing eigen linear algebra backend" << std::endl;
    using EigenBackend = Dumux::EigenLinearAlgebraBackend<Scalar, numEq>;
    testConvergenceRates<Scalar, EigenBackend, EigenSolverTrait>();
#endif

    std::cout << "\nAll tests passed" << std::endl;
    return 0;
}
