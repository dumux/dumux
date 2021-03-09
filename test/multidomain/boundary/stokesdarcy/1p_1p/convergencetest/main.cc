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
 * \brief A test problem for the coupled FreeFlow/Darcy problem (1p).
 */

#include <config.h>

#include <ctime>
#include <iostream>

#include <dune/common/parallel/mpihelper.hh>
#include <dune/common/timer.hh>

#include <dumux/common/properties.hh>
#include <dumux/common/parameters.hh>
#include <dumux/common/partial.hh>
#include <dumux/common/dumuxmessage.hh>
#include <dumux/linear/seqsolverbackend.hh>
#include <dumux/assembly/fvassembler.hh>
#include <dumux/assembly/diffmethod.hh>
#include <dumux/discretization/method.hh>
#include <dumux/io/vtkoutputmodule.hh>
#include <dumux/io/staggeredvtkoutputmodule.hh>
#include <dumux/io/grid/gridmanager_yasp.hh>

#include <dumux/multidomain/staggeredtraits.hh>
#include <dumux/multidomain/fvassembler.hh>
#include <dumux/multidomain/newtonsolver.hh>
#include <test/freeflow/navierstokes/l2error.hh>

#include <dumux/multidomain/boundary/stokesdarcy/couplingmanager.hh>

#include "testcase.hh"
#include "problem_darcy.hh"
#include "problem_stokes.hh"

namespace Dumux {
namespace Properties {

template<class TypeTag>
struct CouplingManager<TypeTag, TTag::FreeFlowOneP>
{
    using Traits = StaggeredMultiDomainTraits<TypeTag, TypeTag, Properties::TTag::DarcyOneP>;
    using type = Dumux::StokesDarcyCouplingManager<Traits>;
};

template<class TypeTag>
struct CouplingManager<TypeTag, TTag::DarcyOneP>
{
    using Traits = StaggeredMultiDomainTraits<Properties::TTag::FreeFlowOneP, Properties::TTag::FreeFlowOneP, TypeTag>;
    using type = Dumux::StokesDarcyCouplingManager<Traits>;
};

} // end namespace Properties
} // end namespace Dumux

/*!
* \brief Creates analytical solution.
* Returns a tuple of the analytical solution for the pressure, the velocity and the velocity at the faces
* \param problem the problem for which to evaluate the analytical solution
*/
template<class Scalar, class Problem>
auto createFreeFlowAnalyticalSolution(const Problem& problem)
{
    const auto& gridGeometry = problem.gridGeometry();
    using GridView = typename std::decay_t<decltype(gridGeometry)>::GridView;

    static constexpr auto dim = GridView::dimension;
    static constexpr auto dimWorld = GridView::dimensionworld;

    using VelocityVector = Dune::FieldVector<Scalar, dimWorld>;

    std::vector<Scalar> analyticalPressure;
    std::vector<VelocityVector> analyticalVelocity;
    std::vector<Scalar> analyticalVelocityOnFace;

    analyticalPressure.resize(gridGeometry.numCellCenterDofs());
    analyticalVelocity.resize(gridGeometry.numCellCenterDofs());
    analyticalVelocityOnFace.resize(gridGeometry.numFaceDofs());

    using Indices = typename Problem::Indices;
    for (const auto& element : elements(gridGeometry.gridView()))
    {
        auto fvGeometry = localView(gridGeometry);
        fvGeometry.bindElement(element);
        for (auto&& scv : scvs(fvGeometry))
        {
            auto ccDofIdx = scv.dofIndex();
            auto ccDofPosition = scv.dofPosition();
            auto analyticalSolutionAtCc = problem.analyticalSolution(ccDofPosition);

            // velocities on faces
            for (auto&& scvf : scvfs(fvGeometry))
            {
                const auto faceDofIdx = scvf.dofIndex();
                const auto faceDofPosition = scvf.center();
                const auto dirIdx = scvf.directionIndex();
                const auto analyticalSolutionAtFace = problem.analyticalSolution(faceDofPosition);
                analyticalVelocityOnFace[faceDofIdx] = analyticalSolutionAtFace[Indices::velocity(dirIdx)];
            }

            analyticalPressure[ccDofIdx] = analyticalSolutionAtCc[Indices::pressureIdx];

            for (int dirIdx = 0; dirIdx < dim; ++dirIdx)
                analyticalVelocity[ccDofIdx][dirIdx] = analyticalSolutionAtCc[Indices::velocity(dirIdx)];
        }
    }

    return std::make_tuple(analyticalPressure, analyticalVelocity, analyticalVelocityOnFace);
}

/*!
* \brief Creates analytical solution.
* Returns a tuple of the analytical solution for the pressure, the velocity and the velocity at the faces
* \param problem the problem for which to evaluate the analytical solution
*/
template<class Scalar, class Problem>
auto createDarcyAnalyticalSolution(const Problem& problem)
{
    const auto& gridGeometry = problem.gridGeometry();
    using GridView = typename std::decay_t<decltype(gridGeometry)>::GridView;

    static constexpr auto dim = GridView::dimension;
    static constexpr auto dimWorld = GridView::dimensionworld;

    using VelocityVector = Dune::FieldVector<Scalar, dimWorld>;

    std::vector<Scalar> analyticalPressure;
    std::vector<VelocityVector> analyticalVelocity;

    analyticalPressure.resize(gridGeometry.numDofs());
    analyticalVelocity.resize(gridGeometry.numDofs());

    for (const auto& element : elements(gridGeometry.gridView()))
    {
        auto fvGeometry = localView(gridGeometry);
        fvGeometry.bindElement(element);
        for (auto&& scv : scvs(fvGeometry))
        {
            const auto ccDofIdx = scv.dofIndex();
            const auto ccDofPosition = scv.dofPosition();
            const auto analyticalSolutionAtCc = problem.analyticalSolution(ccDofPosition);
            analyticalPressure[ccDofIdx] = analyticalSolutionAtCc[dim];

            for (int dirIdx = 0; dirIdx < dim; ++dirIdx)
                analyticalVelocity[ccDofIdx][dirIdx] = analyticalSolutionAtCc[dirIdx];
        }
    }

    return std::make_tuple(analyticalPressure, analyticalVelocity);
}

template<class Problem, class SolutionVector>
void printFreeFlowL2Error(const Problem& problem, const SolutionVector& x)
{
    using namespace Dumux;
    using Scalar = double;
    using TypeTag = Properties::TTag::FreeFlowOneP;
    using Indices = typename GetPropType<TypeTag, Properties::ModelTraits>::Indices;
    using ModelTraits = GetPropType<TypeTag, Properties::ModelTraits>;
    using PrimaryVariables = GetPropType<TypeTag, Properties::PrimaryVariables>;

    using L2Error = NavierStokesTestL2Error<Scalar, ModelTraits, PrimaryVariables>;
    const auto l2error = L2Error::calculateL2Error(problem, x);
    const int numCellCenterDofs = problem.gridGeometry().numCellCenterDofs();
    const int numFaceDofs = problem.gridGeometry().numFaceDofs();
    std::ostream tmpOutputObject(std::cout.rdbuf()); // create temporary output with fixed formatting without affecting std::cout
    tmpOutputObject << std::setprecision(8) << "** L2 error (abs/rel) for "
                    << std::setw(6) << numCellCenterDofs << " cc dofs and " << numFaceDofs << " face dofs (total: " << numCellCenterDofs + numFaceDofs << "): "
                    << std::scientific
                    << "L2(p) = " << l2error.first[Indices::pressureIdx] << " / " << l2error.second[Indices::pressureIdx]
                    << " , L2(vx) = " << l2error.first[Indices::velocityXIdx] << " / " << l2error.second[Indices::velocityXIdx]
                    << " , L2(vy) = " << l2error.first[Indices::velocityYIdx] << " / " << l2error.second[Indices::velocityYIdx]
                    << std::endl;

    // write the norm into a log file
    std::ofstream logFile;
    logFile.open(problem.name() + ".log", std::ios::app);
    logFile << "[ConvergenceTest] L2(p) = " << l2error.first[Indices::pressureIdx] << " L2(vx) = " << l2error.first[Indices::velocityXIdx] << " L2(vy) = " << l2error.first[Indices::velocityYIdx] << std::endl;
    logFile.close();
}

template<class Problem, class SolutionVector>
void printDarcyL2Error(const Problem& problem, const SolutionVector& x)
{
    using namespace Dumux;
    using Scalar = double;

    Scalar l2error = 0.0;

    for (const auto& element : elements(problem.gridGeometry().gridView()))
    {
        auto fvGeometry = localView(problem.gridGeometry());
        fvGeometry.bindElement(element);

        for (auto&& scv : scvs(fvGeometry))
        {
            const auto dofIdx = scv.dofIndex();
            const Scalar delta = x[dofIdx] - problem.analyticalSolution(scv.center())[2/*pressureIdx*/];
            l2error += scv.volume()*(delta*delta);
        }
    }
    using std::sqrt;
    l2error = sqrt(l2error);

    const auto numDofs = problem.gridGeometry().numDofs();
    std::ostream tmp(std::cout.rdbuf());
    tmp << std::setprecision(8) << "** L2 error (abs) for "
            << std::setw(6) << numDofs << " cc dofs "
            << std::scientific
            << "L2 error = " << l2error
            << std::endl;

    // write the norm into a log file
    std::ofstream logFile;
    logFile.open(problem.name() + ".log", std::ios::app);
    logFile << "[ConvergenceTest] L2(p) = " << l2error << std::endl;
    logFile.close();
}

int main(int argc, char** argv)
{
    using namespace Dumux;

    // initialize MPI, finalize is done automatically on exit
    const auto& mpiHelper = Dune::MPIHelper::instance(argc, argv);

    // print dumux start message
    if (mpiHelper.rank() == 0)
        DumuxMessage::print(/*firstCall=*/true);

    // parse command line arguments and input file
    Parameters::init(argc, argv);

    // Define the sub problem type tags
    using FreeFlowTypeTag = Properties::TTag::FreeFlowOneP;
    using DarcyTypeTag = Properties::TTag::DarcyOneP;

    // try to create a grid (from the given grid file or the input file)
    // for both sub-domains
    using DarcyGridManager = Dumux::GridManager<GetPropType<DarcyTypeTag, Properties::Grid>>;
    DarcyGridManager darcyGridManager;
    darcyGridManager.init("Darcy"); // pass parameter group

    using FreeFlowGridManager = Dumux::GridManager<GetPropType<FreeFlowTypeTag, Properties::Grid>>;
    FreeFlowGridManager freeFlowGridManager;
    freeFlowGridManager.init("FreeFlow"); // pass parameter group

    // we compute on the leaf grid view
    const auto& darcyGridView = darcyGridManager.grid().leafGridView();
    const auto& freeFlowGridView = freeFlowGridManager.grid().leafGridView();

    // create the finite volume grid geometry
    using FreeFlowGridGeometry = GetPropType<FreeFlowTypeTag, Properties::GridGeometry>;
    auto freeFlowGridGeometry = std::make_shared<FreeFlowGridGeometry>(freeFlowGridView);
    freeFlowGridGeometry->update();
    using DarcyGridGeometry = GetPropType<DarcyTypeTag, Properties::GridGeometry>;
    auto darcyGridGeometry = std::make_shared<DarcyGridGeometry>(darcyGridView);
    darcyGridGeometry->update();

    using Traits = StaggeredMultiDomainTraits<FreeFlowTypeTag, FreeFlowTypeTag, DarcyTypeTag>;

    // the coupling manager
    using CouplingManager = StokesDarcyCouplingManager<Traits>;
    auto couplingManager = std::make_shared<CouplingManager>(freeFlowGridGeometry, darcyGridGeometry);

    // the indices
    constexpr auto freeFlowCellCenterIdx = CouplingManager::stokesCellCenterIdx;
    constexpr auto freeFlowFaceIdx = CouplingManager::stokesFaceIdx;
    constexpr auto darcyIdx = CouplingManager::darcyIdx;

    // the problem (initial and boundary conditions)
    const auto testCaseName = getParam<std::string>("Problem.TestCase");
    const auto testCase = [&testCaseName]()
    {
        if (testCaseName == "ShiueExampleOne")
            return TestCase::ShiueExampleOne;
        else if (testCaseName == "ShiueExampleTwo")
            return TestCase::ShiueExampleTwo;
        else if (testCaseName == "Rybak")
            return TestCase::Rybak;
        else if (testCaseName == "Schneider")
            return TestCase::Schneider;
        else
            DUNE_THROW(Dune::InvalidStateException, testCaseName + " is not a valid test case");
    }();

    using FreeFlowProblem = GetPropType<FreeFlowTypeTag, Properties::Problem>;
    auto freeFlowProblem = std::make_shared<FreeFlowProblem>(freeFlowGridGeometry, couplingManager, testCase, testCaseName);
    using DarcyProblem = GetPropType<DarcyTypeTag, Properties::Problem>;
    auto spatialParams = std::make_shared<typename DarcyProblem::SpatialParams>(darcyGridGeometry, testCase);
    auto darcyProblem = std::make_shared<DarcyProblem>(darcyGridGeometry, couplingManager, spatialParams, testCase, testCaseName);

    // the solution vector
    Traits::SolutionVector sol;
    sol[freeFlowCellCenterIdx].resize(freeFlowGridGeometry->numCellCenterDofs());
    sol[freeFlowFaceIdx].resize(freeFlowGridGeometry->numFaceDofs());
    sol[darcyIdx].resize(darcyGridGeometry->numDofs());

    // get a solution vector storing references to the two FreeFlow solution vectors
    auto freeFlowSol = partial(sol, freeFlowFaceIdx, freeFlowCellCenterIdx);

    couplingManager->init(freeFlowProblem, darcyProblem, sol);

    // the grid variables
    using FreeFlowGridVariables = GetPropType<FreeFlowTypeTag, Properties::GridVariables>;
    auto freeFlowGridVariables = std::make_shared<FreeFlowGridVariables>(freeFlowProblem, freeFlowGridGeometry);
    freeFlowGridVariables->init(freeFlowSol);
    using DarcyGridVariables = GetPropType<DarcyTypeTag, Properties::GridVariables>;
    auto darcyGridVariables = std::make_shared<DarcyGridVariables>(darcyProblem, darcyGridGeometry);
    darcyGridVariables->init(sol[darcyIdx]);

    // intialize the vtk output module
    using Scalar = typename Traits::Scalar;
    StaggeredVtkOutputModule<FreeFlowGridVariables, decltype(freeFlowSol)> freeFlowVtkWriter(*freeFlowGridVariables, freeFlowSol, freeFlowProblem->name());
    GetPropType<FreeFlowTypeTag, Properties::IOFields>::initOutputModule(freeFlowVtkWriter);
    const auto freeFlowAnalyticalSolution = createFreeFlowAnalyticalSolution<Scalar>(*freeFlowProblem);
    freeFlowVtkWriter.addField(std::get<0>(freeFlowAnalyticalSolution), "pressureExact");
    freeFlowVtkWriter.addField(std::get<1>(freeFlowAnalyticalSolution), "velocityExact");
    freeFlowVtkWriter.addFaceField(std::get<2>(freeFlowAnalyticalSolution), "faceVelocityExact");
    freeFlowVtkWriter.write(0.0);

    VtkOutputModule<DarcyGridVariables, GetPropType<DarcyTypeTag, Properties::SolutionVector>> darcyVtkWriter(*darcyGridVariables, sol[darcyIdx],  darcyProblem->name());
    using DarcyVelocityOutput = GetPropType<DarcyTypeTag, Properties::VelocityOutput>;
    darcyVtkWriter.addVelocityOutput(std::make_shared<DarcyVelocityOutput>(*darcyGridVariables));
    GetPropType<DarcyTypeTag, Properties::IOFields>::initOutputModule(darcyVtkWriter);
    const auto darcyAnalyticalSolution = createDarcyAnalyticalSolution<Scalar>(*darcyProblem);
    darcyVtkWriter.addField(std::get<0>(darcyAnalyticalSolution), "pressureExact");
    darcyVtkWriter.addField(std::get<1>(darcyAnalyticalSolution), "velocityExact");
    darcyVtkWriter.write(0.0);

    // the assembler for a stationary problem
    using Assembler = MultiDomainFVAssembler<Traits, CouplingManager, DiffMethod::numeric>;
    auto assembler = std::make_shared<Assembler>(std::make_tuple(freeFlowProblem, freeFlowProblem, darcyProblem),
                                                 std::make_tuple(freeFlowGridGeometry->faceFVGridGeometryPtr(),
                                                                 freeFlowGridGeometry->cellCenterFVGridGeometryPtr(),
                                                                 darcyGridGeometry),
                                                 std::make_tuple(freeFlowGridVariables->faceGridVariablesPtr(),
                                                                 freeFlowGridVariables->cellCenterGridVariablesPtr(),
                                                                 darcyGridVariables),
                                                 couplingManager);

    // the linear solver
    using LinearSolver = UMFPackBackend;
    auto linearSolver = std::make_shared<LinearSolver>();

    // the non-linear solver
    using NewtonSolver = MultiDomainNewtonSolver<Assembler, LinearSolver, CouplingManager>;
    NewtonSolver nonLinearSolver(assembler, linearSolver, couplingManager);

    // solve the non-linear system
    nonLinearSolver.solve(sol);

    // write vtk output
    freeFlowVtkWriter.write(1.0);
    darcyVtkWriter.write(1.0);

    printFreeFlowL2Error(*freeFlowProblem, freeFlowSol);
    printDarcyL2Error(*darcyProblem, sol[darcyIdx]);

    ////////////////////////////////////////////////////////////
    // finalize, print dumux message to say goodbye
    ////////////////////////////////////////////////////////////

    // print dumux end message
    if (mpiHelper.rank() == 0)
    {
        Parameters::print();
        DumuxMessage::print(/*firstCall=*/false);
    }

    return 0;
} // end main
