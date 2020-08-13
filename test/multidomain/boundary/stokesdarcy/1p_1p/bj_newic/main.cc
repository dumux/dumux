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
 * \brief A test problem for the coupled Stokes/Darcy problem (1p).
 */

#include <config.h>

#include <ctime>
#include <iostream>

#include <dune/common/parallel/mpihelper.hh>
#include <dune/common/timer.hh>
#include <dune/istl/io.hh>

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

#include "problem_darcy.hh"
#include "problem_stokes.hh"

namespace Dumux {
namespace Properties {

template<class TypeTag>
struct CouplingManager<TypeTag, TTag::StokesOneP>
{
    using Traits = StaggeredMultiDomainTraits<TypeTag, TypeTag, Properties::TTag::DarcyOneP>;
    using type = Dumux::StokesDarcyCouplingManager<Traits>;
};

template<class TypeTag>
struct CouplingManager<TypeTag, TTag::DarcyOneP>
{
    using Traits = StaggeredMultiDomainTraits<Properties::TTag::StokesOneP, Properties::TTag::StokesOneP, TypeTag>;
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
auto createStokesAnalyticalSolution(const Problem& problem)
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
void printStokesL2Error(const Problem& problem, const SolutionVector& x)
{
    using namespace Dumux;
    using Scalar = double;
    using TypeTag = Properties::TTag::StokesOneP;
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
            const Scalar delta = x[dofIdx] - problem.analyticalSolution(scv.dofPosition())[2/*pressureIdx*/];
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

int main(int argc, char** argv) try
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
    using StokesTypeTag = Properties::TTag::StokesOneP;
    using DarcyTypeTag = Properties::TTag::DarcyOneP;

    // try to create a grid (from the given grid file or the input file)
    // for both sub-domains
    using DarcyGridManager = Dumux::GridManager<GetPropType<DarcyTypeTag, Properties::Grid>>;
    DarcyGridManager darcyGridManager;
    darcyGridManager.init("Darcy"); // pass parameter group

    using StokesGridManager = Dumux::GridManager<GetPropType<StokesTypeTag, Properties::Grid>>;
    StokesGridManager stokesGridManager;
    stokesGridManager.init("Stokes"); // pass parameter group

    // we compute on the leaf grid view
    const auto& darcyGridView = darcyGridManager.grid().leafGridView();
    const auto& stokesGridView = stokesGridManager.grid().leafGridView();

    // create the finite volume grid geometry
    using StokesGridGeometry = GetPropType<StokesTypeTag, Properties::GridGeometry>;
    auto stokesGridGeometry = std::make_shared<StokesGridGeometry>(stokesGridView);
    stokesGridGeometry->update();
    using DarcyGridGeometry = GetPropType<DarcyTypeTag, Properties::GridGeometry>;
    auto darcyGridGeometry = std::make_shared<DarcyGridGeometry>(darcyGridView);
    darcyGridGeometry->update();

    using Traits = StaggeredMultiDomainTraits<StokesTypeTag, StokesTypeTag, DarcyTypeTag>;

    // the coupling manager
    using CouplingManager = StokesDarcyCouplingManager<Traits>;
    auto couplingManager = std::make_shared<CouplingManager>(stokesGridGeometry, darcyGridGeometry);

    // the indices
    constexpr auto freeFlowCellCenterIdx = CouplingManager::freeFlowCellCenterIdx;
    constexpr auto freeFlowFaceIdx = CouplingManager::freeFlowFaceIdx;
    constexpr auto porousMediumIdx = CouplingManager::porousMediumIdx;

    // the problem (initial and boundary conditions)
    using StokesProblem = GetPropType<StokesTypeTag, Properties::Problem>;
    auto stokesProblem = std::make_shared<StokesProblem>(stokesGridGeometry, couplingManager);
    using DarcyProblem = GetPropType<DarcyTypeTag, Properties::Problem>;
    auto darcyProblem = std::make_shared<DarcyProblem>(darcyGridGeometry, couplingManager);

    // the solution vector
    Traits::SolutionVector sol;
    sol[freeFlowCellCenterIdx].resize(stokesGridGeometry->numCellCenterDofs());
    sol[freeFlowFaceIdx].resize(stokesGridGeometry->numFaceDofs());
    sol[porousMediumIdx].resize(darcyGridGeometry->numDofs());

    // get a solution vector storing references to the two Stokes solution vectors
    auto stokesSol = partial(sol, freeFlowFaceIdx, freeFlowCellCenterIdx);

    couplingManager->init(stokesProblem, darcyProblem, sol);

    // the grid variables
    using StokesGridVariables = GetPropType<StokesTypeTag, Properties::GridVariables>;
    auto stokesGridVariables = std::make_shared<StokesGridVariables>(stokesProblem, stokesGridGeometry);
    stokesGridVariables->init(stokesSol);
    using DarcyGridVariables = GetPropType<DarcyTypeTag, Properties::GridVariables>;
    auto darcyGridVariables = std::make_shared<DarcyGridVariables>(darcyProblem, darcyGridGeometry);
    darcyGridVariables->init(sol[porousMediumIdx]);

    // intialize the vtk output module
    using Scalar = typename Traits::Scalar;
    StaggeredVtkOutputModule<StokesGridVariables, decltype(stokesSol)> stokesVtkWriter(*stokesGridVariables, stokesSol, stokesProblem->name());
    GetPropType<StokesTypeTag, Properties::IOFields>::initOutputModule(stokesVtkWriter);
    const auto stokesAnalyticalSolution = createStokesAnalyticalSolution<Scalar>(*stokesProblem);
    stokesVtkWriter.addField(std::get<0>(stokesAnalyticalSolution), "pressureExact");
    stokesVtkWriter.addField(std::get<1>(stokesAnalyticalSolution), "velocityExact");
    stokesVtkWriter.addFaceField(std::get<2>(stokesAnalyticalSolution), "faceVelocityExact");
    stokesVtkWriter.write(0.0);

    VtkOutputModule<DarcyGridVariables, GetPropType<DarcyTypeTag, Properties::SolutionVector>> darcyVtkWriter(*darcyGridVariables, sol[porousMediumIdx],  darcyProblem->name());
    using DarcyVelocityOutput = GetPropType<DarcyTypeTag, Properties::VelocityOutput>;
    darcyVtkWriter.addVelocityOutput(std::make_shared<DarcyVelocityOutput>(*darcyGridVariables));
    GetPropType<DarcyTypeTag, Properties::IOFields>::initOutputModule(darcyVtkWriter);
    const auto darcyAnalyticalSolution = createDarcyAnalyticalSolution<Scalar>(*darcyProblem);
    darcyVtkWriter.addField(std::get<0>(darcyAnalyticalSolution), "pressureExact");
    darcyVtkWriter.addField(std::get<1>(darcyAnalyticalSolution), "velocityExact");
    darcyVtkWriter.write(0.0);

    // the assembler for a stationary problem
    using Assembler = MultiDomainFVAssembler<Traits, CouplingManager, DiffMethod::numeric>;
    auto assembler = std::make_shared<Assembler>(std::make_tuple(stokesProblem, stokesProblem, darcyProblem),
                                                 std::make_tuple(stokesGridGeometry->faceFVGridGeometryPtr(),
                                                                 stokesGridGeometry->cellCenterFVGridGeometryPtr(),
                                                                 darcyGridGeometry),
                                                 std::make_tuple(stokesGridVariables->faceGridVariablesPtr(),
                                                                 stokesGridVariables->cellCenterGridVariablesPtr(),
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
    stokesVtkWriter.write(1.0);
    darcyVtkWriter.write(1.0);

    printStokesL2Error(*stokesProblem, stokesSol);
    printDarcyL2Error(*darcyProblem, sol[porousMediumIdx]);

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
catch (Dumux::ParameterException &e)
{
    std::cerr << std::endl << e << " ---> Abort!" << std::endl;
    return 1;
}
catch (Dune::DGFException & e)
{
    std::cerr << "DGF exception thrown (" << e <<
                 "). Most likely, the DGF file name is wrong "
                 "or the DGF file is corrupted, "
                 "e.g. missing hash at end of file or wrong number (dimensions) of entries."
                 << " ---> Abort!" << std::endl;
    return 2;
}
catch (Dune::Exception &e)
{
    std::cerr << "Dune reported error: " << e << " ---> Abort!" << std::endl;
    return 3;
}
catch (...)
{
    std::cerr << "Unknown exception thrown! ---> Abort!" << std::endl;
    return 4;
}
