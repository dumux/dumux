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
 * \ingroup NavierStokesTests
 * \brief Test for the instationary staggered grid Navier-Stokes model
 *        with analytical solution.
 */
#include <config.h>

#include <ctime>
#include <iostream>
#include <type_traits>
#include <tuple>

#include <dune/common/parallel/mpihelper.hh>
#include <dune/common/timer.hh>
#include <dune/grid/io/file/dgfparser/dgfexception.hh>
#include <dune/grid/io/file/vtk.hh>

#include <dumux/assembly/staggeredfvassembler.hh>
#include <dumux/assembly/diffmethod.hh>
#include <dumux/common/dumuxmessage.hh>
#include <dumux/common/parameters.hh>
#include <dumux/common/properties.hh>
#include <dumux/common/valgrind.hh>
#include <dumux/io/grid/gridmanager.hh>
#include <dumux/io/staggeredvtkoutputmodule.hh>
#include <dumux/linear/seqsolverbackend.hh>
#include <dumux/nonlinear/newtonsolver.hh>

#include "problem.hh"

/*!
* \brief Creates analytical solution.
* Returns a tuple of the analytical solution for the pressure, the velocity and the velocity at the faces
* \param problem the problem for which to evaluate the analytical solution
*/
template<class Problem>
auto createAnalyticalSolution(const Problem& problem)
{
    const auto& fvGridGeometry = problem.fvGridGeometry();
    using GridView = typename std::decay_t<decltype(fvGridGeometry)>::GridView;

    static constexpr auto dim = GridView::dimension;
    static constexpr auto dimWorld = GridView::dimensionworld;

    using Scalar = double;
    using VelocityVector = Dune::FieldVector<Scalar, dimWorld>;

    std::vector<Scalar> analyticalPressure;
    std::vector<VelocityVector> analyticalVelocity;
    std::vector<VelocityVector> analyticalVelocityOnFace;

    analyticalPressure.resize(fvGridGeometry.numCellCenterDofs());
    analyticalVelocity.resize(fvGridGeometry.numCellCenterDofs());
    analyticalVelocityOnFace.resize(fvGridGeometry.numFaceDofs());

    using Indices = typename Problem::Indices;
    for (const auto& element : elements(fvGridGeometry.gridView()))
    {
        auto fvGeometry = localView(fvGridGeometry);
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
                analyticalVelocityOnFace[faceDofIdx][dirIdx] = analyticalSolutionAtFace[Indices::velocity(dirIdx)];
            }

            analyticalPressure[ccDofIdx] = analyticalSolutionAtCc[Indices::pressureIdx];

            for(int dirIdx = 0; dirIdx < dim; ++dirIdx)
                analyticalVelocity[ccDofIdx][dirIdx] = analyticalSolutionAtCc[Indices::velocity(dirIdx)];
        }
    }

    return std::make_tuple(analyticalPressure, analyticalVelocity, analyticalVelocityOnFace);
}

int main(int argc, char** argv) try
{
    using namespace Dumux;

    // define the type tag for this problem
    using TypeTag = Properties::TTag::SincosSteadyTest;

    // initialize MPI, finalize is done automatically on exit
    const auto& mpiHelper = Dune::MPIHelper::instance(argc, argv);

    // print dumux start message
    if (mpiHelper.rank() == 0)
        DumuxMessage::print(/*firstCall=*/true);

    // parse command line arguments and input file
    Parameters::init(argc, argv);

    // try to create a grid (from the given grid file or the input file)
    GridManager<GetPropType<TypeTag, Properties::Grid>> gridManager;
    gridManager.init();

    ////////////////////////////////////////////////////////////
    // run instationary non-linear problem on this grid
    ////////////////////////////////////////////////////////////

    // we compute on the leaf grid view
    const auto& leafGridView = gridManager.grid().leafGridView();

    // create the finite volume grid geometry
    using FVGridGeometry = GetPropType<TypeTag, Properties::FVGridGeometry>;
    auto fvGridGeometry = std::make_shared<FVGridGeometry>(leafGridView);
    fvGridGeometry->update();

    // the problem (initial and boundary conditions)
    using Problem = GetPropType<TypeTag, Properties::Problem>;
    auto problem = std::make_shared<Problem>(fvGridGeometry);

    // the solution vector
    using SolutionVector = GetPropType<TypeTag, Properties::SolutionVector>;
    SolutionVector x;
    x[FVGridGeometry::cellCenterIdx()].resize(fvGridGeometry->numCellCenterDofs());
    x[FVGridGeometry::faceIdx()].resize(fvGridGeometry->numFaceDofs());

    // the grid variables
    using GridVariables = GetPropType<TypeTag, Properties::GridVariables>;
    auto gridVariables = std::make_shared<GridVariables>(problem, fvGridGeometry);
    gridVariables->init(x);

    // initialize the vtk output module
    StaggeredVtkOutputModule<GridVariables, SolutionVector> vtkWriter(*gridVariables, x, problem->name());
    using IOFields = GetPropType<TypeTag, Properties::IOFields>;
    IOFields::initOutputModule(vtkWriter); // Add model specific output fields

    auto analyticalSolution = createAnalyticalSolution(*problem);
    vtkWriter.addField(std::get<0>(analyticalSolution), "pressureExact");
    vtkWriter.addField(std::get<1>(analyticalSolution), "velocityExact");
    vtkWriter.addFaceField(std::get<2>(analyticalSolution), "faceVelocityExact");
    vtkWriter.write(0.0);

    // the assembler with time loop for instationary problem
    using Assembler = StaggeredFVAssembler<TypeTag, DiffMethod::numeric>;
    auto assembler = std::make_shared<Assembler>(problem, fvGridGeometry, gridVariables);

    // the linear solver
    using LinearSolver = Dumux::UMFPackBackend;
    auto linearSolver = std::make_shared<LinearSolver>();

    // the non-linear solver
    using NewtonSolver = Dumux::NewtonSolver<Assembler, LinearSolver>;
    NewtonSolver nonLinearSolver(assembler, linearSolver);

    const bool printL2Error = getParam<bool>("Problem.PrintL2Error");

    // linearize & solve
    Dune::Timer timer;
    nonLinearSolver.solve(x);

    if (printL2Error)
    {
        using Indices = typename GetPropType<TypeTag, Properties::ModelTraits>::Indices;
        using ModelTraits = GetPropType<TypeTag, Properties::ModelTraits>;
        using PrimaryVariables = GetPropType<TypeTag, Properties::PrimaryVariables>;

        using Scalar = GetPropType<TypeTag, Properties::Scalar>;
        using L2Error = NavierStokesTestL2Error<Scalar, ModelTraits, PrimaryVariables>;
        const auto l2error = L2Error::calculateL2Error(*problem, x);
        const int numCellCenterDofs = fvGridGeometry->numCellCenterDofs();
        const int numFaceDofs = fvGridGeometry->numFaceDofs();
        std::cout << std::setprecision(8) << "** L2 error (abs/rel) for "
                  << std::setw(6) << numCellCenterDofs << " cc dofs and " << numFaceDofs << " face dofs (total: " << numCellCenterDofs + numFaceDofs << "): "
                  << std::scientific
                  << "L2(p) = " << l2error.first[Indices::pressureIdx] << " / " << l2error.second[Indices::pressureIdx]
                  << " , L2(vx) = " << l2error.first[Indices::velocityXIdx] << " / " << l2error.second[Indices::velocityXIdx]
                  << " , L2(vy) = " << l2error.first[Indices::velocityYIdx] << " / " << l2error.second[Indices::velocityYIdx]
                  << std::endl;
    }

    // write vtk output
    analyticalSolution = createAnalyticalSolution(*problem);
    vtkWriter.write(1.0);

    timer.stop();

    const auto& comm = Dune::MPIHelper::getCollectiveCommunication();
    std::cout << "Simulation took " << timer.elapsed() << " seconds on "
              << comm.size() << " processes.\n"
              << "The cumulative CPU time was " << timer.elapsed()*comm.size() << " seconds.\n";

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
