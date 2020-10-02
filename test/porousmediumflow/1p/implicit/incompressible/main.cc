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
 * \ingroup OnePTests
 * \brief Test for the one-phase CC model
 */

#include <config.h>

#include <ctime>
#include <iostream>

#include <dune/common/float_cmp.hh>
#include <dune/common/parallel/mpihelper.hh>
#include <dune/common/timer.hh>
#include <dune/grid/io/file/vtk.hh>

#include <dumux/linear/seqsolverbackend.hh>
#include <dumux/linear/pdesolver.hh>

#include <dumux/common/properties.hh>
#include <dumux/common/parameters.hh>
#include <dumux/common/dumuxmessage.hh>

#include <dumux/io/vtkoutputmodule.hh>
#include <dumux/io/grid/gridmanager.hh>
#include <dumux/io/pointcloudvtkwriter.hh>

#include <dumux/discretization/method.hh>
#include <dumux/discretization/cellcentered/mpfa/scvgradients.hh>

#include <dumux/assembly/fvassembler.hh>

#include "problem.hh"
#include "../internaldirichlet/problem.hh"

//! Function to write out the scv-wise velocities (overload for mpfa)
template<class GridGeometry, class GridVariables, class Sol,
         std::enable_if_t<GridGeometry::discMethod == Dumux::DiscretizationMethod::ccmpfa, int> = 0>
void writeMpfaVelocities(const GridGeometry& gridGeometry,
                         const GridVariables& gridVariables,
                         const Sol& x)
{
    using Scalar = typename GridVariables::Scalar;
    using GlobalPos = typename GridGeometry::SubControlVolume::GlobalPosition;

    const auto velocities = Dumux::CCMpfaScvGradients::computeVelocities(gridGeometry, gridVariables, x, /*phaseIdx*/0);
    Dumux::PointCloudVtkWriter<Scalar, GlobalPos> writer(velocities.first);
    writer.addPointData(velocities.second, "velocity (m/s)");
    writer.write("mpfa_scv_velocities");
}

//! Function to write out the scv-wise velocities (overload for NOT mpfa)
template<class GridGeometry, class GridVariables, class Sol,
         std::enable_if_t<GridGeometry::discMethod != Dumux::DiscretizationMethod::ccmpfa, int> = 0>
void writeMpfaVelocities(const GridGeometry& gridGeometry,
                         const GridVariables& gridVariables,
                         const Sol& x)
{}

int main(int argc, char** argv)
{
    using namespace Dumux;

    using TypeTag = Properties::TTag::TYPETAG;

    // initialize MPI, finalize is done automatically on exit
    const auto& mpiHelper = Dune::MPIHelper::instance(argc, argv);

    // print dumux start message
    if (mpiHelper.rank() == 0)
        DumuxMessage::print(/*firstCall=*/true);

    ////////////////////////////////////////////////////////////
    // parse the command line arguments and input file
    ////////////////////////////////////////////////////////////

    // parse command line arguments
    Parameters::init(argc, argv);

    //////////////////////////////////////////////////////////////////////
    // try to create a grid (from the given grid file or the input file)
    /////////////////////////////////////////////////////////////////////
    using Grid = GetPropType<TypeTag, Properties::Grid>;
    GridManager<Grid> gridManager;
    gridManager.init();

    // we compute on the leaf grid view
    const auto& leafGridView = gridManager.grid().leafGridView();

    // start timer
    Dune::Timer timer;

    // create the finite volume grid geometry
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    auto gridGeometry = std::make_shared<GridGeometry>(leafGridView);
    gridGeometry->update();

    // the problem (boundary conditions)
    using Problem = GetPropType<TypeTag, Properties::Problem>;
    auto problem = std::make_shared<Problem>(gridGeometry);

    // the solution vector
    using SolutionVector = GetPropType<TypeTag, Properties::SolutionVector>;
    SolutionVector x(gridGeometry->numDofs());

    // the grid variables
    using GridVariables = GetPropType<TypeTag, Properties::GridVariables>;
    auto gridVariables = std::make_shared<GridVariables>(problem, gridGeometry);
    gridVariables->init(x);

    // intialize the vtk output module
    VtkOutputModule<GridVariables, SolutionVector> vtkWriter(*gridVariables, x, problem->name());
    using VelocityOutput = GetPropType<TypeTag, Properties::VelocityOutput>;
    vtkWriter.addVelocityOutput(std::make_shared<VelocityOutput>(*gridVariables));
    using IOFields = GetPropType<TypeTag, Properties::IOFields>;
    IOFields::initOutputModule(vtkWriter); // Add model specific output fields
    vtkWriter.write(0.0);

    // create assembler & linear solver
    using Assembler = FVAssembler<TypeTag, NUMDIFFMETHOD>;
    auto assembler = std::make_shared<Assembler>(problem, gridGeometry, gridVariables);

    using LinearSolver = SSORCGBackend;
    auto linearSolver = std::make_shared<LinearSolver>();

    // solver the linear problem
    LinearPDESolver solver(assembler, linearSolver);
    solver.solve(x);

    // output result to vtk
    vtkWriter.write(1.0);

    timer.stop();

    const bool checkIsConstantVelocity = getParam<bool>("Problem.CheckIsConstantVelocity", false);
    if(checkIsConstantVelocity)
    {
        // instantiate the velocity output
        VelocityOutput velocityOutput(*gridVariables);
        using VelocityVector = typename VelocityOutput::VelocityVector;
        VelocityVector velocity;

        constexpr bool isBox = GridGeometry::discMethod == Dumux::DiscretizationMethod::box;
        constexpr int dimWorld = GridGeometry::GridView::dimensionworld;
        const auto numCells = leafGridView.size(0);
        const auto numDofs = gridGeometry->numDofs();
        auto numVelocities = (isBox && dimWorld == 1) ? numCells : numDofs;

        velocity.resize(numVelocities);

        const auto exactVel = problem->velocity();

        for (const auto& element : elements(leafGridView, Dune::Partitions::interior))
        {
            const auto eIdx = gridGeometry->elementMapper().index(element);

            auto fvGeometry = localView(*gridGeometry);
            auto elemVolVars = localView(gridVariables->curGridVolVars());
            auto elemFluxVarsCache = localView(gridVariables->gridFluxVarsCache());

            fvGeometry.bind(element);
            elemVolVars.bind(element, fvGeometry, x);
            elemFluxVarsCache.bind(element, fvGeometry, elemVolVars);

            velocityOutput.calculateVelocity(velocity, element, fvGeometry, elemVolVars, elemFluxVarsCache, 0);

            using Scalar = Grid::ctype;
            // the y-component of the velocity should be exactly reproduced
            // the x-component should be zero
            // use a relative comparison for the y-component and an absolute one for the x-component
            if(Dune::FloatCmp::ne(velocity[eIdx][dimWorld-1], exactVel[dimWorld-1], /*eps*/1e-8) ||
               Dune::FloatCmp::ne<Scalar, Dune::FloatCmp::absolute>(velocity[eIdx][0], exactVel[0], /*eps*/1e-10))
                    DUNE_THROW(Dune::InvalidStateException, "Velocity is not exactly reproduced");
        }
    }

    const auto& comm = Dune::MPIHelper::getCollectiveCommunication();
    if (mpiHelper.rank() == 0)
        std::cout << "Simulation took " << timer.elapsed() << " seconds on "
                  << comm.size() << " processes.\n"
                  << "The cumulative CPU time was " << timer.elapsed()*comm.size() << " seconds.\n";

    // For the mpfa test, write out the gradients in the scv centers
    if (getParam<bool>("IO.WriteMpfaVelocities", false))
        writeMpfaVelocities(*gridGeometry, *gridVariables, x);

    if (mpiHelper.rank() == 0)
        Parameters::print();

    return 0;

}
