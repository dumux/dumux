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
 * \brief 3D Channel flow test for the staggered grid (Navier-)Stokes model
 */

#include <config.h>

#include <ctime>
#include <iostream>

#include <dune/common/parallel/mpihelper.hh>
#include <dune/common/timer.hh>
#include <dune/grid/io/file/dgfparser/dgfexception.hh>
#include <dune/grid/io/file/vtk.hh>
#include <dune/istl/io.hh>

#include <dumux/assembly/staggeredfvassembler.hh>
#include <dumux/assembly/diffmethod.hh>
#include <dumux/common/dumuxmessage.hh>
#include <dumux/common/parameters.hh>
#include <dumux/common/properties.hh>
// #include <dumux/freeflow/navierstokes/staggered/fluxoversurface.hh>
#include <dumux/io/grid/gridmanager_sub.hh>
#include <dumux/io/grid/gridmanager_yasp.hh>
#include <dumux/linear/seqsolverbackend.hh>

#include <dumux/io/vtkoutputmodule.hh>
#include <dumux/freeflow/navierstokes/velocityoutput.hh>
#include <dumux/multidomain/fvassembler.hh>
#include <dumux/multidomain/traits.hh>
#include <dumux/multidomain/staggeredfreeflow/couplingmanager.hh>
#include <dumux/multidomain/newtonsolver.hh>

#include <dumux/freeflow/navierstokes/fluxoveraxisalignedplane.hh>

#include "problem_new.hh"

namespace Dumux::Properties{

// Set the problem property
template<class TypeTag>
struct CouplingManager<TypeTag, TTag::ThreeDChannelTest>
{
private:
    using Traits = MultiDomainTraits<TTag::ThreeDChannelTestMomentum, TTag::ThreeDChannelTestMass>;
public:
    using type = StaggeredFreeFlowCouplingManager<Traits>;
};

}

int main(int argc, char** argv) try
{
    using namespace Dumux;

    // define the type tag for this problem
    using MomentumTypeTag = Properties::TTag::ThreeDChannelTestMomentum;
    using MassTypeTag = Properties::TTag::ThreeDChannelTestMass;

    // initialize MPI, finalize is done automatically on exit
    const auto& mpiHelper = Dune::MPIHelper::instance(argc, argv);

    // print dumux start message
    if (mpiHelper.rank() == 0)
        DumuxMessage::print(/*firstCall=*/true);

    // parse command line arguments and input file
    Parameters::init(argc, argv);

    // create a grid
    using Grid = GetPropType<MomentumTypeTag, Properties::Grid>;
    Dumux::GridManager<Grid> gridManager;

#if HAVE_DUNE_SUBGRID
    using Scalar = GetPropType<MomentumTypeTag, Properties::Scalar>;
    const bool isStaircaseGeometry = getParam<bool>("Problem.IsStaircaseGeometry", false);

    auto selector = [&](const auto& element)
    {
        if (!isStaircaseGeometry)
            return true;

        const Scalar deltaX = 0.003;
        const Scalar deltaZ = 0.000075;
        const Scalar deltaY = 0.0003;

        const Scalar eps = 1e-8;
        const auto globalPos = element.geometry().center();

        return globalPos[2] > (deltaZ/deltaX * globalPos[0] + deltaZ/deltaY * globalPos[1] - deltaZ + eps);
    };

    gridManager.init(selector, "Internal");
#else
    gridManager.init();
#endif

    ////////////////////////////////////////////////////////////
    // run instationary non-linear problem on this grid
    ////////////////////////////////////////////////////////////

    // we compute on the leaf grid view
    const auto& leafGridView = gridManager.grid().leafGridView();

    // create the finite volume grid geometry
    using MomentumGridGeometry = GetPropType<MomentumTypeTag, Properties::GridGeometry>;
    auto momentumGridGeometry = std::make_shared<MomentumGridGeometry>(leafGridView);
    momentumGridGeometry->update();

    using MassGridGeometry = GetPropType<MassTypeTag, Properties::GridGeometry>;
    auto massGridGeometry = std::make_shared<MassGridGeometry>(leafGridView);
    massGridGeometry->update();

    // the coupling manager
    using Traits = MultiDomainTraits<MomentumTypeTag, MassTypeTag>;
    using CouplingManager = StaggeredFreeFlowCouplingManager<Traits>;

    auto couplingManager = std::make_shared<CouplingManager>();

    // the problem (boundary conditions)
    using MomentumProblem = GetPropType<MomentumTypeTag, Properties::Problem>;
    auto momentumProblem = std::make_shared<MomentumProblem>(momentumGridGeometry, couplingManager);

    using MassProblem = GetPropType<MassTypeTag, Properties::Problem>;
    auto massProblem = std::make_shared<MassProblem>(massGridGeometry, couplingManager);

    // the solution vector
    constexpr auto momentumIdx = Dune::index_constant<0>();
    constexpr auto massIdx = Dune::index_constant<1>();
    using SolutionVector = typename Traits::SolutionVector;
    SolutionVector x;
    x[momentumIdx].resize(momentumGridGeometry->numDofs());
    x[massIdx].resize(massGridGeometry->numDofs());

    // the grid variables
    using MomentumGridVariables = GetPropType<MomentumTypeTag, Properties::GridVariables>;
    auto momentumGridVariables = std::make_shared<MomentumGridVariables>(momentumProblem, momentumGridGeometry);

    using MassGridVariables = GetPropType<MassTypeTag, Properties::GridVariables>;
    auto massGridVariables = std::make_shared<MassGridVariables>(massProblem, massGridGeometry);

    couplingManager->init(momentumProblem, massProblem, std::make_tuple(momentumGridVariables, massGridVariables), x);
    massGridVariables->init(x[massIdx]);
    momentumGridVariables->init(x[momentumIdx]);

    using Assembler = MultiDomainFVAssembler<Traits, CouplingManager, DiffMethod::numeric>;
    auto assembler = std::make_shared<Assembler>(std::make_tuple(momentumProblem, massProblem),
                                                 std::make_tuple(momentumGridGeometry, massGridGeometry),
                                                 std::make_tuple(momentumGridVariables, massGridVariables),
                                                 couplingManager);

    // intialize the vtk output module
    using IOFields = GetPropType<MassTypeTag, Properties::IOFields>;
    VtkOutputModule vtkWriter(*massGridVariables, x[massIdx], massProblem->name());
    IOFields::initOutputModule(vtkWriter); // Add model specific output fields
    vtkWriter.addVelocityOutput(std::make_shared<NavierStokesVelocityOutput<MassGridVariables>>());

    // the linear solver
    using LinearSolver = Dumux::UMFPackBackend;
    auto linearSolver = std::make_shared<LinearSolver>();

    // the non-linear solver
    using NewtonSolver = MultiDomainNewtonSolver<Assembler, LinearSolver, CouplingManager>;
    NewtonSolver nonLinearSolver(assembler, linearSolver, couplingManager);

    // set up three planes over which fluxes are calculated
    FluxOverAxisAlignedPlane<MassGridVariables,
                             std::decay_t<decltype(x[massIdx])>,
                             GetPropType<MassTypeTag, Properties::LocalResidual>> flux(*massGridVariables, x[massIdx]);

    using GridView = typename GetPropType<MassTypeTag, Properties::GridGeometry>::GridView;
    using GlobalPosition = Dune::FieldVector<Scalar, GridView::dimensionworld>;

    const Scalar xMin = massGridGeometry->bBoxMin()[0];
    const Scalar xMax = massGridGeometry->bBoxMax()[0];
    const Scalar yMin = massGridGeometry->bBoxMin()[1];
    const Scalar yMax = massGridGeometry->bBoxMax()[1];

#if GRID_DIM == 3
    const Scalar zMin = massGridGeometry->bBoxMin()[2];
    const Scalar zMax = massGridGeometry->bBoxMax()[2];
#endif

    // the first plane is at the inlet
#if GRID_DIM == 3
    const auto inletLowerLeft = GlobalPosition{xMin, yMin, zMin};
    const auto inletUpperRight = GlobalPosition{xMin, yMax, zMax};
    flux.addPlane("inlet", inletLowerLeft, inletUpperRight, 0);
#else
    const auto inletLowerLeft = GlobalPosition{xMin, yMin};
    const auto inletUpperRight = GlobalPosition{xMin, yMax};
    flux.addPlane("inlet", inletLowerLeft, inletUpperRight, 0);
#endif

    // the second plane is at the middle of the channel
    const Scalar planePosMiddleX = xMin + 0.5*(xMax - xMin);
#if GRID_DIM == 3
    const auto middleLowerLeft = GlobalPosition{planePosMiddleX, yMin, zMin};
    const auto middleUpperRight = GlobalPosition{planePosMiddleX, yMax, zMax};
    flux.addPlane("middle", middleLowerLeft, middleUpperRight, 0);
#else
    const auto middleLowerLeft = GlobalPosition{planePosMiddleX, yMin};
    const auto middleUpperRight = GlobalPosition{planePosMiddleX, yMax};
    flux.addPlane("middle", middleLowerLeft, middleUpperRight, 0);
#endif

    // The last plane is placed at the outlet of the channel.
#if GRID_DIM == 3
    const auto outletLowerLeft = GlobalPosition{xMax, yMin, zMin};
    const auto outletUpperRight = GlobalPosition{xMax, yMax, zMax};
    flux.addPlane("outlet", outletLowerLeft, outletUpperRight, 0);
#else
    const auto outletLowerLeft = GlobalPosition{xMax, yMin};
    const auto outletUpperRight = GlobalPosition{xMax, yMax};
    flux.addPlane("outlet", outletLowerLeft, outletUpperRight, 0);
#endif

    // linearize & solve
    Dune::Timer timer;
    nonLinearSolver.solve(x);

    // write vtk output
    vtkWriter.write(1.0);

    // calculate and print mass fluxes over the planes
    flux.calculateAllScalarFluxes();

    std::cout << "mass flux at inlet is: " << flux.netFlux("inlet") << std::endl;
    std::cout << "mass flux at middle is: " << flux.netFlux("middle") << std::endl;
    std::cout << "mass flux at outlet is: " << flux.netFlux("outlet") << std::endl;


    timer.stop();

    std::cout << "analytical mass flux: " << massProblem->analyticalFlux()*1e3 << std::endl;

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
