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

#include <config.h>

#include <iostream>

#include <dune/common/parallel/mpihelper.hh>

#include <dumux/common/properties.hh>
#include <dumux/common/parameters.hh>
#include <dumux/common/partial.hh>
#include <dumux/linear/seqsolverbackend.hh>
#include <dumux/assembly/diffmethod.hh>
#include <dumux/io/vtkoutputmodule.hh>
#include <dumux/io/grid/gridmanager_yasp.hh>

#include <dumux/multidomain/fvassembler.hh>
#include <dumux/multidomain/traits.hh>
#include <dumux/multidomain/staggeredfreeflow/couplingmanager.hh>
#include <dumux/multidomain/newtonsolver.hh>

#include <dumux/io/vtk/intersectionwriter.hh>
#include <dumux/io/vtkoutputmodule.hh>
#include <dumux/freeflow/navierstokes/velocityoutput.hh>

#include "problem_new.hh"
#include "../../analyticalsolution.hh"
#include "../../l2error.hh"
#include <dumux/freeflow/navierstokes/fluxoveraxisalignedplane.hh>

template<class MomentumProblem, class MassProblem, class SolutionVector>
void printL2Error(const MomentumProblem& momentumProblem,
                  const MassProblem& massProblem,
                  const SolutionVector& x)
{
    const auto velocityL2error = calculateL2Error(momentumProblem, x[Dune::index_constant<0>()]);
    const auto pressureL2error = calculateL2Error(massProblem, x[Dune::index_constant<1>()]);
    static constexpr auto pressureIdx = MassProblem::Indices::pressureIdx;
    static constexpr auto velocityXIdx = MomentumProblem::Indices::velocityXIdx;
    static constexpr auto velocityYIdx = MomentumProblem::Indices::velocityYIdx;
    const auto numCellCenterDofs = massProblem.gridGeometry().numDofs();
    const auto numFaceDofs = momentumProblem.gridGeometry().numDofs();

    std::cout << std::setprecision(8) << "** L2 error (abs/rel) for "
                << std::setw(6) << numCellCenterDofs << " cc dofs and " << numFaceDofs << " face dofs (total: " << numCellCenterDofs + numFaceDofs << "): "
                << std::scientific
                << "L2(p) = " << pressureL2error.absolute[pressureIdx] << " / " << pressureL2error.relative[pressureIdx]
                << " , L2(vx) = " << velocityL2error.absolute[velocityXIdx] << " / " << velocityL2error.relative[velocityXIdx]
                << " , L2(vy) = " << velocityL2error.absolute[velocityYIdx] << " / " << velocityL2error.relative[velocityYIdx]
                << std::endl;

    // write the norm into a log file
    std::ofstream logFile;
    logFile.open(momentumProblem.name() + ".log", std::ios::app);
    logFile << "[ConvergenceTest] L2(p) = " << pressureL2error.absolute[pressureIdx] << " L2(vx) = "
            << velocityL2error.absolute[velocityXIdx] << " L2(vy) = " << velocityL2error.absolute[velocityYIdx] << std::endl;
    logFile.close();
}

namespace Dumux::Properties{

// Set the problem property
template<class TypeTag>
struct CouplingManager<TypeTag, TTag::PipeFlow>
{
private:
    using Traits = MultiDomainTraits<TTag::PipeFlowMomentum, TTag::PipeFlowMass>;
public:
    using type = StaggeredFreeFlowCouplingManager<Traits>;
};

}

int main(int argc, char** argv) try
{
    using namespace Dumux;

    // initialize MPI, finalize is done automatically on exit
    const auto& mpiHelper = Dune::MPIHelper::instance(argc, argv);

    // parse command line arguments and input file
    Parameters::init(argc, argv);

    // Define the sub problem type tags
    using MomentumTypeTag = Properties::TTag::PipeFlowMomentum;
    using MassTypeTag = Properties::TTag::PipeFlowMass;

    // try to create a grid (from the given grid file or the input file)
    Dumux::GridManager<GetPropType<MomentumTypeTag, Properties::Grid>> gridManager;
    gridManager.init();

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

    // momentumProblem->applyInitialSolution(x[momentumIdx]);
    // massProblem->applyInitialSolution(x[massIdx]);

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
    const auto exactPressure = getScalarAnalyticalSolution(*massProblem)[GetPropType<MassTypeTag, Properties::ModelTraits>::Indices::pressureIdx];
    const auto exactVelocity = getVelocityAnalyticalSolution(*momentumProblem);
    vtkWriter.addField(exactPressure, "pressureExact");
    vtkWriter.addField(exactVelocity, "velocityExact");

    vtkWriter.write(0.0);

    // the linear solver
    using LinearSolver = Dumux::UMFPackBackend;
    auto linearSolver = std::make_shared<LinearSolver>();

    // the non-linear solver
    using NewtonSolver = MultiDomainNewtonSolver<Assembler, LinearSolver, CouplingManager>;
    NewtonSolver nonLinearSolver(assembler, linearSolver, couplingManager);

    nonLinearSolver.solve(x);
    vtkWriter.write(1.0);

    // set up three planes over which fluxes are calculated
    FluxOverAxisAlignedPlane<MassGridVariables,
                             std::decay_t<decltype(x[massIdx])>,
                             GetPropType<MassTypeTag, Properties::LocalResidual>> flux(*massGridVariables, x[massIdx]);

    using GridView = typename GetPropType<MassTypeTag, Properties::GridGeometry>::GridView;
    using Scalar = double;
    using GlobalPosition = Dune::FieldVector<Scalar, GridView::dimensionworld>;

    const Scalar xMin = massGridGeometry->bBoxMin()[0];
    const Scalar xMax = massGridGeometry->bBoxMax()[0];
    const Scalar yMin = massGridGeometry->bBoxMin()[1];
    const Scalar yMax = massGridGeometry->bBoxMax()[1];

    const auto inletLowerLeft = GlobalPosition{xMin, yMin};
    const auto inletUpperRight = GlobalPosition{xMax, yMin};
    flux.addPlane("inlet", inletLowerLeft, inletUpperRight, 1);

    const Scalar planePosMiddleY = yMin + 0.5*(yMax - yMin);
    const auto middleLowerLeft = GlobalPosition{xMin, planePosMiddleY};
    const auto middleUpperRight = GlobalPosition{xMax, planePosMiddleY};
    flux.addPlane("middle", middleLowerLeft, middleUpperRight, 1);

    const auto outletLowerLeft = GlobalPosition{xMin, yMax};
    const auto outletUpperRight = GlobalPosition{xMax, yMax};
    flux.addPlane("outlet", outletLowerLeft, outletUpperRight, 1);

    // calculate and print mass fluxes over the planes
    flux.calculateAllScalarFluxes();

    std::cout << "mass flux at inlet is: " << flux.netFlux("inlet") << std::endl;
    std::cout << "mass flux at middle is: " << flux.netFlux("middle") << std::endl;
    std::cout << "mass flux at outlet is: " << flux.netFlux("outlet") << std::endl;

    // calculate and print L2 error
    printL2Error(*momentumProblem, *massProblem, x);

    // print dumux end message
    if (mpiHelper.rank() == 0)
        Parameters::print();

    return 0;

} // end main
catch (const Dumux::ParameterException &e)
{
    std::cerr << std::endl << e << " ---> Abort!" << std::endl;
    return 1;
}
catch (const Dune::DGFException & e)
{
    std::cerr << "DGF exception thrown (" << e <<
                 "). Most likely, the DGF file name is wrong "
                 "or the DGF file is corrupted, "
                 "e.g. missing hash at end of file or wrong number (dimensions) of entries."
                 << " ---> Abort!" << std::endl;
    return 2;
}
catch (const Dune::Exception &e)
{
    std::cerr << "Dune reported error: " << e << " ---> Abort!" << std::endl;
    return 3;
}
