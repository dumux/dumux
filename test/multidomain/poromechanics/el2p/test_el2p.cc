// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   See the file COPYING for full copying permissions.                      *
 *                                                                           *
 *   This program is free software: you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation, either version 2 of the License, or       *
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
 * \brief test for a single-phase elastic coupled model.
 */
#include <config.h>
#include <ctime>
#include <iostream>

#include <dune/common/parallel/mpihelper.hh>
#include <dune/common/timer.hh>

#include "2pproblem.hh"
#include "poroelasticproblem.hh"

#include <dumux/common/properties.hh>
#include <dumux/common/parameters.hh>
#include <dumux/common/dumuxmessage.hh>
#include <dumux/common/defaultusagemessage.hh>

#include <dumux/assembly/diffmethod.hh>

#include <dumux/linear/seqsolverbackend.hh>
#include <dumux/multidomain/newtonsolver.hh>
#include <dumux/multidomain/fvassembler.hh>
#include <dumux/multidomain/traits.hh>

#include <dumux/geomechanics/poroelastic/couplingmanager.hh>

#include <dumux/io/vtkoutputmodule.hh>
#include <dumux/io/grid/gridmanager.hh>

// set the coupling manager property in the sub-problems
namespace Dumux {
namespace Properties {

SET_PROP(TwoPSub, CouplingManager)
{
private:
    // define traits etc. as below in main
    using Traits = MultiDomainTraits<TTAG(TwoPSub), TTAG(PoroElasticSub)>;
public:
    using type = PoroMechanicsCouplingManager< Traits >;
};

SET_PROP(PoroElasticSub, CouplingManager)
{
private:
    // define traits etc. as below in main
    using Traits = MultiDomainTraits<TTAG(TwoPSub), TTAG(PoroElasticSub)>;
public:
    using type = PoroMechanicsCouplingManager< Traits >;
};
} // end namespace Properties
} // end namespace Dumux

int main(int argc, char** argv) try
{
    using namespace Dumux;

    ////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////

    // initialize MPI, finalize is done automatically on exit
    const auto& mpiHelper = Dune::MPIHelper::instance(argc, argv);

    // print dumux start message
    if (mpiHelper.rank() == 0)
        DumuxMessage::print(/*firstCall=*/true);

    // initialize parameter tree
    Parameters::init(argc, argv);

    //////////////////////////////////////////////////////////////////////
    // try to create a grid (from the given grid file or the input file)
    /////////////////////////////////////////////////////////////////////
    using TwoPTypeTag = TTAG(TwoPSub);
    using PoroMechTypeTag = TTAG(PoroElasticSub);

    // we simply extract the grid creator from one of the type tags
    using GridManager = Dumux::GridManager<typename GET_PROP_TYPE(TwoPTypeTag, Grid)>;
    GridManager gridManager;
    gridManager.init();

    ////////////////////////////////////////////////////////////
    // run stationary non-linear problem on this grid
    ////////////////////////////////////////////////////////////

    // we compute on the leaf grid view
    const auto& leafGridView = gridManager.grid().leafGridView();

    // create the finite volume grid geometries
    using TwoPFVGridGeometry = typename GET_PROP_TYPE(TwoPTypeTag, FVGridGeometry);
    using PoroMechFVGridGeometry = typename GET_PROP_TYPE(PoroMechTypeTag, FVGridGeometry);
    auto twoPFvGridGeometry = std::make_shared<TwoPFVGridGeometry>(leafGridView);
    auto poroMechFvGridGeometry = std::make_shared<PoroMechFVGridGeometry>(leafGridView);
    twoPFvGridGeometry->update();
    poroMechFvGridGeometry->update();

    // the problems (boundary conditions)
    using TwoPProblem = typename GET_PROP_TYPE(TwoPTypeTag, Problem);
    using PoroMechProblem = typename GET_PROP_TYPE(PoroMechTypeTag, Problem);
    auto twoPProblem = std::make_shared<TwoPProblem>(twoPFvGridGeometry, "TwoP");
    auto poroMechProblem = std::make_shared<PoroMechProblem>(poroMechFvGridGeometry, "PoroElastic");

    // the solution vectors
    using Traits = MultiDomainTraits<TwoPTypeTag, PoroMechTypeTag>;
    using SolutionVector = typename Traits::SolutionVector;
    SolutionVector x;

    static const auto twoPId = Traits::template DomainIdx<0>();
    static const auto poroMechId = Traits::template DomainIdx<1>();
    x[twoPId].resize(twoPFvGridGeometry->numDofs());
    x[poroMechId].resize(poroMechFvGridGeometry->numDofs());
    twoPProblem->applyInitialSolution(x[twoPId]);
    poroMechProblem->applyInitialSolution(x[poroMechId]);
    SolutionVector xOld = x;

    // the coupling manager
    using CouplingManager = PoroMechanicsCouplingManager<Traits>;
    auto couplingManager = std::make_shared<CouplingManager>();
    couplingManager->init(twoPProblem, poroMechProblem, x);

    // set coupling manager pointer in sub-problems
    // (in 1p problem the coupling enters in the parameters)
    twoPProblem->spatialParams().setCouplingManager(couplingManager);
    poroMechProblem->setCouplingManager(couplingManager);

    // the grid variables
    using TwoPGridVariables = typename GET_PROP_TYPE(TwoPTypeTag, GridVariables);
    using PoroMechGridVariables = typename GET_PROP_TYPE(PoroMechTypeTag, GridVariables);
    auto twoPGridVariables = std::make_shared<TwoPGridVariables>(twoPProblem, twoPFvGridGeometry);
    auto poroMechGridVariables = std::make_shared<PoroMechGridVariables>(poroMechProblem, poroMechFvGridGeometry);
    twoPGridVariables->init(x[twoPId]);
    poroMechGridVariables->init(x[poroMechId]);

    // get some time loop parameters
    using Scalar = typename GET_PROP_TYPE(TwoPTypeTag, Scalar);
    const auto tEnd = getParam<Scalar>("TimeLoop.TEnd");
    const auto maxDT = getParam<Scalar>("TimeLoop.MaxTimeStepSize");
    auto dt = getParam<Scalar>("TimeLoop.DtInitial");

    // intialize the dune vtk writers
    using TwoPVtkOutputModule = Dumux::VtkOutputModule<TwoPGridVariables, typename GET_PROP_TYPE(TwoPTypeTag, SolutionVector)>;
    using PoroMechVtkOutputModule = Dumux::VtkOutputModule<PoroMechGridVariables, typename GET_PROP_TYPE(PoroMechTypeTag, SolutionVector)>;
    TwoPVtkOutputModule twoPVtkWriter(*twoPGridVariables, x[twoPId], twoPProblem->name());
    PoroMechVtkOutputModule poroMechVtkWriter(*poroMechGridVariables, x[poroMechId], poroMechProblem->name());

    // add output fields to writers
    using TwoPOutputFields = typename GET_PROP_TYPE(TwoPTypeTag, IOFields);
    using PoroMechOutputFields = typename GET_PROP_TYPE(PoroMechTypeTag, IOFields);
    TwoPOutputFields::init(twoPVtkWriter);
    PoroMechOutputFields::init(poroMechVtkWriter);

    // write initial solution
    twoPVtkWriter.write(0.0);
    poroMechVtkWriter.write(0.0);

    //instantiate time loop
    auto timeLoop = std::make_shared<TimeLoop<Scalar>>(0.0, dt, tEnd);
    timeLoop->setMaxTimeStepSize(maxDT);

    // the assembler
    using Assembler = MultiDomainFVAssembler<Traits, CouplingManager, DiffMethod::numeric, /*implicit?*/true>;
    auto assembler = std::make_shared<Assembler>( std::make_tuple(twoPProblem, poroMechProblem),
                                                  std::make_tuple(twoPFvGridGeometry, poroMechFvGridGeometry),
                                                  std::make_tuple(twoPGridVariables, poroMechGridVariables),
                                                  couplingManager, timeLoop );

    // the linear solver
    using LinearSolver = UMFPackBackend;
    auto linearSolver = std::make_shared<LinearSolver>();

    // the non-linear solver
    using NewtonSolver = Dumux::MultiDomainNewtonSolver<Assembler, LinearSolver, CouplingManager>;
    NewtonSolver nonLinearSolver(assembler, linearSolver, couplingManager);

    // time loop
    timeLoop->start(); do
    {
        // set previous solution for storage evaluations
        assembler->setPreviousSolution(xOld);

        // solve the non-linear system with time step control
        nonLinearSolver.solve(x, *timeLoop);

        // make the new solution the old solution
        xOld = x;

        // advance to the time loop to the next step
        timeLoop->advanceTimeStep();

        // write vtk output
        twoPVtkWriter.write(timeLoop->time());
        poroMechVtkWriter.write(timeLoop->time());

        // report statistics of this time step
        timeLoop->reportTimeStep();

        using TwoPPrimaryVariables = typename GET_PROP_TYPE(TwoPTypeTag, PrimaryVariables);
        TwoPPrimaryVariables storage(0);
        const auto& twoPLocalResidual = assembler->localResidual(twoPId);
        for (const auto& element : elements(leafGridView, Dune::Partitions::interior))
        {
            auto storageVec = twoPLocalResidual.evalStorage(*twoPProblem, element, *twoPFvGridGeometry, *twoPGridVariables, x[twoPId]);
            storage += storageVec[0];
        }
        std::cout << "time, mass CO2 (kg), mass brine (kg):" << std::endl;
        std::cout << timeLoop->time() << " , " << storage[1] << " , " << storage[0] << std::endl;
        std::cout << "***************************************" << std::endl;

        // set new dt as suggested by the Newton solver
        timeLoop->setTimeStepSize(nonLinearSolver.suggestTimeStepSize(timeLoop->timeStepSize()));

    } while (!timeLoop->finished());


    // output some Newton statistics
    nonLinearSolver.report();

    timeLoop->finalize(leafGridView.comm());

    ////////////////////////////////////////////////////////////
    // finalize, print dumux message to say goodbye
    ////////////////////////////////////////////////////////////
    if (mpiHelper.rank() == 0)
    {
        Parameters::print();
        DumuxMessage::print(/*firstCall=*/false);
    }

    return 0;

}
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
