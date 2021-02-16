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
 * \ingroup PoromechanicsTests
 * \brief Test for a single-phase elastic coupled model.
 */

#include <config.h>
#include <ctime>
#include <iostream>

#include <dune/common/parallel/mpihelper.hh>
#include <dune/common/timer.hh>

#include "problem_2p.hh"
#include "problem_poroelastic.hh"

#include <dumux/common/properties.hh>

namespace Dumux::Properties {
    // property forward declaration
    // TODO: We need this in MDTraits such that the MDAssembler
    //       can extract the local operator from the sub-typetags.
    //       Thus, we need to introduce this property for all models?
    //       Should this be circumvented somehow? (In standard models one can
    //       select the operator in the main file so far...)
    template<class TypeTag, class MyTypeTag>
    struct LocalOperator { using type = UndefinedProperty; };
} // end namespace Dumux::Properties

#include <dumux/common/parameters.hh>
#include <dumux/common/dumuxmessage.hh>

#include <dumux/assembly/diffmethod.hh>

#include <dumux/linear/seqsolverbackend.hh>
#include <dumux/multidomain/newtonsolver.hh>
#include <dumux/multidomain/assembler.hh>
#include <dumux/multidomain/fvgridvariables.hh>
#include <dumux/multidomain/traits.hh>

#include <dumux/assembly/fv/localoperator.hh>
#include <dumux/porousmediumflow/immiscible/operators.hh>
#include <dumux/geomechanics/fvoperators.hh>
#include <dumux/timestepping/multistagemethods.hh>

#include <dumux/geomechanics/poroelastic/couplingmanager.hh>

#include <dumux/io/vtkoutputmodule.hh>
#include <dumux/io/grid/gridmanager.hh>

// set the coupling manager property in the sub-problems
namespace Dumux {
namespace Properties {

template<class TypeTag>
struct CouplingManager<TypeTag, TTag::TwoPSub>
{
private:
    // define traits etc. as below in main
    using Traits = MultiDomainTraits<Properties::TTag::TwoPSub, Properties::TTag::PoroElasticSub>;
public:
    using type = PoroMechanicsCouplingManager< Traits >;
};

template<class TypeTag>
struct CouplingManager<TypeTag, TTag::PoroElasticSub>
{
private:
    // define traits etc. as below in main
    using Traits = MultiDomainTraits<Properties::TTag::TwoPSub, Properties::TTag::PoroElasticSub>;
public:
    using type = PoroMechanicsCouplingManager< Traits >;
};

template<class TypeTag>
struct LocalOperator<TypeTag, TTag::TwoPSub>
{
private:
    using ModelTraits = GetPropType<TypeTag, Properties::ModelTraits>;
    using FluxVariables = GetPropType<TypeTag, Properties::FluxVariables>;
    using ElemVariables = typename GetPropType<TypeTag, Properties::GridVariables>::LocalView;
    using ImmiscibleOperators = FVImmiscibleOperators<ModelTraits, FluxVariables, ElemVariables>;

public:
    using type = FVLocalOperator<ElemVariables, ImmiscibleOperators>;
};

template<class TypeTag>
struct LocalOperator<TypeTag, TTag::PoroElasticSub>
{
private:
    using ModelTraits = GetPropType<TypeTag, Properties::ModelTraits>;
    using ElemVariables = typename GetPropType<TypeTag, Properties::GridVariables>::LocalView;
    using StressType = GetPropType<TypeTag, Properties::StressType>;
    using GeomechanicsOperators = FVGeomechanicsOperators<ModelTraits, ElemVariables, StressType>;

public:
    using type = FVLocalOperator<ElemVariables, GeomechanicsOperators>;
};

} // end namespace Properties
} // end namespace Dumux

int main(int argc, char** argv)
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
    using TwoPTypeTag = Properties::TTag::TwoPSub;
    using PoroMechTypeTag = Properties::TTag::PoroElasticSub;

    // we simply extract the grid creator from one of the type tags
    using GridManager = Dumux::GridManager<GetPropType<TwoPTypeTag, Properties::Grid>>;
    GridManager gridManager;
    gridManager.init();

    ////////////////////////////////////////////////////////////
    // run stationary non-linear problem on this grid
    ////////////////////////////////////////////////////////////

    // we compute on the leaf grid view
    const auto& leafGridView = gridManager.grid().leafGridView();

    // create the finite volume grid geometries
    using TwoPFVGridGeometry = GetPropType<TwoPTypeTag, Properties::GridGeometry>;
    using PoroMechFVGridGeometry = GetPropType<PoroMechTypeTag, Properties::GridGeometry>;
    auto twoPFvGridGeometry = std::make_shared<TwoPFVGridGeometry>(leafGridView);
    auto poroMechFvGridGeometry = std::make_shared<PoroMechFVGridGeometry>(leafGridView);
    twoPFvGridGeometry->update();
    poroMechFvGridGeometry->update();

    // the coupling manager
    using Traits = MultiDomainTraits<TwoPTypeTag, PoroMechTypeTag>;
    using CouplingManager = PoroMechanicsCouplingManager<Traits>;
    auto couplingManager = std::make_shared<CouplingManager>();

    // the problems (boundary conditions)
    using TwoPProblem = GetPropType<TwoPTypeTag, Properties::Problem>;
    using PoroMechProblem = GetPropType<PoroMechTypeTag, Properties::Problem>;
    auto twoPSpatialParams = std::make_shared<typename TwoPProblem::SpatialParams>(twoPFvGridGeometry, couplingManager);
    auto twoPProblem = std::make_shared<TwoPProblem>(twoPFvGridGeometry, twoPSpatialParams, "TwoP");
    auto poroMechProblem = std::make_shared<PoroMechProblem>(poroMechFvGridGeometry, couplingManager, "PoroElastic");

    // the solution vectors
    using SolutionVector = typename Traits::SolutionVector;
    SolutionVector x;

    static const auto twoPId = Traits::template SubDomain<0>::Index();
    static const auto poroMechId = Traits::template SubDomain<1>::Index();
    x[twoPId].resize(twoPFvGridGeometry->numDofs());
    x[poroMechId].resize(poroMechFvGridGeometry->numDofs());
    twoPProblem->applyInitialSolution(x[twoPId]);
    poroMechProblem->applyInitialSolution(x[poroMechId]);

    // initialize the coupling manager
    couplingManager->init(twoPProblem, poroMechProblem, x);

    // the grid variables
    using TwoPGridVariables = GetPropType<TwoPTypeTag, Properties::GridVariables>;
    using PoroMechGridVariables = GetPropType<PoroMechTypeTag, Properties::GridVariables>;
    auto twoPGridVariables = std::make_shared<TwoPGridVariables>(twoPProblem, twoPFvGridGeometry);
    auto poroMechGridVariables = std::make_shared<PoroMechGridVariables>(poroMechProblem, poroMechFvGridGeometry);
    twoPGridVariables->init(x[twoPId]);
    poroMechGridVariables->init(x[poroMechId]);

    // get some time loop parameters
    using Scalar = GetPropType<TwoPTypeTag, Properties::Scalar>;
    const auto tEnd = getParam<Scalar>("TimeLoop.TEnd");
    const auto maxDT = getParam<Scalar>("TimeLoop.MaxTimeStepSize");
    auto dt = getParam<Scalar>("TimeLoop.Dt");

    // intialize the vtk output module
    using TwoPVtkOutputModule = Dumux::VtkOutputModule<TwoPGridVariables, GetPropType<TwoPTypeTag, Properties::SolutionVector>>;
    using PoroMechVtkOutputModule = Dumux::VtkOutputModule<PoroMechGridVariables, GetPropType<PoroMechTypeTag, Properties::SolutionVector>>;
    TwoPVtkOutputModule twoPVtkWriter(*twoPGridVariables, x[twoPId], twoPProblem->name());
    PoroMechVtkOutputModule poroMechVtkWriter(*poroMechGridVariables, x[poroMechId], poroMechProblem->name());

    // add output fields to writers
    using TwoPOutputFields = GetPropType<TwoPTypeTag, Properties::IOFields>;
    using PoroMechOutputFields = GetPropType<PoroMechTypeTag, Properties::IOFields>;
    TwoPOutputFields::initOutputModule(twoPVtkWriter);
    PoroMechOutputFields::initOutputModule(poroMechVtkWriter);

    // write initial solution
    // twoPVtkWriter.write(0.0);
    // poroMechVtkWriter.write(0.0);

    //instantiate time loop
    auto timeLoop = std::make_shared<TimeLoop<Scalar>>(0.0, dt, tEnd);
    timeLoop->setMaxTimeStepSize(maxDT);

    // use implicit Euler for time integration
    auto timeMethod = std::make_shared<MultiStage::ImplicitEuler<Scalar>>();

    // the assembler
    using Assembler = MultiDomainAssembler<Traits, CouplingManager, DiffMethod::numeric>;
    auto assembler = std::make_shared<Assembler>( std::make_tuple(twoPFvGridGeometry, poroMechFvGridGeometry),
                                                  couplingManager, *timeMethod);

    auto mdGridVars = MultiDomainFVGridVariables<Traits>(std::make_tuple(twoPGridVariables, poroMechGridVariables));
    mdGridVars.update(x);

    // the linear solver
    using LinearSolver = UMFPackBackend;
    auto linearSolver = std::make_shared<LinearSolver>();

    // the non-linear solver
    using NewtonSolver = Dumux::MultiDomainNewtonSolver<Assembler, LinearSolver, CouplingManager>;
    auto nonLinearSolver = std::make_shared<NewtonSolver>(assembler, linearSolver, couplingManager);

    // the time stepper for time integration
    using TimeStepper = MultiStageTimeStepper<NewtonSolver>;
    TimeStepper timeStepper(nonLinearSolver, timeMethod);

    // time loop
    timeLoop->start(); do
    {
        // do time integraiton
        timeStepper.step(mdGridVars, timeLoop->time(), timeLoop->timeStepSize());

        // advance to the time loop to the next step
        timeLoop->advanceTimeStep();
        x[twoPId] = mdGridVars[twoPId].dofs();
        x[poroMechId] = mdGridVars[poroMechId].dofs();
        // twoPGridVariables->advanceTimeStep();
        // poroMechGridVariables->advanceTimeStep();

        // write vtk output
        twoPVtkWriter.write(timeLoop->time());
        poroMechVtkWriter.write(timeLoop->time());

        // report statistics of this time step
        timeLoop->reportTimeStep();

        // using TwoPPrimaryVariables = GetPropType<TwoPTypeTag, Properties::PrimaryVariables>;
        // TwoPPrimaryVariables storage(0);
        // const auto& twoPLocalResidual = assembler->localResidual(twoPId);
        // for (const auto& element : elements(leafGridView, Dune::Partitions::interior))
        // {
        //     auto storageVec = twoPLocalResidual.evalStorage(*twoPProblem, element, *twoPFvGridGeometry, *twoPGridVariables, x[twoPId]);
        //     storage += storageVec[0];
        // }
        // std::cout << "time, mass CO2 (kg), mass brine (kg):" << std::endl;
        // std::cout << timeLoop->time() << " , " << storage[1] << " , " << storage[0] << std::endl;
        // std::cout << "***************************************" << std::endl;

    } while (!timeLoop->finished());

    // // output some Newton statistics
    // nonLinearSolver.report();

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
