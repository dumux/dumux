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
#include <iostream>

#include <dune/common/parallel/mpihelper.hh>

#include "problem_1p.hh"
#include "problem_poroelastic.hh"

#include <dumux/common/properties.hh>
#include <dumux/common/parameters.hh>
#include <dumux/common/dumuxmessage.hh>

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

template<class TypeTag>
struct CouplingManager<TypeTag, TTag::OnePSub>
{
private:
    // define traits etc. as below in main
    using Traits = MultiDomainTraits<TTag::OnePSub, TTag::PoroElasticSub>;
public:
    using type = PoroMechanicsCouplingManager< Traits >;
};

template<class TypeTag>
struct CouplingManager<TypeTag, TTag::PoroElasticSub>
{
private:
    // define traits etc. as below in main
    using Traits = MultiDomainTraits<TTag::OnePSub, TTag::PoroElasticSub>;
public:
    using type = PoroMechanicsCouplingManager< Traits >;
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
    using OnePTypeTag = Properties::TTag::OnePSub;
    using PoroMechTypeTag = Properties::TTag::PoroElasticSub;

    // we simply extract the grid creator from one of the type tags
    using GridManager = Dumux::GridManager<GetPropType<OnePTypeTag, Properties::Grid>>;
    GridManager gridManager;
    gridManager.init();

    ////////////////////////////////////////////////////////////
    // run stationary non-linear problem on this grid
    ////////////////////////////////////////////////////////////

    // we compute on the leaf grid view
    const auto& leafGridView = gridManager.grid().leafGridView();

    // create the finite volume grid geometries
    using OnePFVGridGeometry = GetPropType<OnePTypeTag, Properties::GridGeometry>;
    using PoroMechFVGridGeometry = GetPropType<PoroMechTypeTag, Properties::GridGeometry>;
    auto onePFvGridGeometry = std::make_shared<OnePFVGridGeometry>(leafGridView);
    auto poroMechFvGridGeometry = std::make_shared<PoroMechFVGridGeometry>(leafGridView);
    onePFvGridGeometry->update();
    poroMechFvGridGeometry->update();

    // the coupling manager
    using Traits = MultiDomainTraits<OnePTypeTag, PoroMechTypeTag>;
    using CouplingManager = PoroMechanicsCouplingManager<Traits>;
    auto couplingManager = std::make_shared<CouplingManager>();

    // the problems (boundary conditions)
    using OnePProblem = GetPropType<OnePTypeTag, Properties::Problem>;
    using PoroMechProblem = GetPropType<PoroMechTypeTag, Properties::Problem>;
    auto onePSpatialParams = std::make_shared<typename OnePProblem::SpatialParams>(onePFvGridGeometry, couplingManager);
    auto onePProblem = std::make_shared<OnePProblem>(onePFvGridGeometry, onePSpatialParams, "OneP");
    auto poroMechProblem = std::make_shared<PoroMechProblem>(poroMechFvGridGeometry, couplingManager, "PoroElastic");

    // the solution vectors
    using SolutionVector = typename Traits::SolutionVector;
    SolutionVector x;

    static const auto onePId = Traits::template SubDomain<0>::Index();
    static const auto poroMechId = Traits::template SubDomain<1>::Index();
    x[onePId].resize(onePFvGridGeometry->numDofs());
    x[poroMechId].resize(poroMechFvGridGeometry->numDofs());
    onePProblem->applyInitialSolution(x[onePId]);
    poroMechProblem->applyInitialSolution(x[poroMechId]);
    SolutionVector xOld = x;

    // initialize the coupling manager
    couplingManager->init(onePProblem, poroMechProblem, x);

    // the grid variables
    using OnePGridVariables = GetPropType<OnePTypeTag, Properties::GridVariables>;
    using PoroMechGridVariables = GetPropType<PoroMechTypeTag, Properties::GridVariables>;
    auto onePGridVariables = std::make_shared<OnePGridVariables>(onePProblem, onePFvGridGeometry);
    auto poroMechGridVariables = std::make_shared<PoroMechGridVariables>(poroMechProblem, poroMechFvGridGeometry);
    onePGridVariables->init(x[onePId]);
    poroMechGridVariables->init(x[poroMechId]);

    // intialize the vtk output module
    using OnePVtkOutputModule = Dumux::VtkOutputModule<OnePGridVariables, GetPropType<OnePTypeTag, Properties::SolutionVector>>;
    using PoroMechVtkOutputModule = Dumux::VtkOutputModule<PoroMechGridVariables, GetPropType<PoroMechTypeTag, Properties::SolutionVector>>;
    OnePVtkOutputModule onePVtkWriter(*onePGridVariables, x[onePId], onePProblem->name());
    PoroMechVtkOutputModule poroMechVtkWriter(*poroMechGridVariables, x[poroMechId], poroMechProblem->name());

    // add output fields to writers
    using OnePOutputFields = GetPropType<OnePTypeTag, Properties::IOFields>;
    using PoroMechOutputFields = GetPropType<PoroMechTypeTag, Properties::IOFields>;
    OnePOutputFields::initOutputModule(onePVtkWriter);
    PoroMechOutputFields::initOutputModule(poroMechVtkWriter);

    // write initial solution
    onePVtkWriter.write(0.0);
    poroMechVtkWriter.write(0.0);

    // the assembler
    using Assembler = MultiDomainFVAssembler<Traits, CouplingManager, DiffMethod::numeric, /*implicit?*/true>;
    auto assembler = std::make_shared<Assembler>( std::make_tuple(onePProblem, poroMechProblem),
                                                  std::make_tuple(onePFvGridGeometry, poroMechFvGridGeometry),
                                                  std::make_tuple(onePGridVariables, poroMechGridVariables),
                                                  couplingManager);

    // the linear solver
    using LinearSolver = ILU0BiCGSTABBackend;
    auto linearSolver = std::make_shared<LinearSolver>();

    // the non-linear solver
    using NewtonSolver = Dumux::MultiDomainNewtonSolver<Assembler, LinearSolver, CouplingManager>;
    auto newtonSolver = std::make_shared<NewtonSolver>(assembler, linearSolver, couplingManager);

    // linearize & solve
    newtonSolver->solve(x);

    // update grid variables for output
    onePGridVariables->update(x[onePId]);
    poroMechGridVariables->update(x[poroMechId]);

    // write vtk output
    onePVtkWriter.write(1.0);
    poroMechVtkWriter.write(1.0);

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
