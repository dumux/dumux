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
#include <iostream>

#include <dune/common/parallel/mpihelper.hh>

#include "1pproblem.hh"
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

#include <dumux/multidomain/poromechanics/poromechanicscouplingmanager.hh>

#include <dumux/io/vtkfunction.hh>
#include <dune/grid/io/file/vtk/vtkwriter.hh>

// set the coupling manager property in the sub-problems
namespace Dumux {
namespace Properties {
NEW_PROP_TAG(CouplingManager);

SET_PROP(OnePSubTypeTag, CouplingManager)
{
private:
    // define traits etc. as below in main
    using Traits = MultiDomainTraits<TTAG(OnePSubTypeTag), TTAG(PoroElasticSubTypeTag)>;
public:
    using type = PoroMechanicsCouplingManager< Traits >;
};

SET_PROP(PoroElasticSubTypeTag, CouplingManager)
{
private:
    // define traits etc. as below in main
    using Traits = MultiDomainTraits<TTAG(OnePSubTypeTag), TTAG(PoroElasticSubTypeTag)>;
public:
    using type = PoroMechanicsCouplingManager< Traits >;
};
}
}

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
    using OnePTypeTag = TTAG(OnePSubTypeTag);
    using PoroMechTypeTag = TTAG(PoroElasticSubTypeTag);

    // we simply extract the grid creator from one of the type tags
    using GridCreator = typename GET_PROP_TYPE(OnePTypeTag, GridCreator);
    GridCreator::makeGrid();
    GridCreator::loadBalance();

    ////////////////////////////////////////////////////////////
    // run stationary non-linear problem on this grid
    ////////////////////////////////////////////////////////////

    // we compute on the leaf grid view
    const auto& leafGridView = GridCreator::grid().leafGridView();

    // create the finite volume grid geometries
    using OnePFVGridGeometry = typename GET_PROP_TYPE(OnePTypeTag, FVGridGeometry);
    using PoroMechFVGridGeometry = typename GET_PROP_TYPE(PoroMechTypeTag, FVGridGeometry);
    auto onePFvGridGeometry = std::make_shared<OnePFVGridGeometry>(leafGridView);
    auto poroMechFvGridGeometry = std::make_shared<PoroMechFVGridGeometry>(leafGridView);
    onePFvGridGeometry->update();
    poroMechFvGridGeometry->update();

    // the problems (boundary conditions)
    using OnePProblem = typename GET_PROP_TYPE(OnePTypeTag, Problem);
    using PoroMechProblem = typename GET_PROP_TYPE(PoroMechTypeTag, Problem);
    auto onePProblem = std::make_shared<OnePProblem>(onePFvGridGeometry, "OneP");
    auto poroMechProblem = std::make_shared<PoroMechProblem>(poroMechFvGridGeometry, "PoroElastic");

    // the solution vectors
    using Traits = MultiDomainTraits<OnePTypeTag, PoroMechTypeTag>;
    using SolutionVector = typename Traits::SolutionVector;
    SolutionVector x;

    static const auto onePId = Traits::template DomainIdx<0>();
    static const auto poroMechId = Traits::template DomainIdx<1>();
    x[onePId].resize(onePFvGridGeometry->numDofs());
    x[poroMechId].resize(poroMechFvGridGeometry->numDofs());
    onePProblem->applyInitialSolution(x[onePId]);
    poroMechProblem->applyInitialSolution(x[poroMechId]);
    SolutionVector xOld = x;

    // the coupling manager
    using CouplingManager = PoroMechanicsCouplingManager<Traits>;
    auto couplingManager = std::make_shared<CouplingManager>();
    couplingManager->init(onePProblem, poroMechProblem, x);

    // set coupling manager pointer in sub-problems
    // (in 1p problem the coupling enters in the parameters)
    onePProblem->spatialParams().setCouplingManager(couplingManager);
    poroMechProblem->setCouplingManager(couplingManager);

    // the grid variables
    using OnePGridVariables = typename GET_PROP_TYPE(OnePTypeTag, GridVariables);
    using PoroMechGridVariables = typename GET_PROP_TYPE(PoroMechTypeTag, GridVariables);
    auto onePGridVariables = std::make_shared<OnePGridVariables>(onePProblem, onePFvGridGeometry);
    auto poroMechGridVariables = std::make_shared<PoroMechGridVariables>(poroMechProblem, poroMechFvGridGeometry);
    onePGridVariables->init(x[onePId]);
    poroMechGridVariables->init(x[poroMechId]);

    // intialize the vtk output modules
    Dune::VTKWriter<typename OnePFVGridGeometry::GridView> onePVtkWriter(onePFvGridGeometry->gridView());
    Dune::VTKWriter<typename PoroMechFVGridGeometry::GridView> poroMechVtkWriter(poroMechFvGridGeometry->gridView());

    using VTKFunction = Dumux::Vtk::VectorP1VTKFunction< typename PoroMechFVGridGeometry::GridView, decltype(x[poroMechId]) >;
    auto displacementFunction = std::make_shared<VTKFunction>( poroMechFvGridGeometry->gridView(), poroMechFvGridGeometry->vertexMapper(),
                                                               x[poroMechId], "u", PoroMechFVGridGeometry::GridView::dimensionworld );

    // add displacement/pressure to poro-elastic vtk writer and write initial solution
    onePVtkWriter.addCellData(x[onePId], "p");
    poroMechVtkWriter.addVertexData(displacementFunction);

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
    onePVtkWriter.write("onep");
    poroMechVtkWriter.write("poromech");

    ////////////////////////////////////////////////////////////
    // finalize, print dumux message to say goodbye
    ////////////////////////////////////////////////////////////
    if (mpiHelper.rank() == 0)
        DumuxMessage::print(/*firstCall=*/false);

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
