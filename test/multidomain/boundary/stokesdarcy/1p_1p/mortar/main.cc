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
 * \brief TODO doc me
 */
#include <config.h>
#include <iostream>

#include <dune/common/timer.hh>
#include <dune/common/parallel/mpihelper.hh>

#include <dune/grid/io/file/vtk/vtkwriter.hh>
#include <dune/foamgrid/foamgrid.hh>

#include <dune/functions/functionspacebases/lagrangebasis.hh>

#include <dumux/common/properties.hh>
#include <dumux/common/parameters.hh>
#include <dumux/common/partial.hh>
#include <dumux/common/dumuxmessage.hh>

#include <dumux/linear/seqsolverbackend.hh>

#include "interfaceoperator.hh"
#include "preconditioner.hh"

#include "subdomainsolvers.hh"
#include "problem_darcy.hh"
#include "problem_stokes.hh"

#include "fluxprojector.hh"

// The type tags of the sub-domains
using DarcyTypeTag = Dumux::Properties::TTag::DarcyOneP;
using StokesTypeTag = Dumux::Properties::TTag::StokesOneP;

// The solution vectors of the sub-domains
using DarcySolutionVector = Dumux::GetPropType<DarcyTypeTag, Dumux::Properties::SolutionVector>;
using StokesSolutionVector = Dumux::GetPropType<StokesTypeTag, Dumux::Properties::SolutionVector>;

// The grid geometries of the sub-domains
using DarcyGridGeometry = Dumux::GetPropType<DarcyTypeTag, Dumux::Properties::FVGridGeometry>;
using StokesGridGeometry = Dumux::GetPropType<StokesTypeTag, Dumux::Properties::FVGridGeometry>;

// Define grid an basis for mortar domain
using MortarScalar = double;
using MortarGrid = Dune::FoamGrid<1, 2>;
using MortarGridView = typename MortarGrid::LeafGridView;
using MortarSolutionVector = Dune::BlockVector<Dune::FieldVector<MortarScalar, 1>>;
using MortarFEBasis = Dune::Functions::LagrangeBasis<MortarGridView, 0>;

// Projection operators
using TheMortarDarcyProjector = Dumux::MortarFluxProjector<MortarFEBasis, MortarSolutionVector,
                                                           DarcyGridGeometry, DarcySolutionVector>;

// Set Projection property in Sub-problems
namespace Dumux {
namespace Properties {

template<class TypeTag>
struct MortarProjector<TypeTag, DarcyTypeTag> { using type = TheMortarDarcyProjector; };

// TODO: STOKES PROJECTOR

} // end namespace Properties
} // end namespace Dumux

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

    // create sub-domain solvers (two times Darcy for now)
    using DarcySolver1 = DarcySolver<DarcyTypeTag>;
    using DarcySolver2 = DarcySolver<DarcyTypeTag>;
    auto darcySolver = std::make_shared< DarcySolver1 >();
    auto darcy2Solver = std::make_shared< DarcySolver2 >();

    darcySolver->init("Darcy");
    darcy2Solver->init("Darcy2");

    // make mortar grid, function space basis and solution
    GridManager<MortarGrid> mortarGridManager;
    mortarGridManager.init("Mortar");

    const auto& mortarGridView = mortarGridManager.grid().leafGridView();
    auto feBasis = std::make_shared<MortarFEBasis>(mortarGridView);

    auto mortarSolution = std::make_shared<MortarSolutionVector>();
    mortarSolution->resize(feBasis->size());
    *mortarSolution = 0.0;

    // create the projectors between mortar and sub-domains
    auto darcyProjector = std::make_shared<TheMortarDarcyProjector>(feBasis, darcySolver->gridGeometryPointer(), "Mortar");
    auto darcyProjector2 = std::make_shared<TheMortarDarcyProjector>(feBasis, darcy2Solver->gridGeometryPointer(), "Mortar");

    // let problem and projectors know about each other
    darcySolver->problemPointer()->setMortarProjector(darcyProjector);
    darcy2Solver->problemPointer()->setMortarProjector(darcyProjector2);

    darcyProjector->setSubDomainSolutionPointer(darcySolver->solutionPointer());
    darcyProjector2->setSubDomainSolutionPointer(darcy2Solver->solutionPointer());

    darcyProjector->setMortarSolutionPointer(mortarSolution);
    darcyProjector2->setMortarSolutionPointer(mortarSolution);

    // write out initial solution
    Dune::VTKWriter<MortarGridView> mortarWriter(mortarGridView);
    mortarWriter.addCellData(*mortarSolution, "flux");
    mortarWriter.write("mortar");
    darcySolver->write(1.0);
    darcy2Solver->write(1.0);

    // compute pressure jump
    using Operator = InterfaceOperator<DarcySolver1, TheMortarDarcyProjector,
                                       DarcySolver2, TheMortarDarcyProjector>;
    Operator op(darcySolver, darcyProjector, darcy2Solver, darcyProjector2);

    darcySolver->problemPointer()->setUseHomogeneousSetup(false);
    darcy2Solver->problemPointer()->setUseHomogeneousSetup(false);

    MortarSolutionVector deltaP;
    op.apply(*mortarSolution, deltaP);

    //
    const double reduction = getParam<double>("InterfaceSolver.ResidualReduction");
    const std::size_t maxIt = getParam<double>("InterfaceSolver.MaxIterations");
    const bool verbose = getParam<double>("InterfaceSolver.Verbosity");

    darcySolver->problemPointer()->setUseHomogeneousSetup(true);
    darcy2Solver->problemPointer()->setUseHomogeneousSetup(true);

    MortarStokesDarcyPreconditioner<MortarSolutionVector> prec;
    Dune::CGSolver<MortarSolutionVector> cgSolver(op, prec, reduction, maxIt, verbose);

    deltaP *= -1.0;
    Dune::InverseOperatorResult result;
    cgSolver.apply(*mortarSolution, deltaP, result);

    if (!result.converged)
        DUNE_THROW(Dune::InvalidStateException, "Linear solver did not converge");

    // solve the sub-domains
    darcySolver->problemPointer()->setUseHomogeneousSetup(false);
    darcy2Solver->problemPointer()->setUseHomogeneousSetup(false);
    op.apply(*mortarSolution, deltaP);
    std::cout << "DeltaP:\n" << deltaP << std::endl;

    // write solutions
    darcySolver->write(1.0);
    darcy2Solver->write(1.0);
    mortarWriter.write("mortar");

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
