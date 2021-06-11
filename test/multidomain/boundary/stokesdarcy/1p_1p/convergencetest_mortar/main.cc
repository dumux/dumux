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
#include <type_traits>

#include <dune/common/timer.hh>
#include <dune/common/exceptions.hh>
#include <dune/common/parallel/mpihelper.hh>

#include <dune/foamgrid/foamgrid.hh>
#include <dune/functions/functionspacebases/lagrangebasis.hh>
#include <dune/functions/functionspacebases/lagrangedgbasis.hh>
#include <dune/functions/gridfunctions/discreteglobalbasisfunction.hh>
#include <dune/functions/gridfunctions/analyticgridviewfunction.hh>

#include <dune/grid/io/file/vtk/vtkwriter.hh>
#include <dune/grid/io/file/vtk/vtksequencewriter.hh>

#include <dumux/common/properties.hh>
#include <dumux/common/parameters.hh>
#include <dumux/common/dumuxmessage.hh>
#include <dumux/common/integrate.hh>
#include <dumux/discretization/fem/fegridgeometry.hh>

#include "normalfluxbasis.hh"
#include "problem_darcy.hh"
#include "problem_stokes.hh"

#include "subdomainsolvers.hh"
#include "mortarvariabletype.hh"

#include "projector.hh"
#include "projectorcreator.hh"
#include "reconstructor.hh"

#include "preconditioner.hh"
#include "interfaceoperator.hh"

template<class Problem, class SolutionVector>
void printFreeFlowL2Error(const Problem& problem, const SolutionVector& x)
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
            const Scalar delta = x[dofIdx] - problem.exactPressure(scv.center());
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

////////////////////////////////////////////
// Some aliases etc to be used in solve() //
////////////////////////////////////////////

// Traits class for a sub-domain
template<class Solver>
struct SubDomainTraits
{
    using SolutionVector = typename Solver::SolutionVector;
    using GridGeometry = typename Solver::FVGridGeometry;
    using GridVariables = typename Solver::GridVariables;
    using FluxVariables = typename Solver::FluxVariables;
};

// Define grid and basis for mortar domain
struct MortarTraits
{
    static constexpr int basisOrder = 0;

    using Scalar = double;
    using Grid = Dune::FoamGrid<1, 2>;
    using GridView = typename Grid::LeafGridView;
    using BlockType = Dune::FieldVector<Scalar, 1>;
    using SolutionVector = Dune::BlockVector<BlockType>;
    using FEBasis = Dune::Functions::LagrangeBasis<GridView, basisOrder>;
    using GridGeometry = Dumux::FEGridGeometry<FEBasis>;
};

template<class SubDomainTypeTag>
using DarcySolverType = Dumux::DarcySolver<SubDomainTypeTag>;

template<class SubDomainTypeTag>
using StokesSolverType = Dumux::StokesSolver<SubDomainTypeTag>;

// translate mortar variable into variable name for output
std::string getMortarVariableName(Dumux::OnePMortarVariableType mv)
{ return mv == Dumux::OnePMortarVariableType::pressure ? "p" : "flux"; }



/////////////////////////////////
// The iterative solve routine //
/////////////////////////////////
template< class Solver1,
          class Solver2,
          class ProjectorCreator = Dumux::DefaultMortarProjectorCreator >
void solveMortar(Dumux::OnePMortarVariableType mv)
{
    using namespace Dumux;

    Dune::Timer watch;

    // create sub-domain solvers
    auto solver1 = std::make_shared< Solver1 >();
    auto solver2 = std::make_shared< Solver2 >();

    solver1->init("Domain1");
    solver2->init("Domain2");

    // make mortar grid, function space basis and solution
    using MortarGrid = typename MortarTraits::Grid;
    using MortarGridView = typename MortarTraits::GridView;
    GridManager<MortarGrid> mortarGridManager;
    mortarGridManager.init("Mortar");

    const auto& mortarGridView = mortarGridManager.grid().leafGridView();
    auto feBasis = std::make_shared<typename MortarTraits::FEBasis>(mortarGridView);
    auto mortarGridGeometry = std::make_shared<typename MortarTraits::GridGeometry>(feBasis);

    using MortarSolution = typename MortarTraits::SolutionVector;
    auto mortarSolution = std::make_shared<MortarSolution>();
    mortarSolution->resize(feBasis->size());
    *mortarSolution = 0.0;

    // create the projectors between mortar and sub-domains
    const auto projectors = ProjectorCreator::template makeProjectors<MortarSolution>(*solver1, *solver2, *mortarGridGeometry, mv);
    auto projector1 = projectors.first;
    auto projector2 = projectors.second;

    // create vtk writer for mortar grid
    auto mortarWriter = std::make_shared<Dune::VTKWriter<MortarGridView>>(mortarGridView, Dune::VTK::nonconforming);
    Dune::VTKSequenceWriter<MortarGridView> mortarSequenceWriter(mortarWriter, "mortar");

    auto mortarGridFunction = Dune::Functions::template makeDiscreteGlobalBasisFunction<typename MortarSolution::block_type>(*feBasis, *mortarSolution);
    const auto fieldInfoMortar = Dune::VTK::FieldInfo({getMortarVariableName(mv), Dune::VTK::FieldInfo::Type::scalar, 1});

    if (MortarTraits::basisOrder == 0)
        mortarWriter->addCellData(mortarGridFunction, fieldInfoMortar);
    else
        mortarWriter->addVertexData(mortarGridFunction, fieldInfoMortar);

    // project initial mortar solution into sub-domains
    solver1->problemPointer()->setMortarProjection( projector1->projectMortarToSubDomain(*mortarSolution) );
    solver2->problemPointer()->setMortarProjection( projector2->projectMortarToSubDomain(*mortarSolution) );

    // write out initial solution
    mortarSequenceWriter.write(0.0);
    solver1->write(0.0);
    solver2->write(0.0);

    // create interface operator
    using Reconstructor1 = MortarReconstructor< SubDomainTraits<Solver1> >;
    using Reconstructor2 = MortarReconstructor< SubDomainTraits<Solver2> >;
    using Operator = OnePMortarInterfaceOperator<Solver1, Reconstructor1,
                                                 Solver2, Reconstructor2, MortarSolution>;
    Operator op(solver1, projector1, solver2, projector2, *mortarGridGeometry, mv);

    // first compute the jump in mortar variable
    solver1->problemPointer()->setUseHomogeneousSetup(false);
    solver2->problemPointer()->setUseHomogeneousSetup(false);

    MortarSolution deltaMortarVariable;
    op.apply(*mortarSolution, deltaMortarVariable);

    // Solve the homogeneous problem with CG solver
    const double reduction = getParam<double>("InterfaceSolver.ResidualReduction");
    const std::size_t maxIt = getParam<double>("InterfaceSolver.MaxIterations");
    const int verbosity = getParam<int>("InterfaceSolver.Verbosity");

    solver1->problemPointer()->setUseHomogeneousSetup(true);
    solver2->problemPointer()->setUseHomogeneousSetup(true);

    // create preconditioner
    using Prec = Dumux::OnePMortarPreconditioner<Solver1, Reconstructor1,
                                                 Solver2, Reconstructor2, MortarSolution>;
    Prec prec(solver1, projector1, solver2, projector2, *mortarGridGeometry, mv);

    // apply linear solver using our linear operator
    deltaMortarVariable *= -1.0;
    Dune::InverseOperatorResult result;

    const auto lsType = getParam<std::string>("InterfaceSolver.LinearSolverType");
    if (lsType == "CG")
    {
        Dune::CGSolver<MortarSolution> cgSolver(op, prec, reduction, maxIt, verbosity);
        cgSolver.apply(*mortarSolution, deltaMortarVariable, result);
    }
    else if (lsType == "GMRes")
    {
        Dune::RestartedGMResSolver<MortarSolution> gmresSolver(op, prec, reduction, maxIt, maxIt, verbosity);
        gmresSolver.apply(*mortarSolution, deltaMortarVariable, result);
    }

    if (!result.converged)
        DUNE_THROW(Dune::InvalidStateException, "CG solver did not converge with given maximum number of iterations");

    // solve the sub-domains again to get the right output
    solver1->problemPointer()->setUseHomogeneousSetup(false);
    solver2->problemPointer()->setUseHomogeneousSetup(false);
    op.apply(*mortarSolution, deltaMortarVariable);

    // write solutions
    mortarSequenceWriter.write(1.0);
    solver1->write(1.0);
    solver2->write(1.0);

    // compute L2 error
    printFreeFlowL2Error(*solver2->problemPointer(), *solver2->solutionPointer());
    printDarcyL2Error((*solver1->problemPointer()), *solver1->solutionPointer());

    // const auto l2Error1 = solver1->problemPointer()->calculateL2Error(*solver1->solutionPointer());
    // const auto l2Error2 = solver2->problemPointer()->calculateL2Error(*solver2->solutionPointer());

    // // write into file
    // std::ofstream errorFile(getParam<std::string>("L2Error.OutputFile"), std::ios::app);
    // errorFile << l2Error1 << ","
    //           << l2Error2[0] << ","
    //           << l2Error2[1] << ","
    //           << l2Error2[2] << std::endl;
    // errorFile.close();

    // print time necessary for solve
    std::cout << "\n#####################################################\n\n"
              << "Iterative scheme took " << watch.elapsed() << " seconds\n"
              << "\n#####################################################\n\n";
}

///////////////////////////////////////////////////////////////////
// Main Program. Selects the solvers etc to be passed to solve() //
///////////////////////////////////////////////////////////////////
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

    // get solver types of the two subdomains
    const auto solver1Type = getParam<std::string>("Domain1.SolverType");
    const auto solver2Type = getParam<std::string>("Domain2.SolverType");

    // discretization scheme used in the sub-domains
    const auto discScheme1 = getParam<std::string>("Domain1.DiscretizationScheme");
    const auto discScheme2 = getParam<std::string>("Domain2.DiscretizationScheme");

    // determine what the mortar variable is
    const auto mortarVariableType = getParam<std::string>("Mortar.VariableType");

    //////////////////////////////////////////
    // Check validity of the specifications //
    //////////////////////////////////////////
    for (unsigned int i = 0; i < 2; ++i)
    {
        const auto& s = i == 0 ? solver1Type : solver2Type;
        const auto& d = i == 0 ? discScheme1 : discScheme2;

        if (s != "Darcy" && s != "Stokes")
            DUNE_THROW(Dune::InvalidStateException, "Invalid solver type -" << s << "- provided!");

        if (s == "Darcy" && d != "Tpfa" && d != "Mpfa" && d != "Box")
            DUNE_THROW(Dune::InvalidStateException, "Invalid Darcy discretization scheme -" << d << "- provided!");
        else if (s == "Stokes" && d != "Staggered")
            DUNE_THROW(Dune::InvalidStateException, "Invalid Stokes discretization scheme -" << d << "- provided!");
    }

    // // check validity of mortar variable specification
    // if (mortarVariableType != "Pressure" && mortarVariableType != "Flux")
    //     DUNE_THROW(Dune::InvalidStateException, "Invalid mortar variable type -" << mortarVariableType << "- provided!");


    ///////////////////////////////////////////////////////////////////////
    // Select the classes depending on input file setup and call solve() //
    ///////////////////////////////////////////////////////////////////////
    OnePMortarVariableType mvType = mortarVariableType == "Pressure" ? OnePMortarVariableType::pressure
                                                                     : OnePMortarVariableType::flux;

    using TTDarcyTpfa = Properties::TTag::DarcyOnePTpfa;
    using TTDarcyBox = Properties::TTag::DarcyOnePBox;
    using TTStokesStaggered = Properties::TTag::StokesOneP;

    // stokes-darcy coupling
    if (solver1Type == "Darcy" && solver2Type == "Stokes")
    {
        if (discScheme1 == "Tpfa")
            solveMortar< DarcySolverType<TTDarcyTpfa>, StokesSolverType<TTStokesStaggered> >(mvType);
        else if (discScheme1 == "Box")
            solveMortar< DarcySolverType<TTDarcyBox>, StokesSolverType<TTStokesStaggered> >(mvType);
    }

    else
        DUNE_THROW(Dune::InvalidStateException, "Solver combination not implemented!");


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
