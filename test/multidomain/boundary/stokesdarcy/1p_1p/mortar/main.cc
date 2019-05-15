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
#include <dune/grid/io/file/vtk/vtkwriter.hh>
#include <dune/grid/io/file/vtk/vtksequencewriter.hh>

#include <dumux/common/properties.hh>
#include <dumux/common/parameters.hh>
#include <dumux/common/dumuxmessage.hh>

#include "typetags.hh"
#include "subdomainsolvers.hh"

#include "fluxinterfaceoperator.hh"
#include "pressureinterfaceoperator.hh"
#include "fluxpreconditioner.hh"
#include "pressurepreconditioner.hh"

////////////////////////////////////////////
// Some aliases etc to be used in solve() //
////////////////////////////////////////////
template<class TypeTag>
using DarcySolverType = Dumux::DarcySolver<TypeTag>;
template<class TypeTag>
using StokesSolverType = Dumux::StokesSolver<TypeTag>;

template<class TypeTag1, class Solver1, class TypeTag2, class Solver2>
using FluxInterfaceOperator = Dumux::FluxInterfaceOperator< Solver1, Dumux::GetPropType<TypeTag1, Dumux::Properties::MortarProjector>,
                                                            Solver2, Dumux::GetPropType<TypeTag2, Dumux::Properties::MortarProjector> >;

template<class TypeTag1, class Solver1, class TypeTag2, class Solver2>
using PressureInterfaceOperator = Dumux::PressureInterfaceOperator< Solver1, Dumux::GetPropType<TypeTag1, Dumux::Properties::MortarProjector>,
                                                                    Solver2, Dumux::GetPropType<TypeTag2, Dumux::Properties::MortarProjector> >;

// enum class for mortar variable types
enum class MortarVariableType { pressure, flux };

// translate mortar variable into variable name for output
std::string getMortarVariableName(MortarVariableType mv)
{ return mv == MortarVariableType::pressure ? "p" : "flux"; }

// alias that selects the preconditioner depending on mortar variable
template<MortarVariableType mv>
using Preconditioner = std::conditional_t< mv == MortarVariableType::pressure,
                                           Dumux::MortarPressurePreconditioner<typename Dumux::MortarSpaceTraits::SolutionVector>,
                                           Dumux::MortarFluxPreconditioner<typename Dumux::MortarSpaceTraits::SolutionVector> >;

// alias that selects the right operator depending on mortar variable
template<class TT1, class S1, class TT2, class S2, MortarVariableType mv>
using InterfaceOperator = std::conditional_t< mv == MortarVariableType::pressure,
                                              PressureInterfaceOperator<TT1, S1, TT2, S2>,
                                              FluxInterfaceOperator<TT1, S1, TT2, S2> >;



/////////////////////////////////
// The iterative solve routine //
/////////////////////////////////
template< class TypeTag1, class Solver1,
          class TypeTag2, class Solver2, MortarVariableType mv >
void solveMortar()
{
    using namespace Dumux;

    // create sub-domain solvers
    auto solver1 = std::make_shared< Solver1 >();
    auto solver2 = std::make_shared< Solver2 >();

    solver1->init("Domain1");
    solver2->init("Domain2");

    // make mortar grid, function space basis and solution
    using MortarGrid = typename MortarSpaceTraits::Grid;
    using MortarGridView = typename MortarSpaceTraits::GridView;
    GridManager<MortarGrid> mortarGridManager;
    mortarGridManager.init("Mortar");

    const auto& mortarGridView = mortarGridManager.grid().leafGridView();
    auto feBasis = std::make_shared<typename MortarSpaceTraits::FEBasis>(mortarGridView);

    using MortarSolution = typename MortarSpaceTraits::SolutionVector;
    auto mortarSolution = std::make_shared<MortarSolution>();
    mortarSolution->resize(feBasis->size());
    *mortarSolution = 0.0;

    // create the projectors between mortar and sub-domains
    using Projector1 = Dumux::GetPropType<TypeTag1, Properties::MortarProjector>;
    using Projector2 = Dumux::GetPropType<TypeTag2, Properties::MortarProjector>;

    auto projector1 = std::make_shared<Projector1>(feBasis, solver1->gridGeometryPointer(), solver1->gridVariablesPointer(), "Mortar");
    auto projector2 = std::make_shared<Projector2>(feBasis, solver2->gridGeometryPointer(), solver2->gridVariablesPointer(), "Mortar");

    // let problem and projectors know about each other
    solver1->problemPointer()->setMortarProjector(projector1);
    solver2->problemPointer()->setMortarProjector(projector2);

    projector1->setSubDomainSolutionPointer(solver1->solutionPointer());
    projector2->setSubDomainSolutionPointer(solver2->solutionPointer());

    projector1->setMortarSolutionPointer(mortarSolution);
    projector2->setMortarSolutionPointer(mortarSolution);

    // create vtk writer for mortar grid
    auto mortarWriter = std::make_shared<Dune::VTKWriter<MortarGridView>>(mortarGridView);
    Dune::VTKSequenceWriter<MortarGridView> mortarSequenceWriter(mortarWriter, "mortar");
    mortarSequenceWriter.addCellData(*mortarSolution, getMortarVariableName(mv));

    // write out initial solution
    mortarSequenceWriter.write(0.0);
    solver1->write(0.0);
    solver2->write(0.0);

    // create interface operator
    using Operator = InterfaceOperator<TypeTag1, Solver1, TypeTag2, Solver2, mv>;
    Operator op(solver1, projector1, solver2, projector2);

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

    Preconditioner<mv> prec;
    Dune::CGSolver<MortarSolution> cgSolver(op, prec, reduction, maxIt, verbosity);

    deltaMortarVariable *= -1.0;
    Dune::InverseOperatorResult result;
    cgSolver.apply(*mortarSolution, deltaMortarVariable, result);

    if (!result.converged)
        DUNE_THROW(Dune::InvalidStateException, "CG solver did not converge");

    // solve the sub-domains again to get the right output
    solver1->problemPointer()->setUseHomogeneousSetup(false);
    solver2->problemPointer()->setUseHomogeneousSetup(false);
    op.apply(*mortarSolution, deltaMortarVariable);

    // write solutions
    mortarSequenceWriter.write(1.0);
    solver1->write(1.0);
    solver2->write(1.0);
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

    // check validity of mortar variable specification
    if (mortarVariableType != "Pressure" && mortarVariableType != "Flux")
        DUNE_THROW(Dune::InvalidStateException, "Invalid mortar variable type -" << mortarVariableType << "- provided!");


    ///////////////////////////////////////////////////////////////////////
    // Select the classes depending on input file setup and call solve() //
    ///////////////////////////////////////////////////////////////////////

    static constexpr MortarVariableType pmv = MortarVariableType::pressure;
    static constexpr MortarVariableType fmv = MortarVariableType::flux;

    // darcy-darcy type coupling
    if (solver1Type == "Darcy" && solver2Type == "Darcy")
    {
        if (discScheme1 == "Tpfa" && discScheme2 == "Tpfa")
        {
            if (mortarVariableType == "Pressure")
                solveMortar< Properties::TTag::DarcyOnePTpfaPressure, DarcySolverType<Properties::TTag::DarcyOnePTpfaPressure>,
                             Properties::TTag::DarcyOnePTpfaPressure, DarcySolverType<Properties::TTag::DarcyOnePTpfaPressure>, pmv >();
            else
                solveMortar< Properties::TTag::DarcyOnePTpfaFlux, DarcySolverType<Properties::TTag::DarcyOnePTpfaFlux>,
                             Properties::TTag::DarcyOnePTpfaFlux, DarcySolverType<Properties::TTag::DarcyOnePTpfaFlux>, fmv >();
        }

        else if (discScheme1 == "Box" && discScheme2 == "Box")
        {
            if (mortarVariableType == "Pressure")
                solveMortar< Properties::TTag::DarcyOnePBoxPressure, DarcySolverType<Properties::TTag::DarcyOnePBoxPressure>,
                             Properties::TTag::DarcyOnePBoxPressure, DarcySolverType<Properties::TTag::DarcyOnePBoxPressure>, pmv >();
            else
                solveMortar< Properties::TTag::DarcyOnePBoxFlux, DarcySolverType<Properties::TTag::DarcyOnePBoxFlux>,
                             Properties::TTag::DarcyOnePBoxFlux, DarcySolverType<Properties::TTag::DarcyOnePBoxFlux>, fmv >();
        }

        else if (discScheme1 == "Box" && discScheme2 == "Tpfa")
        {
            if (mortarVariableType == "Pressure")
                solveMortar< Properties::TTag::DarcyOnePBoxPressure, DarcySolverType<Properties::TTag::DarcyOnePBoxPressure>,
                             Properties::TTag::DarcyOnePTpfaPressure, DarcySolverType<Properties::TTag::DarcyOnePTpfaPressure>, pmv >();
            else
                solveMortar< Properties::TTag::DarcyOnePBoxFlux, DarcySolverType<Properties::TTag::DarcyOnePBoxFlux>,
                             Properties::TTag::DarcyOnePTpfaFlux, DarcySolverType<Properties::TTag::DarcyOnePTpfaFlux>, fmv >();
        }

        else if (discScheme1 == "Tpfa" && discScheme2 == "Box")
        {
            if (mortarVariableType == "Pressure")
                solveMortar< Properties::TTag::DarcyOnePTpfaPressure, DarcySolverType<Properties::TTag::DarcyOnePTpfaPressure>,
                             Properties::TTag::DarcyOnePBoxPressure, DarcySolverType<Properties::TTag::DarcyOnePBoxPressure>, pmv >();
            else
                solveMortar< Properties::TTag::DarcyOnePTpfaFlux, DarcySolverType<Properties::TTag::DarcyOnePTpfaFlux>,
                             Properties::TTag::DarcyOnePBoxFlux, DarcySolverType<Properties::TTag::DarcyOnePBoxFlux>, fmv >();
        }
    }
    else if (solver1Type == "Darcy" && solver2Type == "Stokes")
    {
        if (discScheme1 == "Tpfa")
        {
            if (mortarVariableType == "Pressure")
                solveMortar< Properties::TTag::DarcyOnePTpfaPressure, DarcySolverType<Properties::TTag::DarcyOnePTpfaPressure>,
                             Properties::TTag::StokesOnePPressure, StokesSolverType<Properties::TTag::StokesOnePPressure>, pmv >();
            else
                solveMortar< Properties::TTag::DarcyOnePTpfaFlux, DarcySolverType<Properties::TTag::DarcyOnePTpfaFlux>,
                             Properties::TTag::StokesOnePFlux, StokesSolverType<Properties::TTag::StokesOnePFlux>, fmv >();
        }

        else if (discScheme1 == "Box")
        {
            if (mortarVariableType == "Pressure")
                solveMortar< Properties::TTag::DarcyOnePBoxPressure, DarcySolverType<Properties::TTag::DarcyOnePBoxPressure>,
                             Properties::TTag::StokesOnePPressure, StokesSolverType<Properties::TTag::StokesOnePPressure>, pmv >();
            else
                solveMortar< Properties::TTag::DarcyOnePBoxFlux, DarcySolverType<Properties::TTag::DarcyOnePBoxFlux>,
                             Properties::TTag::StokesOnePFlux, StokesSolverType<Properties::TTag::StokesOnePFlux>, fmv >();
        }

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
