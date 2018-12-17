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
 * \brief test for the one-phase facet coupling model
 */
#include <config.h>

#include <iostream>

#include <dune/common/parallel/mpihelper.hh>

#include "problem_bulk.hh"
#include "problem_lowdim.hh"

#include <dumux/common/properties.hh>
#include <dumux/common/parameters.hh>
#include <dumux/common/dumuxmessage.hh>
#include <dumux/common/defaultusagemessage.hh>

#include <dumux/assembly/diffmethod.hh>
#include <dumux/linear/seqsolverbackend.hh>

#include <dumux/multidomain/newtonsolver.hh>
#include <dumux/multidomain/fvassembler.hh>
#include <dumux/multidomain/traits.hh>

#include <dumux/multidomain/facet/gridmanager.hh>
#include <dumux/multidomain/facet/couplingmapper.hh>
#include <dumux/multidomain/facet/couplingmanager.hh>

#include <dumux/io/vtkoutputmodule.hh>

// obtain/define some types to be used below in the property definitions and in main
template< class BulkTypeTag, class LowDimTypeTag >
class TestTraits
{
    using BulkFVGridGeometry = typename GET_PROP_TYPE(BulkTypeTag, FVGridGeometry);
    using LowDimFVGridGeometry = typename GET_PROP_TYPE(LowDimTypeTag, FVGridGeometry);
public:
    using MDTraits = Dumux::MultiDomainTraits<BulkTypeTag, LowDimTypeTag>;
    using CouplingMapper = Dumux::FacetCouplingMapper<BulkFVGridGeometry, LowDimFVGridGeometry>;
    using CouplingManager = Dumux::FacetCouplingManager<MDTraits, CouplingMapper>;
};

namespace Dumux {
namespace Properties {

// set cm property in the sub-problems
using TpfaTraits = TestTraits<TTAG(OnePBulkTpfa), TTAG(OnePLowDimTpfa)>;
SET_TYPE_PROP(OnePBulkTpfa, CouplingManager, typename TpfaTraits::CouplingManager);
SET_TYPE_PROP(OnePLowDimTpfa, CouplingManager, typename TpfaTraits::CouplingManager);

} // end namespace Properties
} // end namespace Dumux

// main program
int main(int argc, char** argv) try
{
    using namespace Dumux;

    //////////////////////////////////////////////////////
    //////////////////////////////////////////////////////

    // initialize MPI, finalize is done automatically on exit
    const auto& mpiHelper = Dune::MPIHelper::instance(argc, argv);

    // print dumux start message
    if (mpiHelper.rank() == 0)
        DumuxMessage::print(/*firstCall=*/true);

    // initialize parameter tree
    Parameters::init(argc, argv);

    //////////////////////////////////////////////////////
    // try to create the grids (from the given grid file)
    //////////////////////////////////////////////////////
    using BulkProblemTypeTag = TTAG(BULKTYPETAG);
    using LowDimProblemTypeTag = TTAG(LOWDIMTYPETAG);
    using BulkGrid = typename GET_PROP_TYPE(BulkProblemTypeTag, Grid);
    using LowDimGrid = typename GET_PROP_TYPE(LowDimProblemTypeTag, Grid);

    using GridManager = FacetCouplingGridManager<BulkGrid, LowDimGrid>;
    GridManager gridManager;
    gridManager.init();
    gridManager.loadBalance();

    ////////////////////////////////////////////////////////////
    // run stationary, non-linear problem on this grid
    ////////////////////////////////////////////////////////////

    // we compute on the leaf grid views
    const auto& bulkGridView = gridManager.template grid<0>().leafGridView();
    const auto& lowDimGridView = gridManager.template grid<1>().leafGridView();

    // create the finite volume grid geometries
    using BulkFVGridGeometry = typename GET_PROP_TYPE(BulkProblemTypeTag, FVGridGeometry);
    using LowDimFVGridGeometry = typename GET_PROP_TYPE(LowDimProblemTypeTag, FVGridGeometry);
    auto bulkFvGridGeometry = std::make_shared<BulkFVGridGeometry>(bulkGridView);
    auto lowDimFvGridGeometry = std::make_shared<LowDimFVGridGeometry>(lowDimGridView);
    bulkFvGridGeometry->update();
    lowDimFvGridGeometry->update();

    // the coupling mapper
    using TestTraits = TestTraits<BulkProblemTypeTag, LowDimProblemTypeTag>;
    auto couplingMapper = std::make_shared<typename TestTraits::CouplingMapper>();
    couplingMapper->update(*bulkFvGridGeometry, *lowDimFvGridGeometry, gridManager.getEmbeddings());

    // the coupling manager
    using CouplingManager = typename TestTraits::CouplingManager;
    auto couplingManager = std::make_shared<CouplingManager>();

    // the problems (boundary conditions)
    using BulkProblem = typename GET_PROP_TYPE(BulkProblemTypeTag, Problem);
    using LowDimProblem = typename GET_PROP_TYPE(LowDimProblemTypeTag, Problem);
    auto bulkSpatialParams = std::make_shared<typename BulkProblem::SpatialParams>(bulkFvGridGeometry, "Bulk");
    auto bulkProblem = std::make_shared<BulkProblem>(bulkFvGridGeometry, bulkSpatialParams, couplingManager, "Bulk");
    auto lowDimSpatialParams = std::make_shared<typename LowDimProblem::SpatialParams>(lowDimFvGridGeometry, "LowDim");
    auto lowDimProblem = std::make_shared<LowDimProblem>(lowDimFvGridGeometry, lowDimSpatialParams, couplingManager, "LowDim");

    // the solution vector
    using MDTraits = typename TestTraits::MDTraits;
    using SolutionVector = typename MDTraits::SolutionVector;
    SolutionVector x;

    static const auto bulkId = typename MDTraits::template SubDomain<0>::Index();
    static const auto lowDimId = typename MDTraits::template SubDomain<1>::Index();
    x[bulkId].resize(bulkFvGridGeometry->numDofs());
    x[lowDimId].resize(lowDimFvGridGeometry->numDofs());
    bulkProblem->applyInitialSolution(x[bulkId]);
    lowDimProblem->applyInitialSolution(x[lowDimId]);

    // initialize coupling manager
    couplingManager->init(bulkProblem, lowDimProblem, couplingMapper, x);

    // the grid variables
    using BulkGridVariables = typename GET_PROP_TYPE(BulkProblemTypeTag, GridVariables);
    using LowDimGridVariables = typename GET_PROP_TYPE(LowDimProblemTypeTag, GridVariables);
    auto bulkGridVariables = std::make_shared<BulkGridVariables>(bulkProblem, bulkFvGridGeometry);
    auto lowDimGridVariables = std::make_shared<LowDimGridVariables>(lowDimProblem, lowDimFvGridGeometry);
    bulkGridVariables->init(x[bulkId]);
    lowDimGridVariables->init(x[lowDimId]);

    // intialize the vtk output module
    using BulkSolutionVector = std::decay_t<decltype(x[bulkId])>;
    using LowDimSolutionVector = std::decay_t<decltype(x[lowDimId])>;
    VtkOutputModule<BulkGridVariables, BulkSolutionVector> bulkVtkWriter(*bulkGridVariables, x[bulkId], bulkProblem->name(), "Bulk");
    VtkOutputModule<LowDimGridVariables, LowDimSolutionVector> lowDimVtkWriter(*lowDimGridVariables, x[lowDimId], lowDimProblem->name(), "LowDim");

    // container for the output of the exact solutions
    std::vector<typename GET_PROP_TYPE(BulkProblemTypeTag, Scalar)> bulkExact;
    std::vector<typename GET_PROP_TYPE(LowDimProblemTypeTag, Scalar)> lowDimExact;

    // Add model specific output fields
    using BulkIOFields = typename GET_PROP_TYPE(BulkProblemTypeTag, IOFields);
    using LowDimIOFields = typename GET_PROP_TYPE(LowDimProblemTypeTag, IOFields);
    BulkIOFields::initOutputModule(bulkVtkWriter);
    LowDimIOFields::initOutputModule(lowDimVtkWriter);

    // write initial solution
    bulkVtkWriter.write(0.0);
    lowDimVtkWriter.write(0.0);

    // the assembler
    using Assembler = MultiDomainFVAssembler<MDTraits, CouplingManager, DiffMethod::numeric, /*implicit?*/true>;
    auto assembler = std::make_shared<Assembler>( std::make_tuple(bulkProblem, lowDimProblem),
                                                  std::make_tuple(bulkFvGridGeometry, lowDimFvGridGeometry),
                                                  std::make_tuple(bulkGridVariables, lowDimGridVariables),
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
    bulkGridVariables->update(x[bulkId]);
    lowDimGridVariables->update(x[lowDimId]);

    // write vtk output
    bulkVtkWriter.write(1.0);
    lowDimVtkWriter.write(1.0);

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
