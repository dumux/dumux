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
 *
 * \brief test for the one-phase facet coupling model
 */
#include <config.h>

#include <iostream>

#include <dune/common/parallel/mpihelper.hh>

#include "bulkproblem.hh"
#include "lowdimproblem.hh"

#include <dumux/common/properties.hh>
#include <dumux/common/parameters.hh>
#include <dumux/common/dumuxmessage.hh>
#include <dumux/common/defaultusagemessage.hh>

#include <dumux/assembly/diffmethod.hh>

#include <dumux/linear/seqsolverbackend.hh>
#include <dumux/multidomain/newtonsolver.hh>
#include <dumux/multidomain/fvassembler.hh>
#include <dumux/multidomain/traits.hh>

#include <dumux/multidomain/facet/gridcreator.hh>
#include <dumux/multidomain/facet/cellcentered/tpfa/couplingmapper.hh>
#include <dumux/multidomain/facet/cellcentered/tpfa/couplingmanager.hh>

#include <dumux/io/vtkoutputmodule.hh>

// set the coupling manager property in the sub-problems
namespace Dumux {
namespace Properties {

NEW_PROP_TAG(CouplingManager);

SET_PROP(OnePLowDimTpfa, CouplingManager)
{
private:
    // define traits etc. as below in main
    using Traits = MultiDomainTraits<TTAG(OnePBulkTpfa), TTAG(OnePLowDimTpfa)>;
    using BulkFVG = typename GET_PROP_TYPE(TTAG(OnePBulkTpfa), FVGridGeometry);
    using LowDimFVG = typename GET_PROP_TYPE(TTAG(OnePLowDimTpfa), FVGridGeometry);
    using CouplingMapper = CCTpfaFacetCouplingMapper<BulkFVG, LowDimFVG>;
public:
    using type = CCTpfaFacetCouplingManager<Traits, CouplingMapper>;
};

SET_PROP(OnePBulkTpfa, CouplingManager)
{
private:
    // define traits etc. as below in main
    using Traits = MultiDomainTraits<TTAG(OnePBulkTpfa), TTAG(OnePLowDimTpfa)>;
    using BulkFVG = typename GET_PROP_TYPE(TTAG(OnePBulkTpfa), FVGridGeometry);
    using LowDimFVG = typename GET_PROP_TYPE(TTAG(OnePLowDimTpfa), FVGridGeometry);
    using CouplingMapper = CCTpfaFacetCouplingMapper<BulkFVG, LowDimFVG>;
public:
    using type = CCTpfaFacetCouplingManager<Traits, CouplingMapper>;
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
    using BulkProblemTypeTag = TTAG(OnePBulkTpfa);
    using LowDimProblemTypeTag = TTAG(OnePLowDimTpfa);
    using BulkGrid = typename GET_PROP_TYPE(BulkProblemTypeTag, Grid);
    using LowDimGrid = typename GET_PROP_TYPE(LowDimProblemTypeTag, Grid);

    using GridCreator = FacetCouplingGridCreator<BulkGrid, LowDimGrid>;
    GridCreator gridCreator;
    gridCreator.makeGrids(getParam<std::string>("Grid.File"));
    gridCreator.loadBalance();

    ////////////////////////////////////////////////////////////
    // run stationary non-linear problem on this grid
    ////////////////////////////////////////////////////////////

    // we compute on the leaf grid views
    const auto& bulkGridView = gridCreator.template grid<0>().leafGridView();
    const auto& lowDimGridView = gridCreator.template grid<1>().leafGridView();

    // create the finite volume grid geometries
    using BulkFVGridGeometry = typename GET_PROP_TYPE(BulkProblemTypeTag, FVGridGeometry);
    using LowDimFVGridGeometry = typename GET_PROP_TYPE(LowDimProblemTypeTag, FVGridGeometry);
    auto bulkFvGridGeometry = std::make_shared<BulkFVGridGeometry>(bulkGridView);
    auto lowDimFvGridGeometry = std::make_shared<LowDimFVGridGeometry>(lowDimGridView);
    bulkFvGridGeometry->update();
    lowDimFvGridGeometry->update();

    // the problems (boundary conditions)
    using BulkProblem = typename GET_PROP_TYPE(BulkProblemTypeTag, Problem);
    using LowDimProblem = typename GET_PROP_TYPE(LowDimProblemTypeTag, Problem);
    auto bulkSpatialParams = std::make_shared<typename BulkProblem::SpatialParams>(bulkFvGridGeometry, "Bulk");
    auto bulkProblem = std::make_shared<BulkProblem>(bulkFvGridGeometry, bulkSpatialParams, "Bulk");
    auto lowDimSpatialParams = std::make_shared<typename LowDimProblem::SpatialParams>(lowDimFvGridGeometry, "LowDim");
    auto lowDimProblem = std::make_shared<LowDimProblem>(lowDimFvGridGeometry, lowDimSpatialParams, "LowDim");

    // the solution vector
    using Traits = MultiDomainTraits<BulkProblemTypeTag, LowDimProblemTypeTag>;
    using SolutionVector = typename Traits::SolutionVector;
    SolutionVector x;

    static const auto bulkId = Traits::template DomainIdx<0>();
    static const auto lowDimId = Traits::template DomainIdx<1>();
    x[bulkId].resize(bulkFvGridGeometry->numDofs());
    x[lowDimId].resize(lowDimFvGridGeometry->numDofs());
    bulkProblem->applyInitialSolution(x[bulkId]);
    lowDimProblem->applyInitialSolution(x[lowDimId]);

    // the coupling mapper
    using CouplingMapper = CCTpfaFacetCouplingMapper<BulkFVGridGeometry, LowDimFVGridGeometry>;
    auto couplingMapper = std::make_shared<CouplingMapper>();
    couplingMapper->update(*bulkFvGridGeometry, *lowDimFvGridGeometry, gridCreator);

    // the coupling manager
    using CouplingManager = CCTpfaFacetCouplingManager<Traits, CouplingMapper>;
    auto couplingManager = std::make_shared<CouplingManager>();
    couplingManager->init(bulkProblem, lowDimProblem, couplingMapper, x);

    // set coupling manager pointer in sub-problems
    bulkProblem->setCouplingManager(couplingManager);
    lowDimProblem->setCouplingManager(couplingManager);

    // the grid variables
    using BulkGridVariables = typename GET_PROP_TYPE(BulkProblemTypeTag, GridVariables);
    using LowDimGridVariables = typename GET_PROP_TYPE(LowDimProblemTypeTag, GridVariables);
    auto bulkGridVariables = std::make_shared<BulkGridVariables>(bulkProblem, bulkFvGridGeometry);
    auto lowDimGridVariables = std::make_shared<LowDimGridVariables>(lowDimProblem, lowDimFvGridGeometry);
    bulkGridVariables->init(x[bulkId]);
    lowDimGridVariables->init(x[lowDimId]);

    // intialize the vtk output module
    VtkOutputModule<BulkProblemTypeTag> bulkVtkWriter(*bulkProblem, *bulkFvGridGeometry, *bulkGridVariables, x[bulkId], bulkProblem->name());
    VtkOutputModule<LowDimProblemTypeTag> lowDimVtkWriter(*lowDimProblem, *lowDimFvGridGeometry, *lowDimGridVariables, x[lowDimId], lowDimProblem->name());

    // Add model specific output fields
    using BulkVtkOutputFields = typename GET_PROP_TYPE(BulkProblemTypeTag, VtkOutputFields);
    using LowDimVtkOutputFields = typename GET_PROP_TYPE(LowDimProblemTypeTag, VtkOutputFields);
    BulkVtkOutputFields::init(bulkVtkWriter);
    LowDimVtkOutputFields::init(lowDimVtkWriter);
    bulkVtkWriter.write(0.0);
    lowDimVtkWriter.write(0.0);

    // the assembler
    using Assembler = MultiDomainFVAssembler<Traits, CouplingManager, DiffMethod::numeric, /*implicit?*/true>;
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
