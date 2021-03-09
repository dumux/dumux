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
 * \ingroup FacetTests
 * \brief Test for the one-phase facet coupling model.
 */

#include <config.h>

#include <iostream>

#include <dune/common/parallel/mpihelper.hh>

#include <dumux/common/properties.hh>
#include <dumux/common/parameters.hh>
#include <dumux/common/dumuxmessage.hh>

#include <dumux/assembly/diffmethod.hh>
#include <dumux/linear/seqsolverbackend.hh>

#include <dumux/multidomain/newtonsolver.hh>
#include <dumux/multidomain/fvassembler.hh>
#include <dumux/multidomain/traits.hh>

#include <dumux/multidomain/facet/gridmanager.hh>
#include <dumux/multidomain/facet/couplingmapper.hh>
#include <dumux/multidomain/facet/couplingmanager.hh>

#include <dumux/io/vtkoutputmodule.hh>

#include "properties_bulk.hh"
#include "properties_lowdim.hh"

namespace Dumux {

// obtain/define some types to be used below in the property definitions and in main
template< class BulkTypeTag, class LowDimTypeTag >
class TestTraits
{
    using BulkFVGridGeometry = GetPropType<BulkTypeTag, Properties::GridGeometry>;
    using LowDimFVGridGeometry = GetPropType<LowDimTypeTag, Properties::GridGeometry>;
public:
    using MDTraits = Dumux::MultiDomainTraits<BulkTypeTag, LowDimTypeTag>;
    using CouplingMapper = Dumux::FacetCouplingMapper<BulkFVGridGeometry, LowDimFVGridGeometry>;
    using CouplingManager = Dumux::FacetCouplingManager<MDTraits, CouplingMapper>;
};

namespace Properties {

// set cm property in the sub-problems
using TpfaTraits = TestTraits<TTag::OnePBulkTpfa, TTag::OnePLowDimTpfa>;
using MpfaTraits = TestTraits<TTag::OnePBulkMpfa, TTag::OnePLowDimMpfa>;
template<class TypeTag> struct CouplingManager<TypeTag, TTag::OnePBulkTpfa> { using type = typename TpfaTraits::CouplingManager; };
template<class TypeTag> struct CouplingManager<TypeTag, TTag::OnePLowDimTpfa> { using type = typename TpfaTraits::CouplingManager; };
template<class TypeTag> struct CouplingManager<TypeTag, TTag::OnePBulkMpfa> { using type = typename MpfaTraits::CouplingManager; };
template<class TypeTag> struct CouplingManager<TypeTag, TTag::OnePLowDimMpfa> { using type = typename MpfaTraits::CouplingManager; };
} // end namespace Properties
} // end namespace Dumux

// main program
int main(int argc, char** argv)
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
    using BulkProblemTypeTag = Properties::TTag::BULKTYPETAG;
    using LowDimProblemTypeTag = Properties::TTag::LOWDIMTYPETAG;
    using BulkGrid = GetPropType<BulkProblemTypeTag, Properties::Grid>;
    using LowDimGrid = GetPropType<LowDimProblemTypeTag, Properties::Grid>;

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
    using BulkFVGridGeometry = GetPropType<BulkProblemTypeTag, Properties::GridGeometry>;
    using LowDimFVGridGeometry = GetPropType<LowDimProblemTypeTag, Properties::GridGeometry>;
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
    using BulkProblem = GetPropType<BulkProblemTypeTag, Properties::Problem>;
    using LowDimProblem = GetPropType<LowDimProblemTypeTag, Properties::Problem>;
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
    using BulkGridVariables = GetPropType<BulkProblemTypeTag, Properties::GridVariables>;
    using LowDimGridVariables = GetPropType<LowDimProblemTypeTag, Properties::GridVariables>;
    auto bulkGridVariables = std::make_shared<BulkGridVariables>(bulkProblem, bulkFvGridGeometry);
    auto lowDimGridVariables = std::make_shared<LowDimGridVariables>(lowDimProblem, lowDimFvGridGeometry);
    bulkGridVariables->init(x[bulkId]);
    lowDimGridVariables->init(x[lowDimId]);

    // intialize the vtk output module
    using BulkSolutionVector = std::decay_t<decltype(x[bulkId])>;
    using LowDimSolutionVector = std::decay_t<decltype(x[lowDimId])>;
    VtkOutputModule<BulkGridVariables, BulkSolutionVector> bulkVtkWriter(*bulkGridVariables, x[bulkId], bulkProblem->name(), "Bulk");
    VtkOutputModule<LowDimGridVariables, LowDimSolutionVector> lowDimVtkWriter(*lowDimGridVariables, x[lowDimId], lowDimProblem->name(), "LowDim");

    // Add model specific output fields
    using BulkIOFields = GetPropType<BulkProblemTypeTag, Properties::IOFields>;
    using LowIOFields = GetPropType<LowDimProblemTypeTag, Properties::IOFields>;
    BulkIOFields::initOutputModule(bulkVtkWriter);
    LowIOFields::initOutputModule(lowDimVtkWriter);

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
