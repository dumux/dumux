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
 * \ingroup NavierStokesTests
 * \brief Test for the staggered grid Navier-Stokes model.
 */

#include <config.h>

#include <iostream>
#include <dune/common/parallel/mpihelper.hh>
#include <dune/common/fvector.hh>

#include <dumux/common/dumuxmessage.hh>
#include <dumux/common/parameters.hh>
#include <dumux/common/properties.hh>

#include <dumux/assembly/fvassembler.hh>

#include <dumux/io/grid/gridmanager.hh>
#include <dumux/io/vtkoutputmodule.hh>
#include <dumux/io/vtk/function.hh>
#include <dumux/io/vtk/intersectionwriter.hh>

#include <dumux/linear/seqsolverbackend.hh>
#include <dumux/nonlinear/newtonsolver.hh>

#include <dumux/linear/linearsolvertraits.hh>
#include <dumux/linear/istlsolverfactorybackend.hh>

#include <test/freeflow/navierstokes/analyticalsolutionvectors.hh>
#include <test/freeflow/navierstokes/errors.hh>

#include "properties_momentum.hh"

namespace Dumux {

template<class Error>
void writeError_(std::ofstream& logFile, const Error& error, const std::string& format = "{:.5e}")
{
    for (const auto& e : error)
        logFile << Fmt::format(", " + format, e);
}

template<class Problem, class GridVariables, class SolutionVector>
void printErrors(std::shared_ptr<Problem> problem,
                 const GridVariables& gridVariables,
                 const SolutionVector& x)
{
    using GridGeometry = std::decay_t<decltype(std::declval<Problem>().gridGeometry())>;
    static constexpr int dim = GridGeometry::GridView::dimension;
    const bool printErrors = getParam<bool>("Problem.PrintErrors", false);

    if (printErrors)
    {
        NavierStokesTest::ErrorsSubProblem errors(problem, x);

        std::ofstream logFile("errors.csv", std::ios::app);
        auto totalVolume = errors.totalVolume();
        // For the staggered scheme, the control volumes are overlapping
        if constexpr (GridGeometry::discMethod == Dumux::DiscretizationMethods::fcstaggered)
            totalVolume /= dim;

        logFile << Fmt::format("{:.5e}", errors.time()) << ", ";
        logFile << problem->gridGeometry().numDofs() << ", ";
        logFile << std::pow(totalVolume / problem->gridGeometry().numDofs(), 1.0/dim);
        const auto& componentErrors = errors.l2Absolute();
        // Calculate L2-error for velocity field
        Dune::FieldVector<double, 1> velError(0.0);
        velError[0] = std::sqrt(componentErrors * componentErrors);

        writeError_(logFile, velError);
        //writeError_(logFile, errors.l2Relative());
        //writeError_(logFile, errors.lInfAbsolute());
        //writeError_(logFile, errors.lInfRelative());

        logFile << "\n";
    }
}

template<class GridGeometry, class GridVariables, class SolutionVector>
void updateVelocities(
    std::vector<Dune::FieldVector<double, 3>>& velocity,
    std::vector<Dune::FieldVector<double, 3>>& faceVelocity,
    const GridGeometry& gridGeometry,
    const GridVariables& gridVariables,
    const SolutionVector& x
){
    auto fvGeometry = localView(gridGeometry);
    auto elemVolVars = localView(gridVariables.curGridVolVars());

    using ShapeValue = typename Dune::FieldVector<double, 1>;
    std::vector<ShapeValue> shapeValues;

    for (const auto& element : elements(gridGeometry.gridView()))
    {
        fvGeometry.bind(element);
        elemVolVars.bind(element, fvGeometry, x);
        const auto eIdx = gridGeometry.elementMapper().index(element);

        if constexpr (GridGeometry::discMethod == Dumux::DiscretizationMethods::fcstaggered)
        {
            for (const auto& scv : scvs(fvGeometry))
            {
                const auto& vars = elemVolVars[scv];
                velocity[eIdx][scv.dofAxis()] += 0.5*vars.velocity();
                faceVelocity[scv.dofIndex()][scv.dofAxis()] = vars.velocity();
            }
        }
        else if constexpr (GridGeometry::discMethod == Dumux::DiscretizationMethods::fcdiamond)
        {
            auto elemSol = elementSolution(element, x, gridGeometry);
            const auto elemGeom = element.geometry();

            // TODO use evalSolution (would be good to have eval solution for a list of coordinates too)
            const auto& localBasis = fvGeometry.feLocalBasis();
            const auto centerLocal = elemGeom.local(elemGeom.center());
            localBasis.evaluateFunction(centerLocal, shapeValues);

            velocity[eIdx] = 0.0;
            for (const auto& scv : scvs(fvGeometry))
            {
                const auto& volVars = elemVolVars[scv];
                velocity[eIdx].axpy(shapeValues[scv.indexInElement()][0], volVars.velocity());
                faceVelocity[scv.dofIndex()] = volVars.velocity();
            }
        }
        else
            DUNE_THROW(Dune::Exception, "Unknown discretization type: " << GridGeometry::discMethod);
    }
}

template<class GridGeometry>
void updateRank(
    std::vector<int>& rank,
    const GridGeometry& gridGeometry
){
    for (const auto& element : elements(gridGeometry.gridView(), Dune::Partitions::interior))
    {
        const auto eIdxGlobal = gridGeometry.elementMapper().index(element);
        rank[eIdxGlobal] = gridGeometry.gridView().comm().rank();
    }
}

} // end namespace Dumux


int main(int argc, char** argv)
{
    using namespace Dumux;

    using TypeTag = Properties::TTag::TrigonometricTestMomentum;

    const auto& mpiHelper = Dune::MPIHelper::instance(argc, argv);

    if (mpiHelper.rank() == 0)
        DumuxMessage::print(/*firstCall=*/true);

    Parameters::init(argc, argv);

    using GridManager = Dumux::GridManager<GetPropType<TypeTag, Properties::Grid>>;
    GridManager gridManager;
    gridManager.init();

    ////////////////////////////////////////////////////////////
    // solve Stokes problem on this grid
    ////////////////////////////////////////////////////////////

    const auto& leafGridView = gridManager.grid().leafGridView();

    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    auto gridGeometry = std::make_shared<GridGeometry>(leafGridView);

    using Problem = GetPropType<TypeTag, Properties::Problem>;
    auto problem = std::make_shared<Problem>(gridGeometry);

    using SolutionVector = GetPropType<TypeTag, Properties::SolutionVector>;
    SolutionVector x(gridGeometry->numDofs());

    using GridVariables = GetPropType<TypeTag, Properties::GridVariables>;
    auto gridVariables = std::make_shared<GridVariables>(problem, gridGeometry);
    gridVariables->init(x);

    std::vector<Dune::FieldVector<double, 3>> velocity(gridGeometry->gridView().size(0));
    std::vector<Dune::FieldVector<double, 3>> faceVelocity(x.size());
    Dumux::updateVelocities(velocity, faceVelocity, *gridGeometry, *gridVariables, x);

    std::vector<int> rank(gridGeometry->gridView().size(0));
    Dumux::updateRank(rank, *gridGeometry);

    std::vector<std::size_t> dofIdx(x.size());
    for (const auto& facet : facets(gridGeometry->gridView()))
    {
        const auto idx = gridGeometry->gridView().indexSet().index(facet);
        dofIdx[idx] = idx;
    }

    std::string discSuffix = std::string("_") + GridGeometry::DiscretizationMethod::name();
    std::string rankSuffix = std::string("_") + std::to_string(gridGeometry->gridView().comm().rank());
    Dune::VTKWriter<typename GridGeometry::GridView> writer(gridGeometry->gridView());
    using Field = Vtk::template Field<typename GridGeometry::GridView>;
    writer.addCellData(Field(
        gridGeometry->gridView(), gridGeometry->elementMapper(), velocity,
        "velocity", /*numComp*/3, /*codim*/0
    ).get());
    writer.addCellData(rank, "rank");
    writer.write("trigonometric_momentum" + discSuffix + "_0");

    ConformingIntersectionWriter faceVtk(gridGeometry->gridView());
    faceVtk.addField(dofIdx, "dofIdx");
    faceVtk.addField(faceVelocity, "velocityVector");
    faceVtk.addField([&](const auto& is, const auto idx) {
        const auto& facet = is.inside().template subEntity <1> (is.indexInInside());
        return facet.partitionType();
    }, "partitionType");
    faceVtk.write("trigonometric_momentum_face" + discSuffix + rankSuffix + "_0", Dune::VTK::ascii);

    using Assembler = FVAssembler<TypeTag, DiffMethod::numeric>;
    auto assembler = std::make_shared<Assembler>(problem, gridGeometry, gridVariables);

    using LinearSolver = IstlSolverFactoryBackend<LinearSolverTraits<GridGeometry>>;
    const auto dofMapper = LinearSolverTraits<GridGeometry>::DofMapper(gridGeometry->gridView());
    auto linearSolver = std::make_shared<LinearSolver>(gridGeometry->gridView(), dofMapper);

    using NewtonSolver = Dumux::NewtonSolver<Assembler, LinearSolver>;
    NewtonSolver nonLinearSolver(assembler, linearSolver);

    nonLinearSolver.solve(x);

    ////////////////////////////////////////////////////////////
    // write VTK output
    ////////////////////////////////////////////////////////////
    Dumux::updateVelocities(velocity, faceVelocity, *gridGeometry, *gridVariables, x);
    writer.write("trigonometric_momentum" + discSuffix + "_1");
    faceVtk.write("trigonometric_momentum_face" + discSuffix + rankSuffix + "_1", Dune::VTK::ascii);

    Dumux::printErrors(problem, *gridVariables, x);

    ////////////////////////////////////////////////////////////
    // finalize, print parameters and Dumux message to say goodbye
    ////////////////////////////////////////////////////////////

    // print dumux end message
    if (mpiHelper.rank() == 0)
    {
        Parameters::print();
        DumuxMessage::print(/*firstCall=*/false);
    }

    return 0;
}
