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
 * \brief Test for the staggered grid Navier-Stokes model (Donea 2003, \cite Donea2003).
 */

#include <config.h>

#include <iostream>
#include <dune/common/parallel/mpihelper.hh>

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
#include <dumux/linear/algebratraits.hh>
#include <dumux/linear/istlsolverfactorybackend.hh>

#include "properties_momentum.hh"

template<class GridGeometry, class GridVariables, class SolutionVector>
void updateVelocities(
    std::vector<Dune::FieldVector<double, 2>>& velocity,
    std::vector<Dune::FieldVector<double, 2>>& faceVelocity,
    const GridGeometry& gridGeometry,
    const GridVariables& gridVariables,
    const SolutionVector& x
){
    auto fvGeometry = localView(gridGeometry);
    auto elemVolVars = localView(gridVariables.curGridVolVars());
    for (const auto& element : elements(gridGeometry.gridView()))
    {
        fvGeometry.bind(element);
        elemVolVars.bind(element, fvGeometry, x);

        for (const auto& scv : scvs(fvGeometry))
        {
            const auto& vars = elemVolVars[scv];
            velocity[scv.elementIndex()][scv.dofAxis()] += 0.5*vars.velocity();
            faceVelocity[scv.dofIndex()][scv.dofAxis()] = vars.velocity();
        }
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

int main(int argc, char** argv)
{
    using namespace Dumux;

    using TypeTag = Properties::TTag::DoneaTestMomentum;

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

    std::vector<Dune::FieldVector<double, 2>> velocity(gridGeometry->gridView().size(0));
    std::vector<Dune::FieldVector<double, 2>> faceVelocity(x.size());
    updateVelocities(velocity, faceVelocity, *gridGeometry, *gridVariables, x);

    std::vector<int> rank(gridGeometry->gridView().size(0));
    updateRank(rank, *gridGeometry);

    std::vector<std::size_t> dofIdx(x.size());
    for (const auto& facet : facets(gridGeometry->gridView()))
    {
        const auto idx = gridGeometry->gridView().indexSet().index(facet);
        dofIdx[idx] = idx;
    }

    Dune::VTKWriter<typename GridGeometry::GridView> writer(gridGeometry->gridView());
    using Field = Vtk::template Field<typename GridGeometry::GridView>;
    writer.addCellData(Field(
        gridGeometry->gridView(), gridGeometry->elementMapper(), velocity,
        "velocity", /*numComp*/2, /*codim*/0
    ).get());
    writer.addCellData(rank, "rank");
    writer.write("donea_momentum_0");

    ConformingIntersectionWriter faceVtk(gridGeometry->gridView());
    faceVtk.addField(dofIdx, "dofIdx");
    faceVtk.addField(faceVelocity, "velocityVector");
    faceVtk.addField([&](const auto& is, const auto idx) {
        const auto& facet = is.inside().template subEntity <1> (is.indexInInside());
        return facet.partitionType();
    }, "partitionType");
    faceVtk.write("donea_momentum_face_0_" + std::to_string(gridGeometry->gridView().comm().rank()), Dune::VTK::ascii);

    using Assembler = FVAssembler<TypeTag, DiffMethod::numeric>;
    auto assembler = std::make_shared<Assembler>(problem, gridGeometry, gridVariables);

    using LinearSolver = IstlSolverFactoryBackend<LinearSolverTraits<GridGeometry>,
                                                  LinearAlgebraTraitsFromAssembler<Assembler>>;
    const auto dofMapper = LinearSolverTraits<GridGeometry>::DofMapper(gridGeometry->gridView());
    auto linearSolver = std::make_shared<LinearSolver>(gridGeometry->gridView(), dofMapper);

    using NewtonSolver = Dumux::NewtonSolver<Assembler, LinearSolver>;
    NewtonSolver nonLinearSolver(assembler, linearSolver);

    nonLinearSolver.solve(x);

    ////////////////////////////////////////////////////////////
    // write VTK output
    ////////////////////////////////////////////////////////////
    updateVelocities(velocity, faceVelocity, *gridGeometry, *gridVariables, x);
    writer.write("donea_momentum_1");
    faceVtk.write("donea_momentum_face_1_" + std::to_string(gridGeometry->gridView().comm().rank()), Dune::VTK::ascii);

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
