// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup Discretization
 * \brief Test for the Dirichlet constraints of the local dofs on the boundary faces
 *
 * On the unit square or cube, equation 0 is constrained (by primary variable 1) on the faces at x = 0
 * and equation 1 (by primary variable 0) on the faces at y = 0. The constraints are compared with those
 * expected from the position of every local dof of the scheme, on cube and simplex grids.
 */
#include <config.h>

#include <array>
#include <cstddef>
#include <iostream>
#include <memory>
#include <string>
#include <vector>

#include <dune/common/exceptions.hh>
#include <dune/common/fvector.hh>
#include <dune/grid/utility/structuredgridfactory.hh>
#include <dune/grid/yaspgrid.hh>
#if HAVE_DUNE_ALUGRID
#include <dune/alugrid/grid.hh>
#endif

#include <dumux/common/constraintinfo.hh>
#include <dumux/common/initialize.hh>
#include <dumux/discretization/dirichletconstraints.hh>
#include <dumux/discretization/box/fvgridgeometry.hh>
#include <dumux/discretization/pq1bubble/fvgridgeometry.hh>
#include <dumux/discretization/pq2/fvgridgeometry.hh>
#include <dumux/discretization/pq3/fvgridgeometry.hh>
#include <dumux/discretization/facecentered/diamond/fvgridgeometry.hh>
#include <dumux/discretization/pq1/fegriddiscretization.hh>
#include <dumux/discretization/pq2/fegriddiscretization.hh>
#include <dumux/discretization/pq3/fegriddiscretization.hh>

namespace Dumux::Test {

static constexpr int numEq = 2;
using ConstraintInfo = DirichletConstraintInfo<numEq>;
using ConstraintValues = Dune::FieldVector<double, numEq>;
using Constraint = DirichletConstraintData<ConstraintInfo, ConstraintValues, std::size_t>;

constexpr double eps = 1e-10;

template<class GlobalPosition>
ConstraintValues exactValues(const GlobalPosition& globalPos)
{ return {globalPos[0] + 10.0*globalPos[1], 100.0 + globalPos[0]}; }

template<class GridDiscretization>
void testBoundaryDirichletConstraints(const GridDiscretization& gridDiscretization, const std::string& name)
{
    ConstraintInfo existingInfo;
    existingInfo.set(/*pvIdx=*/0, /*eqIdx=*/0);
    const Constraint existing{existingInfo, ConstraintValues(-1.0), gridDiscretization.numDofs() - 1};
    std::vector<Constraint> constraints{existing};

    std::vector<int> numValueEvaluations(gridDiscretization.numDofs(), 0);
    appendBoundaryDirichletConstraints(
        constraints, gridDiscretization,
        [] (const auto& elemDisc, const auto& boundaryFace, const auto& localDof)
        {
            ConstraintInfo info;
            if (boundaryFace.center()[0] < eps)
                info.set(/*pvIdx=*/1, /*eqIdx=*/0);
            if (boundaryFace.center()[1] < eps)
                info.set(/*pvIdx=*/0, /*eqIdx=*/1);
            return info;
        },
        [&] (const auto& elemDisc, const auto& localDof)
        {
            ++numValueEvaluations[localDof.dofIndex()];
            return exactValues(ipData(elemDisc, localDof).global());
        }
    );

    // the equations to constrain at every dof, from its position
    std::vector<std::array<bool, numEq>> expected(gridDiscretization.numDofs(), {false, false});
    std::vector<ConstraintValues> expectedValues(gridDiscretization.numDofs());
    auto elemDisc = localView(gridDiscretization);
    for (const auto& element : elements(gridDiscretization.gridView()))
    {
        elemDisc.bindElement(element);
        for (const auto& localDof : localDofs(elemDisc))
        {
            const auto globalPos = ipData(elemDisc, localDof).global();
            expected[localDof.dofIndex()] = {globalPos[0] < eps, globalPos[1] < eps};
            expectedValues[localDof.dofIndex()] = exactValues(globalPos);
        }
    }

    const auto& kept = constraints.front();
    if (kept.dofIndex() != existing.dofIndex() || kept.values() != existing.values()
        || !kept.constraintInfo().isConstraintEquation(0) || kept.constraintInfo().eqToPriVarIndex(0) != 0
        || kept.constraintInfo().isConstraintEquation(1))
        DUNE_THROW(Dune::Exception, name << ": the constraint already in the vector was changed");

    std::vector<int> numConstraints(gridDiscretization.numDofs(), 0);
    for (std::size_t i = 1; i < constraints.size(); ++i)
    {
        const auto& constraint = constraints[i];
        const auto dofIdx = constraint.dofIndex();
        ++numConstraints[dofIdx];

        const auto& info = constraint.constraintInfo();
        for (int eqIdx = 0; eqIdx < numEq; ++eqIdx)
            if (info.isConstraintEquation(eqIdx) != expected[dofIdx][eqIdx])
                DUNE_THROW(Dune::Exception, name << ": dof " << dofIdx << ", equation " << eqIdx
                                            << " constrained: " << info.isConstraintEquation(eqIdx)
                                            << ", expected: " << expected[dofIdx][eqIdx]);

        if (info.isConstraintEquation(0) && info.eqToPriVarIndex(0) != 1)
            DUNE_THROW(Dune::Exception, name << ": dof " << dofIdx << ", equation 0 constrains primary variable "
                                        << info.eqToPriVarIndex(0) << " instead of 1");
        if (info.isConstraintEquation(1) && info.eqToPriVarIndex(1) != 0)
            DUNE_THROW(Dune::Exception, name << ": dof " << dofIdx << ", equation 1 constrains primary variable "
                                        << info.eqToPriVarIndex(1) << " instead of 0");

        auto difference = constraint.values();
        difference -= expectedValues[dofIdx];
        if (difference.infinity_norm() > eps)
            DUNE_THROW(Dune::Exception, name << ": dof " << dofIdx << " has values " << constraint.values()
                                        << ", expected " << expectedValues[dofIdx]);
    }

    std::size_t numExpected = 0;
    for (std::size_t dofIdx = 0; dofIdx < gridDiscretization.numDofs(); ++dofIdx)
    {
        const bool isConstrained = expected[dofIdx][0] || expected[dofIdx][1];
        numExpected += isConstrained;
        if (numConstraints[dofIdx] != isConstrained)
            DUNE_THROW(Dune::Exception, name << ": dof " << dofIdx << " has " << numConstraints[dofIdx]
                                        << " constraints, expected " << isConstrained);
        if (numValueEvaluations[dofIdx] != isConstrained)
            DUNE_THROW(Dune::Exception, name << ": the values of dof " << dofIdx << " were evaluated "
                                        << numValueEvaluations[dofIdx] << " times, expected " << isConstrained);
    }

    // dofs off the constrained faces get no constraint at all
    std::vector<Constraint> none;
    appendBoundaryDirichletConstraints(
        none, gridDiscretization,
        [] (const auto&, const auto&, const auto&) { return ConstraintInfo{}; },
        [] (const auto& elemDisc, const auto& localDof) { return ConstraintValues(0.0); }
    );
    if (!none.empty())
        DUNE_THROW(Dune::Exception, name << ": " << none.size() << " constraints without constrained equations");

    std::cout << name << ": " << constraints.size() - 1 << " constrained dofs of "
              << gridDiscretization.numDofs() << ", as expected " << numExpected << std::endl;
}

template<class GridView>
void testAllSchemes(const GridView& gridView, const std::string& prefix)
{
    testBoundaryDirichletConstraints(BoxFVGridGeometry<double, GridView, true>(gridView), prefix + "box (cached)");
    testBoundaryDirichletConstraints(BoxFVGridGeometry<double, GridView, false>(gridView), prefix + "box");
    testBoundaryDirichletConstraints(PQ1BubbleFVGridGeometry<double, GridView, true>(gridView), prefix + "pq1bubble");
    testBoundaryDirichletConstraints(PQ2FVGridGeometry<double, GridView, true>(gridView), prefix + "pq2");
    testBoundaryDirichletConstraints(PQ3FVGridGeometry<double, GridView, true>(gridView), prefix + "pq3");
    testBoundaryDirichletConstraints(FaceCenteredDiamondFVGridGeometry<GridView, true>(gridView), prefix + "diamond");
    testBoundaryDirichletConstraints(Experimental::PQ1FEGridDiscretization<double, GridView>(gridView), prefix + "pq1 fe");
    testBoundaryDirichletConstraints(Experimental::PQ2FEGridDiscretization<double, GridView>(gridView), prefix + "pq2 fe");
    testBoundaryDirichletConstraints(Experimental::PQ3FEGridDiscretization<double, GridView>(gridView), prefix + "pq3 fe");
}

template<int dim>
void testCubeGrid()
{
    using Grid = Dune::YaspGrid<dim>;
    Dune::FieldVector<double, dim> upperRight(1.0);
    std::array<int, dim> cells; cells.fill(3);
    cells[0] = 4;
    const Grid grid(upperRight, cells);
    testAllSchemes(grid.leafGridView(), std::to_string(dim) + "d cube ");
}

#if HAVE_DUNE_ALUGRID
template<int dim>
void testSimplexGrid()
{
    using Grid = Dune::ALUGrid<dim, dim, Dune::simplex, Dune::conforming>;
    const Dune::FieldVector<double, dim> lowerLeft(0.0);
    const Dune::FieldVector<double, dim> upperRight(1.0);
    std::array<unsigned int, dim> cells; cells.fill(3);
    cells[0] = 4;
    const auto grid = Dune::StructuredGridFactory<Grid>::createSimplexGrid(lowerLeft, upperRight, cells);
    testAllSchemes(grid->leafGridView(), std::to_string(dim) + "d simplex ");
}
#endif

} // end namespace Dumux::Test

int main(int argc, char* argv[])
{
    Dumux::initialize(argc, argv);
    Dumux::Test::testCubeGrid<2>();
    Dumux::Test::testCubeGrid<3>();
#if HAVE_DUNE_ALUGRID
    Dumux::Test::testSimplexGrid<2>();
    Dumux::Test::testSimplexGrid<3>();
#endif
    return 0;
}
