// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \brief Test of the dof indices and the basis functions of the elements used by the embedded coupling
 *
 * For control-volume finite element and finite element discretizations
 * with dofs on vertices, edges and in elements, the test checks that
 * - the dof indices of the elements are in the range of the dofs and cover all dofs,
 * - the dof of an interpolation point has the same global index in all elements sharing it,
 * - for nodal basis functions, the basis functions of an element sum to one, and interpolating an
 *   affine function from its values at the dof positions with the dof indices and basis functions
 *   reproduces the function at random points in the elements.
 */
#include <config.h>

#include <algorithm>
#include <cmath>
#include <iostream>
#include <map>
#include <random>
#include <string>
#include <vector>

#include <dune/common/exceptions.hh>
#include <dune/common/fvector.hh>
#include <dune/grid/yaspgrid.hh>

#include <dumux/common/initialize.hh>
#include <dumux/discretization/box.hh>
#include <dumux/discretization/pq1.hh>
#include <dumux/discretization/pq1bubble.hh>
#include <dumux/discretization/pq2.hh>
#include <dumux/multidomain/embedded/elementdofs.hh>

namespace Dumux {

constexpr int dim = 3;
using GlobalPosition = Dune::FieldVector<double, dim>;
using Grid = Dune::YaspGrid<dim>;
using GridView = Grid::LeafGridView;

double affineFunction(const GlobalPosition& x)
{ return 1.0 + 2.0*x[0] + 3.0*x[1] - 0.5*x[2]; }

// nodalBasis: the basis functions form a partition of unity and the coefficient of a function
// is its value at the dof position
template<class GridGeometry>
void check(const std::string& name, const GridGeometry& gridGeometry, bool nodalBasis = true)
{
    std::cout << "Checking " << name << ", " << gridGeometry.numDofs() << " dofs" << std::endl;

    std::vector<double> dofValues(gridGeometry.numDofs(), 0.0);
    std::vector<bool> dofSet(gridGeometry.numDofs(), false);
    std::mt19937 generator(1);
    std::uniform_real_distribution<double> distribution(0.0, 1.0);

    auto elemDisc = localView(gridGeometry);

    // the dof values of the affine function at the positions of the dofs
    for (const auto& element : elements(gridGeometry.gridView()))
    {
        elemDisc.bindElement(element);
        const auto dofIndices = EmbeddedCoupling::elementDofIndices(elemDisc);
        if (dofIndices.size() != gridGeometry.feCache().get(element.type()).localBasis().size())
            DUNE_THROW(Dune::InvalidStateException, name << ": number of dof indices differs from number of basis functions");

        for (const auto& localDof : localDofs(elemDisc))
        {
            if (dofIndices[localDof.index()] != localDof.dofIndex())
                DUNE_THROW(Dune::InvalidStateException, name << ": dof indices are not in the order of the local dof indices");
            if (localDof.dofIndex() >= gridGeometry.numDofs())
                DUNE_THROW(Dune::InvalidStateException, name << ": dof index " << localDof.dofIndex() << " out of range");

            const auto value = affineFunction(ipData(elemDisc, localDof).global());
            if (dofSet[localDof.dofIndex()] && std::abs(dofValues[localDof.dofIndex()] - value) > 1e-13)
                DUNE_THROW(Dune::InvalidStateException, name << ": dof " << localDof.dofIndex() << " has different positions in different elements");
            dofValues[localDof.dofIndex()] = value;
            dofSet[localDof.dofIndex()] = true;
        }
    }

    if (!std::all_of(dofSet.begin(), dofSet.end(), [](bool b){ return b; }))
        DUNE_THROW(Dune::InvalidStateException, name << ": not all dofs belong to an element");

    if (!nodalBasis)
        return;

    // interpolation at random points
    double maxUnityError = 0.0, maxInterpolationError = 0.0;
    for (const auto& element : elements(gridGeometry.gridView()))
    {
        elemDisc.bindElement(element);
        const auto dofIndices = EmbeddedCoupling::elementDofIndices(elemDisc);
        const auto geometry = element.geometry();
        for (int i = 0; i < 20; ++i)
        {
            GlobalPosition local;
            for (auto& x : local)
                x = distribution(generator);
            const auto globalPos = geometry.global(local);

            std::vector<Dune::FieldVector<double, 1>> shapeValues;
            EmbeddedCoupling::basisFunctionValues(gridGeometry, geometry, globalPos, shapeValues);

            double sum = 0.0, interpolated = 0.0;
            for (std::size_t j = 0; j < dofIndices.size(); ++j)
            {
                sum += shapeValues[j][0];
                interpolated += shapeValues[j][0]*dofValues[dofIndices[j]];
            }
            maxUnityError = std::max(maxUnityError, std::abs(sum - 1.0));
            maxInterpolationError = std::max(maxInterpolationError, std::abs(interpolated - affineFunction(globalPos)));
        }
    }

    std::cout << "  largest error of the sum of the basis functions: " << maxUnityError
              << ", of the affine interpolation: " << maxInterpolationError << std::endl;
    if (maxUnityError > 1e-12 || maxInterpolationError > 1e-12)
        DUNE_THROW(Dune::InvalidStateException, name << ": basis function error");
}

} // end namespace Dumux

int main(int argc, char** argv)
{
    using namespace Dumux;

    Dumux::initialize(argc, argv);

    Grid grid({1.0, 1.0, 1.0}, {3, 3, 3});
    const auto gridView = grid.leafGridView();

    check("box", BoxFVGridGeometry<double, GridView, true>(gridView));
    check("pq1bubble", PQ1BubbleFVGridGeometry<double, GridView, true>(gridView));
    {
        using QuadTraits = PQ1BubbleQuadratureTraits<GridView>;
        using Traits = HybridPQ1BubbleCVFEGridGeometryTraits<
            PQ1BubbleDefaultGridGeometryTraits<GridView, PQ1BubbleMapperTraits<GridView, 2>, QuadTraits>
        >;
        // two bubble functions per cube: the second one is not a nodal function
        check("pq1bubble hybrid", PQ1BubbleFVGridGeometry<double, GridView, false, Traits>(gridView), false);
    }
    check("pq2 hybrid", PQ2FVGridGeometry<double, GridView, true>(gridView));
    check("pq1 finite element", Experimental::PQ1FEGridDiscretization<double, GridView, true>(gridView));
    check("pq2 finite element", Experimental::PQ2FEGridDiscretization<double, GridView, true>(gridView));

    return 0;
}
