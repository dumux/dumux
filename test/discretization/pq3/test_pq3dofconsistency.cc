// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \brief Checks that neighbouring elements map the two dofs of a shared edge to the same global indices,
 *        independent of the orientation in which they traverse the edge.
 */
#include <config.h>

#include <cmath>
#include <iostream>
#include <map>
#include <vector>

#include <dune/common/fvector.hh>
#include <dune/common/exceptions.hh>
#include <dune/localfunctions/common/localkey.hh>
#include <dune/grid/common/gridenums.hh>
#include <dune/grid/utility/structuredgridfactory.hh>
#include <dune/alugrid/grid.hh>

#include <dumux/common/initialize.hh>
#include <dumux/discretization/pq3/dofhelper.hh>
#include <dumux/discretization/pq3/fegriddiscretization.hh>

template<class Grid>
int check(const std::string& name, const Grid& gridRef)
{
    using namespace Dumux;
    using GridView = typename Grid::LeafGridView;
    using Scalar = typename GridView::ctype;
    using GridDiscretization = Experimental::PQ3FEGridDiscretization<Scalar, GridView>;
    using DofHelper = PQ3LagrangeDofHelper<GridView>;

    GridDiscretization gridDisc(gridRef.leafGridView());
    auto elemDisc = localView(gridDisc);

    // the element solution gathers the dofs through dofIndex(element, localKey) if it exists and
    // otherwise through subIndex + localKey.index(), which ignores the orientation of the edge
    static_assert(requires (const GridDiscretization& gd,
                            const typename GridView::template Codim<0>::Entity& e,
                            const Dune::LocalKey& lk) { gd.dofIndex(e, lk); },
                  "The PQ3 grid discretization must expose a flip-aware dofIndex(element, localKey)");

    const auto& idSet = gridDisc.gridView().grid().globalIdSet();
    std::map<std::size_t, Dune::FieldVector<Scalar, GridView::dimensionworld>> position;
    std::map<std::size_t, int> codim;
    int conflicts = 0;

    // check if DofHelper::localDofPosition provides the correct positions:
    // the basis function of a local dof is one at the position of that dof and zero at all others
    int nodalPropertyViolations = 0;
    {
        elemDisc.bind(*Dune::elements(gridDisc.gridView()).begin());
        const auto& localBasis = elemDisc.feLocalBasis();
        const auto& localCoefficients = elemDisc.feLocalCoefficients();
        const auto type = elemDisc.element().type();
        std::vector<typename std::decay_t<decltype(localBasis)>::Traits::RangeType> shapeValues;
        for (std::size_t nodeIdx = 0; nodeIdx < localCoefficients.size(); ++nodeIdx)
        {
            localBasis.evaluateFunction(DofHelper::localDofPosition(type, localCoefficients.localKey(nodeIdx)), shapeValues);
            for (std::size_t basisIdx = 0; basisIdx < localCoefficients.size(); ++basisIdx)
            {
                const Scalar expectedValue = basisIdx == nodeIdx ? 1.0 : 0.0;
                if (std::abs(shapeValues[basisIdx][0] - expectedValue) > 1e-10)
                {
                    std::cout << name << ": basis function " << basisIdx << " is " << shapeValues[basisIdx][0]
                              << " at the position of local dof " << nodeIdx << std::endl;
                    ++nodalPropertyViolations;
                }
            }
        }
    }

    for (const auto& element : Dune::elements(gridDisc.gridView()))
    {
        elemDisc.bind(element);
        const auto& localCoefficients = elemDisc.feLocalCoefficients();
        const auto geometry = element.geometry();

        for (std::size_t i = 0; i < elemDisc.numLocalDofs(); ++i)
        {
            const auto& localKey = localCoefficients.localKey(i);
            const auto dofIdx = DofHelper::dofIndex(gridDisc.dofMapper(), element, localKey, idSet);
            const auto pos = DofHelper::dofPosition(geometry, localKey);

            const auto it = position.find(dofIdx);
            if (it == position.end())
            {
                position[dofIdx] = pos;
                codim[dofIdx] = localKey.codim();
            }
            else if ((it->second - pos).two_norm() > 1e-10)
            {
                std::cout << name << ": global dof " << dofIdx << " is at " << it->second
                          << " in one element and at " << pos << " in another" << std::endl;
                ++conflicts;
            }
        }
    }

    std::size_t unvisited = 0, outOfRange = 0;
    for (std::size_t i = 0; i < gridDisc.numDofs(); ++i)
        if (!position.count(i)) ++unvisited;
    for (const auto& [idx, pos] : position)
        if (idx >= gridDisc.numDofs()) ++outOfRange;

    std::cout << name << ": dof classification  ";
    for (std::size_t i = 0; i < gridDisc.numDofs(); ++i)
        std::cout << i << ":c" << codim[i] << " ";
    std::cout << std::endl;

    std::cout << name << ": " << gridDisc.numDofs() << " dofs, "
              << position.size() << " visited, " << unvisited << " never claimed, "
              << outOfRange << " out of range, "
              << conflicts << " position conflicts" << std::endl;
    return conflicts + nodalPropertyViolations + static_cast<int>(unvisited) + static_cast<int>(outOfRange);
}

int main (int argc, char *argv[])
{
    using namespace Dumux;
    initialize(argc, argv);

    using SimplexGrid = Dune::ALUGrid<2, 2, Dune::simplex, Dune::nonconforming>;
    using CubeGrid = Dune::ALUGrid<2, 2, Dune::cube, Dune::nonconforming>;

    Dune::FieldVector<double, 2> lower(0.0), upper(1.0);
    std::array<unsigned int, 2> els{{2, 2}};

    int failed = 0;
    auto simplexGrid = Dune::StructuredGridFactory<SimplexGrid>::createSimplexGrid(lower, upper, els);
    auto cubeGrid = Dune::StructuredGridFactory<CubeGrid>::createCubeGrid(lower, upper, els);
    failed += check("simplex 2x2", *simplexGrid);
    failed += check("cube 2x2", *cubeGrid);

    if (failed)
        DUNE_THROW(Dune::Exception, "PQ3 dof map is not consistent between neighbouring elements");

    return 0;
}
