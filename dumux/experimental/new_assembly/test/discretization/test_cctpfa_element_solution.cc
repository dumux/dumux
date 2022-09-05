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
 * \brief Tests for the `elementSolution` function for TPFA grid geometries.
 */
#include <iostream>
#include <concepts>
#include <ranges>
#include <cmath>

#include <dune/common/exceptions.hh>
#include <dune/grid/yaspgrid.hh>

#include <dumux/experimental/new_assembly/dumux/discretization/cctpfa/gridgeometry.hh>
#include <dumux/experimental/new_assembly/dumux/discretization/elementsolution.hh>
#include <dumux/experimental/new_assembly/dumux/linear/systemtraits.hh>
#include <dumux/experimental/new_assembly/dumux/common/indexstrategies.hh>
#include <dumux/experimental/new_assembly/dumux/common/concepts.hh>

struct SolutionAtDof
{
    double priVar0 = 0.0;
    double priVar1 = 0.0;
};

struct SolutionVector
{
public:
    explicit SolutionVector(std::size_t size)
    : x_(size)
    {
        std::ranges::for_each(std::views::iota(std::size_t{0}, size), [&] (auto idx) {
            x_[idx].priVar0 = static_cast<double>(idx);
            x_[idx].priVar1 = static_cast<double>(idx) + 0.5;
        });
    }

    const SolutionAtDof& operator[](std::size_t i) const
    { return x_[i]; }

private:
    std::vector<SolutionAtDof> x_;
};

namespace Dumux::LinearSystem::Traits {

template<>
struct Scalar<SolutionVector>
: public std::type_identity<double>
{};

template<std::integral I0, std::integral I1>
struct VectorAccess<SolutionVector, Dumux::MultiIndex<I0, I1>>
{
    using MultiIndex = Dumux::MultiIndex<I0, I1>;
    static double get(const SolutionVector& x, const MultiIndex& i)
    {
        if (i.template get<1>() > 1)
            DUNE_THROW(Dune::InvalidStateException, "Unexpected privar index");
        return i.template get<1>() == 0 ? x[i.template get<0>()].priVar0
                                        : x[i.template get<0>()].priVar1;
    }
};

} // namespace Dumux::LinearSystem::Traits


int main (int argc, char *argv[])
{
    const unsigned int numCells = 10;
    Dune::YaspGrid<2> grid{
        {1.0, 1.0},
        {numCells, numCells}
    };

    Dumux::CCTpfaGridGeometry gridGeometry{grid.leafGridView()};
    Dumux::BlockedIndexStrategy<2> indexStrategy{};
    std::integral_constant<int, 2> numEq;
    SolutionVector dofs(gridGeometry.numDofs());

    for (const auto& element : elements(gridGeometry.gridView()))
    {
        const auto eIdx = gridGeometry.elementMapper().index(element);
        const auto localGeometry = localView(gridGeometry).bind(element);
        const auto elemSol = Dumux::elementSolution(localGeometry,
                                                    indexStrategy,
                                                    dofs,
                                                    numEq);
        if (elemSol.size() != 1)
            DUNE_THROW(Dune::InvalidStateException, "Unexpected size");
        if (std::abs(elemSol[0][0] - dofs[eIdx].priVar0) > 1e-7)
            DUNE_THROW(Dune::InvalidStateException, "Unexpected first privar");
        if (std::abs(elemSol[0][1] - dofs[eIdx].priVar1) > 1e-7)
            DUNE_THROW(Dune::InvalidStateException, "Unexpected second privar");
    }

    std::cout << "\nAll tests passed" << std::endl;
    return 0;
}
