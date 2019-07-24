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
 * \ingroup FEDiscretization
 * \brief The local element solution class for the finite element schemes
 */
#ifndef DUMUX_FEM_ELEMENT_SOLUTION_HH
#define DUMUX_FEM_ELEMENT_SOLUTION_HH

#include <type_traits>
#include <dune/istl/bvector.hh>
#include <dumux/discretization/method.hh>

namespace Dumux {

/*!
 * \ingroup FEDiscretization
 * \brief The element solution vector
 */
template<class FEElementGeometry, class PV>
class FEElementSolution
{
    using GridGeometry = typename FEElementGeometry::GridGeometry;
    using GridView = typename GridGeometry::GridView;
    using Element = typename GridView::template Codim<0>::Entity;

public:
    //! export the primary variables type
    using PrimaryVariables = PV;

    //! Default constructor
    FEElementSolution() = default;

    //! Constructor with element, solution and function space basis
    template<class SolutionVector, class FunctionSpaceBasis>
    FEElementSolution(const Element& element,
                      const SolutionVector& sol,
                      const FunctionSpaceBasis& feBasis)
    {
        update(element, sol, feBasis);
    }

    //! Constructor with element, solution and grid geometry
    template<class SolutionVector>
    FEElementSolution(const Element& element,
                      const SolutionVector& sol,
                      const GridGeometry& gridGeometry)
    {
        update(element, sol, gridGeometry);
    }

    //! extract element solution from element, solution and function space basis
    template<class SolutionVector, class FunctionSpaceBasis>
    void update(const Element& element,
                const SolutionVector& sol,
                const FunctionSpaceBasis& feBasis)
    {
        auto localView = feBasis.localView();
        localView.bind(element);

        priVars_.resize(localView.size());
        for (int i = 0; i < localView.size(); i++)
            priVars_[i] = sol[localView.index(i)];
    }

    //! extract element solution from element, solution and grid geometry
    template<class SolutionVector>
    void update(const Element& element,
                const SolutionVector& sol,
                const GridGeometry& gridGeometry)
    {
        update(element, sol, gridGeometry.feBasis());
    }

    //! return the size of the element solution
    std::size_t size() const
    { return priVars_.size(); }

    //! bracket operator const access
    template<typename IndexType>
    const PrimaryVariables& operator [](IndexType i) const
    { return priVars_[i]; }

    //! bracket operator access
    template<typename IndexType>
    PrimaryVariables& operator [](IndexType i)
    { return priVars_[i]; }

private:
    Dune::BlockVector<PrimaryVariables> priVars_;
};

/*!
 * \ingroup FemDiscretization
 * \brief  Make an element solution for finite element schemes
 */
template<class Element, class SolutionVector, class GridGeometry>
auto elementSolution(const Element& element, const SolutionVector& sol, const GridGeometry& gg)
-> std::enable_if_t< GridGeometry::discMethod == DiscretizationMethod::fem,
                     FEElementSolution<typename GridGeometry::LocalView, std::decay_t<decltype(std::declval<SolutionVector>()[0])>> >
{
    using PrimaryVariables = std::decay_t<decltype(std::declval<SolutionVector>()[0])>;
    return FEElementSolution<typename GridGeometry::LocalView, PrimaryVariables>(element, sol, gg);
}

} // end namespace Dumux

#endif
