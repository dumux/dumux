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
 * \ingroup Discretization
 * \brief The grid variable class for finite element schemes.
 * \note This default implementation does not store any additional data on the grid.
 */
#ifndef DUMUX_FE_GRID_VARIABLES_HH
#define DUMUX_FE_GRID_VARIABLES_HH

#include <cassert>
#include <dumux/discretization/localview.hh>

#include "ipvariablesbase.hh"

namespace Dumux {

/*!
 * \ingroup Discretization
 * \brief Local view on the grid variables for finite element schemes.
 * \note This default implementation only stores a pointer on the grid variables.
 * \tparam GV The grid variables class
 */
template<class GV>
class FEGridVariablesLocalView
{
    using GridGeometry = typename GV::GridGeometry;
    using FEElementGeometry = typename GridGeometry::LocalView;

    using GridView = typename GridGeometry::GridView;
    using Element = typename GridView::template Codim<0>::Entity;

public:
    using GridVariables = GV;

    //! The constructor
    FEGridVariablesLocalView(const GridVariables& gridVariables)
    : gridVariables_(&gridVariables)
    {}

    /*!
     * \brief Bind this local view to a grid element.
     * \param element The grid element
     * \param feGeometry Local view on the grid geometry
     * \todo TODO: In fv schemes here the current solution enters
     *             (in the vol vars/flux vars cache) in order to bind
     *             the local views to a potentially deflected solution
     *             Do we require the grid variables to carry the notion
     *             of the deflected solution or should this be injected here!?
     */
    void bind(const Element& element,
              const FEElementGeometry& feGeometry)
    {}

    /*!
     * \brief Return a reference to the grid variables.
     */
    const GridVariables& gridVariables() const
    { return *gridVariables_; }

private:
    const GridVariables* gridVariables_;
};

//! TODO: have these defined at a central place
template<class Problem>
struct ProblemTraits;

/*!
 * \ingroup Discretization
 * \brief The grid variable class for finite element schemes.
 * \note This default implementation does not store any additional data on the grid.
 * \tparam P the problem to be solved
 * \tparam X the type used to represent solution vectors on the grid
 * \tparam IPV the variable type to be constructed at integration points
 */
template<class P, class X, class IPV = IntegrationPointVariablesBase<typename X::value_type>>
class FEGridVariables
{
    using ThisType = FEGridVariables<P, X, IPV>;

public:
    //! export the underlying problem
    using Problem = P;

    //! export the type used for solution vectors
    using SolutionVector = X;

    //! export type of the finite volume grid geometry
    using GridGeometry = typename ProblemTraits<Problem>::GridGeometry;

    //! export primary variable type
    using PrimaryVariables = typename SolutionVector::value_type;

    //! export scalar type
    using Scalar = typename PrimaryVariables::value_type;

    //! export type of variables constructed at integration points
    using IntegrationPointVariables = IPV;

    //! The local view on these grid variables
    using LocalView = FEGridVariablesLocalView<ThisType>;

    //! constructor
    FEGridVariables(std::shared_ptr<const Problem> problem)
    : problem_(problem)
    {}

    //! return the underlying problem
    const Problem& problem() const
    { return *problem_; }

    //! return the finite volume grid geometry
    const GridGeometry& gridGeometry() const
    { return problem_->gridGeometry(); }

    //! return the solution for which the grid variables were updated
    const SolutionVector& solutionVector() const
    { assert(x_ && "No solution vector set!"); return *x_; }

    ///////////////
    // TODO: The following interfaces are defined in FVGridVariables.
    //       DO WE NEED ALL OF THEM?
    //////////////

    //! initialize all variables
    void init(const SolutionVector& x)
    {
        x_ = &x;
    }

    //! update all variables
    void update(const SolutionVector& x, bool forceUpdate = false)
    {
        x_ = &x;
    }

    //! update all variables after grid adaption
    void updateAfterGridAdaption(const SolutionVector& x)
    {
        x_ = &x;
    }

    //! Sets the current state as the previous for next time step
    void advanceTimeStep()
    {}

    //! resets state to the one before time integration
    void resetTimeStep(const SolutionVector& x)
    {
        x_ = &x;
    }

private:
    std::shared_ptr<const Problem> problem_; //!< pointer to the problem to be solved
    const SolutionVector* x_ = nullptr; //!< the solution corresponding to this class' current state
};

} // end namespace Dumux

#endif
