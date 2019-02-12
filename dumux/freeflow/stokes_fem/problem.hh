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
 * \brief Base class for all stokes problems which use the box scheme.
 */
#ifndef DUMUX_STOKES_PROBLEM_HH
#define DUMUX_STOKES_PROBLEM_HH

#include <dumux/implicit/problem.hh>

#include "properties.hh"

namespace Dumux
{
/*!
 * \ingroup ImplicitBaseProblems
 * \ingroup BoxStokesModel
 * \brief Base class for all problems which use the Stokes box model.
 *
 * This implements gravity (if desired) and a function returning the temperature.
 */
template<class TypeTag>
class StokesProblem : public ImplicitProblem<TypeTag>
{
    using ParentType = ImplicitProblem<TypeTag>;
    using Implementation = typename GET_PROP_TYPE(TypeTag, Problem);

    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
    using TimeManager = typename GET_PROP_TYPE(TypeTag, TimeManager);
    using Grid = typename GridView::Grid;
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);

    using Element = typename GridView::template Codim<0>::Entity;
    using Intersection = typename GridView::Intersection;


    enum {
        dim = Grid::dimension,
        dimWorld = Grid::dimensionworld
    };

    using GlobalPosition = Dune::FieldVector<Scalar, dimWorld>;

public:
    StokesProblem(TimeManager &timeManager, const GridView &gridView)
        : ParentType(timeManager, gridView),
          gravity_(0)
    {
        if (GET_PARAM_FROM_GROUP(TypeTag, bool, Problem, EnableGravity))
            gravity_[dim-1]  = -9.81;
    }


    /*!
     * \brief Returns the temperature within the domain.
     *
     * This method MUST be overwritten by the actual problem.
     */
    Scalar temperature() const
    { DUNE_THROW(Dune::NotImplemented, "temperature() method not implemented by the actual problem"); }

    /*!
     * \brief Returns the acceleration due to gravity.
     *
     * If the <tt>EnableGravity</tt> property is true, this means
     * \f$\boldsymbol{g} = ( 0,\dots,\ -9.81)^T \f$, else \f$\boldsymbol{g} = ( 0,\dots, 0)^T \f$
     */
    const GlobalPosition &gravity() const
    { return gravity_; }

    // \}

private:
    //! Returns the implementation of the problem (i.e. static polymorphism)
    Implementation &asImp_()
    { return *static_cast<Implementation *>(this); }

    //! \copydoc asImp_()
    const Implementation &asImp_() const
    { return *static_cast<const Implementation *>(this); }

    GlobalPosition gravity_;
};

}

#endif
