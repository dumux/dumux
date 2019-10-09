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
 * \brief free functions for the evaluation of primary variable gradients inside elements.
 */
#ifndef DUMUX_DISCRETIZATION_EVAL_GRADIENTS_HH
#define DUMUX_DISCRETIZATION_EVAL_GRADIENTS_HH

#include <dune/localfunctions/lagrange/pqkfactory.hh>
#include <dune/common/exceptions.hh>

#include <dumux/common/typetraits/state.hh>
#include <dumux/common/typetraits/isvalid.hh>
#include <dumux/discretization/box/elementsolution.hh>
#include <dumux/discretization/cellcentered/elementsolution.hh>

#include "evalsolution.hh"

namespace Dumux {

/*!
 * \brief Evaluates the gradient of a given box element solution to a given global position.
 * \ingroup Discretization
 *
 * \param element The element
 * \param geometry The element geometry
 * \param gridGeometry The finite volume grid geometry
 * \param elemSol The primary variables at the dofs of the element
 * \param globalPos The global position
 * \param ignoreState If true, the state of primary variables is ignored
 *
 * \return Dune::FieldVector with as many entries as dimension of
 *         the PrimaryVariables object (i.e. numEq). Each entry is
 *         a GlobalCoordinate object holding the priVar gradient.
 */
template<class Element, class FVElementGeometry, class PrimaryVariables>
auto evalGradients(const Element& element,
                   const typename Element::Geometry& geometry,
                   const typename FVElementGeometry::GridGeometry& gridGeometry,
                   const BoxElementSolution<FVElementGeometry, PrimaryVariables>& elemSol,
                   const typename Element::Geometry::GlobalCoordinate& globalPos,
                   bool ignoreState = false)
{
    // determine if all states are the same at all vertices
    using HasState = decltype(isValid(Detail::hasState())(elemSol[0]));
    bool allStatesEqual = ignoreState || Detail::allStatesEqual(elemSol, HasState{});

    if (allStatesEqual)
    {
        using GlobalPosition = typename Element::Geometry::GlobalCoordinate;

        // evaluate gradients using the local finite element basis
        const auto& localBasis = gridGeometry.feCache().get(geometry.type()).localBasis();

        // evaluate the shape function gradients at the scv center
        using ShapeJacobian = typename std::decay_t< decltype(localBasis) >::Traits::JacobianType;
        const auto localPos = geometry.local(globalPos);
        std::vector< ShapeJacobian > shapeJacobian;
        localBasis.evaluateJacobian(localPos, shapeJacobian);

        // the inverse transposed of the jacobian matrix
        const auto jacInvT = geometry.jacobianInverseTransposed(localPos);

        // interpolate the gradients
        Dune::FieldVector<GlobalPosition, PrimaryVariables::dimension> result( GlobalPosition(0.0) );
        for (int i = 0; i < element.subEntities(Element::Geometry::mydimension); ++i)
        {
            // the global shape function gradient
            GlobalPosition gradN;
            jacInvT.mv(shapeJacobian[i][0], gradN);

            // add gradient to global privar gradients
            for (unsigned int pvIdx = 0; pvIdx < PrimaryVariables::dimension; ++pvIdx)
            {
                GlobalPosition tmp(gradN);
                tmp *= elemSol[i][pvIdx];
                result[pvIdx] += tmp;
            }
        }

        return result;
    }
    else
    {
        DUNE_THROW(Dune::NotImplemented, "Element vertices have different phase states. Enforce calculation by setting ignoreState to true.");
    }
}

/*!
 * \ingroup Discretization
 * \brief Evaluates the gradient of a given CCElementSolution to a given global position.
 *        This function is only here for (compilation) compatibility reasons with the box scheme.
 *        The solution within the control volumes is constant and thus gradients are zero.
 *        One can compute gradients towards the sub-control volume faces after reconstructing
 *        the solution on the faces.
 *
 * \param element The element
 * \param geometry The element geometry
 * \param gridGeometry The finite volume grid geometry
 * \param elemSol The primary variables at the dofs of the element
 * \param globalPos The global position
 * \throws Dune::NotImplemented
 *
 * \return Dune::FieldVector with as many entries as dimension of
 *         the PrimaryVariables object (i.e. numEq). Each entry is
 *         a GlobalCoordinate object holding the priVar gradient.
 */
template<class Element, class FVElementGeometry, class PrimaryVariables>
Dune::FieldVector<typename Element::Geometry::GlobalCoordinate, PrimaryVariables::dimension>
evalGradients(const Element& element,
              const typename Element::Geometry& geometry,
              const typename FVElementGeometry::GridGeometry& gridGeometry,
              const CCElementSolution<FVElementGeometry, PrimaryVariables>& elemSol,
              const typename Element::Geometry::GlobalCoordinate& globalPos)
{ DUNE_THROW(Dune::NotImplemented, "General gradient evaluation for cell-centered methods"); }

} // namespace Dumux

#endif
