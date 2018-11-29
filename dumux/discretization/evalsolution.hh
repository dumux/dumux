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
 * \ingroup Discretization
 * \brief free functions for the evaluation of primary variables inside elements.
 */
#ifndef DUMUX_DISCRETIZATION_EVAL_SOLUTION_HH
#define DUMUX_DISCRETIZATION_EVAL_SOLUTION_HH

#include <iterator>
#include <algorithm>

#include <dune/localfunctions/lagrange/pqkfactory.hh>

#include <dumux/common/typetraits/state.hh>
#include <dumux/discretization/box/elementsolution.hh>
#include <dumux/discretization/cellcentered/elementsolution.hh>

namespace Dumux {

// The following function declarations are required because of the default
// parameters together with the friend declarations further below.
template<class Element, class FVElementGeometry, class PrimaryVariables>
PrimaryVariables evalSolution(const Element& element,
                              const typename Element::Geometry& geometry,
                              const typename FVElementGeometry::FVGridGeometry& fvGridGeometry,
                              const BoxElementSolution<FVElementGeometry, PrimaryVariables>& elemSol,
                              const typename Element::Geometry::GlobalCoordinate& globalPos,
                              bool ignoreState = false);

template<class Element, class FVElementGeometry, class PrimaryVariables>
PrimaryVariables evalSolution(const Element& element,
                              const typename Element::Geometry& geometry,
                              const BoxElementSolution<FVElementGeometry, PrimaryVariables>& elemSol,
                              const typename Element::Geometry::GlobalCoordinate& globalPos,
                              bool ignoreState = false);

class SolutionEvaluator
{
    // declare the free functions `evalSolution` as friends
    template<class Element, class FVElementGeometry, class PrimaryVariables>
    friend PrimaryVariables evalSolution(const Element&,
                                         const typename Element::Geometry&,
                                         const typename FVElementGeometry::FVGridGeometry&,
                                         const BoxElementSolution<FVElementGeometry, PrimaryVariables>&,
                                         const typename Element::Geometry::GlobalCoordinate&,
                                         bool ignoreState);

    template<class Element, class FVElementGeometry, class PrimaryVariables>
    friend PrimaryVariables evalSolution(const Element&,
                                         const typename Element::Geometry&,
                                         const BoxElementSolution<FVElementGeometry, PrimaryVariables>&,
                                         const typename Element::Geometry::GlobalCoordinate&,
                                         bool ignoreState);

    /*!
     * \brief Interpolates a given box element solution at a given global position.
     *        Uses the finite element cache of the grid geometry.
     *
     * Specialization for primary variables without state. Always does an
     * interpolation using the shape functions.
     *
     * \return the interpolated primary variables
     * \param element The element
     * \param geometry The element geometry
     * \param fvGridGeometry The finite volume grid geometry
     * \param elemSol The primary variables at the dofs of the element
     * \param globalPos The global position
     * \param ignoreState If true, the state of primary variables is ignored
     */
    template<class Element, class FVElementGeometry, class PrimaryVariables>
    static auto eval_(const Element& element,
                      const typename Element::Geometry& geometry,
                      const typename FVElementGeometry::FVGridGeometry& fvGridGeometry,
                      const BoxElementSolution<FVElementGeometry, PrimaryVariables>& elemSol,
                      const typename Element::Geometry::GlobalCoordinate& globalPos,
                      bool ignoreState)
    -> typename std::enable_if_t<!decltype(isValid(Detail::hasState())(elemSol[0]))::value, PrimaryVariables>
    {
        return evalIgnoringState_(element, geometry, fvGridGeometry, elemSol, globalPos);
    }

    /*!
     * \brief Interpolates a given box element solution at a given global position.
     *        Uses the finite element cache of the grid geometry.
     *
     * Specialization for primary variables with state. Does an interpolation using
     * the shape functions if states at all vertices are equal. Otherwise, it takes
     * the state and value from the vertex closest to the evaluation point.
     *
     * \return the interpolated primary variables
     * \param element The element
     * \param geometry The element geometry
     * \param fvGridGeometry The finite volume grid geometry
     * \param elemSol The primary variables at the dofs of the element
     * \param globalPos The global position
     * \param ignoreState If true, the state of primary variables is ignored
     */
    template<class Element, class FVElementGeometry, class PrimaryVariables>
    static auto eval_(const Element& element,
                      const typename Element::Geometry& geometry,
                      const typename FVElementGeometry::FVGridGeometry& fvGridGeometry,
                      const BoxElementSolution<FVElementGeometry, PrimaryVariables>& elemSol,
                      const typename Element::Geometry::GlobalCoordinate& globalPos,
                      bool ignoreState)
    -> typename std::enable_if_t<decltype(isValid(Detail::hasState())(elemSol[0]))::value, PrimaryVariables>
    {
        // determine if all states are the same at all vertices
        bool allStatesEqual = true;
        if (!ignoreState)
        {
            auto firstState = elemSol[0].state();
            for (int i = 1; i < geometry.corners(); ++i)
            {
                if (elemSol[i].state() != firstState)
                {
                    allStatesEqual = false;
                    break;
                }
            }
        }

        if (allStatesEqual)
        {
            auto result = evalIgnoringState_(element, geometry, fvGridGeometry, elemSol, globalPos);

            if (!ignoreState)
                result.setState(elemSol[0].state());

            return result;
        }
        else
        {
            if (!warnedAboutUsingMinDist)
            {
                std::cout << "Warning: Using nearest-neighbor interpolation in evalSolution"
                          << "\nbecause not all states are equal and ignoreState is false!"
                          << std::endl;
                warnedAboutUsingMinDist = true;
            }

            auto minDistanceIdx = indexOfMinDistVertex_(geometry, globalPos);

            return elemSol[minDistanceIdx];
        }
    }

    /*!
     * \brief Interpolates a given box element solution at a given global position.
     *
     * Overload of the above evalSolution() function without a given fvGridGeometry.
     * The local basis is computed on the fly. Specialization for primary variables
     * without state. Always does an interpolation using the shape functions.
     *
     * \return the interpolated primary variables
     * \param element The element
     * \param geometry The element geometry
     * \param elemSol The primary variables at the dofs of the element
     * \param globalPos The global position
     * \param ignoreState If true, the state of primary variables is ignored
     */
    template<class Element, class FVElementGeometry, class PrimaryVariables>
    static auto eval_(const Element& element,
                      const typename Element::Geometry& geometry,
                      const BoxElementSolution<FVElementGeometry, PrimaryVariables>& elemSol,
                      const typename Element::Geometry::GlobalCoordinate& globalPos,
                      bool ignoreState)
    -> typename std::enable_if_t<!decltype(isValid(Detail::hasState())(elemSol[0]))::value, PrimaryVariables>
    {
        return evalIgnoringState_(element, geometry, elemSol, globalPos);
    }

    /*!
     * \brief Interpolates a given box element solution at a given global position.
     *
     * Overload of the above evalSolution() function without a given fvGridGeometry.
     * The local basis is computed on the fly. Specialization for primary variables
     * with state. Does an interpolation using the shape functions if states at all
     * vertices are equal. Otherwise, it takes the state and value from the vertex
     * closest to the evaluation point.
     *
     * \return the interpolated primary variables
     * \param element The element
     * \param geometry The element geometry
     * \param elemSol The primary variables at the dofs of the element
     * \param globalPos The global position
     * \param ignoreState If true, the state of primary variables is ignored
     */
    template<class Element, class FVElementGeometry, class PrimaryVariables>
    static auto eval_(const Element& element,
                      const typename Element::Geometry& geometry,
                      const BoxElementSolution<FVElementGeometry, PrimaryVariables>& elemSol,
                      const typename Element::Geometry::GlobalCoordinate& globalPos,
                      bool ignoreState)
    -> typename std::enable_if_t<decltype(isValid(Detail::hasState())(elemSol[0]))::value, PrimaryVariables>
    {
        // determine if all states are the same at all vertices
        bool allStatesEqual = true;
        if (!ignoreState)
        {
            auto firstState = elemSol[0].state();
            for (int i = 1; i < geometry.corners(); ++i)
            {
                if (elemSol[i].state() != firstState)
                {
                    allStatesEqual = false;
                    break;
                }
            }
        }

        if (allStatesEqual)
        {
            auto result = evalIgnoringState_(element, geometry, elemSol, globalPos);

            if (!ignoreState)
                result.setState(elemSol[0].state());

            return result;
        }
        else
        {
            if (!warnedAboutUsingMinDist)
            {
                std::cout << "Warning: Using nearest-neighbor interpolation in evalSolution"
                          << "\nbecause not all states are equal and ignoreState is false!"
                          << std::endl;
                warnedAboutUsingMinDist = true;
            }

            auto minDistanceIdx = indexOfMinDistVertex_(geometry, globalPos);

            return elemSol[minDistanceIdx];
        }
    }

    template<class Element, class FVElementGeometry, class PrimaryVariables>
    static PrimaryVariables evalIgnoringState_(const Element& element,
                                               const typename Element::Geometry& geometry,
                                               const typename FVElementGeometry::FVGridGeometry& fvGridGeometry,
                                               const BoxElementSolution<FVElementGeometry, PrimaryVariables>& elemSol,
                                               const typename Element::Geometry::GlobalCoordinate& globalPos)
    {
        using Scalar = typename PrimaryVariables::value_type;

        // interpolate the solution
        const auto& localBasis = fvGridGeometry.feCache().get(geometry.type()).localBasis();

        // evaluate the shape functions at the scv center
        const auto localPos = geometry.local(globalPos);
        std::vector< Dune::FieldVector<Scalar, 1> > shapeValues;
        localBasis.evaluateFunction(localPos, shapeValues);

        PrimaryVariables result(0.0);
        for (int i = 0; i < geometry.corners(); ++i)
        {
            auto value = elemSol[i];
            value *= shapeValues[i][0];
            result += value;
        }

        return result;
    }

    template<class Element, class FVElementGeometry, class PrimaryVariables>
    static PrimaryVariables evalIgnoringState_(const Element& element,
                                               const typename Element::Geometry& geometry,
                                               const BoxElementSolution<FVElementGeometry, PrimaryVariables>& elemSol,
                                               const typename Element::Geometry::GlobalCoordinate& globalPos)
    {
        using Scalar = typename PrimaryVariables::value_type;
        using CoordScalar = typename Element::Geometry::GlobalCoordinate::value_type;
        static constexpr int dim = Element::Geometry::mydimension;

        //! The box scheme always uses linear Ansatz functions
        using FeCache = Dune::PQkLocalFiniteElementCache<CoordScalar, Scalar, dim, 1>;
        using ShapeValue = typename FeCache::FiniteElementType::Traits::LocalBasisType::Traits::RangeType;

        // obtain local finite element basis
        FeCache feCache;
        const auto& localBasis = feCache.get(geometry.type()).localBasis();

        // evaluate the shape functions at the scv center
        const auto localPos = geometry.local(globalPos);
        std::vector< ShapeValue > shapeValues;
        localBasis.evaluateFunction(localPos, shapeValues);

        PrimaryVariables result(0.0);
        for (int i = 0; i < geometry.corners(); ++i)
        {
            auto value = elemSol[i];
            value *= shapeValues[i][0];
            result += value;
        }

        return result;
    }

    template<class Geometry>
    static auto indexOfMinDistVertex_(const Geometry& geometry,
                                      const typename Geometry::GlobalCoordinate& globalPos)
    {
        // calculate the distances from the evaluation point to the vertices
        std::vector<typename Geometry::ctype> distances(geometry.corners());
        for (int i = 0; i < geometry.corners(); ++i)
            distances[i] = (geometry.corner(i) - globalPos).two_norm2();

        // determine the minimum distance and corresponding vertex
        auto minDistanceIt = std::min_element(distances.begin(), distances.end());
        return std::distance(distances.begin(), minDistanceIt);
    }

public:
    static bool warnedAboutUsingMinDist;
};

bool SolutionEvaluator::warnedAboutUsingMinDist = false;

/*!
 * \brief Interpolates a given box element solution at a given global position.
 *        Uses the finite element cache of the grid geometry.
 * \ingroup Discretization
 *
 * \return the interpolated primary variables
 * \param element The element
 * \param geometry The element geometry
 * \param fvGridGeometry The finite volume grid geometry
 * \param elemSol The primary variables at the dofs of the element
 * \param globalPos The global position
 * \param ignoreState If true, the state of primary variables is ignored
 */
template<class Element, class FVElementGeometry, class PrimaryVariables>
PrimaryVariables evalSolution(const Element& element,
                              const typename Element::Geometry& geometry,
                              const typename FVElementGeometry::FVGridGeometry& fvGridGeometry,
                              const BoxElementSolution<FVElementGeometry, PrimaryVariables>& elemSol,
                              const typename Element::Geometry::GlobalCoordinate& globalPos,
                              bool ignoreState)
{
    return SolutionEvaluator::eval_(element, geometry, fvGridGeometry, elemSol, globalPos, ignoreState);
}

/*!
 * \ingroup Discretization
 * \brief Interpolates a given box element solution at a given global position.
 *
 * Overload of the above evalSolution() function without a given fvGridGeometry.
 * The local basis is computed on the fly.
 *
 * \return the interpolated primary variables
 * \param element The element
 * \param geometry The element geometry
 * \param elemSol The primary variables at the dofs of the element
 * \param globalPos The global position
 * \param ignoreState If true, the state of primary variables is ignored
 */
template<class Element, class FVElementGeometry, class PrimaryVariables>
PrimaryVariables evalSolution(const Element& element,
                              const typename Element::Geometry& geometry,
                              const BoxElementSolution<FVElementGeometry, PrimaryVariables>& elemSol,
                              const typename Element::Geometry::GlobalCoordinate& globalPos,
                              bool ignoreState)
{
    return SolutionEvaluator::eval_(element, geometry, elemSol, globalPos, ignoreState);
}

/*!
 * \brief Interpolates a given cell-centered element solution at a given global position.
 * \ingroup Discretization
 *
 * \return the primary variables (constant over the element)
 * \param element The element
 * \param geometry The element geometry
 * \param fvGridGeometry The finite volume grid geometry
 * \param elemSol The primary variables at the dofs of the element
 * \param globalPos The global position
 * \param ignoreState If true, the state of primary variables is ignored
 */
template<class Element, class FVElementGeometry, class PrimaryVariables>
PrimaryVariables evalSolution(const Element& element,
                              const typename Element::Geometry& geometry,
                              const typename FVElementGeometry::FVGridGeometry& fvGridGeometry,
                              const CCElementSolution<FVElementGeometry, PrimaryVariables>& elemSol,
                              const typename Element::Geometry::GlobalCoordinate& globalPos,
                              bool ignoreState = false)
{
    return elemSol[0];
}

/*!
 * \brief Interpolates a given cell-centered element solution at a given global position.
 *        Overload of the above evalSolution() function without a given fvGridGeometry.
 *        For compatibility reasons with the box scheme.
 * \ingroup Discretization
 *
 * \return the primary variables (constant over the element)
 * \param element The element
 * \param geometry The element geometry
 * \param elemSol The primary variables at the dofs of the element
 * \param globalPos The global position
 * \param ignoreState If true, the state of primary variables is ignored
 */
template< class Element, class FVElementGeometry, class PrimaryVariables>
PrimaryVariables evalSolution(const Element& element,
                              const typename Element::Geometry& geometry,
                              const CCElementSolution<FVElementGeometry, PrimaryVariables>& elemSol,
                              const typename Element::Geometry::GlobalCoordinate& globalPos,
                              bool ignoreState = false)
{
    return elemSol[0];
}

} // namespace Dumux

#endif
