// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup Discretization
 * \brief free functions for the evaluation of primary variables inside elements.
 */
#ifndef DUMUX_DISCRETIZATION_EVAL_SOLUTION_HH
#define DUMUX_DISCRETIZATION_EVAL_SOLUTION_HH

#include <iterator>
#include <algorithm>
#include <type_traits>
#include <memory>

#include <dune/localfunctions/lagrange/pqkfactory.hh>

#include <dumux/common/typetraits/state.hh>
#include <dumux/common/typetraits/isvalid.hh>
#include <dumux/discretization/method.hh>
#include <dumux/discretization/cvfe/elementsolution.hh>
#include <dumux/discretization/cellcentered/elementsolution.hh>

namespace Dumux {

// some implementation details
namespace Detail {

//! returns true if all states in an element solution are the same
template<class ElementSolution>
bool allStatesEqual(const ElementSolution& elemSol, std::true_type hasState)
{
    // determine if all states are the same at all vertices
    auto firstState = elemSol[0].state();
    for (std::size_t i = 1; i < elemSol.size(); ++i)
        if (elemSol[i].state() != firstState)
            return false;

    return true;
}

//! overload if the solution is stateless
template<class ElementSolution>
bool allStatesEqual(const ElementSolution& elemSol, std::false_type hasState)
{
    return true;
}

//! return the solution at the closest dof
template<class Geometry, class ElementSolution>
auto minDistVertexSol(const Geometry& geometry, const typename Geometry::GlobalCoordinate& globalPos,
                      const ElementSolution& elemSol)
{
    // calculate the distances from the evaluation point to the vertices
    std::vector<typename Geometry::ctype> distances(geometry.corners());
    for (int i = 0; i < geometry.corners(); ++i)
        distances[i] = (geometry.corner(i) - globalPos).two_norm2();

    // determine the minimum distance and corresponding vertex
    auto minDistanceIt = std::min_element(distances.begin(), distances.end());
    return elemSol[std::distance(distances.begin(), minDistanceIt)];
}

/*!
 * \brief Interpolates a given control-volume finite element solution at a given global position.
 *        Uses the finite element cache of the grid geometry.
 * \ingroup Discretization
 *
 * \return the interpolated primary variables
 * \param element The element
 * \param geometry The element geometry
 * \param gridGeometry The finite volume grid geometry
 * \param elemSol The primary variables at the dofs of the element
 * \param globalPos The global position
 * \param ignoreState If true, the state of primary variables is ignored
 */
template<class Element, class GridGeometry, class CVFEElemSol>
typename CVFEElemSol::PrimaryVariables
evalCVFESolution(const Element& element,
                 const typename Element::Geometry& geometry,
                 const GridGeometry& gridGeometry,
                 const CVFEElemSol& elemSol,
                 const typename Element::Geometry::GlobalCoordinate& globalPos,
                 bool ignoreState = false)
{
    // determine if all states are the same at all vertices
    using HasState = decltype(isValid(Detail::hasState())(elemSol[0]));
    bool allStatesEqual = ignoreState || Detail::allStatesEqual(elemSol, HasState{});

    if (allStatesEqual)
    {
        using Scalar = typename CVFEElemSol::PrimaryVariables::value_type;

        // interpolate the solution
        const auto& localBasis = gridGeometry.feCache().get(geometry.type()).localBasis();

        // evaluate the shape functions at the scv center
        const auto localPos = geometry.local(globalPos);
        std::vector< Dune::FieldVector<Scalar, 1> > shapeValues;
        localBasis.evaluateFunction(localPos, shapeValues);

        typename CVFEElemSol::PrimaryVariables result(0.0);
        for (int i = 0; i < shapeValues.size(); ++i)
        {
            auto value = elemSol[i];
            value *= shapeValues[i][0];
            result += value;
        }

        // set an arbitrary state if the model requires a state (models constexpr if)
        if constexpr (HasState{})
            if (!ignoreState)
                result.setState(elemSol[0].state());

        return result;
    }
    else
    {
        static bool warnedAboutUsingMinDist = false;
        if (!warnedAboutUsingMinDist)
        {
            std::cout << "Warning: Using nearest-neighbor interpolation in evalSolution"
            << "\nbecause not all states are equal and ignoreState is false!"
            << std::endl;
            warnedAboutUsingMinDist = true;
        }

        return Detail::minDistVertexSol(geometry, globalPos, elemSol);
    }
}

} // end namespace Detail

/*!
 * \brief Interpolates a given box element solution at a given global position.
 *        Uses the finite element cache of the grid geometry.
 * \ingroup Discretization
 *
 * \return the interpolated primary variables
 * \param element The element
 * \param geometry The element geometry
 * \param gridGeometry The finite volume grid geometry
 * \param elemSol The primary variables at the dofs of the element
 * \param globalPos The global position
 * \param ignoreState If true, the state of primary variables is ignored
 */
template<class Element, class FVElementGeometry, class PrimaryVariables>
PrimaryVariables evalSolution(const Element& element,
                              const typename Element::Geometry& geometry,
                              const typename FVElementGeometry::GridGeometry& gridGeometry,
                              const CVFEElementSolution<FVElementGeometry, PrimaryVariables>& elemSol,
                              const typename Element::Geometry::GlobalCoordinate& globalPos,
                              bool ignoreState = false)
{
    return Detail::evalCVFESolution(element, geometry, gridGeometry, elemSol, globalPos, ignoreState);
}

/*!
 * \ingroup Discretization
 * \brief Interpolates a given box element solution at a given global position.
 *
 * Overload of the above evalSolution() function without a given gridGeometry.
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
                              const CVFEElementSolution<FVElementGeometry, PrimaryVariables>& elemSol,
                              const typename Element::Geometry::GlobalCoordinate& globalPos,
                              bool ignoreState = false)
{
    static_assert(FVElementGeometry::GridGeometry::discMethod == DiscretizationMethods::box, "Only implemented for the Box method");

    // determine if all states are the same at all vertices
    using HasState = decltype(isValid(Detail::hasState())(elemSol[0]));
    bool allStatesEqual = ignoreState || Detail::allStatesEqual(elemSol, HasState{});

    if (allStatesEqual)
    {
        using Scalar = typename PrimaryVariables::value_type;
        using CoordScalar = typename Element::Geometry::GlobalCoordinate::value_type;
        static constexpr int dim = Element::Geometry::mydimension;

        // The box scheme always uses linear Ansatz functions
        using FEFactory = Dune::PQkLocalFiniteElementFactory<CoordScalar, Scalar, dim, 1>;
        using FiniteElement = typename FEFactory::FiniteElementType;
        std::unique_ptr<FiniteElement> fe(FEFactory::create(geometry.type()));
        const auto& localBasis = fe->localBasis();

        // evaluate the shape functions at the scv center
        const auto localPos = geometry.local(globalPos);
        using ShapeValue = typename FiniteElement::Traits::LocalBasisType::Traits::RangeType;
        std::vector< ShapeValue > shapeValues;
        localBasis.evaluateFunction(localPos, shapeValues);

        PrimaryVariables result(0.0);
        for (int i = 0; i < geometry.corners(); ++i)
        {
            auto value = elemSol[i];
            value *= shapeValues[i][0];
            result += value;
        }

        // set an arbitrary state if the model requires a state (models constexpr if)
        if constexpr (HasState{})
            if (!ignoreState)
                result.setState(elemSol[0].state());

        return result;
    }
    else
    {
        static bool warnedAboutUsingMinDist = false;
        if (!warnedAboutUsingMinDist)
        {
            std::cout << "Warning: Using nearest-neighbor interpolation in evalSolution"
            << "\nbecause not all states are equal and ignoreState is false!"
            << std::endl;
            warnedAboutUsingMinDist = true;
        }

        return Detail::minDistVertexSol(geometry, globalPos, elemSol);
    }
}

/*!
 * \ingroup Discretization
 * \brief Interpolates a given cell-centered element solution at a given global position.
 *
 * \return the primary variables (constant over the element)
 * \param element The element
 * \param geometry The element geometry
 * \param gridGeometry The finite volume grid geometry
 * \param elemSol The primary variables at the dofs of the element
 * \param globalPos The global position
 * \param ignoreState If true, the state of primary variables is ignored
 */
template<class Element, class FVElementGeometry, class PrimaryVariables>
PrimaryVariables evalSolution(const Element& element,
                              const typename Element::Geometry& geometry,
                              const typename FVElementGeometry::GridGeometry& gridGeometry,
                              const CCElementSolution<FVElementGeometry, PrimaryVariables>& elemSol,
                              const typename Element::Geometry::GlobalCoordinate& globalPos,
                              bool ignoreState = false)
{
    return elemSol[0];
}

/*!
 * \brief Interpolates a given cell-centered element solution at a given global position.
 *        Overload of the above evalSolution() function without a given gridGeometry.
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
