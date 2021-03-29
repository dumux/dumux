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
 * \brief Discretization-related concepts.
 */
#ifndef DUMUX_DISCRETIZATION_CONCEPTS_HH
#define DUMUX_DISCRETIZATION_CONCEPTS_HH

#include <utility>
#include <dune/common/concept.hh>

namespace Dumux::Experimental::Concept {

/*!
 * \ingroup Discretization
 * \brief Concept for (still) experimental variables of a numerical model.
 */
 struct Variables
 {
     template<class T>
     auto require(const T& t) -> decltype(
         Dune::Concept::requireType<typename T::SolutionVector>(),
         Dune::Concept::requireType<typename T::TimeLevel>(),
         Dune::Concept::requireConvertible<typename T::SolutionVector>(t.dofs()),
         Dune::Concept::requireConvertible<typename T::TimeLevel>(t.timeLevel()),
         const_cast<T&>(t).updateTime(std::declval<typename T::TimeLevel>()),
         const_cast<T&>(t).update(std::declval<typename T::SolutionVector>(),
                                  std::declval<typename T::TimeLevel>()),
         const_cast<T&>(t).update(std::declval<typename T::SolutionVector>())
     );
 };

/*!
 * \ingroup Discretization
 * \brief Concept for (still) experimental grid variables.
 */
 struct GridVariables : Dune::Concept::Refines<Variables>
 {
     template<class T>
     auto require(const T& t) -> decltype(
         Dune::Concept::requireType<typename T::GridGeometry>(),
         Dune::Concept::requireConvertible<typename T::GridGeometry>(t.gridGeometry())
     );
 };

 /*!
  * \ingroup Discretization
  * \brief Concept for (still) experimental finite volume grid variables.
  */
  struct FVGridVariables : Dune::Concept::Refines<GridVariables>
  {
      template<class T>
      auto require(const T& t) -> decltype(
          Dune::Concept::requireType<typename T::GridVolumeVariables>(),
          Dune::Concept::requireType<typename T::GridFluxVariablesCache>(),
          Dune::Concept::requireConvertible<typename T::GridVolumeVariables>(t.gridVolVars()),
          Dune::Concept::requireConvertible<typename T::GridFluxVariablesCache>(t.gridFluxVarsCache())
      );
  };

} // end namespace Dumux::Experimental::Concept

#endif
