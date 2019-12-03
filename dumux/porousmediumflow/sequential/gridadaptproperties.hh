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
 * \ingroup IMPETProperties
 * \ingroup IMPET
 * \file
 *
 * \brief Defines a type tag and some fundamental properties for
 *        linear solvers
 */
#ifndef DUMUX_GRIDADAPT_PROPERTIES_HH
#define DUMUX_GRIDADAPT_PROPERTIES_HH

#include <dumux/common/properties.hh>
#include <dumux/common/parameters.hh>

namespace Dumux
{
// forward declarations
template<class TypeTag>
class GridAdaptInitializationIndicatorDefault;

template<class TypeTag, bool>
class GridAdapt;


namespace Properties
{
//! Grid adaption type tag for all sequential models.
namespace TTag {
struct GridAdapt {};
}

//! Defines if the grid is h-adaptive
template<class TypeTag, class MyTypeTag>
struct  AdaptiveGrid { using type = UndefinedProperty; };

//! The type of grid adaptation
template<class TypeTag, class MyTypeTag>
struct  GridAdaptModel  { using type = UndefinedProperty; };

//! Class defining the refinement/coarsening indicator
template<class TypeTag, class MyTypeTag>
struct AdaptionIndicator { using type = UndefinedProperty; };

//! Class defining the refinement/coarsening indicator for grid initialization
template<class TypeTag, class MyTypeTag>
struct AdaptionInitializationIndicator { using type = UndefinedProperty; };

//! Tolerance for refinement
template<class TypeTag, class MyTypeTag>
struct GridAdaptRefineThreshold { using type = UndefinedProperty; };

//! Tolerance for coarsening
template<class TypeTag, class MyTypeTag>
struct GridAdaptCoarsenThreshold { using type = UndefinedProperty; };

//no adaptive grid
template<class TypeTag>
struct AdaptiveGrid<TypeTag, TTag::GridAdapt> { static constexpr bool value = false; };

//Set default class for adaptation initialization indicator
template<class TypeTag>
struct AdaptionInitializationIndicator<TypeTag, TTag::GridAdapt> { using type = GridAdaptInitializationIndicatorDefault<TypeTag>; };
//Set default class for adaptation
template<class TypeTag>
struct GridAdaptModel<TypeTag, TTag::GridAdapt> { using type = GridAdapt<TypeTag, getPropValue<TypeTag, Properties::AdaptiveGrid>()>; };

//standard setting
template<class TypeTag>
struct GridAdaptRefineThreshold<TypeTag, TTag::GridAdapt> { static constexpr auto value  = 0.0; };
template<class TypeTag>
struct GridAdaptCoarsenThreshold<TypeTag, TTag::GridAdapt> { static constexpr auto value  = 0.0; };
} // namespace Properties
} // namespace Dumux


#endif
