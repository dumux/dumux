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
NEW_TYPE_TAG(GridAdapt);

//! Defines if the grid is h-adaptive
NEW_PROP_TAG( AdaptiveGrid);

//! The type of grid adaptation
NEW_PROP_TAG( GridAdaptModel );

//! Class defining the refinement/coarsening indicator
NEW_PROP_TAG(AdaptionIndicator);

//! Class defining the refinement/coarsening indicator for grid initialization
NEW_PROP_TAG(AdaptionInitializationIndicator);

//! Tolerance for refinement
NEW_PROP_TAG(GridAdaptRefineThreshold);

//! Tolerance for coarsening
NEW_PROP_TAG(GridAdaptCoarsenThreshold);

//no adaptive grid
SET_BOOL_PROP(GridAdapt, AdaptiveGrid, false);

//Set default class for adaptation initialization indicator
SET_TYPE_PROP(GridAdapt,  AdaptionInitializationIndicator, GridAdaptInitializationIndicatorDefault<TypeTag>);
//Set default class for adaptation
SET_TYPE_PROP(GridAdapt,  GridAdaptModel, GridAdapt<TypeTag, GET_PROP_VALUE(TypeTag, AdaptiveGrid)>);


//standard setting
SET_SCALAR_PROP(GridAdapt, GridAdaptRefineThreshold, 0.0);
SET_SCALAR_PROP(GridAdapt, GridAdaptCoarsenThreshold, 0.0);
} // namespace Properties
} // namespace Dumux


#endif
