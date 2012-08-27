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
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the           *
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

#include <dumux/common/basicproperties.hh>

namespace Dumux
{
namespace Properties
{
//! Grid adaption type tag for all decoupled models.
NEW_TYPE_TAG(GridAdaptTypeTag);

//! Defines if the grid is h-adaptive
NEW_PROP_TAG( AdaptiveGrid);

//! Class defining the refinement/coarsening indicator
NEW_PROP_TAG(AdaptionIndicator);

//! Class defining the refinement/coarsening indicator for grid initialization
NEW_PROP_TAG(AdaptionInitializationIndicator);

//! Switch the use of initial grid adaption on/of
NEW_PROP_TAG(GridAdaptEnableInitializationIndicator);
NEW_PROP_TAG(EnableInitializationIndicator);//DEPRECATED

//! Mimimum allowed level
NEW_PROP_TAG(GridAdaptMinLevel);
NEW_PROP_TAG(MinLevel);//DEPRECATED

//! Maximum allowed level
NEW_PROP_TAG(GridAdaptMaxLevel);
NEW_PROP_TAG(MaxLevel);//DEPRECATED

//! Tolerance for refinement
NEW_PROP_TAG(GridAdaptRefineTolerance);
NEW_PROP_TAG(RefineTolerance);//DEPRECATED

//! Tolerance for coarsening
NEW_PROP_TAG(GridAdaptCoarsenTolerance);
NEW_PROP_TAG(CoarsenTolerance);//DEPRECATED

//! Tolerance for refinement
NEW_PROP_TAG(GridAdaptRefineThreshold);
NEW_PROP_TAG(RefineThreshold);//DEPRECATED

//! Tolerance for coarsening
NEW_PROP_TAG(GridAdaptCoarsenThreshold);
NEW_PROP_TAG(CoarsenThreshold);//DEPRECATED

//! Time step interval for adaption
NEW_PROP_TAG(GridAdaptAdaptionInterval);
NEW_PROP_TAG(AdaptionInterval);//DEPRECATED

//! Switch for refinement at dirichlet BC's -> not used by all indicators!
NEW_PROP_TAG(GridAdaptRefineAtDirichletBC);
NEW_PROP_TAG(RefineAtDirichletBC);//DEPRECATED

//! Switch for refinement at neumann BC's -> not used by all indicators!
NEW_PROP_TAG(GridAdaptRefineAtFluxBC);
NEW_PROP_TAG(RefineAtFluxBC);//DEPRECATED

//! Switch for refinement at sources -> not used by all indicators!
NEW_PROP_TAG(GridAdaptRefineAtSource);
NEW_PROP_TAG(RefineAtSource);//DEPRECATED

//no adaptive grid
SET_BOOL_PROP(GridAdaptTypeTag, AdaptiveGrid, false);

//standard setting
SET_INT_PROP(GridAdaptTypeTag, GridAdaptMinLevel, GET_PROP_VALUE(TypeTag, MinLevel));
SET_INT_PROP(GridAdaptTypeTag, MinLevel, 0);//DEPRECATED
SET_INT_PROP(GridAdaptTypeTag, GridAdaptMaxLevel, GET_PROP_VALUE(TypeTag, MaxLevel));
SET_INT_PROP(GridAdaptTypeTag, MaxLevel, 1);//DEPRECATED
SET_SCALAR_PROP(GridAdaptTypeTag, GridAdaptRefineTolerance, GET_PROP_VALUE(TypeTag, RefineTolerance));
SET_SCALAR_PROP(GridAdaptTypeTag, RefineTolerance, 0.05);//DEPRECATED
SET_SCALAR_PROP(GridAdaptTypeTag, GridAdaptCoarsenTolerance, GET_PROP_VALUE(TypeTag, CoarsenTolerance));
SET_SCALAR_PROP(GridAdaptTypeTag, CoarsenTolerance, 0.001);//DEPRECATED
SET_SCALAR_PROP(GridAdaptTypeTag, GridAdaptRefineThreshold, GET_PROP_VALUE(TypeTag, RefineThreshold));
SET_SCALAR_PROP(GridAdaptTypeTag, RefineThreshold, 0.0);//DEPRECATED
SET_SCALAR_PROP(GridAdaptTypeTag, GridAdaptCoarsenThreshold, GET_PROP_VALUE(TypeTag, CoarsenThreshold));
SET_SCALAR_PROP(GridAdaptTypeTag, CoarsenThreshold, 0.0);//DEPRECATED
SET_INT_PROP(GridAdaptTypeTag, GridAdaptAdaptionInterval, GET_PROP_VALUE(TypeTag, AdaptionInterval));
SET_INT_PROP(GridAdaptTypeTag, AdaptionInterval, 1);//DEPRECATED
//Switch initial grid adaption off per default
SET_BOOL_PROP(GridAdaptTypeTag, GridAdaptEnableInitializationIndicator, GET_PROP_VALUE(TypeTag, EnableInitializationIndicator));
SET_BOOL_PROP(GridAdaptTypeTag, EnableInitializationIndicator, false);//DEPRECATED

// Switch of extra refinement strategy at boundaries/sources
SET_BOOL_PROP(GridAdaptTypeTag, GridAdaptRefineAtDirichletBC, GET_PROP_VALUE(TypeTag, RefineAtDirichletBC));
SET_BOOL_PROP(GridAdaptTypeTag, RefineAtDirichletBC, false);//DEPRECATED
SET_BOOL_PROP(GridAdaptTypeTag, GridAdaptRefineAtFluxBC, GET_PROP_VALUE(TypeTag, RefineAtFluxBC));
SET_BOOL_PROP(GridAdaptTypeTag, RefineAtFluxBC, false);//DEPRECATED
SET_BOOL_PROP(GridAdaptTypeTag, GridAdaptRefineAtSource, GET_PROP_VALUE(TypeTag, RefineAtSource));
SET_BOOL_PROP(GridAdaptTypeTag, RefineAtSource, false);//DEPRECATED
} // namespace Properties
} // namespace Dumux


#endif
