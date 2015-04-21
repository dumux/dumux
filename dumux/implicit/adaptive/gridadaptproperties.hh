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
 * \ingroup ImplicitGridAdaptProperties
 * \ingroup ImplicitGridAdapt
 * \file
 *
 * \brief Defines a type tag and some fundamental properties for
 *        linear solvers
 */
#ifndef DUMUX_IMPLICIT_GRIDADAPT_PROPERTIES_HH
#define DUMUX_IMPLICIT_GRIDADAPT_PROPERTIES_HH

#include <dumux/common/basicproperties.hh>

namespace Dumux
{
namespace Properties
{
//! Grid adaptation type tag for all decoupled models.
NEW_TYPE_TAG(GridAdapt);

//! Defines if the grid is h-adaptive
NEW_PROP_TAG(AdaptiveGrid);

//! Class defining the refinement/coarsening indicator
NEW_PROP_TAG(AdaptationIndicator);

//! Class defining the refinement/coarsening indicator for grid initialization
NEW_PROP_TAG(AdaptationInitializationIndicator);

//! Switch the use of initial grid adaptation on/off
NEW_PROP_TAG(GridAdaptEnableInitializationIndicator);

//! Mimimum allowed level
NEW_PROP_TAG(GridAdaptMinLevel);

//! Maximum allowed level
NEW_PROP_TAG(GridAdaptMaxLevel);

//! Tolerance for refinement
NEW_PROP_TAG(GridAdaptRefineTolerance);

//! Tolerance for coarsening
NEW_PROP_TAG(GridAdaptCoarsenTolerance);

//! Tolerance for refinement
NEW_PROP_TAG(GridAdaptRefineThreshold);

//! Tolerance for coarsening
NEW_PROP_TAG(GridAdaptCoarsenThreshold);

//! Time step interval for adaptation
NEW_PROP_TAG(GridAdaptAdaptationInterval);

//! Switch for refinement at Dirichlet BC's -> not used by all indicators!
NEW_PROP_TAG(GridAdaptRefineAtDirichletBC);

//! Switch for refinement at Neumann BC's -> not used by all indicators!
NEW_PROP_TAG(GridAdaptRefineAtFluxBC);

//! Switch for refinement at sources -> not used by all indicators!
NEW_PROP_TAG(GridAdaptRefineAtSource);

} // namespace Properties
} // namespace Dumux

#endif
