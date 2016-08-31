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
 * \ingroup ImplicitGridAdaptPropertyDefaults
 * \ingroup ImplicitGridAdapt
 * \file
 *
 * \brief Defines a type tag and some fundamental properties for
 *        linear solvers
 */
#ifndef DUMUX_IMPLICIT_GRIDADAPT_PROPERTY_DEFAULTS_HH
#define DUMUX_IMPLICIT_GRIDADAPT_PROPERTY_DEFAULTS_HH

#include <dumux/common/basicproperties.hh>
#include "gridadaptproperties.hh"
#include "gridadaptindicatordefault.hh"
#include "gridadaptinitializationindicator.hh"
#include "adaptionhelper.hh"

namespace Dumux
{
namespace Properties
{

//no adaptive grid
SET_BOOL_PROP(GridAdapt, AdaptiveGrid, false);

//standard setting
SET_INT_PROP(GridAdapt, GridAdaptMinLevel, 0);
SET_INT_PROP(GridAdapt, GridAdaptMaxLevel, 1);
SET_SCALAR_PROP(GridAdapt, GridAdaptRefineTolerance, 0.05);
SET_SCALAR_PROP(GridAdapt, GridAdaptCoarsenTolerance, 0.001);
SET_SCALAR_PROP(GridAdapt, GridAdaptRefineThreshold, 0.0);
SET_SCALAR_PROP(GridAdapt, GridAdaptCoarsenThreshold, 0.0);
SET_INT_PROP(GridAdapt, GridAdaptAdaptionInterval, 1);
//Switch initial grid adaption off per default
SET_BOOL_PROP(GridAdapt, GridAdaptEnableInitializationIndicator, false);

// Switch of extra refinement strategy at boundaries/sources
SET_BOOL_PROP(GridAdapt, GridAdaptRefineAtDirichletBC, false);
SET_BOOL_PROP(GridAdapt, GridAdaptRefineAtFluxBC, false);
SET_BOOL_PROP(GridAdapt, GridAdaptRefineAtSource, false);

//! Set the default indicator class models for adaption or coarsening
SET_TYPE_PROP(GridAdapt, AdaptionIndicator, ImplicitGridAdaptIndicatorDefault<TypeTag>);
//!Set default class for adaption initialization indicator
SET_TYPE_PROP(GridAdapt,  AdaptionInitializationIndicator, ImplicitGridAdaptInitializationIndicatorDefault<TypeTag>);
//!Set default class for adaption helper
SET_TYPE_PROP(GridAdapt, AdaptionHelper, ImplicitAdaptionHelper<TypeTag>);

} // namespace Properties
} // namespace Dumux

#endif
