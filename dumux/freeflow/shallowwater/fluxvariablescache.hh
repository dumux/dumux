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
 * \ingroup ShallowWaterModel
 * \brief Base class for the flux variables
 */
#ifndef DUMUX_FREEFLOW_SHALLOW_WATER_FLUXVARIABLESCACHE_HH
#define DUMUX_FREEFLOW_SHALLOW_WATER_FLUXVARIABLESCACHE_HH

#include <dumux/common/properties.hh>
#include <dumux/discretization/method.hh>
#include <dumux/discretization/box/fluxvariablescache.hh>
#include <dumux/flux/fluxvariablescaching.hh>

namespace Dumux {

// forward declaration
template<class TypeTag, DiscretizationMethod discMethod>
class ShallowWaterFluxVariablesCacheImplementation;

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
//! The cache is dependent on the active physical processes (advection, diffusion, heat conduction)
//! For each type of process there is a base cache storing the data required to compute the respective fluxes
//! Specializations of the overall cache are provided for combinations of processes
//!
//! Caches are not used for the Shallow Water Model since we use a Riemann Solver for all kind of fluxes
///////////////////////////////////////////////////////////////////////////////////////////////////////////////

/*!
 * \ingroup ImplicitModel
 * \brief The flux variables cache classes.
 *        Store data required for flux calculation. For each type of physical process (advection, diffusion, heat conduction)
 *        there is a base cache storing the data required to compute the respective fluxes. Specializations of the overall
 *        cache class are provided for different combinations of processes.
 */
template<class TypeTag>
using ShallowWaterFluxVariablesCache = ShallowWaterFluxVariablesCacheImplementation<TypeTag, GetPropType<TypeTag, Properties::FVGridGeometry>::discMethod>;


// the following classes choose the cache type: empty if the law disabled and the law's cache if it's enabled
// if advections is disabled the advection type is still instatiated if we use std::conditional_t and has to be a full type
// in order to prevent that instead of std::conditional_t we use this helper type is only dependent on the advection type
// if advection is enabled otherwise its an empty cache type
template<class TypeTag, bool EnableAdvection> class AdvectionCacheChooser : public FluxVariablesCaching::EmptyAdvectionCache {};
template<class TypeTag> class AdvectionCacheChooser<TypeTag, true> : public GetPropType<TypeTag, Properties::AdvectionType>::Cache {};
template<class TypeTag, bool EnableDiffusion> class DiffusionCacheChooser : public FluxVariablesCaching::EmptyDiffusionCache {};
template<class TypeTag> class DiffusionCacheChooser<TypeTag, true> : public GetPropType<TypeTag, Properties::AdvectionType>::Cache {};


// specialization for the cell centered tpfa method
template<class TypeTag>
class ShallowWaterFluxVariablesCacheImplementation<TypeTag, DiscretizationMethod::cctpfa>
{};


} // end namespace Dumux

#endif
