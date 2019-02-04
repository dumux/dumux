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
 * \brief Base class for the flux variables
 */
#ifndef DUMUX_GODUNOV_FLUXVARIABLESCACHE_HH
#define DUMUX_GODUNOV_FLUXVARIABLESCACHE_HH

#include <dumux/common/properties.hh>
#include <dumux/discretization/methods.hh>
#include <dumux/discretization/fluxvariablescaching.hh>

namespace Dumux {

// forward declaration
template<class TypeTag, DiscretizationMethod discMethod>
class SweFluxVariablesCacheImplementation;

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
//! The cache is dependent on the active physical processes (advection, diffusion, heat conduction)
//! For each type of process there is a base cache storing the data required to compute the respective fluxes
//! Specializations of the overall cache are provided for combinations of processes
///////////////////////////////////////////////////////////////////////////////////////////////////////////////

/*!
 * \ingroup ImplicitModel
 * \brief The flux variables cache classes for porous media.
 *        Store data required for flux calculation. For each type of physical process (advection, diffusion, heat conduction)
 *        there is a base cache storing the data required to compute the respective fluxes. Specializations of the overall
 *        cache class are provided for different combinations of processes.
 */
template<class TypeTag>
using SweFluxVariablesCache = SweFluxVariablesCacheImplementation<TypeTag, GET_PROP_TYPE(TypeTag, FVGridGeometry)::discMethod>;


// the following classes choose the cache type: empty if the law disabled and the law's cache if it's enabled
// if advections is disabled the advection type is still instatiated if we use std::conditional_t and has to be a full type
// in order to prevent that instead of std::conditional_t we use this helper type is only dependent on the advection type
// if advection is enabled otherwise its an empty cache type
//template<class TypeTag, bool EnableAdvection> class AdvectionCacheChooser : public FluxVariablesCaching::EmptyAdvectionCache {};
//template<class TypeTag> class AdvectionCacheChooser<TypeTag, true> : public GET_PROP_TYPE(TypeTag, AdvectionType)::Cache {};
//template<class TypeTag, bool EnableMolecularDiffusion> class DiffusionCacheChooser : public FluxVariablesCaching::EmptyDiffusionCache {};
//template<class TypeTag> class DiffusionCacheChooser<TypeTag, true> : public GET_PROP_TYPE(TypeTag, MolecularDiffusionType)::Cache {};
//template<class TypeTag, bool EnableEnergyBalance> class EnergyCacheChooser : public FluxVariablesCaching::EmptyHeatConductionCache {};
//template<class TypeTag> class EnergyCacheChooser<TypeTag, true> : public GET_PROP_TYPE(TypeTag, HeatConductionType)::Cache {};


// specialization for the cell centered tpfa method
template<class TypeTag>
class SweFluxVariablesCacheImplementation<TypeTag, DiscretizationMethod::godunov>
{};


} // end namespace Dumux

#endif
