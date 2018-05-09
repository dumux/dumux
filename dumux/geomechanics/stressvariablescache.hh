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
 * \ingroup Geomechanics
 * \brief Base class for the stress variables cache
 */
#ifndef DUMUX_GEOMECHANICS_STRESSVARIABLESCACHE_HH
#define DUMUX_GEOMECHANICS_STRESSVARIABLESCACHE_HH

#include <dune/common/exceptions.hh>

#include <dumux/discretization/methods.hh>
#include <dumux/discretization/fluxvariablescaching.hh>
#include <dumux/discretization/box/fluxvariablescache.hh>

namespace Dumux {

/*!
 * \ingroup Geomechanics
 * \brief The stress variables cache classes for models involving geomechanics.
 *        Store data required for stress calculation.
 */
template< class Scalar, class FVGridGeometry, DiscretizationMethod dm = FVGridGeometry::discMethod >
class StressVariablesCache;

//! We only store discretization-related quantities for the box method.
template< class Scalar, class FVGridGeometry >
class StressVariablesCache<Scalar, FVGridGeometry, DiscretizationMethod::box>
: public BoxFluxVariablesCache< Scalar, FVGridGeometry >
{};

// specialization for the cell centered tpfa method
template< class Scalar, class FVGridGeometry >
class StressVariablesCache<Scalar, FVGridGeometry, DiscretizationMethod::cctpfa> : public FluxVariablesCaching::_EmptyCache
{
public:
    template<typename... Args>
    void update(Args&&... args) { DUNE_THROW(Dune::NotImplemented, "Geomechanics with cell-centered schemes"); }
};

// specialization for the cell centered mpfa method
template< class Scalar, class FVGridGeometry >
class StressVariablesCache<Scalar, FVGridGeometry, DiscretizationMethod::ccmpfa>
: public StressVariablesCache<Scalar, FVGridGeometry, DiscretizationMethod::cctpfa>
{};

} // end namespace Dumux

#endif
