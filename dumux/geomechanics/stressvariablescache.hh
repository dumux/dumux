// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup GeomechanicsModels
 * \brief Base class for the stress variables cache
 */
#ifndef DUMUX_GEOMECHANICS_STRESSVARIABLESCACHE_HH
#define DUMUX_GEOMECHANICS_STRESSVARIABLESCACHE_HH

#include <dune/common/exceptions.hh>

#include <dumux/discretization/method.hh>
#include <dumux/flux/fluxvariablescaching.hh>
#include <dumux/discretization/cvfe/fluxvariablescache.hh>

namespace Dumux {

/*!
 * \ingroup GeomechanicsModels
 * \brief The stress variables cache classes for models involving geomechanics.
 *        Store data required for stress calculation.
 */
template< class Scalar, class GridGeometry, class DiscretizationMethod = typename GridGeometry::DiscretizationMethod >
class StressVariablesCache;

//! We only store discretization-related quantities for the box method.
template< class Scalar, class GridGeometry >
class StressVariablesCache<Scalar, GridGeometry, DiscretizationMethods::Box>
: public CVFEFluxVariablesCache< Scalar, GridGeometry >
{};

// specialization for the cell centered tpfa method
template< class Scalar, class GridGeometry >
class StressVariablesCache<Scalar, GridGeometry, DiscretizationMethods::CCTpfa>
: public FluxVariablesCaching::_EmptyCache
{
public:
    /*!
     * \brief Currently, we do not consider cell-centered schemes for geomechanics.
     *        In case this is to be integrated, one would have to rethink the structure
     *        of e.g. the elastic volume variables and what quantities they store that
     *        are necessary for flux/stress evaluation. In the porous medium flow context we
     *        define permeabilities/thermal conductivities etc. in the volume variables
     *        and then do the harmonic average at scvf integration points during flux calculation.
     *        However, for the box scheme one could evaluate the parameters required for flux
     *        computations directly at the integration point without averaging. For compatibility
     *        reasons with cell-centered schemes, and because one can derive an expression for the
     *        harmonic average, we do not do so in the porous medium framework. For now, we choose
     *        a more FEM-like mentality in the geomechanics framework and call the parameters in the
     *        spatial parameters, required for stress calculations, for a position in the element
     *        when assembling the stress tensors. We do not store them in the volume variables! This
     *        means that in case cell-centered geomechanical models are considered in the future, both
     *        the volume variables as well as the stress tensor assembly laws have to be restructured!
     */
    template<typename... Args>
    void update(Args&&... args)
    { DUNE_THROW(Dune::NotImplemented, "Geomechanics with cell-centered schemes"); }
};

// specialization for the cell centered mpfa method
template< class Scalar, class GridGeometry >
class StressVariablesCache<Scalar, GridGeometry, DiscretizationMethods::CCMpfa>
: public StressVariablesCache<Scalar, GridGeometry, DiscretizationMethods::CCTpfa>
{};

} // end namespace Dumux

#endif
