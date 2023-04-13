// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup Flux
 * \brief Constant velocity advective law for transport models.
 *        This file contains the data which is required to calculate
 *        volume and mass fluxes of fluid phases over a face of a finite volume.
 *        A stationary velocity field is given by the user for use in tracer models.
 */
#ifndef DUMUX_DISCRETIZATION_STATIONARY_VELOCITY_FIELD_HH
#define DUMUX_DISCRETIZATION_STATIONARY_VELOCITY_FIELD_HH

#include <type_traits>
#include <dumux/flux/traits.hh>

#include <dumux/discretization/method.hh>
#include <dumux/flux/fluxvariablescaching.hh>

namespace Dumux {

/*!
 * \ingroup Flux
 * \brief Evaluates a user given velocity field
 */
template <class Scalar>
class StationaryVelocityField
{
public:
    using DiscretizationMethod = DiscretizationMethods::None;
    //! state the discretization method this implementation belongs to
    static constexpr DiscretizationMethod discMethod{};

    //! state the type for the corresponding cache
    using Cache = FluxVariablesCaching::EmptyAdvectionCache;

    //! returns the volume flux given in the spatial params
    template<class Problem, class Element,
             class FVElementGeometry,
             class ElementVolumeVariables,
             class ElementFluxVarsCache>
    static Scalar flux(const Problem& problem,
                       const Element& element,
                       const FVElementGeometry& fvGeometry,
                       const ElementVolumeVariables& elemVolVars,
                       const typename FVElementGeometry::SubControlVolumeFace& scvf,
                       int phaseIdx,
                       const ElementFluxVarsCache& elemFluxVarsCache)
    {
        //! Obtain the volume flux from the user, specified in the spatial params in m^3/s
        return problem.spatialParams().volumeFlux(element, fvGeometry, elemVolVars, scvf);
    }
};

//! Set stationary velocity field to true in the FluxTraits
template<class Scalar>
struct HasStationaryVelocityField<StationaryVelocityField<Scalar>>
: public std::true_type {};

} // end namespace Dumux

#endif
