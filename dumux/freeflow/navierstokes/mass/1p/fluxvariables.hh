// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup NavierStokesModel
 * \copydoc Dumux::NavierStokesMassOnePFluxVariables
 */
#ifndef DUMUX_FREEFLOW_NAVIERSTOKES_MASS_1P_FLUXVARIABLES_HH
#define DUMUX_FREEFLOW_NAVIERSTOKES_MASS_1P_FLUXVARIABLES_HH

#include <dumux/flux/upwindscheme.hh>
#include <dumux/freeflow/navierstokes/scalarfluxvariables.hh>

#include "advectiveflux.hh"

namespace Dumux {

/*!
 * \ingroup NavierStokesModel
 * \brief The flux variables class for the single-phase flow Navier-Stokes model.
 */
template<class Problem,
         class ModelTraits,
         class FluxTs,
         class ElementVolumeVariables,
         class ElementFluxVariablesCache,
         class UpwindScheme = UpwindScheme<typename ProblemTraits<Problem>::GridGeometry>>
class NavierStokesMassOnePFluxVariables
: public NavierStokesScalarConservationModelFluxVariables<Problem,
                                                          ModelTraits,
                                                          FluxTs,
                                                          ElementVolumeVariables,
                                                          ElementFluxVariablesCache,
                                                          UpwindScheme>
{
    using ParentType = NavierStokesScalarConservationModelFluxVariables<Problem,
                                                                        ModelTraits,
                                                                        FluxTs,
                                                                        ElementVolumeVariables,
                                                                        ElementFluxVariablesCache,
                                                                        UpwindScheme>;

    using VolumeVariables = typename ElementVolumeVariables::VolumeVariables;
    using NumEqVector = typename VolumeVariables::PrimaryVariables;

public:

    /*!
     * \brief Returns the advective mass flux in kg/s.
     */
    NumEqVector advectiveFlux(int phaseIdx = 0) const
    {
        NumEqVector result(0.0);
        // g++ requires to capture 'this' by value
        const auto upwinding = [this](const auto& term) { return this->getAdvectiveFlux(term); };
        AdvectiveFlux<ModelTraits>::addAdvectiveFlux(result, upwinding);
        return result;
    }

    /*!
     * \brief Returns all fluxes for the single-phase flow Navier-Stokes model: the
     *        advective mass flux in kg/s and the energy flux in J/s (for nonisothermal models).
     */
    NumEqVector flux(int phaseIdx = 0) const
    {
        NumEqVector flux = advectiveFlux(phaseIdx);
        ParentType::addHeatFlux(flux);
        return flux;
    }
};

} // end namespace Dumux

#endif
