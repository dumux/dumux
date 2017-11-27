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
 *
 * \brief Element-wise calculation of the local residual for non-isothermal
 *        fully implicit models.
 */
#ifndef DUMUX_ENERGY_LOCAL_RESIDUAL_HH
#define DUMUX_ENERGY_LOCAL_RESIDUAL_HH

namespace Dumux
{

// property forward declarations
namespace Properties
{
NEW_PROP_TAG(Indices);
}

// forward declaration
template<class TypeTag, bool enableEneryBalance>
class EnergyLocalResidualImplementation;

/*!
 * \ingroup NIModel
 * \ingroup ImplicitLocalResidual
 * \brief Element-wise calculation of the energy residual for non-isothermal problems.
 */
template<class TypeTag>
using EnergyLocalResidual = EnergyLocalResidualImplementation<TypeTag, GET_PROP_VALUE(TypeTag, EnableEnergyBalance)>;

template<class TypeTag>
class EnergyLocalResidualImplementation<TypeTag, false>
{
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using ResidualVector = typename GET_PROP_TYPE(TypeTag, NumEqVector);
    using VolumeVariables = typename GET_PROP_TYPE(TypeTag, VolumeVariables);
    using SubControlVolume = typename GET_PROP_TYPE(TypeTag, SubControlVolume);
    using FluxVariables = typename GET_PROP_TYPE(TypeTag, FluxVariables);

public:
    //! The energy storage in the fluid phase with index phaseIdx
    static void fluidPhaseStorage(ResidualVector& storage,
                                  const SubControlVolume& scv,
                                  const VolumeVariables& volVars,
                                  int phaseIdx)
    {}

    //! The energy storage in the solid matrix
    static void solidPhaseStorage(ResidualVector& storage,
                                  const SubControlVolume& scv,
                                  const VolumeVariables& volVars)
    {}

    //! The advective phase energy fluxes
    static void heatConvectionFlux(ResidualVector& flux,
                                   FluxVariables& fluxVars,
                                   int phaseIdx)
    {}

    //! The diffusive energy fluxes
    static void heatConductionFlux(ResidualVector& flux,
                                   FluxVariables& fluxVars)
    {}
};

template<class TypeTag>
class EnergyLocalResidualImplementation<TypeTag, true>
{
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using ResidualVector = typename GET_PROP_TYPE(TypeTag, NumEqVector);
    using VolumeVariables = typename GET_PROP_TYPE(TypeTag, VolumeVariables);
    using SubControlVolume = typename GET_PROP_TYPE(TypeTag, SubControlVolume);
    using FluxVariables = typename GET_PROP_TYPE(TypeTag, FluxVariables);
    using Indices = typename GET_PROP_TYPE(TypeTag, Indices);

    enum { energyEqIdx = Indices::energyEqIdx };

public:

    //! The energy storage in the fluid phase with index phaseIdx
    static void fluidPhaseStorage(ResidualVector& storage,
                                  const SubControlVolume& scv,
                                  const VolumeVariables& volVars,
                                  int phaseIdx)
    {
        storage[energyEqIdx] += volVars.porosity()
                                * volVars.density(phaseIdx)
                                * volVars.internalEnergy(phaseIdx)
                                * volVars.saturation(phaseIdx);
    }

    //! The energy storage in the solid matrix
    static void solidPhaseStorage(ResidualVector& storage,
                                  const SubControlVolume& scv,
                                  const VolumeVariables& volVars)
    {
        storage[energyEqIdx] += volVars.temperature()
                                * volVars.solidHeatCapacity()
                                * volVars.solidDensity()
                                * (1.0 - volVars.porosity());
    }

    //! The advective phase energy fluxes
    static void heatConvectionFlux(ResidualVector& flux,
                                   FluxVariables& fluxVars,
                                   int phaseIdx)
    {
        auto upwindTerm = [phaseIdx](const auto& volVars)
        { return volVars.density(phaseIdx)*volVars.mobility(phaseIdx)*volVars.enthalpy(phaseIdx); };

        flux[energyEqIdx] += fluxVars.advectiveFlux(phaseIdx, upwindTerm);
    }

    //! The diffusive energy fluxes
    static void heatConductionFlux(ResidualVector& flux,
                                   FluxVariables& fluxVars)
    {
        flux[energyEqIdx] += fluxVars.heatConductionFlux();
    }
};

} // end namespace Dumux

#endif
