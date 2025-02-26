// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup NIModel
 * \brief Element-wise calculation of the local residual for non-isothermal
 *        fully implicit models. See NIModel for the detailed description.
 */

#ifndef DUMUX_ENERGY_LOCAL_RESIDUAL_HH
#define DUMUX_ENERGY_LOCAL_RESIDUAL_HH

#include <dumux/common/properties.hh>
#include <dumux/common/numeqvector.hh>

namespace Dumux {

// forward declaration
template<class TypeTag, bool enableEneryBalance>
class EnergyLocalResidualImplementation;

template<class TypeTag>
using EnergyLocalResidual = EnergyLocalResidualImplementation<TypeTag, GetPropType<TypeTag, Properties::ModelTraits>::enableEnergyBalance()>;

/*!
 * \ingroup NIModel
 * \brief Element-wise calculation of the energy residual for isothermal problems.
 */
template<class TypeTag>
class EnergyLocalResidualImplementation<TypeTag, false>
{
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using NumEqVector = Dumux::NumEqVector<GetPropType<TypeTag, Properties::PrimaryVariables>>;
    using Problem = GetPropType<TypeTag, Properties::Problem>;
    using VolumeVariables = GetPropType<TypeTag, Properties::VolumeVariables>;
    using FVElementGeometry = typename GetPropType<TypeTag, Properties::GridGeometry>::LocalView;
    using SubControlVolume = typename FVElementGeometry::SubControlVolume;
    using FluxVariables = GetPropType<TypeTag, Properties::FluxVariables>;

public:
    template <typename T = void>
    static void fluidPhaseStorage(NumEqVector& storage,
                                  const SubControlVolume& scv,
                                  const VolumeVariables& volVars,
                                  int phaseIdx)
    {
        static_assert(AlwaysFalse<T>::value, "Deprecated interface that has been removed! Use new interface with additional argument problem instead.");
    }

    /*!
     * \brief The energy storage in the fluid phase with index phaseIdx.
     */
    static void fluidPhaseStorage(NumEqVector& storage,
                                  const Problem& problem,
                                  const SubControlVolume& scv,
                                  const VolumeVariables& volVars,
                                  int phaseIdx)
    {}

    /*!
     * \brief The energy storage in the solid matrix.
     *
     * \param storage The mass of the component within the sub-control volume
     * \param scv The sub-control volume
     * \param volVars The volume variables
     */
    static void solidPhaseStorage(NumEqVector& storage,
                                  const SubControlVolume& scv,
                                  const VolumeVariables& volVars)
    {}

    /*!
     * \brief The advective phase energy fluxes.
     *
     * \param flux The flux
     * \param fluxVars The flux variables.
     * \param phaseIdx The phase index
     */
    static void heatConvectionFlux(NumEqVector& flux,
                                   FluxVariables& fluxVars,
                                   int phaseIdx)
    {}

    /*!
     * \brief The diffusive energy fluxes
     *
     * \param flux The flux
     * \param fluxVars The flux variables.
     */
    static void heatConductionFlux(NumEqVector& flux,
                                   FluxVariables& fluxVars)
    {}

    /*!
     * \brief The energy fluxes related to molecular diffusion
     *
     * \param flux The flux
     * \param fluxVars The flux variables.
     */
    static void heatMolecularDiffusionFlux(NumEqVector& flux,
                                           FluxVariables& fluxVars,
                                           int phaseIdx)
    {}

    /*!
     * \brief The dispersive energy fluxes
     *
     * \param flux The flux
     * \param fluxVars The flux variables.
     */
    static void heatDispersionFlux(NumEqVector& flux,
                                   FluxVariables& fluxVars)
    {}
};

/*!
 * \ingroup NIModel
 * \brief Element-wise calculation of the energy residual for non-isothermal problems.
 */
template<class TypeTag>
class EnergyLocalResidualImplementation<TypeTag, true>
{
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using NumEqVector = Dumux::NumEqVector<GetPropType<TypeTag, Properties::PrimaryVariables>>;
    using Problem = GetPropType<TypeTag, Properties::Problem>;
    using VolumeVariables = GetPropType<TypeTag, Properties::VolumeVariables>;
    using FVElementGeometry = typename GetPropType<TypeTag, Properties::GridGeometry>::LocalView;
    using SubControlVolume = typename FVElementGeometry::SubControlVolume;
    using FluxVariables = GetPropType<TypeTag, Properties::FluxVariables>;
    using GridView = typename GetPropType<TypeTag, Properties::GridGeometry>::GridView;
    using Element = typename GridView::template Codim<0>::Entity;
    using ElementVolumeVariables = typename GetPropType<TypeTag, Properties::GridVolumeVariables>::LocalView;
    using ModelTraits = GetPropType<TypeTag, Properties::ModelTraits>;
    using Indices = typename ModelTraits::Indices;

    static constexpr int numPhases = ModelTraits::numFluidPhases();
    enum { energyEqIdx = Indices::energyEqIdx };

public:

    /*!
     * \brief The energy storage in the fluid phase with index phaseIdx.
     */
    static void fluidPhaseStorage(NumEqVector& storage,
                                  const Problem& problem,
                                  const SubControlVolume& scv,
                                  const VolumeVariables& volVars,
                                  int phaseIdx)
    {
        // this implementation of the potential energy contribution is only correct
        // if gravity vector is constant in space and time
        const auto& x = scv.dofPosition();
        const auto gravityPotential = x*problem.spatialParams().gravity(x);

        storage[energyEqIdx] += volVars.porosity()
                                * volVars.density(phaseIdx)
                                * volVars.saturation(phaseIdx)
                                * (volVars.internalEnergy(phaseIdx) - gravityPotential);
    }

    template <typename T = void>
    static void fluidPhaseStorage(NumEqVector& storage,
                                  const SubControlVolume& scv,
                                  const VolumeVariables& volVars,
                                  int phaseIdx)
    {
        static_assert(AlwaysFalse<T>::value, "Deprecated interface that has been removed! Use new interface with additional argument problem instead.");
    }

    /*!
     * \brief The energy storage in the solid matrix
     *
     * \param storage The mass of the component within the sub-control volume
     * \param scv The sub-control volume
     * \param volVars The volume variables
     */
    static void solidPhaseStorage(NumEqVector& storage,
                                  const SubControlVolume& scv,
                                  const VolumeVariables& volVars)
    {
        storage[energyEqIdx] += volVars.temperature()
                                * volVars.solidHeatCapacity()
                                * volVars.solidDensity()
                                * (1.0 - volVars.porosity());
    }

    /*!
     * \brief The advective phase energy fluxes
     *
     * \param flux The flux
     * \param fluxVars The flux variables.
     * \param phaseIdx The phase index
     */
    static void heatConvectionFlux(NumEqVector& flux,
                                   FluxVariables& fluxVars,
                                   int phaseIdx)
    {
        // this implementation of the potential energy contribution is only correct
        // if gravity vector is constant in space and time
        const auto& x = fluxVars.scvFace().ipGlobal();
        const auto gravityPotential = x*fluxVars.problem().spatialParams().gravity(x);

        auto upwindTerm = [=](const auto& volVars){
            return volVars.density(phaseIdx)*volVars.mobility(phaseIdx)
                * (volVars.enthalpy(phaseIdx) - gravityPotential);
        };

        flux[energyEqIdx] += fluxVars.advectiveFlux(phaseIdx, upwindTerm);
    }

    /*!
     * \brief The diffusive energy fluxes
     *
     * \param flux The flux
     * \param fluxVars The flux variables.
     */
    static void heatConductionFlux(NumEqVector& flux,
                                   FluxVariables& fluxVars)
    {
        flux[energyEqIdx] += fluxVars.heatConductionFlux();
    }

    /*!
     * \brief The energy fluxes related to molecular diffusion
     *
     * \param flux The flux
     * \param fluxVars The flux variables.
     */
    static void heatMolecularDiffusionFlux(NumEqVector& flux,
                                           FluxVariables& fluxVars,
                                           int phaseIdx)
    {
        flux[energyEqIdx] += fluxVars.heatMolecularDiffusionFlux(phaseIdx);
    }

    /*!
     * \brief The dispersive energy fluxes
     *
     * \param flux The flux
     * \param fluxVars The flux variables.
     */
    static void heatDispersionFlux(NumEqVector& flux,
                                   FluxVariables& fluxVars)
    {

        if constexpr (ModelTraits::enableThermalDispersion())
        {
            flux[energyEqIdx] += fluxVars.thermalDispersionFlux();
        }
    }


    /*!
     * \brief heat transfer between the phases for nonequilibrium models
     *
     * \param source The source which ought to be simulated
     * \param element An element which contains part of the control volume
     * \param fvGeometry The finite-volume geometry
     * \param elemVolVars The volume variables of the current element
     * \param scv The sub-control volume over which we integrate the source term
     */
    static void computeSourceEnergy(NumEqVector& source,
                                    const Element& element,
                                    const FVElementGeometry& fvGeometry,
                                    const ElementVolumeVariables& elemVolVars,
                                    const SubControlVolume &scv)
    {}
};

} // end namespace Dumux

#endif
