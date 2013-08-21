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
 * \brief This file contains the parts of the local residual to
 *        calculate the heat flux in the fully coupled MpNc model.
 */
#ifndef DUMUX_MPNC_LOCAL_RESIDUAL_ENERGY_HH
#define DUMUX_MPNC_LOCAL_RESIDUAL_ENERGY_HH

#include <dumux/implicit/mpnc/mpncproperties.hh>

namespace Dumux {

/*!
 * \brief Specialization of the energy module for the isothermal case.
 *
 * This class just does nothing.
 */
template <class TypeTag, bool enableEnergy/*=false*/, bool kineticEnergyTransfer /*=false*/>
class MPNCLocalResidualEnergy
{
    static_assert(!(kineticEnergyTransfer && !enableEnergy),
                  "No kinetic energy transfer may only be enabled "
                  "if energy is enabled in general.");
    static_assert(!kineticEnergyTransfer,
                  "No kinetic energy transfer module included, "
                  "but kinetic energy transfer enabled.");

    typedef typename GET_PROP_TYPE(TypeTag, PrimaryVariables) PrimaryVariables;
    typedef typename GET_PROP_TYPE(TypeTag, VolumeVariables) VolumeVariables;
    typedef typename GET_PROP_TYPE(TypeTag, FluxVariables) FluxVariables;
    typedef typename GET_PROP_TYPE(TypeTag, ElementVolumeVariables) ElementVolumeVariables;
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;

    enum { numPhases        = GET_PROP_VALUE(TypeTag, NumPhases) };
    enum { numComponents    = GET_PROP_VALUE(TypeTag, NumComponents) };

    typedef typename Dune::FieldVector<Scalar, numComponents>  ComponentVector;

public:
    /*!
     * \brief Evaluate the amount all conservation quantities
     *        (e.g. phase mass) within a sub-control volume.
     *
     * The result should be averaged over the volume (e.g. phase mass
     * inside a sub-control volume divided by the volume)
     *
     *  \param storage The mass of the component within the sub-control volume
     *  \param volVars the Volume Variables
     */
    static void computeStorage(PrimaryVariables &storage,
                               const VolumeVariables &volVars)
    {
        // do nothing, we're isothermal!
    }

    /*!
     * \brief Calculate the storage for all mass balance equations
     *        within a single fluid phase
     *
     *  \param storage The mass of the component within the sub-control volume
     *  \param volVars the Volume Variables
     *  \param phaseIdx The local index of the phases
     */
    static void addPhaseStorage(PrimaryVariables &storage,
                                const VolumeVariables &volVars,
                                const unsigned int phaseIdx)
    {
        // do nothing, we're isothermal!
    }
    /*!
        * \brief the advective Flux of the enthalpy
        *
        *  \param enthalpy Flux advective Flux of the enthalpy
        *  \param phaseIdx The local index of the phases
        *  \param compMolFlux
        *  \param volVars the Volume Variables
        *  \param fluxVars the flux Variables
        */
    static void phaseEnthalpyFlux(PrimaryVariables &enthalpyFlux,
                                  const unsigned int phaseIdx,
                                  const PrimaryVariables &compMolFlux,
                                  const ElementVolumeVariables &volVars,
                                  const FluxVariables &fluxVars)
    {
        // do nothing, we're isothermal!
    }
    /*!
        * \brief The heat conduction in the phase
        *
        *  \param heatConduction
        *  \param volVars the Volume Variables
        *  \param fluxVars the flux Variables
        */
    static void heatConduction(PrimaryVariables &heatConduction,
                               const ElementVolumeVariables &volVars,
                               const FluxVariables &fluxVars)
    {
        // do nothing, we're isothermal!
    }

    /*!
     * \brief Evaluates the total flux of all conservation quantities
     *        over a face of a sub-control volume.
     *
     * \param flux The flux over the SCV (sub-control-volume) face for each component
     * \param volVars The volume variables
     * \param molarPhaseComponentValuesMassTransport[numPhases]
     */
    static void computeFlux(PrimaryVariables & flux,
                                const FluxVariables & fluxVars,
                                const ElementVolumeVariables & volVars,
                                const ComponentVector molarPhaseComponentValuesMassTransport[numPhases])
    {
        // do nothing, we're isothermal!
    }
    /*!
     * \brief Calculate the source term of the equation
     *
     * \param source The source/sink in the sub-control volume for each component
     * \param componentIntoPhaseMassTransfer[numPhases]
     */
    static void computeSource(PrimaryVariables &source,
                              const VolumeVariables &volVars,
                              const ComponentVector componentIntoPhaseMassTransfer[numPhases])
    {
        // do nothing, we're isothermal!
    }
};


template <class TypeTag>
class MPNCLocalResidualEnergy<TypeTag, /*enableEnergy=*/true, /*kineticenergyTransfer=*/false>
{
    typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;
    typedef typename GET_PROP_TYPE(TypeTag, FluidState) FluidState;
    typedef typename FluidSystem::ParameterCache ParameterCache;
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
    typedef typename GET_PROP_TYPE(TypeTag, PrimaryVariables) PrimaryVariables;
    typedef typename GET_PROP_TYPE(TypeTag, FluxVariables) FluxVariables;
    typedef typename GET_PROP_TYPE(TypeTag, VolumeVariables) VolumeVariables;
    typedef typename GET_PROP_TYPE(TypeTag, ElementVolumeVariables) ElementVolumeVariables;
    typedef typename GET_PROP_TYPE(TypeTag, Indices) Indices;

    enum { dim = GridView::dimension };
    enum { numPhases = GET_PROP_VALUE(TypeTag, NumPhases) };
    enum { numComponents = GET_PROP_VALUE(TypeTag, NumComponents) };
    enum { conti0EqIdx = Indices::conti0EqIdx };
    enum { energyEqIdx = Indices::energyEqIdx };

    typedef typename Dune::FieldVector<Scalar, numComponents> ComponentVector;

public:
    /*!
     * \brief Evaluate the amount all conservation quantities
     *        (e.g. phase mass) within a sub-control volume.
     *
     * The result should be averaged over the volume (e.g. phase mass
     * inside a sub-control volume divided by the volume)
     *
     *  \param storage The mass of the component within the sub-control volume
     *  \param volVars the Volume Variables
     */
    static void computeStorage(PrimaryVariables &storage,
                               const VolumeVariables &volVars)
    {
        storage[energyEqIdx] = 0;

        // energy of the fluids
        for (int phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx) {
            addPhaseStorage(storage, volVars, phaseIdx);
        }

        // heat stored in the rock matrix
        storage[energyEqIdx] +=
            volVars.fluidState().temperature(/*phaseIdx=*/0)
            * volVars.soilDensity()
            * (1.0 - volVars.porosity())
            * volVars.heatCapacity();
    }
    /*!
     * \brief Calculate the storage for all mass balance equations
     *        within a single fluid phase
     *
     *  \param storage The mass of the component within the sub-control volume
     *  \param volVars the Volume Variables
     *  \param phaseIdx The local index of the phases
     */
    static void addPhaseStorage(PrimaryVariables &storage,
                                const VolumeVariables &volVars,
                                const unsigned int phaseIdx)
    {
        const FluidState &fs = volVars.fluidState();

        // energy of the fluid
        storage[energyEqIdx] +=
            fs.density(phaseIdx)
            * fs.internalEnergy(phaseIdx)
            * fs.saturation(phaseIdx)
            * volVars.porosity();
    }
    /*!
      * \brief Evaluates the total flux of all conservation quantities
      *        over a face of a sub-control volume.
      *
      * \param flux The flux over the SCV (sub-control-volume) face for each component
      * \param fluxVars The flux Variables
      * \param volVars The volume variables
      * \param molarPhaseComponentValuesMassTransport[numPhases]
      */
    static void computeFlux(PrimaryVariables & flux,
                            const FluxVariables & fluxVars,
                            const ElementVolumeVariables & elemVolVars,
                            const ComponentVector molarPhaseComponentValuesMassTransport[numPhases])
    {
        flux[energyEqIdx] = 0.0;

        // fluid phases transport enthalpy individually
        for(int phaseIdx=0; phaseIdx<numPhases; ++phaseIdx)
            computePhaseEnthalpyFlux(flux,
                                fluxVars,
                                elemVolVars,
                                phaseIdx,
                                molarPhaseComponentValuesMassTransport[phaseIdx]);

        //conduction is treated lumped in this model
        computeHeatConduction(flux,
                             fluxVars,
                             elemVolVars);
    }
    /*!
        * \brief the advective Flux of the enthalpy
        *        \param flux The flux over the SCV (sub-control-volume) face for each component
        *        \param fluxVars The flux Variables
        *        \param volVars The volume variables
        *        \param phaseIdx The local index of the phases
        */
    static void computePhaseEnthalpyFlux(PrimaryVariables & flux,
                                         const FluxVariables & fluxVars,
                                         const ElementVolumeVariables & elemVolVars,
                                         const unsigned int phaseIdx,
                                         const ComponentVector & molarComponentValuesMassTransport)
    {
        Scalar massFlux = 0;

        // calculate the mass flux in the phase i.e. make mass flux out of mole flux and add up the fluxes of a phase
        for (int compIdx = 0; compIdx < numComponents; ++compIdx)
            massFlux += molarComponentValuesMassTransport[compIdx]
                                    * FluidSystem::molarMass(compIdx);

        unsigned int upIdx = fluxVars.face().i;
        if (massFlux < 0) upIdx = fluxVars.face().j;

        // use the phase enthalpy of the upstream vertex to calculate
        // the enthalpy transport
        const VolumeVariables &up = elemVolVars[upIdx];
        flux[energyEqIdx] += up.fluidState().enthalpy(phaseIdx) * massFlux;
    }
    /*!
        * \brief The heat conduction in the phase
        *
        *        \param flux The flux over the SCV (sub-control-volume) face for each component
        *        \param fluxVars The flux Variables
        *        \param elemVolVars The volume variables of the current element
        *
        */
    static void computeHeatConduction(PrimaryVariables & flux,
                                    const FluxVariables & fluxVars,
                                    const ElementVolumeVariables & elemVolVars)
    {
        //lumped heat conduction of the rock matrix and the fluid phases
        Scalar lumpedConductivity   = fluxVars.fluxVarsEnergy().lambdaPm() ;
        Scalar temperatureGradientNormal  = fluxVars.fluxVarsEnergy().temperatureGradientNormal() ;
        Scalar lumpedHeatConduction = - lumpedConductivity * temperatureGradientNormal ;
        flux[energyEqIdx] += lumpedHeatConduction ;
    }

    /*!
     * \brief Calculate the source term of the equation
     *
     * \param source The source/sink in the sub-control volume for each component
     * \param volVars The volume variables
     * \param componentIntoPhaseMassTransfer[numPhases]
     */

    static void computeSource(PrimaryVariables &source,
                              const VolumeVariables &volVars,
                              const ComponentVector componentIntoPhaseMassTransfer[numPhases])
    {
        source[energyEqIdx] = 0.0;
    }
};





}

#endif // DUMUX_MPNC_ENERGY_HH
