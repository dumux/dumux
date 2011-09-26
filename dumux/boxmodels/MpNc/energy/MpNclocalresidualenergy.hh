/*****************************************************************************
 *   Copyright (C) 2009-2011 by Andreas Lauser                               *
 *   Institute of Hydraulic Engineering                                      *
 *   University of Stuttgart, Germany                                        *
 *   email: <givenname>.<name>@iws.uni-stuttgart.de                          *
 *                                                                           *
 *   This program is free software: you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation, either version 2 of the License, or       *
 *   (at your option) any later version.                                     *
 *                                                                           *
 *   This program is distributed in the hope that it will be useful,         *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of          *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the           *
 *   GNU General Public License for more details.                            *
 *                                                                           *
 *   You should have received a copy of the GNU General Public License       *
 *   along with this program.  If not, see <http://www.gnu.org/licenses/>.   *
 *****************************************************************************/
/*!
 * \file
 *
 * \brief This file contains the parts of the local residual to
 *        calculate the heat flux in the fully coupled two-phase
 *        N-component model
 */
#ifndef DUMUX_MPNC_LOCAL_RESIDUAL_ENERGY_HH
#define DUMUX_MPNC_LOCAL_RESIDUAL_ENERGY_HH

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

    typedef typename GET_PROP_TYPE(TypeTag, PTAG(PrimaryVariables)) PrimaryVariables;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(VolumeVariables)) VolumeVariables;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(FluxVariables)) FluxVariables;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(ElementVolumeVariables)) ElementVolumeVariables;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(Scalar)) Scalar;

    enum { numPhases        = GET_PROP_VALUE(TypeTag, PTAG(NumPhases)) };
    enum { numComponents    = GET_PROP_VALUE(TypeTag, PTAG(NumComponents)) };
    typedef typename Dune::FieldVector<Scalar, numComponents>  ComponentVector;

public:
    static void computeStorage(PrimaryVariables &result,
                               const VolumeVariables &volVars)
    {
        // do nothing, we're isothermal!
    };


    static void addPhaseStorage(PrimaryVariables &storage,
                                const VolumeVariables &volVars,
                                int phaseIdx)
    {
        // do nothing, we're isothermal!
    };

    static void phaseEnthalpyFlux(PrimaryVariables &result,
                                  int phaseIdx,
                                  const PrimaryVariables &compMolFlux,
                                  const ElementVolumeVariables &volVars,
                                  const FluxVariables &fluxVars)
    {
        // do nothing, we're isothermal!
    };

    static void heatConduction(PrimaryVariables &result,
                               const ElementVolumeVariables &volVars,
                               const FluxVariables &fluxVars)
    {
        // do nothing, we're isothermal!
    };


    static void computeFlux(PrimaryVariables & flux,
                                const FluxVariables & fluxVars,
                                const ElementVolumeVariables & volVars,
                                const ComponentVector molarPhaseComponentValuesMassTransport[numPhases])
    {
        // do nothing, we're isothermal!
    }

    static void computeSource(PrimaryVariables &result,
                              const VolumeVariables &volVars)
    {
        // do nothing, we're isothermal!
    }
};

template <class TypeTag>
class MPNCLocalResidualEnergy<TypeTag, /*enableEnergy=*/true, /*kineticenergyTransfer=*/false>
{
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(FluidSystem)) FluidSystem;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(Scalar)) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(GridView)) GridView;

    typedef typename GET_PROP_TYPE(TypeTag, PTAG(PrimaryVariables)) PrimaryVariables;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(FluxVariables)) FluxVariables;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(VolumeVariables)) VolumeVariables;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(ElementVolumeVariables)) ElementVolumeVariables;

    typedef typename GET_PROP_TYPE(TypeTag, PTAG(MPNCIndices)) Indices;

    enum { dim = GridView::dimension };
    enum { numPhases = GET_PROP_VALUE(TypeTag, PTAG(NumPhases)) };
    enum { numComponents = GET_PROP_VALUE(TypeTag, PTAG(NumComponents)) };
    enum { conti0EqIdx = Indices::conti0EqIdx };
    enum { energyEqIdx = Indices::energyEqIdx };

    typedef Dune::FieldVector<Scalar, dim> Vector;
    typedef typename Dune::FieldVector<Scalar, numComponents> ComponentVector;
    typedef typename Dune::FieldMatrix<Scalar, numPhases, numComponents>    PhaseComponentMatrix;



public:
    static void computeStorage(PrimaryVariables &result,
                               const VolumeVariables &volVars)
    {
        result[energyEqIdx] = 0;

        // energy of the fluids
        for (int phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx) {
            addPhaseStorage(result, volVars, phaseIdx);
        }

        // heat stored in the rock matrix
        result[energyEqIdx] +=
            volVars.temperature()
            * volVars.soilDensity()
            * (1.0 - volVars.porosity())
            * volVars.heatCapacity();
    }

    static void addPhaseStorage(PrimaryVariables &storage,
                                const VolumeVariables &volVars,
                                int phaseIdx)
    {
        const typename VolumeVariables::FluidState &fs =
            volVars.fluidState();

        // energy of the fluid
        storage[energyEqIdx] +=
            fs.density(phaseIdx)
            * fs.internalEnergy(phaseIdx)
            * fs.saturation(phaseIdx)
            * volVars.porosity();
    }

    static void computeFlux(PrimaryVariables & flux,
                            const FluxVariables & fluxVars,
                            const ElementVolumeVariables & volVars,
                            const ComponentVector molarPhaseComponentValuesMassTransport[numPhases])
    {
        flux[energyEqIdx] = 0.0;

        // fluid phases transport enthalpy individually
        for(int phaseIdx=0; phaseIdx<numPhases; ++phaseIdx)
            computePhaseEnthalpyFlux(flux,
                                fluxVars,
                                volVars,
                                phaseIdx,
                                molarPhaseComponentValuesMassTransport[phaseIdx]);

        //conduction is treated lumped in this model
        computeHeatConduction(flux,
                             fluxVars,
                             volVars);
    }

    static void computePhaseEnthalpyFlux(PrimaryVariables & result,
                                         const FluxVariables & fluxVars,
                                         const ElementVolumeVariables & volVars,
                                         const int phaseIdx,
                                         const ComponentVector & molarComponentValuesMassTransport)
    {
        Scalar massFlux = 0;

        // calculate the mass flux in the phase i.e. make mass flux out of mole flux and add up the fluxes of a phase
        for (int compIdx = 0; compIdx < numComponents; ++compIdx)
            massFlux += molarComponentValuesMassTransport[compIdx]
                                    * FluidSystem::molarMass(compIdx);

        int upIdx = fluxVars.face().i;
        if (massFlux < 0) upIdx = fluxVars.face().j;

        // use the phase enthalpy of the upstream vertex to calculate
        // the enthalpy transport
        const VolumeVariables &up = volVars[upIdx];
        result[energyEqIdx] += up.fluidState().enthalpy(phaseIdx) * massFlux;
    }

    static void computeHeatConduction(PrimaryVariables & result,
                                    const FluxVariables & fluxVars,
                                    const ElementVolumeVariables & volVars)
    {
        //lumped heat conduction of the rock matrix and the fluid phases
        Scalar lumpedConductivity   = fluxVars.energyData().lambdaPm() ;
        Scalar temperatureGradientNormal  = fluxVars.energyData().temperatureGradientNormal() ;
        Scalar lumpedHeatConduction = - lumpedConductivity * temperatureGradientNormal ;
        result[energyEqIdx] += lumpedHeatConduction ;
    }



    static void computeSource(PrimaryVariables &result,
                              const VolumeVariables &volVars)
    {
        result[energyEqIdx] = 0.0;
    }
};





};

#endif // DUMUX_MPNC_ENERGY_HH
