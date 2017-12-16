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
 * \brief The mass conservation part of the MpNc model.
 */
#ifndef DUMUX_MPNC_LOCAL_RESIDUAL_MASS_HH
#define DUMUX_MPNC_LOCAL_RESIDUAL_MASS_HH

#include <dune/common/fvector.hh>

#include <dumux/porousmediumflow/mpnc/implicit/properties.hh>
#include <dumux/material/constraintsolvers/compositionfromfugacities.hh>
#include <dumux/common/math.hh>
#include <dumux/common/spline.hh>

#include "../diffusion/diffusion.hh"
#include "../energy/localresidual.hh"

namespace Dumux
{
/*!
 * \brief This class represents methods which are shared amongst all
 *        mass conservation modules.
 */
template<class TypeTag>
class MPNCLocalResidualMassCommon
{
protected:
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using VolumeVariables = typename GET_PROP_TYPE(TypeTag, VolumeVariables);
    using FluxVariables = typename GET_PROP_TYPE(TypeTag, FluxVariables);

    enum { numPhases        = GET_PROP_VALUE(TypeTag, NumPhases) };
    enum { numComponents    = GET_PROP_VALUE(TypeTag, NumComponents) };
    enum { enableDiffusion  = GET_PROP_VALUE(TypeTag, EnableDiffusion) };
    enum { enableEnergy     = GET_PROP_VALUE(TypeTag, EnableEnergy) };
    enum { numEnergyEquations     = GET_PROP_VALUE(TypeTag, NumEnergyEquations) };

    using ComponentVector = Dune::FieldVector<Scalar, numComponents>;
    using Diffusion = MPNCDiffusion<TypeTag, enableDiffusion>;
    using EnergyResid = MPNCLocalResidualEnergy<TypeTag, enableEnergy, numEnergyEquations>;

public:
    /*!
     * \brief Evaluate the amount moles within a sub-control volume in
     *        a phase.
     *
     *    \param storage The mass of the component within the sub-control volume
     *    \param volVars The volume variables
     *    \param phaseIdx phaseIdx The index of the fluid phase
     *
     * The result should be averaged over the volume.
     */
    static void computePhaseStorage(ComponentVector &storage,
                                    const VolumeVariables &volVars,
                                    const unsigned int phaseIdx)
    {
        // compute storage term of all components within all phases
        storage = 0;
        for (int compIdx = 0; compIdx < numComponents; ++ compIdx) {
            storage[compIdx] +=
                volVars.saturation(phaseIdx)*
                volVars.molarity(phaseIdx, compIdx);
#ifndef NDEBUG
using std::isfinite;
if (!isfinite(storage[compIdx]))
    DUNE_THROW(NumericalProblem, "Calculated non-finite storage");
#endif
        }

        storage *= volVars.porosity();

    }

    /*!
     * \brief Evaluates the advective flux of all conservation
     *        quantities over a face of a subcontrol volume via a
     *        fluid phase.
     *
     *        \param flux The flux over the SCV (sub-control-volume) face for each component
     *        \param fluxVars The flux Variables
     *        \param phaseIdx phaseIdx The index of the fluid phase
     */
    static void computeAdvectivePhaseFlux(ComponentVector &flux,
                                          const FluxVariables &fluxVars,
                                          const unsigned int phaseIdx)
    {

        const Scalar volumeFlux =  fluxVars.volumeFlux(phaseIdx) ;


        // retrieve the upwind weight for the mass conservation equations. Use the value
        // specified via the property system as default, and overwrite
        // it by the run-time parameter from the Dune::ParameterTree
        const Scalar massUpwindWeight = GET_PARAM_FROM_GROUP(TypeTag, Scalar, Implicit, MassUpwindWeight);

        static bool enableSmoothUpwinding_ = GET_PARAM_FROM_GROUP(TypeTag, bool, Implicit, EnableSmoothUpwinding);

        // data attached to upstream and the downstream vertices
        // of the current phase
        unsigned int upIdx = fluxVars.upstreamIdx(phaseIdx);
        unsigned int dnIdx = fluxVars.downstreamIdx(phaseIdx);

        const VolumeVariables &up = fluxVars.volVars(upIdx);
        const VolumeVariables &dn = fluxVars.volVars(dnIdx);
using std::isfinite;
if (!isfinite(volumeFlux))
    DUNE_THROW(NumericalProblem, "Calculated non-finite normal flux in phase" << phaseIdx);
        ////////
        // advective fluxes of all components in the phase
        ////////
        for (unsigned int compIdx = 0; compIdx < numComponents; ++compIdx) {
            // add advective flux of current component in current
            // phase. we use full upwinding.
            if (enableSmoothUpwinding_) {
                const Scalar kGradPNormal   = fluxVars.kGradPNormal(phaseIdx);
                const Scalar mobUp          = up.mobility(phaseIdx);
                const Scalar conUp          = up.molarity(phaseIdx, compIdx);

                const Scalar mobDn  = dn.mobility(phaseIdx);
                const Scalar conDn  = dn.molarity(phaseIdx, compIdx);

                const Scalar mobConUp   = mobUp*conUp;
                const Scalar mobConDn   = mobDn*conDn;
                const Scalar meanMobCon = harmonicMean(mobConUp, mobConDn);

                using std::abs;
                const Scalar x      = abs(kGradPNormal);
                const Scalar sign   = (kGradPNormal > 0)?-1:1;

                // approximate the mean viscosity at the face
                const Scalar meanVisc = (up.viscosity(phaseIdx) +
                                   dn.viscosity(phaseIdx))/2;

                // put the mean viscosity and permeanbility in
                // relation to the viscosity of water at
                // approximatly 20 degrees Celsius.
                const Scalar pGradRef   = 10; // [Pa/m]
                const Scalar muRef      = 1e-3; // [Ns/m^2]
                const Scalar Kref       = 1e-12; // [m^2] = approx 1 Darcy

                const Scalar faceArea   = fluxVars.face().normal.two_norm();
                const Scalar eps        = pGradRef * Kref * faceArea * meanVisc/muRef; // * (1e3/18e-3)/meanC;

                Scalar compFlux;
                if (x >= eps) {
                    // we only do tricks if x is below the epsilon
                    // value
                    compFlux = x*mobConUp;
                }
                else {
                    const Scalar xPos[] = { 0, eps };
                    const Scalar yPos[] = { 0, eps*mobConUp };
                    const Spline<Scalar> sp2(xPos, yPos, meanMobCon, mobConUp);
                    compFlux = sp2.eval(x);
                }
                #ifndef NDEBUG
                using std::isfinite;
                if (!isfinite(compFlux))
                    DUNE_THROW(NumericalProblem, "Calculated non-finite normal flux in smooth upwinding");
                #endif

                flux[compIdx] = sign*compFlux;
            }
            else
            {// not use smooth upwinding
                flux[compIdx] =
                        volumeFlux *
                        ((     massUpwindWeight)*up.molarity(phaseIdx, compIdx)
                                +
                        (  1. - massUpwindWeight)*dn.molarity(phaseIdx, compIdx) );
                        using std::isfinite;
                        if (!isfinite(flux[compIdx]))
                            DUNE_THROW(NumericalProblem, "Calculated non-finite normal flux in phase " <<  phaseIdx << " comp " << compIdx << "T: "<<  up.temperature(phaseIdx) << "S "<<up.saturation(phaseIdx)  ) ;
            }
        }
    }


    /*!
     * \brief Evaluates the advective flux of all conservation
     *        quantities over a face of a subcontrol volume via a
     *        fluid phase.
     *
     *        \param flux The flux over the SCV (sub-control-volume) face for each component
     *        \param fluxVars The flux Variables
     *        \param phaseIdx phaseIdx The index of the fluid phase
     */
    static void computeDiffusivePhaseFlux(ComponentVector &flux,
                                          const FluxVariables &fluxVars,
                                          const unsigned int phaseIdx)
    {
        if (!enableDiffusion) {
            flux = 0.0;
            return;
        }


        const VolumeVariables &volVarsI = fluxVars.volVars(fluxVars.face().i);
        const VolumeVariables &volVarsJ = fluxVars.volVars(fluxVars.face().j);
        if (volVarsI.saturation(phaseIdx) < 1e-4 ||
            volVarsJ.saturation(phaseIdx) < 1e-4)
        {
            return; // phase is not present in one of the finite volumes
        }

        // approximate the total concentration of the phase at the
        // integration point by the arithmetic mean of the
        // concentration of the sub-control volumes
        Scalar molarDensityAtIP;
        molarDensityAtIP = volVarsI.molarDensity(phaseIdx);
        molarDensityAtIP += volVarsJ.molarDensity(phaseIdx);
        molarDensityAtIP /= 2;

        Diffusion::flux(flux, phaseIdx, fluxVars, molarDensityAtIP);
    }
};

/*!
 * \brief The mass conservation part of the MpNc model.
 *
 * This is the specialization for the case where kinetic mass transfer
 * is _not_ considered.
 */
template<class TypeTag, bool enableKinetic /*=false*/>
class MPNCLocalResidualMass
{
    static_assert(!enableKinetic,
                  "No kinetic mass transfer module included, "
                  "but kinetic mass transfer enabled.");

    using MassCommon = MPNCLocalResidualMassCommon<TypeTag>;
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using Indices = typename GET_PROP_TYPE(TypeTag, Indices);
    using PrimaryVariables = typename GET_PROP_TYPE(TypeTag, PrimaryVariables);
    using VolumeVariables = typename GET_PROP_TYPE(TypeTag, VolumeVariables);
    using ElementVolumeVariables = typename GET_PROP_TYPE(TypeTag, ElementVolumeVariables);
    using FluxVariables = typename GET_PROP_TYPE(TypeTag, FluxVariables);

    enum { numPhases        = GET_PROP_VALUE(TypeTag, NumPhases) };
    enum { numComponents    = GET_PROP_VALUE(TypeTag, NumComponents) };
    enum { conti0EqIdx      = Indices::conti0EqIdx };
    enum { enableEnergy     = GET_PROP_VALUE(TypeTag, EnableEnergy) };
    enum { numEnergyEquations     = GET_PROP_VALUE(TypeTag, NumEnergyEquations) };

    using ComponentVector = Dune::FieldVector<Scalar, numComponents>;
    using EnergyResid = MPNCLocalResidualEnergy<TypeTag, enableEnergy, numEnergyEquations>;

public:
    /*!
     * \brief Calculate the storage for all mass balance equations
     *
     *    \param storage The mass of the component within the sub-control volume
     *    \param volVars The volume variables
     */
    static void computeStorage(PrimaryVariables &storage,
                               const VolumeVariables &volVars)
    {
        for (int compIdx = 0; compIdx < numComponents; ++compIdx)
            storage[conti0EqIdx + compIdx] = 0.0;

        for (int phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx) {
            addPhaseStorage(storage, volVars, phaseIdx);
        }
    }

    /*!
     * \brief Calculate the storage for all mass balance equations
     *        within a single fluid phase
     *
     *    \param storage The mass of the component within the sub-control volume
     *    \param volVars The volume variables
     *    \param phaseIdx phaseIdx The index of the fluid phase
     */
    static void addPhaseStorage(PrimaryVariables &storage,
                                const VolumeVariables &volVars,
                                const unsigned int phaseIdx)
    {
        // calculate the component-wise mass storage
        ComponentVector phaseComponentValues;
        MassCommon::computePhaseStorage(phaseComponentValues,
                                        volVars,
                                        phaseIdx);

        // copy to the primary variables
        for (int compIdx = 0; compIdx < numComponents; ++compIdx)
            storage[conti0EqIdx + compIdx] += phaseComponentValues[compIdx];
    }

    /*!
     * \brief Calculate the storage for all mass balance equations
     *
     *        \param flux The flux over the SCV (sub-control-volume) face for each component
     *        \param fluxVars The flux Variables
     *        \param elemVolVars The volume variables of the current element
     */
    static void computeFlux(PrimaryVariables &flux,
                            const FluxVariables &fluxVars,
                            const ElementVolumeVariables & elemVolVars)
    {
        for (int compIdx = 0; compIdx < numComponents; ++compIdx)
            flux[conti0EqIdx + compIdx] = 0.0;

        ComponentVector phaseComponentValuesAdvection(0.);
        ComponentVector phaseComponentValuesDiffusion(0.);
        ComponentVector phaseComponentValuesMassTransport[numPhases]; // what goes into the energy module

        for (int phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx) {
            MassCommon::computeAdvectivePhaseFlux(phaseComponentValuesAdvection, fluxVars, phaseIdx);
            Valgrind::CheckDefined(phaseComponentValuesAdvection);
            for (int compIdx = 0; compIdx < numComponents; ++compIdx)
                flux[conti0EqIdx + compIdx] +=
                        phaseComponentValuesAdvection[compIdx];

            MassCommon::computeDiffusivePhaseFlux(phaseComponentValuesDiffusion, fluxVars, phaseIdx);
            Valgrind::CheckDefined(phaseComponentValuesDiffusion);
            for (int compIdx = 0; compIdx < numComponents; ++compIdx){
                flux[conti0EqIdx + compIdx] +=
                        phaseComponentValuesDiffusion[compIdx];
            }



            // Right now I think that adding the two contributions individually into the flux is best for debugging and understanding.
            // The Energy module needs both contributions.
            phaseComponentValuesMassTransport[phaseIdx] = phaseComponentValuesDiffusion + phaseComponentValuesAdvection ;

            Valgrind::CheckDefined(flux);
        }
//        std::cout<< "KPN Mass: flux: " << flux << endl;


        // The computeflux() of the Energy module needs a
        // component-wise flux (for the diffusive enthalpy transport)
        // It makes some sense calling energy from here, because energy
        // is carried by mass. However, it is not really a clean
        // solution.

        // energy transport in fluid phases
        EnergyResid::computeFlux(flux,
                                 fluxVars,
                                 elemVolVars,
                                 phaseComponentValuesMassTransport);
        Valgrind::CheckDefined(flux);
    }



    /*!
     * \brief Calculate the source terms for all mass balance
     *        equations
     *
     *         \param source The source/sink in the sub-control volume for each component
     *         \param volVars the volume variables
     */
    static void computeSource(PrimaryVariables &source,
                              const VolumeVariables &volVars)
    {
//      static_assert(not enableKineticEnergy, // enableKinetic is disabled, in this specialization
//              "In the case of kinetic energy transfer the advective energy transport between the phases has to be considered. "
//              "It is hard (technically) to say how much mass got transfered in the case of chemical equilibrium. "
//              "Therefore, kineticEnergy and no kinetic mass does not fit (yet).");


        // mass transfer is not considered in this mass module
        for (int compIdx = 0; compIdx < numComponents; ++compIdx)
            source[conti0EqIdx + compIdx] = 0.0;


        PrimaryVariables tmpPriVars(0);
        // Similar to the compute Flux, the energy residual needs to be called from the
        // mass residual.
        ComponentVector dummy[numPhases];
            for (int iDummy =0; iDummy <numPhases; ++iDummy)
                dummy[iDummy] = 0.;
        EnergyResid::computeSource(tmpPriVars,
                                   volVars,
                                   dummy);
        source += tmpPriVars;
        Valgrind::CheckDefined(source);
    }
};

} // end namespace

#endif
