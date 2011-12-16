// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   Copyright (C) 2009-2011 by Andreas Lauser                               *
 *   Copyright (C) 2011 by Philipp Nuske                                     *
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
#ifndef DUMUX_MPNC_LOCAL_RESIDUAL_MASS_HH
#define DUMUX_MPNC_LOCAL_RESIDUAL_MASS_HH

#include <dune/common/fvector.hh>
#include <dumux/material/constraintsolvers/compositionfromfugacities.hh>

#include "../diffusion/diffusion.hh"

namespace Dumux
{
/*!
 * \brief The mass conservation part of the Mp-Nc model.
 *
 * This is the class represents methods which are shared amongst all
 * mass conservation modules.
 */
template<class TypeTag>
class MPNCLocalResidualMassCommon
{
protected:
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(GridView)) GridView;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(Scalar)) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(VolumeVariables)) VolumeVariables;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(FluxVariables)) FluxVariables;

    enum { dim = GridView::dimension };
    enum { numPhases        = GET_PROP_VALUE(TypeTag, PTAG(NumPhases)) };
    enum { numComponents    = GET_PROP_VALUE(TypeTag, PTAG(NumComponents)) };
    enum { enableDiffusion  = GET_PROP_VALUE(TypeTag, PTAG(EnableDiffusion)) };
    enum { enableEnergy     = GET_PROP_VALUE(TypeTag, PTAG(EnableEnergy)) };
    enum { enableKineticEnergy     = GET_PROP_VALUE(TypeTag, PTAG(EnableKineticEnergy)) };

    typedef typename Dune::FieldVector<Scalar, dim> Vector;
    typedef typename Dune::FieldVector<Scalar, numComponents> ComponentVector;
    typedef MPNCDiffusion<TypeTag, enableDiffusion> Diffusion;

    typedef MPNCLocalResidualEnergy<TypeTag, enableEnergy, enableKineticEnergy> EnergyResid;

public:
    /*!
     * \brief Evaluate the amount moles within a sub-control volume in
     *        a phase.
     *
     * The result should be averaged over the volume.
     */
    static void computePhaseStorage(ComponentVector &result,
                                    const VolumeVariables &volVars,
                                    int phaseIdx)
    {
        // compute storage term of all components within all phases
        result = 0;
        for (int compIdx = 0; compIdx < numComponents; ++ compIdx) {
            result[compIdx] +=
                volVars.fluidState().saturation(phaseIdx)*
                volVars.fluidState().molarity(phaseIdx, compIdx);
        }

        result *= volVars.porosity();
    }

    /*!
     * \brief Evaluates the advective flux of all conservation
     *        quantities over a face of a subcontrol volume via a
     *        fluid phase.
     */
    static void computeAdvectivePhaseFlux(ComponentVector &phaseComponentValues,
                                          const FluxVariables &fluxVars,
                                          const int phaseIdx)
    {
        static bool enableSmoothUpwinding_ = GET_PARAM(TypeTag, bool, EnableSmoothUpwinding);

        Vector tmpVec;
        fluxVars.intrinsicPermeability().mv(fluxVars.potentialGrad(phaseIdx),
                                            tmpVec);

        // advective part is flux over face, therefore needs to be multiplied by normal vector
        // (length: area of face)
        Scalar normalFlux = - (tmpVec*fluxVars.face().normal);

        // data attached to upstream and the downstream vertices
        // of the current phase
        int upIdx = fluxVars.face().i;
        int dnIdx = fluxVars.face().j;
#ifndef NDEBUG
        if (!std::isfinite(normalFlux))
            DUNE_THROW(NumericalProblem, "Calculated non-finite normal flux");
#endif

        if (normalFlux < 0) std::swap(upIdx, dnIdx);
        const VolumeVariables &up = fluxVars.volVars(upIdx);
        const VolumeVariables &dn = fluxVars.volVars(dnIdx);

        ////////
        // advective fluxes of all components in the phase
        ////////
        for (int compIdx = 0; compIdx < numComponents; ++compIdx) {
            // add advective flux of current component in current
            // phase. we use full upwinding.

            if (enableSmoothUpwinding_) {
                Scalar mobUp = up.mobility(phaseIdx);
                Scalar cUp = up.fluidState().molarity(phaseIdx, compIdx);

                Scalar mobDn = dn.mobility(phaseIdx);
                Scalar cDn = dn.fluidState().molarity(phaseIdx, compIdx);

                Scalar mUp = mobUp*cUp;
                Scalar mDn = mobDn*cDn;
                Scalar m0 = Dumux::harmonicMean(mUp, mDn);

                Scalar x = std::abs(normalFlux);
                Scalar sign = (normalFlux < 0)?-1:1;

                // the direction from upstream to downstream
                //GlobalPosition delta = this->curElement_().geometry().corner(upIdx);
                //delta -= this->curElement_().geometry().corner(dnIdx);

                // approximate the mean viscosity at the face
                Scalar meanVisc = (up.fluidState().viscosity(phaseIdx) +
                                   dn.fluidState().viscosity(phaseIdx))/2;

                // put the mean viscosity and permeanbility in
                // relation to the viscosity of water at
                // approximatly 20 degrees Celsius.
                const Scalar pGradRef = 10; // [Pa/m]
                const Scalar muRef = 1e-3; // [Ns/m^2]
                const Scalar Kref = 1e-12; // [m^2] = approx 1 Darcy

                Scalar faceArea = fluxVars.face().normal.two_norm();
                Scalar eps = pGradRef * Kref * faceArea * meanVisc/muRef; // * (1e3/18e-3)/meanC;

                Scalar compFlux;
                if (x >= eps) {
                    // we only do tricks if x is below the epsilon
                    // value
                    compFlux = x*mUp;
                }
                else {
                    Scalar xPos[] = { 0, eps };
                    Scalar yPos[] = { 0, eps*mUp };
                    Spline<Scalar> sp2(xPos, yPos, m0, mUp);
                    compFlux = sp2.eval(x);
                }

                phaseComponentValues[compIdx] = sign*compFlux;
            }
            else {// !use smooth upwinding
                phaseComponentValues[compIdx] =
                    up.mobility(phaseIdx) *
                    up.fluidState().molarity(phaseIdx, compIdx) *
                    normalFlux;
            }
        }
    }


    /*!
     * \brief Evaluates the advective flux of all conservation
     *        quantities over a face of a subcontrol volume via a
     *        fluid phase.
     */
    static void computeDiffusivePhaseFlux(ComponentVector &flux,
                                          const FluxVariables &fluxVars,
                                          int phaseIdx)
    {
        if (!enableDiffusion) {
            flux = 0.0;
            return;
        }


        const VolumeVariables &volVarsI = fluxVars.volVars(fluxVars.face().i);
        const VolumeVariables &volVarsJ = fluxVars.volVars(fluxVars.face().j);
        if (volVarsI.fluidState().saturation(phaseIdx) < 1e-4 ||
            volVarsJ.fluidState().saturation(phaseIdx) < 1e-4)
        {
            return; // phase is not present in one of the finite volumes
        }

        // approximate the total concentration of the phase at the
        // integration point by the arithmetic mean of the
        // concentration of the sub-control volumes
        Scalar molarDensityAtIP;
        molarDensityAtIP = volVarsI.fluidState().molarDensity(phaseIdx);
        molarDensityAtIP += volVarsJ.fluidState().molarDensity(phaseIdx);
        molarDensityAtIP /= 2;

        Diffusion::flux(flux, phaseIdx, fluxVars, molarDensityAtIP);
    }
};

/*!
 * \brief The mass conservation part of the Mp-Nc model.
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

    typedef MPNCLocalResidualMassCommon<TypeTag> MassCommon;

    typedef typename GET_PROP_TYPE(TypeTag, PTAG(Scalar)) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(MPNCIndices)) Indices;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(PrimaryVariables)) PrimaryVariables;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(VolumeVariables)) VolumeVariables;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(ElementVolumeVariables)) ElementVolumeVariables;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(FluxVariables)) FluxVariables;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(FluidSystem)) FluidSystem;

    enum { numPhases        = GET_PROP_VALUE(TypeTag, PTAG(NumPhases)) };
    enum { numComponents    = GET_PROP_VALUE(TypeTag, PTAG(NumComponents)) };
    enum { conti0EqIdx      = Indices::conti0EqIdx };
    enum { numEnergyEqs     = Indices::NumPrimaryEnergyVars};

    typedef typename Dune::FieldVector<Scalar, numComponents> ComponentVector;

    enum { enableEnergy     = GET_PROP_VALUE(TypeTag, PTAG(EnableEnergy)) };
    enum { enableKineticEnergy     = GET_PROP_VALUE(TypeTag, PTAG(EnableKineticEnergy)) };
    typedef MPNCLocalResidualEnergy<TypeTag, enableEnergy, enableKineticEnergy> EnergyResid;

public:
    /*!
     * \brief Calculate the storage for all mass balance equations
     */
    static void computeStorage(PrimaryVariables &storage,
                               const VolumeVariables &volVars)
    {
        for (int compIdx = 0; compIdx < numComponents; ++compIdx)
            storage[conti0EqIdx + compIdx] = 0.0;

        for (int phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx) {
            addPhaseStorage(storage, volVars, phaseIdx);
        }
    };

    /*!
     * \brief Calculate the storage for all mass balance equations
     *        within a single fluid phase
     */
    static void addPhaseStorage(PrimaryVariables &storage,
                                const VolumeVariables &volVars,
                                int phaseIdx)
    {
        // calculate the component-wise mass storage
        ComponentVector phaseComponentValues;
        MassCommon::computePhaseStorage(phaseComponentValues,
                                        volVars,
                                        phaseIdx);

        // copy to the primary variables
        for (int compIdx = 0; compIdx < numComponents; ++compIdx)
            storage[conti0EqIdx + compIdx] += phaseComponentValues[compIdx];
    };

    /*!
     * \brief Calculate the storage for all mass balance equations
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
            for (int compIdx = 0; compIdx < numComponents; ++compIdx)
                flux[conti0EqIdx + compIdx] +=
                        phaseComponentValuesDiffusion[compIdx];

            // Right now I think that adding the two contributions individually into the flux is best for debugging and understanding.
            // The Energy module needs both contributions.
            phaseComponentValuesMassTransport[phaseIdx] = phaseComponentValuesDiffusion + phaseComponentValuesAdvection ;
            Valgrind::CheckDefined(flux);
        }

        // \todo
        //
        // The computeflux() of the Energy module needs a
        // component-wise flux (for the diffusive enthalpy transport)
        // It makes some sense calling energy from here, because energy
        // is carried by mass However, it is not really a clean
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
     */
    static void computeSource(PrimaryVariables &source,
                              const VolumeVariables &volVars)
    {
    static_assert(not enableKineticEnergy,
                      "In the case of kinetic energy transfer the advective energy transport between the phases has to be considered. "
                      "It is hard (technically) to say how much mass got transfered in the case of chemical equilibrium. "
                      "Therefore, kineticEnergy and no kinetic mass does not fit (yet).");

        // mass transfer is not considered in this mass module
        for (int compIdx = 0; compIdx < numComponents; ++compIdx)
            source[conti0EqIdx + compIdx] = 0.0;


        PrimaryVariables tmp(0);
        // Similar to the compute Flux, the energy residual needs to be called from the
        // mass residual.
        ComponentVector dummy[numPhases];
            for (int iDummy =0; iDummy <numPhases; ++iDummy)
                dummy[iDummy] = 0.;
        EnergyResid::computeSource(tmp,
                                   volVars,
                                   dummy);
        source += tmp;
        Valgrind::CheckDefined(source);
    };
};

} // end namepace

#endif
