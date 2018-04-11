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
 * \ingroup OnePNCMinTests
 * \brief Class for the evaluation of the reaction rate of Calciumoxide to Halciumhydroxide
 *
 * It contains simple and advanced reaction kinetics according to Nagel et al. (2014).
 */
#ifndef DUMUX_THERMOCHEM_REACTION_HH
#define DUMUX_THERMOCHEM_REACTION_HH

namespace Dumux {

/*!
 * \ingroup OnePNCMinTests
 * \brief Class for the evaluation of the reaction rate of Calciumoxide to Halciumhydroxide
 *
 * It contains simple and advanced reaction kinetics according to Nagel et al. (2014).
 */
class ThermoChemReaction
{
public:
    /*!
     * \brief evaluates the reaction kinetics (see Nagel et al. 2014)
     */
    template<class VolumeVariables>
    typename VolumeVariables::PrimaryVariables::value_type
    thermoChemReaction(const VolumeVariables &volVars) const
    {
        using FluidSystem = typename VolumeVariables::FluidSystem;
        using Scalar = typename VolumeVariables::PrimaryVariables::value_type;

        // calculate the equilibrium temperature Teq
        Scalar T= volVars.temperature();
        Scalar Teq = 0;

        Scalar moleFractionVapor = 1e-3;

        if(volVars.moleFraction(FluidSystem::gasPhaseIdx, FluidSystem::H2OIdx) > 1e-3)
            moleFractionVapor = volVars.moleFraction(FluidSystem::gasPhaseIdx, FluidSystem::H2OIdx);

        if(volVars.moleFraction(FluidSystem::gasPhaseIdx, FluidSystem::H2OIdx) >= 1.0) moleFractionVapor = 1;

        Scalar pV = volVars.pressure(FluidSystem::gasPhaseIdx) *moleFractionVapor;
        Scalar vaporPressure = pV*1.0e-5;
        Scalar pFactor = log(vaporPressure);

        Teq = -12845;
        Teq /= (pFactor - 16.508); //the equilibrium temperature

        Scalar realSolidDensityAverage = (volVars.precipitateVolumeFraction(FluidSystem::hPhaseIdx)*volVars.density(FluidSystem::hPhaseIdx)
                                        + volVars.precipitateVolumeFraction(FluidSystem::cPhaseIdx)*volVars.density(FluidSystem::cPhaseIdx))
                                        / (volVars.precipitateVolumeFraction(FluidSystem::hPhaseIdx)
                                        + volVars.precipitateVolumeFraction(FluidSystem::cPhaseIdx));

        if(realSolidDensityAverage <= volVars.density(FluidSystem::cPhaseIdx))
        {
            realSolidDensityAverage = volVars.density(FluidSystem::cPhaseIdx);
        }

        if(realSolidDensityAverage >= volVars.density(FluidSystem::hPhaseIdx))
        {
            realSolidDensityAverage = volVars.density(FluidSystem::hPhaseIdx);
        }

        Scalar qMass = 0.0;

        //discharge or hydration
        if (T < Teq){
            Scalar dXH_dt1 = 0.0;
            Scalar dXH_dt2 = 0.0;

            Scalar xH = (realSolidDensityAverage-volVars.density(FluidSystem::cPhaseIdx))/(volVars.density(FluidSystem::hPhaseIdx)- volVars.density(FluidSystem::cPhaseIdx));

            if(xH < 1.0e-5) {xH = 1.0e-5; }
            if(xH >=1 ) {xH = 1 - 1e-5; }

            Scalar R = 8.314 ; // [J/molK]
            Scalar peq = 1e5*exp( (-12845)/T + 16.508);

            if(peq >= pV) {peq=899954;}

            Scalar dXH_dt = 0;

            // reaction kinetics for T-Teq > 50 K
            if(Teq -T > 50.25){

               Scalar A =exp(-8.9486e4/(R*T));
               Scalar B = pow(((pV/peq)-1),0.83);
               Scalar D = 1-xH;
               Scalar C = pow((-log(D)),0.666);

               dXH_dt = 1.3945e4*A *B*3*D*C;

            }

            // reaction kinetics for T-Teq < 50 K
            if(Teq -T < 49.75){

                Scalar E = exp(5.3332e4/T);
                Scalar F = pow((pV*1e-5),6);
                Scalar G = (1-xH);

                dXH_dt = 1.004e-34 *E * F * G;

            }

            // linearization of the point of discontinuity
            if(Teq-T <=50.25 && Teq-T >=49.75){

               Scalar op = ((Teq-T)- 49.75)*2;
               Scalar A =exp(-8.9486e4/(R*T));
               Scalar B = pow(((pV/peq)-1),0.83);
               Scalar D = 1-xH;
               Scalar C = pow((-log(D)),0.666);
               dXH_dt1 = 1.3945e4*A *B*3*D*C;

               Scalar E = exp(5.3332e4/T);
               Scalar F = pow((pV*1e-5),6);
               Scalar G = (1-xH);
               dXH_dt2 = 1.004e-34 *E * F * G;

               dXH_dt = dXH_dt1*op +  dXH_dt2*(1-op);
            }

            // no reaction at equilibrium
            if(Teq -T <= 0)
                dXH_dt = 0;
            Scalar deltaRhoS = volVars.density(FluidSystem::hPhaseIdx) - volVars.density(FluidSystem::cPhaseIdx);
            qMass = dXH_dt*deltaRhoS;
        }

        return qMass;
    }


    /*!
     * \brief evaluates the simple chemical reaction kinetics (see Nagel et al. 2014)
     */
    template<class VolumeVariables>
    typename VolumeVariables::PrimaryVariables::value_type
    thermoChemReactionSimple(const VolumeVariables &volVars) const
    {
        using FluidSystem = typename VolumeVariables::FluidSystem;
        using Scalar = typename VolumeVariables::PrimaryVariables::value_type;

        // calculate the equilibrium temperature Teq
        Scalar T= volVars.temperature();
        Scalar Teq = 0;

        Scalar moleFractionVapor = 1e-3;

        if(volVars.moleFraction(FluidSystem::gasPhaseIdx, FluidSystem::H2OIdx) > 1e-3)
            moleFractionVapor = volVars.moleFraction(FluidSystem::gasPhaseIdx, FluidSystem::H2OIdx);

        if(volVars.moleFraction(FluidSystem::gasPhaseIdx, FluidSystem::H2OIdx) >= 1.0) moleFractionVapor = 1;

        Scalar pV = volVars.pressure(FluidSystem::gasPhaseIdx) *moleFractionVapor;
        Scalar vaporPressure = pV*1.0e-5;
        Scalar pFactor = log(vaporPressure);

        Teq = -12845;
        Teq /= (pFactor - 16.508); //the equilibrium temperature

        Scalar realSolidDensityAverage = (volVars.precipitateVolumeFraction(FluidSystem::hPhaseIdx)*volVars.density(FluidSystem::hPhaseIdx)
                                        + volVars.precipitateVolumeFraction(FluidSystem::cPhaseIdx)*volVars.density(FluidSystem::cPhaseIdx))
                                        / (volVars.precipitateVolumeFraction(FluidSystem::hPhaseIdx)
                                        + volVars.precipitateVolumeFraction(FluidSystem::cPhaseIdx));

        if(realSolidDensityAverage <= volVars.density(FluidSystem::cPhaseIdx))
        {
            realSolidDensityAverage = volVars.density(FluidSystem::cPhaseIdx);
        }

        if(realSolidDensityAverage >= volVars.density(FluidSystem::hPhaseIdx))
        {
            realSolidDensityAverage = volVars.density(FluidSystem::hPhaseIdx);
        }

        Scalar qMass = 0.0;

         //discharge or hydration
        if (T < Teq){
            Scalar massFracH2O_fPhase = volVars.massFraction(FluidSystem::gasPhaseIdx, FluidSystem::H2OIdx);
            Scalar krh = 0.2;

            Scalar rHydration = - massFracH2O_fPhase* (volVars.density(FluidSystem::hPhaseIdx)- realSolidDensityAverage)
                                                     * krh * (T-Teq)/ Teq;

            Scalar qMass =  rHydration;
        }

        //charge or hydration
        else if(T > Teq){

            Scalar krd = 0.05;

            Scalar rDehydration = -(volVars.density(FluidSystem::cPhaseIdx)- realSolidDensityAverage)
                                                     * krd * (Teq-T)/ Teq;

            qMass =  rDehydration;
        }

        if(Teq -T == 0) qMass = 0;

        return qMass;
    }

};

} // namespace Dumux

#endif
