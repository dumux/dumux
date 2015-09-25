/*****************************************************************************
 *   Copyright (C) 2009 by Benjamin Faigle                                   *
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
 * \brief Calculates the 2p2c phase state for compositional models.
 */
#ifndef DUMUX_DEC2P2C_FLUID_STATE_HH
#define DUMUX_DEC2P2C_FLUID_STATE_HH

#include <dumux/material/fluidstate.hh>
#include <dumux/decoupled/2p2c/2p2cproperties.hh>

namespace Dumux
{
/*!
 * \ingroup multiphysics multiphase
 * \brief Calculates the phase state from the primary variables in the
 *        sequential 2p2c model.
 *
 *        This boils down to so-called "flash calculation", in this case isothermal and isobaric.
 *
 *  \tparam TypeTag The property Type Tag
 */
template <class TypeTag>
class DecTwoPTwoCFluidState : public FluidState<typename GET_PROP_TYPE(TypeTag, PTAG(Scalar)),
                                           DecTwoPTwoCFluidState<TypeTag> >
{
    typedef DecTwoPTwoCFluidState<TypeTag> ThisType;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(Scalar))      Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(FluidSystem)) FluidSystem;

    typedef typename GET_PROP_TYPE(TypeTag, PTAG(TwoPIndices)) Indices;
    static const int pressureType = GET_PROP_VALUE(TypeTag, PTAG(PressureFormulation)); //!< gives kind of pressure used (\f$ 0 = p_w\f$, \f$ 1 = p_n\f$, \f$ 2 = p_{global}\f$)


    enum {
        wPhaseIdx = Indices::wPhaseIdx,
        nPhaseIdx = Indices::nPhaseIdx,

        wCompIdx = Indices::wPhaseIdx,
        nCompIdx = Indices::nPhaseIdx,
    };

public:
    enum {     numPhases = GET_PROP_VALUE(TypeTag, PTAG(NumPhases)),
            numComponents = GET_PROP_VALUE(TypeTag, PTAG(NumComponents)),};

public:
    /*!
     * \name flash calculation routines
     * Routines to determine the phase composition after the transport step.
     */
    //@{

    //! the update routine equals the 2p2c - concentration flash for constant p & t.
    /*!
     * Routine goes as follows:
     * - determination of the equilibrium constants from the fluid system
     * - determination of maximum solubilities (mole fractions) according to phase pressures
     * - comparison with Z1 to determine phase presence => phase mass fractions
     * - round off fluid properties
     * \param Z1 Feed mass fraction: Mass of comp1 per total mass \f$\mathrm{[-]}\f$
     * \param pw Pressure of the wetting phase \f$\mathrm{[Pa]}\f$
     * \param poro Porosity \f$\mathrm{[-]}\f$
     * \param temperature Temperature \f$\mathrm{[K]}\f$
     */
    void update(Scalar Z1, Scalar pw, Scalar poro, Scalar temperature)
    {
        if (pressureType == Indices::pressureGlobal)
        {
            DUNE_THROW(Dune::NotImplemented, "Pressure type not supported in fluidState!");
        }
        else
        phasePressure_[wPhaseIdx] = pw;
        phasePressure_[nPhaseIdx] = pw;        // as long as capillary pressure is neglected
        temperature_=temperature;


        //mole equilibrium ratios K for in case wPhase is reference phase
        double k1 = FluidSystem::activityCoeff(wPhaseIdx, wCompIdx, temperature_, phasePressure_[nPhaseIdx], *this)
                    / phasePressure_[nPhaseIdx];    // = p^wComp_vap / p_n
        double k2 = FluidSystem::activityCoeff(wPhaseIdx, nCompIdx, temperature_, phasePressure_[nPhaseIdx], *this)
                    / phasePressure_[nPhaseIdx];    // = H^nComp_w / p_n

        // get mole fraction from equilibrium konstants
        Scalar xw1 = (1. - k2) / (k1 -k2);
        Scalar xn1 = xw1 * k1;


        // transform mole to mass fractions
        massfrac_[wPhaseIdx][wCompIdx] = xw1 * FluidSystem::molarMass(wCompIdx)
            / ( xw1 * FluidSystem::molarMass(wCompIdx) + (1.-xw1) * FluidSystem::molarMass(nCompIdx) );
        massfrac_[nPhaseIdx][wCompIdx] = xn1 * FluidSystem::molarMass(wCompIdx)
            / ( xn1 * FluidSystem::molarMass(wCompIdx) + (1.-xn1) * FluidSystem::molarMass(nCompIdx) );


        //mass equilibrium ratios
        equilRatio_[nPhaseIdx][wCompIdx] = massfrac_[nPhaseIdx][wCompIdx] / massfrac_[wPhaseIdx][wCompIdx];     // = Xn1 / Xw1 = K1
        equilRatio_[nPhaseIdx][nCompIdx] = (1.-massfrac_[nPhaseIdx][wCompIdx])/ (1.-massfrac_[wPhaseIdx][wCompIdx]); // =(1.-Xn1) / (1.-Xw1)     = K2
        equilRatio_[wPhaseIdx][nCompIdx] = equilRatio_[wPhaseIdx][wCompIdx] = 1.;

        // phase fraction of nPhase [mass/totalmass]
        nu_[nPhaseIdx] = 0;

        // check if there is enough of component 1 to form a phase
        if (Z1 > massfrac_[nPhaseIdx][wCompIdx] && Z1 < massfrac_[wPhaseIdx][wCompIdx])
            nu_[nPhaseIdx] = -((equilRatio_[nPhaseIdx][wCompIdx]-1)*Z1 + (equilRatio_[nPhaseIdx][nCompIdx]-1)*(1-Z1)) / (equilRatio_[nPhaseIdx][wCompIdx]-1) / (equilRatio_[nPhaseIdx][nCompIdx] -1);
        else if (Z1 <= massfrac_[nPhaseIdx][wCompIdx]) // too little wComp to form a phase
        {
            nu_[nPhaseIdx] = 1; // only nPhase
            massfrac_[nPhaseIdx][wCompIdx] = Z1; // hence, assign complete mass soluted into nPhase
            massfrac_[wPhaseIdx][wCompIdx] = 1.;
        }
        else     // (Z1 >= Xw1) => no nPhase
        {
            nu_[nPhaseIdx] = 0; // no second phase
            massfrac_[wPhaseIdx][wCompIdx] = Z1;
            massfrac_[nPhaseIdx][wCompIdx] = 0.;
        }

        // complete array of mass fractions
        massfrac_[wPhaseIdx][nCompIdx] = 1. - massfrac_[wPhaseIdx][wCompIdx];
        massfrac_[nPhaseIdx][nCompIdx] = 1. - massfrac_[nPhaseIdx][wCompIdx];

        // complete phase mass fractions
        nu_[wPhaseIdx] = 1. - nu_[nPhaseIdx];

//        // let FluidSystem compute partial pressures
//        FluidSystem::computePartialPressures(temperature_, phasePressure_[nPhaseIdx], *this);

        // get densities with correct composition
        density_[wPhaseIdx] = FluidSystem::phaseDensity(wPhaseIdx,
                                                        temperature,
                                                        phasePressure_[wPhaseIdx],
                                                        *this);
        density_[nPhaseIdx] = FluidSystem::phaseDensity(nPhaseIdx,
                                                        temperature,
                                                        phasePressure_[nPhaseIdx],
                                                        *this);

        Sw_ = (nu_[wPhaseIdx]) / density_[wPhaseIdx];
        Sw_ /= (nu_[wPhaseIdx]/density_[wPhaseIdx] + nu_[nPhaseIdx]/density_[nPhaseIdx]);

        massConcentration_[wCompIdx] =
                poro * (massfrac_[wPhaseIdx][wCompIdx] * Sw_ * density_[wPhaseIdx]
                        + massfrac_[nPhaseIdx][wCompIdx] * (1.-Sw_) * density_[nPhaseIdx]);
        massConcentration_[nCompIdx] =
                poro * (massfrac_[wPhaseIdx][nCompIdx] * Sw_ * density_[wPhaseIdx]
                        + massfrac_[nPhaseIdx][nCompIdx] * (1-Sw_) * density_[nPhaseIdx]);
    }

    //! a flash routine for 2p2c if the saturation instead of total concentration is known.
    /*!
     * Routine goes as follows:
     * - determination of the equilibrium constants from the fluid system
     * - determination of maximum solubilities (mole fractions) according to phase pressures
     * - round off fluid properties
     * \param sat Saturation of phase 1 \f$\mathrm{[-]}\f$
     * \param pw Pressure of the wetting phase \f$\mathrm{[Pa]}\f$
     * \param poro Porosity \f$\mathrm{[-]}\f$
     * \param temperature Temperature \f$\mathrm{[K]}\f$
     */
    void satFlash(Scalar sat, Scalar pw, Scalar poro, Scalar temperature)
    {
        if (pressureType == Indices::pressureGlobal)
        {
            DUNE_THROW(Dune::NotImplemented, "Pressure type not supported in fluidState!");
        }
//        else  if (sat <= 0. || sat >= 1.)
//            Dune::dinfo << "saturation initial and boundary conditions set to zero or one!"
//                << " assuming fully saturated compositional conditions" << std::endl;

        // assign values
        Sw_ = sat;
        phasePressure_[wPhaseIdx] = pw;
        phasePressure_[nPhaseIdx] = pw;        // as long as capillary pressure is neglected
        temperature_=temperature;
        nu_[nPhaseIdx] = nu_[wPhaseIdx] = NAN;    //in contrast to the standard update() method, satflash() does not calculate nu.


        //mole equilibrium ratios K for in case wPhase is reference phase
        double k1 = FluidSystem::activityCoeff(wPhaseIdx, wCompIdx, temperature_, phasePressure_[nPhaseIdx], *this)
                    / phasePressure_[nPhaseIdx];
        double k2 = FluidSystem::activityCoeff(wPhaseIdx, nCompIdx, temperature_, phasePressure_[nPhaseIdx], *this)
                    / phasePressure_[nPhaseIdx];

        // get mole fraction from equilibrium konstants
        Scalar xw1 = (1. - k2) / (k1 -k2);
        Scalar xn1 = xw1 * k1;


        // transform mole to mass fractions
        massfrac_[wPhaseIdx][wCompIdx] = xw1 * FluidSystem::molarMass(wCompIdx)
            / ( xw1 * FluidSystem::molarMass(wCompIdx) + (1.-xw1) * FluidSystem::molarMass(nCompIdx) );
        massfrac_[nPhaseIdx][wCompIdx] = xn1 * FluidSystem::molarMass(wCompIdx)
            / ( xn1 * FluidSystem::molarMass(wCompIdx) + (1.-xn1) * FluidSystem::molarMass(nCompIdx) );

        // complete array of mass fractions
        massfrac_[wPhaseIdx][nCompIdx] = 1. - massfrac_[wPhaseIdx][wCompIdx];
        massfrac_[nPhaseIdx][nCompIdx] = 1. - massfrac_[nPhaseIdx][wCompIdx];

        //mass equilibrium ratios
        equilRatio_[nPhaseIdx][wCompIdx] = massfrac_[nPhaseIdx][wCompIdx] / massfrac_[wPhaseIdx][wCompIdx];     // = Xn1 / Xw1 = K1
        equilRatio_[nPhaseIdx][nCompIdx] = (1.-massfrac_[nPhaseIdx][wCompIdx])/ (1.-massfrac_[wPhaseIdx][wCompIdx]); // =(1.-Xn1) / (1.-Xw1)     = K2
        equilRatio_[wPhaseIdx][nCompIdx] = equilRatio_[wPhaseIdx][wCompIdx] = 1.;

        // get densities with correct composition
        density_[wPhaseIdx] = FluidSystem::phaseDensity(wPhaseIdx,
                                                        temperature,
                                                        phasePressure_[wPhaseIdx],
                                                        *this);
        density_[nPhaseIdx] = FluidSystem::phaseDensity(nPhaseIdx,
                                                        temperature,
                                                        phasePressure_[nPhaseIdx],
                                                        *this);

        massConcentration_[wCompIdx] =
                poro * (massfrac_[wPhaseIdx][wCompIdx] * Sw_ * density_[wPhaseIdx]
                        + massfrac_[nPhaseIdx][wCompIdx] * (1.-Sw_) * density_[nPhaseIdx]);
        massConcentration_[nCompIdx] =
                poro * (massfrac_[wPhaseIdx][nCompIdx] * Sw_ * density_[wPhaseIdx]
                        + massfrac_[nPhaseIdx][nCompIdx] * (1-Sw_) * density_[nPhaseIdx]);
    }
    //@}
    /*!
     * \name acess functions
     */
    //@{
    /*!
     * \brief Returns the saturation of a phase.
     *
     * \param phaseIdx the index of the phase
     */
    Scalar saturation(int phaseIdx) const
    {
        if (phaseIdx == wPhaseIdx)
            return Sw_;
        else
            return Scalar(1.0) - Sw_;
    };

    /*!
     * \brief Returns the mass fraction of a component in a phase.
     *
     * \param phaseIdx the index of the phase
     * \param compIdx the index of the component
     */
    Scalar massFrac(int phaseIdx, int compIdx) const
    {
        return massfrac_[phaseIdx][compIdx];

    }

    /*!
     * \brief Returns the molar fraction of a component in a fluid phase.
     *
     * \param phaseIdx the index of the phase
     * \param compIdx the index of the component
     */
    Scalar moleFrac(int phaseIdx, int compIdx) const
    {
        // as the moass fractions are calculated, it is used to determine the mole fractions
        double moleFrac_ = ( massfrac_[phaseIdx][compIdx] / FluidSystem::molarMass(compIdx) );    // = moles of compIdx
        moleFrac_ /= ( massfrac_[phaseIdx][wCompIdx] / FluidSystem::molarMass(wCompIdx)
                       + massfrac_[phaseIdx][nCompIdx] / FluidSystem::molarMass(nCompIdx) );    // /= total moles in phase
        return moleFrac_;
    }

    /*!
     * \brief Returns the total mass concentration of a component \f$\mathrm{[kg/m^3]}\f$.
     *
     * This is equivalent to the sum of the component concentrations for all
     * phases multiplied with the phase density.
     *
     * \param compIdx the index of the component
     */
    Scalar massConcentration(int compIdx) const
    {
        return massConcentration_[compIdx];
    };


    /*!
     * \brief Returns the density of a phase \f$\mathrm{[kg/m^3]}\f$.
     *
     * \param phaseIdx the index of the phase
     */
    Scalar density(int phaseIdx) const
    { return density_[phaseIdx]; }

    /*!
     * \brief Return the partial pressure of a component in the gas phase.
     *
     * For an ideal gas, this means \f$ R*T*c \f$.
     * Unit: \f$\mathrm{[Pa] = [N/m^2]}\f$
     *
     * \param componentIdx the index of the component
     */
    Scalar partialPressure(int componentIdx) const
    {
        if(componentIdx==nCompIdx)
            return phasePressure_[nPhaseIdx]*moleFrac(nPhaseIdx, nCompIdx);
        if(componentIdx == wCompIdx)
            return phasePressure_[nPhaseIdx]*moleFrac(nPhaseIdx, wCompIdx);
        else
            DUNE_THROW(Dune::NotImplemented, "component not found in fluidState!");
    }

    /*!
     * \brief Returns the pressure of a fluid phase \f$\mathrm{[Pa]}\f$.
     *
     * \param phaseIdx the index of the phase
     */
    Scalar phasePressure(int phaseIdx) const
    { return phasePressure_[phaseIdx]; }

    /*!
     * \brief Returns the capillary pressure \f$\mathrm{[Pa]}\f$
     */
    Scalar capillaryPressure() const
    { return phasePressure_[nPhaseIdx] - phasePressure_[wPhaseIdx]; }

    /*!
     * \brief Returns the temperature of the fluids \f$\mathrm{[K]}\f$.
     *
     * Note that we assume thermodynamic equilibrium, so all fluids
     * and the rock matrix exhibit the same temperature.
     */
    Scalar temperature() const
    { return temperature_; };

    /*!
     * \brief Returns the phase mass fraction. phase mass per total mass \f$\mathrm{[kg/kg]}\f$.
     *
     * \param phaseIdx the index of the phase
     */
    Scalar phaseMassFraction(int phaseIdx)
    {
        if (std::isnan(nu_[phaseIdx]))  //in contrast to the standard update() method, satflash() does not calculate nu.
        {
            nu_[wPhaseIdx] = Sw_ * density_[wPhaseIdx] / (Sw_ * density_[wPhaseIdx] + (1. - Sw_) * density_[nPhaseIdx]);
            nu_[nPhaseIdx] = 1. - nu_[wPhaseIdx];
            return nu_[phaseIdx];
        }
        else
            return nu_[phaseIdx];
    }
    //@}

private:
    Scalar density_[numPhases];
    Scalar massConcentration_[numComponents];
    Scalar massfrac_[numPhases][numComponents];
    Scalar equilRatio_[numPhases][numComponents];
    Scalar phasePressure_[numPhases];
    Scalar temperature_;
    Scalar Sw_;
    Scalar nu_[numPhases]; //phase mass fraction
};

} // end namepace

#endif
