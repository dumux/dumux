// $Id$

#ifndef DUNE_MULTICOMPONENTRELATIONS_HH
#define DUNE_MULTICOMPONENTRELATIONS_HH

#include <dumux/material/phaseproperties/phaseproperties_waterair.hh>

/**
 * \ingroup material
 * \defgroup properties
 * \author Klaus Mosthaf
 */

namespace Dune
{

/** \ingroup properties
 *  \brief base class for the computation of multicomponent relations
 */
class MultiComp
{

public:
    /*! \brief solubility of a component (water) in the non-wetting phase
     *  \param pressureN non-wetting phase pressure \f$ \left[ Pa \right] \f$
     *  \param temperature temperature
     *  \return the mass fraction of water in the non-wetting phase \f$ \left[ kg/kg \right] \f$
     */
    virtual double xWN (const double pressureN, double temperature=283.15) = 0;

    /*! \brief solubility of a component (air) in the wetting phase
     *  \param pressureN non-wetting phase pressure \f$ \left[ Pa \right] \f$
     *  \param temperature temperature [K]
     *  \return mass fraction of gas in the wetting phase \f$ \left[ kg/kg \right] \f$
     */
    virtual double xAW (const double pressureN, double temperature=283.15) = 0;

    /*! \brief solubility of a component (air) in the wetting phase
     *  \param pressureN non-wetting phase pressure \f$ \left[ Pa \right] \f$
     *  \param temperature temperature [K]
     *  \param X_NaCl mass fraction of salt dissolved in the wetting phase
     *  \return mass fraction of gas in the wetting phase \f$ \left[ kg/kg \right] \f$
     */

    virtual double xWNmolar (const double pressureN, double temperature=283.15) = 0;

    /*! \brief solubility of a component (air) in the wetting phase
     *  \param pressureN non-wetting phase pressure \f$ \left[ Pa \right] \f$
     *  \param temperature temperature [K]
     *  \return the mass fraction of water in the non-wetting phase \f$ \left[ mol/mol \right] \f$
     */
    virtual double xAWmolar (const double pressureN, double temperature=283.15) = 0;

    /** @brief Henry coefficient
     * \param pressureN non-wetting phase pressure \f$ \left[ Pa \right] \f$
     * @param temperature Temperature \f$ \left[ K \right] \f$
     * @return Henry coefficient \f$ \left[ 1/Pa \right] \f$
     */
    virtual double henry (double temperature=283.15) const = 0;

    /*! \brief Antoine equation for calculating the vapor pressure
     *  \param temperature temperature [K]
     *  \return vapor pressure [Pa]
     */
    virtual double vaporPressure (double temperature=283.15) const = 0;

    /*! \brief converts mole fractions into mass fractions
     *  \param massfrac mole fraction [mol/mol]
     *  \param phase phase for which a conversion is to be done [-]
     *  \return mass fraction [kg/kg]
     */
    virtual double convertMoleToMassFraction(double massfrac, int phase) const = 0;

    /*! \brief converts mass fractions into mole fractions
     *  \param massfrac mass fraction [kg/kg]
     *  \param phase phase for which a conversion is to be done [-]
     *  \return mole fraction [mol/mol]
     */
    virtual double convertMassToMoleFraction(double massfrac, int phase) const = 0;


    MultiComp(const Liquid_GL& wP = *(new Liq_WaterAir),
              const Gas_GL& nwP = *(new Gas_WaterAir))
    {     }

    virtual ~MultiComp()
    {}

};

/** \todo Please doc me! */

class CWaterAir : public MultiComp
{
public:
    /*! \brief equation for calculating the mass fraction in the nonwetting phase
     *    \param pressureN non-wetting phase pressure \f$ \left[ Pa \right] \f$
     *  \param temperature temperature \f$ \left[ K \right] \f$
     *  \return mass fraction \f$ \left[ kg/kg \right] \f$
     */
    double xWN (const double pressureN, const double temperature=283.15)
    {
        double pWSat;
        double result;

        pWSat = vaporPressure(temperature);
        result = pWSat / pressureN; //xWNmolar

        result = convertMoleToMassFraction(result, 1);

        return(result);
    }


    /*! \brief equation for calculating the mass fraction in the wetting phase
     *    \param pressureN non-wetting phase pressure \f$ \left[ Pa \right] \f$
     *  \param temperature temperature \f$ \left[ K \right] \f$
     *  \return mass fraction \f$ \left[ kg/kg \right] \f$
     */
    double xAW (const double pressureN, const double temperature=283.15)
    {
        double pan;
        double result;
        double hagw;

        pan = pressureN * (1-xWNmolar(pressureN,temperature));
        hagw = henry(temperature);
        result = pan * hagw; //xAWmolar

        result = convertMoleToMassFraction(result, 0);
        return(result);
    }

    /*! \brief equation for calculating the mole fraction in the nonwetting phase
     *    \param pressureN non-wetting phase pressure \f$ \left[ Pa \right] \f$
     *  \param temperature temperature \f$ \left[ K \right] \f$
     *  \return mole fraction \f$ \left[ mol/mol \right] \f$
     */
    double xWNmolar (const double pressureN, const double temperature=283.15)
    {
        double pWSat;
        double result;

        pWSat = vaporPressure(temperature);
        result = pWSat / pressureN;

        return(result);
    }

    /*! \brief equation for calculating the mole fraction in the wetting phase
     *    \param pressureN non-wetting phase pressure \f$ \left[ Pa \right] \f$
     *  \param temperature temperature \f$ \left[ K \right] \f$
     *  \return mole fraction \f$ \left[ mol/mol \right] \f$
     */
    double xAWmolar (const double pressureN, const double temperature=283.15)
    {
        double pan;
        double result;
        double hagw;

        pan = pressureN * (1-xWNmolar(pressureN,temperature));
        hagw = henry(temperature);
        result = pan * hagw;

        return(result);
    }

    /*! \brief equation for calculating the inverse Henry coefficient
     *  \param temperature temperature \f$ \left[ K \right] \f$
     *  \return Henry Coefficient \f$ \left[ 1/Pa \right] \f$
     */
    double henry(double temperature=283.15) const
    {
        double celsius = temperature - 273.15;
        double result = ((0.8942 + 1.47 * exp(-0.04394*celsius) )*1e-10);

        return (result); // [1/Pa]
    }

    /** @brief calculates vapor pressure
     *  @param temperature temperature \f$ \left[ K \right] \f$
     *  @return vapor pressure \f$ \left[ Pa \right] \f$
     */
    double vaporPressure(double temperature=283.15) const
    {
        const double constA = 8.19621;
        const double constB = 1730.63;
        const double constC = 233.426;

        double celsius;
        double exponent, psat;

        celsius = temperature - 273.15;

        exponent = constA - (constB / (celsius + constC));


        psat = pow (10.0, exponent) *100; //1mbar = 100Pa

        return(psat);
    }

    /** @brief converts mole fractions into mass fractions
     */
    double convertMoleToMassFraction(double molefrac, int phase) const
    {
        enum {wPhase = 0, nPhase = 1};

        double result;
        double molarMass1=0, molarMass2=0;

        if (phase == wPhase){
            molarMass1 = wettingPhase.molarMass_a();
            molarMass2 = wettingPhase.molarMass_w();
        }
        else if (phase == nPhase){
            molarMass1 = nonwettingPhase.molarMass_w();
            molarMass2 = nonwettingPhase.molarMass_a();
        }

        result = molefrac * molarMass1 / (molarMass1*molefrac + molarMass2*(1-molefrac));

        return (result);
    }

    /** @brief converts mass fractions into mole fractions
     */
    double convertMassToMoleFraction(double massfrac, int phase) const
    {
        enum {wPhase = 0, nPhase = 1};

        double result;
        double molarMass1 = 0, molarMass2 = 0;

        if (phase == wPhase){
            molarMass1 = wettingPhase.molarMass_a();
            molarMass2 = wettingPhase.molarMass_w();
        }
        else if (phase == nPhase){
            molarMass1 = nonwettingPhase.molarMass_w();
            molarMass2 = nonwettingPhase.molarMass_a();
        }

        result = massfrac * molarMass2 / (molarMass1*(1-massfrac) + molarMass2*massfrac);

        return (result);
    }

    CWaterAir(Liquid_GL& wP = *(new Liq_WaterAir),
              Gas_GL& nwP = *(new Gas_WaterAir))
        : MultiComp(wP, nwP), wettingPhase(wP), nonwettingPhase(nwP)
    {     }

    Liquid_GL& wettingPhase; //!< contains properties of the wetting phase
    Gas_GL& nonwettingPhase; //!< contains properties of the nonwetting phase
};

}
#endif
