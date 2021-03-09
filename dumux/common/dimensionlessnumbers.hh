// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   See the file COPYING for full copying permissions.                      *
 *                                                                           *
 *   This program is free software: you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation, either version 3 of the License, or       *
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
 * \ingroup Common
 * \brief Collection of functions, calculating dimensionless numbers.
 *
 * All the input to the dimensionless numbers has to be provided as function arguments.
 * Rendering this collection generic in the sense that it can be used by any model.
 */
#ifndef DIMENSIONLESS_NUMBERS_HH
#define DIMENSIONLESS_NUMBERS_HH

#include <cmath>
#include <iostream>

#include <dune/common/exceptions.hh>
#include <dune/common/math.hh>

namespace Dumux {

/*!
 * \brief A container for possible values of the property for selecting which nusselt parametrization to choose.
 *        The actual value is set vie the property NusseltFormulation
 */
enum class NusseltFormulation
{
    dittusBoelter, WakaoKaguei, VDI
};

/*!
 * \brief A container for possible values of the property for selecting which sherwood parametrization to choose.
 *        The actual value is set vie the property SherwoodFormulation
 */
enum class SherwoodFormulation
{
    WakaoKaguei
};

/*!
 * \ingroup Common
 * \brief Collection of functions which calculate dimensionless numbers.
 * Each number has it's own function.
 * All the parameters for the calculation have to be handed over.
 * Rendering this collection generic in the sense that it can be used by any model.
 */
template <class Scalar>
class DimensionlessNumbers
{

public:
/*!
 * \brief   Calculate the Reynolds Number [-] (Re).
 *
 * The Reynolds number is a measure for the relation of inertial to viscous forces.
 * The bigger the value, the more important inertial (as compared to viscous) effects become.
 * According to Bear [Dynamics of fluids in porous media (1972)] Darcy's law is valid for Re<1.
 *
 * Source for Reynolds number definition: http://en.wikipedia.org/wiki/Reynolds_number
 *
 * \param darcyMagVelocity      The absolute value of the darcy velocity. In the context of box models this
 *                              leads to a problem: the velocities are defined on the faces while other things (storage, sources, output)
 *                              are defined for the volume/vertex. Therefore, some sort of decision needs to be made which velocity to put
 *                              into this function (e.g.: face-area weighted average). [m/s]
 * \param charcteristicLength   Typically, in the context of porous media flow, the mean grain size is taken as the characteristic length
 *                              for calculation of Re. [m]
 * \param kinematicViscosity    Is defined as the dynamic (or absolute) viscos  ity divided by the density.
 *                              http://en.wikipedia.org/wiki/Viscosity#Dynamic_viscosity. [m^2/s]
 *
 * \return                      The Reynolds Number as calculated from the input parameters
 */
static Scalar reynoldsNumber(const Scalar darcyMagVelocity,
                             const Scalar charcteristicLength,
                             const Scalar kinematicViscosity)
{
    return darcyMagVelocity * charcteristicLength / kinematicViscosity ;
}

/*!
 * \brief   Calculate the Prandtl Number [-] (Pr).
 *
 *          The Prandtl Number is a measure for the relation of viscosity and thermal diffusivity (temperaturleitfaehigkeit).
 *
 *          It is defined as
 *          \f[
 *          \textnormal{Pr}= \frac{\nu}{\alpha} = \frac{c_p \mu}{\lambda}\, ,
 *          \f]
 *          with kinematic viscosity\f$\nu\f$, thermal diffusivity \f$\alpha\f$, heat capacity \f$c_p\f$,
 *          dynamic viscosity \f$\mu\f$ and thermal conductivity \f$\lambda\f$.
 *          Therefore, Pr is a material specific property (i.e.: not a function of flow directly
 *          but only of temperature, pressure and fluid).
 *
 *          source for Prandtl number definition: http://en.wikipedia.org/wiki/Prandtl_number
 *
 * \param dynamicViscosity      Dynamic (absolute) viscosity over density.
 *                              http://en.wikipedia.org/wiki/Viscosity#Dynamic_viscosity [m^2/s]
 * \param heatCapacity          Heat capacity at constant pressure.
 *                              Specifies the energy change for a given temperature change [J / (kg K)]
 * \param thermalConductivity   Conductivity to heat. Specifies how well matter transfers energy without moving. [W/(m K)]
 * \return                      The Prandtl Number as calculated from the input parameters.
 */
static Scalar prandtlNumber(const Scalar dynamicViscosity,
                            const Scalar heatCapacity,
                            const Scalar thermalConductivity)
{
    return dynamicViscosity * heatCapacity / thermalConductivity;
}

/*!
 * \brief   Calculate the Nusselt Number [-] (Nu).
 *
 *          The Nusselt Number is a measure for the relation of convective- to conductive heat exchange.
 *
 *          The Nusselt number is defined as Nu = h d / k,
 *          with h= heat transfer coefficient, d=characteristic length, k=heat conductivity(stagnant).
 *          However, the heat transfer coefficient from one phase to another is typically not known.
 *          Therefore, Nusselt numbers are usually given as *empirical* Nu(Reynolds, Prandtl) for a given flow
 *          field --forced convection-- and *empirical* Nu(Rayleigh, Prandtl) for flow caused by temperature
 *          differences --free convection--. The fluid characteristics enter via the Prandtl number.
 *
 *          This function implements an *empirical* correlation for the case of porous media flow
 *          (packed bed flow as the chemical engineers call it).
 *
 *          source for Nusselt number definition: http://en.wikipedia.org/wiki/Nusselt_number
 *          source for further empirical correlations for Nusselt Numbers:
 *          VDI-Gesellschaft, VDI-Waermeatlas, VDI-Verlag Duesseldorf, 2006
 *
 * \param reynoldsNumber    Dimensionless number relating inertial and viscous forces [-].
 * \param prandtlNumber     Dimensionless number relating viscosity and thermal diffusivity (temperaturleitfaehigkeit) [-].
 * \param porosity          The fraction of the porous medium which is void space.
 * \param formulation       Switch for deciding which parametrization of the Nusselt number is to be used.
 *                          Set via the property NusseltFormulation.
 * \return                  The Nusselt number as calculated from the input parameters [-].
 */
static Scalar nusseltNumberForced(const Scalar reynoldsNumber,
                                  const Scalar prandtlNumber,
                                  const Scalar porosity,
                                  NusseltFormulation formulation)
{
    if (formulation == NusseltFormulation::dittusBoelter){
       /* example: very common and simple case: flow straight circular pipe, only convection (no boiling),
        * 10000<Re<120000, 0.7<Pr<120, far from pipe entrance, smooth surface of pipe ...
        * Dittus, F.W and Boelter, L.M.K, Heat Transfer in Automobile Radiators of the Tubular Type,
        * Publications in Engineering, Vol. 2, pages 443-461, 1930
        */
       using std::pow;
       return 0.023 * pow(reynoldsNumber, 0.8) * pow(prandtlNumber,0.33);
    }

    else if (formulation == NusseltFormulation::WakaoKaguei){
        /* example: flow through porous medium *single phase*, fit to many different data
         * Wakao and Kaguei, Heat and mass Transfer in Packed Beds, Gordon and Breach Science Publishers, page 293
         */
        using std::pow;
        return 2. + 1.1 * pow(prandtlNumber,(1./3.)) * pow(reynoldsNumber, 0.6);
    }

    else if (formulation == NusseltFormulation::VDI){
       /* example: VDI Waermeatlas 10. Auflage 2006, flow in packed beds, page Gj1, see also other sources and limitations therein.
        * valid for 0.1<Re<10000, 0.6<Pr/Sc<10000, packed beds of perfect spheres.
        *
        */
        using std::sqrt;
        using std::pow;
        using Dune::power;
        Scalar numerator    = 0.037 * pow(reynoldsNumber,0.8) * prandtlNumber ;
        Scalar reToMin01    = pow(reynoldsNumber,-0.1);
        Scalar prTo23       = pow(prandtlNumber, (2./3. ) ) ; // MIND THE pts! :-( otherwise the integer exponent version is chosen
        Scalar denominator  = 1+ 2.443 * reToMin01 * (prTo23 -1.) ;

        Scalar nusseltTurbular       = numerator / denominator;
        Scalar nusseltLaminar        = 0.664 * sqrt(reynoldsNumber) * pow(prandtlNumber, (1./3.) );
        Scalar nusseltSingleSphere   = 2 + sqrt( power(nusseltLaminar,2) + power(nusseltTurbular,2));

        Scalar funckyFactor           = 1 + 1.5 * (1.-porosity); // for spheres of same size
        Scalar nusseltNumber          = funckyFactor * nusseltSingleSphere  ;

        return nusseltNumber;
    }

    else {
        DUNE_THROW(Dune::NotImplemented, "wrong index");
    }
}


/*!
 * \brief   Calculate the Schmidt Number [-] (Sc).
 *
 *          The Schmidt Number is a measure for the relation of viscosity and mass diffusivity.
 *
 *          It is defined as
 *          \f[
 *          \textnormal{Sc}= \frac{\nu}{D} = \frac{\mu}{\rho D}\, ,
 *          \f]
 *          with kinematic viscosity\f$\nu\f$, diffusion coefficient \f$D\f$, dynamic viscosity
 *          \f$\mu\f$ and mass density\f$\rho\f$. Therefore, Sc is a material specific property
 *          (i.e.: not a function of flow directly but only of temperature, pressure and fluid).
 *
 *          source for Schmidt number definition: http://en.wikipedia.org/wiki/Schmidt_number
 *
 * \param dynamicViscosity      Dynamic (absolute) viscosity over density.
 *                              http://en.wikipedia.org/wiki/Viscosity#Dynamic_viscosity [m^2/s]
 * \param massDensity           Mass density of the considered phase. [kg / m^3]
 * \param diffusionCoefficient  Measure for how well a component can move through a phase due to a concentration gradient. [m^2/s]
 * \return                      The Schmidt Number as calculated from the input parameters.
 */
static Scalar schmidtNumber(const Scalar dynamicViscosity,
                            const Scalar massDensity,
                            const Scalar diffusionCoefficient)
{
    return dynamicViscosity  / (massDensity * diffusionCoefficient);
}

/*!
 * \brief   Calculate the Sherwood Number [-] (Sh).
 *
 *          The Sherwood Number is a measure for the relation of convective- to diffusive mass exchange.
 *
 *          The Sherwood number is defined as Sh = K L/D,
 *          with K= mass transfer coefficient, L=characteristic length, D=mass diffusivity (stagnant).
 *
 *          However, the mass transfer coefficient from one phase to another is typically not known.
 *          Therefore, Sherwood numbers are usually given as *empirical* Sh(Reynolds, Schmidt) for a given flow
 *          field (and fluid).
 *
 *          Often, even the Sherwood number is not known. By means of the Chilton-Colburn analogy it can be deduced
 *          from the Nusselt number. According to the Chilton-Colburn analogy in a known Nusselt correltion one
 *          basically replaces Pr with Sc and Nu with Sh. For some very special cases this is actually accurate.
 *          (Source: Course Notes, Waerme- und Stoffuebertragung, Prof. Hans Hasse, Uni Stuttgart)
 *
 *          This function implements an *empirical* correlation for the case of porous media flow
 *          (packed bed flow as the chemical engineers call it).
 *
 *          source for Sherwood number definition: http://en.wikipedia.org/wiki/Sherwood_number
 *
 * \param schmidtNumber     Dimensionless number relating viscosity and mass diffusivity [-].
 * \param reynoldsNumber    Dimensionless number relating inertial and viscous forces [-].
 * \param formulation       Switch for deciding which parametrization of the Sherwood number is to be used.
 *                          Set via the property SherwoodFormulation.
 * \return                  The Nusselt number as calculated from the input parameters [-].
 */

static Scalar sherwoodNumber(const Scalar reynoldsNumber,
                             const Scalar schmidtNumber,
                             SherwoodFormulation formulation)
{
    if (formulation == SherwoodFormulation::WakaoKaguei){
        /* example: flow through porous medium *single phase*
         * Wakao and Kaguei, Heat and mass Transfer in Packed Beds, Gordon and Breach Science Publishers, page 156
         */
        using std::cbrt;
        using std::pow;
        return 2. + 1.1 * cbrt(schmidtNumber) * pow(reynoldsNumber, 0.6);
    }

    else {
        DUNE_THROW(Dune::NotImplemented, "wrong index");
    }
}


/*!
 * \brief   Calculate the thermal diffusivity alpha [m^2/s].
 *
 *          The thermal diffusivity is a measure for how fast "temperature (not heat!) spreads".
 *          It is defined as \f$\alpha = \frac{k}{\rho c_p}\f$
 *          with \f$\alpha\f$: \f$k\f$: thermal conductivity [W/mK], \f$\rho\f$: density [kg/m^3],
 *          \f$c_p\f$: cpecific heat capacity at constant pressure [J/kgK].
 *
 *          Source for thermal diffusivity definition: http://en.wikipedia.org/wiki/Thermal_diffusivity
 *
 * \param   thermalConductivity A material property defining how well heat is transported via conduction [W/(mK)].
 * \param   phaseDensity        The density of the phase for which the thermal diffusivity is to be calculated [kg/m^3].
 * \param   heatCapacity        A measure for how a much a material changes temperature for a given change of energy (at p=const.) [J/(kgm^3)].
 * \return  The thermal diffusivity as calculated from the input parameters [m^2/s].
 */
static Scalar thermalDiffusivity(const Scalar & thermalConductivity ,
                                  const Scalar & phaseDensity ,
                                  const Scalar & heatCapacity)
{
    return thermalConductivity / (phaseDensity * heatCapacity);
}

}; // end class DimensionlessNumbers

} // end namespace Dumux

#endif
