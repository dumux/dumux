// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup Binarycoefficients
 * \brief Binary coefficients for water and air.
 */
#ifndef DUMUX_BINARY_COEFF_H2O_AIR_HH
#define DUMUX_BINARY_COEFF_H2O_AIR_HH

#include <cmath>

#include <dumux/material/components/air.hh>
#include <dumux/material/binarycoefficients/h2o_ar.hh>
#include <dumux/material/binarycoefficients/h2o_co2.hh>
#include <dumux/material/binarycoefficients/h2o_n2.hh>
#include <dumux/material/binarycoefficients/h2o_o2.hh>
#include <dumux/common/tag.hh>

namespace Dumux::BinaryCoeff {
namespace Detail::H2O_Air{

//! Finsterle implementation for Henry coefficient
struct FinsterleLaw : public Utility::Tag<FinsterleLaw> {
    static std::string name() { return "Finsterle"; }
};

//! Mixture implementation for Henry coefficient
struct IdealMixtureLaw : public Utility::Tag<IdealMixtureLaw> {
    static std::string name() { return "Mixture"; }
};

//! Default rule (for backwards compatibility, will change to Mixture after 3.11)
using DefaultRule = FinsterleLaw;


/*!
 * \ingroup Finsterle_Implementation
 * \brief Henry coefficient implementation using Finsterle rule
 *
 * This implements the Henry coefficient calculation based on the Finsterle approach.
 * See:
 * Stefan Finsterle (1993, page 33 Formula (2.9)) \cite finsterle1993
 * (fitted to data from Tchobanoglous & Schroeder, 1985 \cite tchobanoglous1985)
 *
 * \param temperature the temperature \f$\mathrm{[K]}\f$
 * \return Henry coefficient \f$\mathrm{[Pa]}\f$
 */
template<typename Scalar>
Scalar henry(Scalar temperature, FinsterleLaw) {
    using std::exp;
    Scalar r = (0.8942 + 1.47 * exp(-0.04394 * (temperature - 273.15))) * 1.E-10;
    return 1. / r;
}

/*!
 * \ingroup TchobanoglousCubic_Implementation
 * \brief Henry coefficient implementation using a cubic function to data from
 * Tchobanoglous & Schroeder, 1985 \cite tchobanoglous1985
 * The data covers the interval from 0°C to 60°C. Outside, the quality of the extrapolation
 * is unknown.
 *
 * \param temperature the temperature \f$\mathrm{[K]}\f$
 * \return Henry coefficient \f$\mathrm{[Pa]}\f$
 */
template <class Scalar>
Scalar henry(Scalar temperature, TchobanoglousCubicImplementation)
{
    Scalar t = temperature - 273.15; // convert from K to °C
    Scalar r = (-5.55556E-6 * t*t*t + 0.000895238 * t*t - 0.0556825 * t + 2.30524) * 1e-5;
    r /= 101325; // convert from atm to Pa

    return 1./r;
}

/*!
 * \ingroup Mixture_Implementation
 * \brief Henry coefficient implementation using Mixture rule
 *
 * This implementation calculates the effective Henry coefficient for air as a mixture
 * based on the individual Henry coefficients of its components (N₂, O₂, Ar, CO₂).
 *
 * For a gas mixture like air in water, the effective Henry coefficient can be derived
 * from the individual Henry coefficients by:
 * \f[K_{H,air} = \left(\sum_i \frac{y_{i,air}}{K_{H,i}}\right)^{-1}\f]
 *
 * Where:
 * - \f$y_{i,air}\f$ is the mole fraction of component \f$i\f$ in air
 * - \f$K_{H,i}\f$ is the Henry coefficient for component \f$i\f$ in water
 *
 * The implementation considers the four main components of dry-air with their respective
 * mole fractions:
 *
 * \copydetails Dumux::Components::Air::composition::dryMoleFraction
 *
 * \param temperature the temperature \f$\mathrm{[K]}\f$
 * \return Henry coefficient for air in water \f$\mathrm{[Pa]}\f$
 */
template<typename Scalar>
Scalar henry(Scalar temperature, IdealMixtureLaw) {
    static constexpr Scalar yO2_air = Dumux::Components::Air<Scalar>::composition::dryMoleFraction::O2;
    static constexpr Scalar yN2_air = Dumux::Components::Air<Scalar>::composition::dryMoleFraction::N2;
    static constexpr Scalar yCO2_air = Dumux::Components::Air<Scalar>::composition::dryMoleFraction::CO2;
    static constexpr Scalar yAr_air = Dumux::Components::Air<Scalar>::composition::dryMoleFraction::Ar;

    Scalar kh_O2 = H2O_O2::henry(temperature);
    Scalar kh_N2 = H2O_N2::henry(temperature);
    Scalar kh_CO2 = H2O_CO2::henry(temperature);
    Scalar kh_Ar = H2O_AR::henry(temperature);

    Scalar r = (yO2_air / kh_O2 + yN2_air / kh_N2 + yCO2_air / kh_CO2 + yAr_air / kh_Ar);
    return 1 / r;
}
} // end namespace Detail::H2O_Air

/*!
 * \ingroup Binarycoefficients
 * \brief Binary coefficients for water and air.
 */
        class H2O_Air
        {
        public:

            /*!
            * \brief Default Henry coefficient \f$\mathrm{[Pa]}\f$ for air in liquid water.
            * \note The default implementation currently uses the Finsterle approach but will change to
            * the mixture approach after version 3.11. Use explicit tag dispatch to choose a specific implementation.
            *
            * \section finsterle_section Finsterle Law Implementation
            * \copydetails Dumux::BinaryCoeff::Detail::H2O_Air::henry(Scalar,FinsterleLaw)
            *
            * \param temperature the temperature \f$\mathrm{[K]}\f$
            * \return Henry coefficient \f$\mathrm{[Pa]}\f$
            */
            template <typename Scalar>
            static Scalar henry(Scalar temperature)
            {
                // Warning about future changes (if needed)
                [[maybe_unused]] static bool _ = []() {
                    std::cout << "Warning: H2O_Air::henry(T) default implementation will change from Finsterle to Mixture implementation after 3.11. Use explicit tag to disable this warning and choose a specific implementation.\n";
                    return true;
                }();


                return henry<Scalar>(temperature, Detail::H2O_Air::DefaultRule{});
            }

            /*!
            * \brief Henry coefficient \f$\mathrm{[Pa]}\f$ for air in liquid water.
            *
            * Two implementations are available:
            *
            * \section finsterle_section Finsterle Law Implementation
            * \copydetails Dumux::BinaryCoeff::Detail::H2O_Air::henry(Scalar,FinsterleLaw)
            *
            * \section mixture_section Mixture Law Implementation
            * \copydetails Dumux::BinaryCoeff::Detail::H2O_Air::henry(Scalar,IdealMixtureLaw)
            *
            * \note The default implementation currently uses the Finsterle approach but will change to
            * the mixture approach after version 3.11. Use explicit tag dispatch to choose a specific implementation.
            *
            * \param temperature the temperature \f$\mathrm{[K]}\f$
            * \return Henry coefficient \f$\mathrm{[Pa]}\f$
            */
            template <typename Scalar, typename Rule>
            static Scalar henry(Scalar temperature, Rule rule)
            {
                return henry<Scalar>(temperature, rule);
            }

            /*!
             * \brief Binary diffusion coefficient \f$\mathrm{[m^2/s]}\f$ for molecular water and air
             *
             * \param temperature the temperature \f$\mathrm{[K]}\f$
             * \param pressure the phase pressure \f$\mathrm{[Pa]}\f$
             * Vargaftik: Tables on the thermophysical properties of liquids and gases.
             * John Wiley & Sons, New York, 1975. \cite vargaftik1975 <BR>
             * Walker, Sabey, Hampton: Studies of heat transfer and water migration in soils.
             * Dep. of Agricultural and Chemical Engineering, Colorado State University,
             * Fort Collins, 1981. \cite walker1981
             */
            template <class Scalar>
            static Scalar gasDiffCoeff(Scalar temperature, Scalar pressure)
            {
                const Scalar Theta = 1.8;
                const Scalar Daw = 2.13e-5;  /* reference value */
                const Scalar pg0 = 1.e5;     /* reference pressure */
                const Scalar T0 = 273.15;    /* reference temperature */
                Scalar Dgaw;

                using std::pow;
                Dgaw = Daw * (pg0 / pressure) * pow((temperature / T0), Theta);

                return Dgaw;
            }

            /*!
             * Lacking better data on water-air diffusion in liquids, we use at the
             * moment the diffusion coefficient of the air's main component nitrogen!!
             * \brief Diffusion coefficient \f$\mathrm{[m^2/s]}\f$ for molecular nitrogen in liquid water.
             *
             * \param temperature the temperature \f$\mathrm{[K]}\f$
             * \param pressure the phase pressure \f$\mathrm{[Pa]}\f$
             *
             * The empirical equations for estimating the diffusion coefficient in
             * infinite solution which are presented in Reid, 1987 all show a
             * linear dependency on temperature. We thus simply scale the
             * experimentally obtained diffusion coefficient of Ferrell and
             * Himmelblau by the temperature.
             *
             * See:
             * R. Reid et al. (1987, pp. 599) \cite reid1987 <BR>
             * R. Ferrell, D. Himmelblau (1967, pp. 111-115) \cite ferrell1967
             */
            template <class Scalar>
            static Scalar liquidDiffCoeff(Scalar temperature, Scalar pressure)
            {
                const Scalar Texp = 273.15 + 25; // [K]
                const Scalar Dexp = 2.01e-9; // [m^2/s]
                return Dexp * temperature / Texp;
            }
        };
} // end namespace Dumux::BinaryCoeff

#endif
