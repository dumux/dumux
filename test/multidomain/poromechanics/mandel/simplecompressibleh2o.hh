// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup Components
 * \brief A simple implementation of pure water
 */
#ifndef DUMUX_SIMPLE_COMPRESSIBLE_H2O_HH
#define DUMUX_SIMPLE_COMPRESSIBLE_H2O_HH

#include <cmath>

#include <dumux/common/parameters.hh>
#include <dumux/material/components/base.hh>
#include <dumux/material/components/liquid.hh>

namespace Dumux::Components {

/*!
 * \ingroup Components
 * \brief A simple implementation of pure compressible water.
 *        The compressibility is defined by the bulk modulus.
 *
 * \tparam Scalar The type used for scalar values
 */
template <class Scalar>
class SimpleCompressibleH2O
: public Components::Base<Scalar, SimpleCompressibleH2O<Scalar> >
, public Components::Liquid<Scalar, SimpleCompressibleH2O<Scalar> >
{
public:
    /*!
     * \brief The density of pure water at a given pressure and temperature \f$\mathrm{kg/m^3}\f$.
     * \param temperature temperature of component in \f$\mathrm{K}\f$
     * \param pressure pressure of component in \f$\mathrm{Pa}\f$
     */
    static Scalar liquidDensity(Scalar temperature, Scalar pressure)
    {
        static const Scalar bulkModulus = getParam<Scalar>("Problem.WaterBulkModulus", 2.0e9);
        static const Scalar pInit = 1.0e5;
        static const Scalar rhoInit = 1000.0;
        Scalar deltaP = pressure - pInit;
        Scalar rho = rhoInit+ rhoInit*deltaP/bulkModulus;
        return rho;
    }

    /*!
     * \brief The dynamic viscosity \f$\mathrm{Pa*s}\f$ of pure water.
     * \param temperature temperature of component in \f$\mathrm{K}\f$
     * \param pressure pressure of component in \f$\mathrm{Pa}\f$
     */
    static Scalar liquidViscosity(Scalar temperature, Scalar pressure)
    {
        static const Scalar mu = getParam<Scalar>("Problem.Viscosity", 1e-3);
        return mu;
    }

};

template <class Scalar>
struct IsAqueous<SimpleCompressibleH2O<Scalar>> : public std::true_type {};

} // end namespace Dumux::Components

#endif
