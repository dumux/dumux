// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup PhaseFieldModels
 * \brief A logarithmic (Flory-Huggins) free-energy density.
 */
#ifndef DUMUX_PHASEFIELD_FREEENERGY_LOGARITHMIC_HH
#define DUMUX_PHASEFIELD_FREEENERGY_LOGARITHMIC_HH

#include <cmath>
#include <algorithm>

namespace Dumux::PhaseField::FreeEnergy {

/*!
 * \ingroup PhaseFieldModels
 * \brief Logarithmic (Flory-Huggins) free-energy density
 *        \f$ f(c) = \theta \left[ c \ln c + (1-c) \ln(1-c) \right] + \theta_c\, c (1-c) \f$.
 *
 * Approaches the double well as \f$ \theta \to 0 \f$ relative to \f$ \theta_c \f$.
 * The logarithm is only defined on \f$ (0,1) \f$; the argument is clamped to
 * \f$ [\epsilon, 1-\epsilon] \f$ to keep the density and its derivatives finite.
 */
template<class Scalar>
class Logarithmic
{
public:
    Logarithmic(Scalar theta, Scalar thetaC)
    : theta_(theta), thetaC_(thetaC)
    {}

    //! the free-energy density f(c)
    Scalar value(Scalar c) const
    {
        c = clamp_(c);
        using std::log;
        return theta_*(c*log(c) + (1.0-c)*log(1.0-c)) + thetaC_*c*(1.0-c);
    }

    //! the derivative df/dc
    Scalar derivative(Scalar c) const
    {
        c = clamp_(c);
        using std::log;
        return theta_*log(c/(1.0-c)) + thetaC_*(1.0-2.0*c);
    }

    //! the second derivative d^2f/dc^2
    Scalar secondDerivative(Scalar c) const
    {
        c = clamp_(c);
        return theta_/(c*(1.0-c)) - 2.0*thetaC_;
    }

private:
    static Scalar clamp_(Scalar c)
    {
        static constexpr Scalar eps = 1e-9;
        using std::clamp;
        return clamp(c, eps, Scalar(1.0)-eps);
    }

    Scalar theta_, thetaC_;
};

} // end namespace Dumux::PhaseField::FreeEnergy

#endif
